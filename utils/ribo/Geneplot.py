#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-11
# Version: 0.2.8-dev.005
# Function: Provide RiboParser core functions for gene-level RPF/RNA density plotting.
# Input: RPF/RNA density file in JSONL or TXT format and one target gene/transcript ID and optional genome annotation.
# Output: Gene-level RPF/RNA density profile, optional bedGraph-like track table, and IGV-like gene plot.

"""Gene-level RPF/RNA density plotting utilities for RiboParser.

This module draws IGV-like RPF density profiles for a target gene or
transcript. Current JSONL density files are read by streaming transcript records
so that plotting one gene does not require expanding the whole density file into
memory. Legacy TXT density files are read in chunks. When an external genePred
or RiboParser norm annotation is provided, TXT profiles can also be projected to
genome coordinates.
"""

from __future__ import annotations

import re
from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, Callable, Iterable

from . import RPFs

import matplotlib

matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FRAME_ORDER = ["0", "1", "2"]
FRAME_COLORS = {
    "0": "#E64B35",
    "1": "#4DBBD5",
    "2": "#00A087",
}
LINE_COLOR = "#3C5488"
RNA_COLOR = "#3C5488"
UTR_COLOR = "#BDBDBD"
CDS_COLOR = "#4DBBD5"
INTRON_COLOR = "#4A4A4A"
RPM_SCALE = 1_000_000.0
CHUNK_SIZE = 500_000
ANNOTATION_CACHE: dict[str, list[dict[str, Any]]] = {}
PLOT_TRANSFORM_LABELS = {
    "none": "",
    "sqrt": "sqrt",
    "log": "log1p",
    "log1p": "log1p",
    "log2": "log2(x + 1)",
    "log10": "log10(x + 1)",
}
BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]


@dataclass
class TargetMeta:
    """Metadata for the selected target transcript."""

    target_id: str
    gene_id: str
    transcript_id: str
    name: str
    chrom: str | None = None
    strand: str | None = None
    candidate_count: int = 1
    selected_by: str = "first"


@dataclass
class GenomeMap:
    """Genome mapping information for one transcript."""

    chrom: str
    strand: str
    exons: list[tuple[int, int]]
    tx_to_genome: list[int]


class Geneplot(object):
    """Draw an IGV-like RPF density plot for one target gene or transcript.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Geneplot``.
    """

    def __init__(self, args):
        # Input and output.
        self.rpf = args.rpf
        self.output = args.output
        self.target = args.target
        self.id_type = args.id_type
        self.sample_arg = args.sample
        self.annotation = getattr(args, "annotation", None)

        # Plot parameters.
        self.coordinate = args.coordinate
        self.mode = args.mode
        self.data_type = str(getattr(args, "data_type", "ribo")).lower()
        self.frame = str(args.frame)
        if self.data_type == "rna":
            self.frame = "all"
        self.norm = bool(args.normal)
        self.plot_transform = args.plot_transform
        self.utr_gray = bool(args.utr_gray)
        self.line_width = float(args.line_width)
        self.bar_width = float(args.bar_width)
        self.spike_clip = bool(args.spike_clip)
        self.clip_quantile = float(args.clip_quantile)
        self.clip_value = args.clip_value
        self.export_track = bool(args.export_track)
        self.figure_width = float(args.width)
        self.per_sample_height = float(args.per_sample_height)
        self.structure_height = float(args.structure_height)
        self.max_height = float(args.max_height)
        self.y_max = args.y_max
        self.x_margin = float(args.x_margin)
        self.intron_scale = float(args.intron_scale)
        self.select_transcript = args.select_transcript
        self.output_format = args.output_format
        self.dpi = int(args.dpi)
        self.title = args.title
        self.font_size = float(args.font_size)
        self.title_size = float(args.title_size) if args.title_size is not None else self.font_size + 2.0
        self.label_size = float(args.label_size) if args.label_size is not None else self.font_size
        self.tick_size = float(args.tick_size) if args.tick_size is not None else max(self.font_size - 1.0, 1.0)
        self.sample_label_size = float(args.sample_label_size) if args.sample_label_size is not None else self.font_size
        self.legend_size = float(args.legend_size) if args.legend_size is not None else max(self.font_size - 1.0, 1.0)

        # Imported data and outputs.
        self.file_format: str | None = None
        self.resolved_coordinate: str | None = None
        self.sample_name: list[str] = []
        self.total_rpf_num: dict[str, float] = {}
        self.target_meta: TargetMeta | None = None
        self.genome_map: GenomeMap | None = None
        self.profile: pd.DataFrame | None = None
        self.feature_segments: pd.DataFrame | None = None
        self.output_files = OrderedDict()
        self.warnings: list[str] = []

    # ------------------------------------------------------------------
    # Input helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_name(name: str) -> str:
        """Return a safe string for output file names."""
        return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name)).strip("_") or "target"

    @staticmethod
    def _parse_sample_names(value: str | None) -> list[str] | None:
        """Parse a comma-separated sample list."""
        if value is None or str(value).strip() == "":
            return None
        return [item.strip() for item in str(value).split(",") if item.strip()]

    @staticmethod
    def _to_int(value: Any, default: int = 0) -> int:
        """Safely convert one value to an integer."""
        try:
            return int(value)
        except (TypeError, ValueError):
            try:
                return int(float(value))
            except (TypeError, ValueError):
                return default

    @staticmethod
    def _to_float(value: Any, default: float = 0.0) -> float:
        """Safely convert one value to a float."""
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

    def _warn(self, message: str) -> None:
        """Print and store one warning message."""
        self.warnings.append(message)
        print("Warning: " + message, flush=True)

    @staticmethod
    def _record_meta(record: dict[str, Any]) -> dict[str, str]:
        """Extract common target identifiers from one JSON transcript record."""
        annotation = record.get("annotation") if isinstance(record.get("annotation"), dict) else {}
        name = str(record.get("name") or record.get("transcript_id") or annotation.get("transcript_id") or "NA")
        transcript_id = str(record.get("transcript_id") or annotation.get("transcript_id") or name)
        gene_id = str(record.get("gene_id") or annotation.get("gene_id") or name)
        return {"name": name, "transcript_id": transcript_id, "gene_id": gene_id}

    def _record_matches_target(self, record: dict[str, Any]) -> bool:
        """Return whether one JSON record matches the requested target."""
        meta = self._record_meta(record)
        target = str(self.target)

        if self.id_type == "transcript":
            return target in {meta["name"], meta["transcript_id"]}
        if self.id_type == "gene":
            return target == meta["gene_id"]
        return target in {meta["name"], meta["transcript_id"], meta["gene_id"]}

    def _record_length_score(self, record: dict[str, Any]) -> int:
        """Return transcript-length score for candidate selection."""
        trim = record.get("trim") if isinstance(record.get("trim"), dict) else {}
        annotation = record.get("annotation") if isinstance(record.get("annotation"), dict) else {}
        return max(
            self._to_int(trim.get("trim_length_nt"), 0),
            self._to_int(trim.get("codon_count"), 0) * 3,
            self._to_int(annotation.get("transcript_length"), 0),
        )

    def _record_density_score(self, record: dict[str, Any], sample_names: list[str]) -> float:
        """Return total RPF count across selected samples for candidate selection."""
        score = 0.0
        for sample in sample_names:
            score += self._density_sum(record, sample)
        return score

    def _select_json_record(self, candidates: list[dict[str, Any]], sample_names: list[str]) -> dict[str, Any]:
        """Select one transcript record when a gene ID matches multiple isoforms."""
        if not candidates:
            raise ValueError("No target gene/transcript was found in the JSON density file.")

        if len(candidates) == 1:
            return candidates[0]

        if self.select_transcript == "first":
            selected = candidates[0]
        elif self.select_transcript == "highest":
            selected = max(candidates, key=lambda record: self._record_density_score(record, sample_names))
        elif self.select_transcript == "longest":
            selected = max(candidates, key=self._record_length_score)
        else:
            raise ValueError("Unsupported transcript selection method: {method}".format(method=self.select_transcript))

        selected_meta = self._record_meta(selected)
        self._warn(
            "Target matched {count} transcript records; selected {transcript} by {method}.".format(
                count=len(candidates),
                transcript=selected_meta["transcript_id"],
                method=self.select_transcript,
            )
        )
        return selected

    @staticmethod
    def _density_sum(record: dict[str, Any], sample_name: str) -> float:
        """Return total density count for one sample in one JSON record."""
        sample_entry = record.get("samples", {}).get(sample_name)
        if not isinstance(sample_entry, dict):
            return 0.0
        density = sample_entry.get("density")
        if not isinstance(density, dict):
            return 0.0

        encoding = str(density.get("encoding", RPFs.SPARSE_ENCODING))
        if encoding in {RPFs.DENSE_ENCODING, "dense"}:
            return float(sum(density.get("f0", [])) + sum(density.get("f1", [])) + sum(density.get("f2", [])))
        return float(sum(density.get("count", [])))

    @staticmethod
    def _parse_int_list(value: Any) -> list[int]:
        """Parse a list-like coordinate value."""
        if value is None:
            return []
        if isinstance(value, str):
            value = value.strip().rstrip(",")
            if not value:
                return []
            return [int(float(item)) for item in re.split(r"[,;\s]+", value) if item != ""]
        if isinstance(value, Iterable):
            return [int(float(item)) for item in value]
        return []

    @classmethod
    def _extract_exons(cls, mapping: dict[str, Any]) -> list[tuple[int, int]]:
        """Extract exon intervals as 0-based half-open coordinates."""
        exon_values = None
        for key in ("exons", "exon", "exon_intervals", "blocks"):
            if key in mapping:
                exon_values = mapping.get(key)
                break

        exons: list[tuple[int, int]] = []
        if isinstance(exon_values, list):
            for item in exon_values:
                if isinstance(item, dict):
                    start = item.get("start", item.get("tx_start", item.get("exon_start")))
                    end = item.get("end", item.get("tx_end", item.get("exon_end")))
                elif isinstance(item, (list, tuple)) and len(item) >= 2:
                    start, end = item[0], item[1]
                else:
                    continue
                start_i = cls._to_int(start, -1)
                end_i = cls._to_int(end, -1)
                if start_i >= 0 and end_i > start_i:
                    exons.append((start_i, end_i))

        if not exons:
            start_keys = ("exon_starts", "exonStarts", "exon_start", "block_starts", "blockStarts")
            end_keys = ("exon_ends", "exonEnds", "exon_end", "block_ends", "blockEnds")
            starts: list[int] = []
            ends: list[int] = []
            for key in start_keys:
                starts = cls._parse_int_list(mapping.get(key))
                if starts:
                    break
            for key in end_keys:
                ends = cls._parse_int_list(mapping.get(key))
                if ends:
                    break
            exons = [(start, end) for start, end in zip(starts, ends) if end > start]

        return sorted(set(exons), key=lambda item: item[0])

    def _extract_genome_map(self, record: dict[str, Any]) -> GenomeMap | None:
        """Extract genome mapping from one JSON transcript record when available."""
        mapping = record.get("genome_mapping")
        annotation = record.get("annotation") if isinstance(record.get("annotation"), dict) else {}

        if not isinstance(mapping, dict):
            mapping = {}
        merged_mapping = dict(annotation)
        merged_mapping.update(mapping)

        chrom = merged_mapping.get("chrom") or merged_mapping.get("chromosome") or merged_mapping.get("seqname")
        strand = merged_mapping.get("strand") or merged_mapping.get("gene_strand") or "+"
        exons = self._extract_exons(merged_mapping)

        if chrom is None or not exons:
            return None

        strand = str(strand)
        if strand not in {"+", "-"}:
            strand = "+"

        tx_to_genome: list[int] = []
        if strand == "+":
            for start, end in exons:
                tx_to_genome.extend(range(start + 1, end + 1))
        else:
            for start, end in reversed(exons):
                tx_to_genome.extend(range(end, start, -1))

        if not tx_to_genome:
            return None

        return GenomeMap(chrom=str(chrom), strand=strand, exons=exons, tx_to_genome=tx_to_genome)


    # ------------------------------------------------------------------
    # External annotation helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _canonical_key(value: Any) -> str:
        """Return a normalized key used for flexible annotation parsing."""
        return re.sub(r"[^0-9a-z]+", "", str(value).strip().lower())

    @classmethod
    def _is_annotation_header(cls, fields: list[str]) -> bool:
        """Return whether the first annotation line looks like a header."""
        known_keys = {
            "name",
            "transcript",
            "transcriptid",
            "transcriptname",
            "gene",
            "geneid",
            "genename",
            "name2",
            "chrom",
            "chromosome",
            "seqname",
            "strand",
            "txstart",
            "txend",
            "cdsstart",
            "cdsend",
            "exonstarts",
            "exonends",
            "blockstarts",
            "blockends",
            "blocksizes",
        }
        normalized = {cls._canonical_key(field.lstrip("#")) for field in fields}
        return bool(normalized & known_keys)

    @staticmethod
    def _split_annotation_line(line: str) -> list[str]:
        """Split one annotation line by tab or repeated whitespace."""
        return re.split(r"\t|\s+", line.strip())

    @staticmethod
    def _clean_annotation_value(value: Any) -> Any:
        """Clean pandas values from annotation rows."""
        if value is None:
            return None
        try:
            if pd.isna(value):
                return None
        except (TypeError, ValueError):
            pass
        if isinstance(value, str):
            value = value.strip()
            return value if value != "" else None
        return value

    def _read_annotation_rows(self) -> list[dict[str, Any]]:
        """Read a genePred or RiboParser norm annotation table as row dictionaries."""
        if not self.annotation:
            return []

        annotation_file = str(self.annotation)
        if annotation_file in ANNOTATION_CACHE:
            return ANNOTATION_CACHE[annotation_file]

        first_fields: list[str] | None = None
        with open(annotation_file, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    maybe_fields = self._split_annotation_line(line.lstrip("#"))
                    if self._is_annotation_header(maybe_fields):
                        first_fields = [field.lstrip("#") for field in maybe_fields]
                        break
                    continue
                first_fields = self._split_annotation_line(line)
                break

        if first_fields is None:
            raise ValueError("Annotation file is empty: {file}".format(file=annotation_file))

        has_header = self._is_annotation_header(first_fields)
        if has_header:
            if first_fields[0].startswith("#"):
                first_fields[0] = first_fields[0].lstrip("#")
            data = pd.read_csv(
                annotation_file,
                sep=r"\t|\s+",
                engine="python",
                comment="#",
                header=0,
                dtype=str,
            )
            # If the header line started with '#', pandas skipped it. Re-read
            # with explicit column names in that common genePred-header case.
            if data.empty or not ({self._canonical_key(col) for col in data.columns} & {"name", "transcriptid", "chrom"}):
                data = pd.read_csv(
                    annotation_file,
                    sep=r"\t|\s+",
                    engine="python",
                    comment="#",
                    header=None,
                    names=first_fields,
                    dtype=str,
                )
        else:
            data = pd.read_csv(
                annotation_file,
                sep=r"\t|\s+",
                engine="python",
                comment="#",
                header=None,
                dtype=str,
            )
            gene_pred_columns = [
                "name",
                "chrom",
                "strand",
                "txStart",
                "txEnd",
                "cdsStart",
                "cdsEnd",
                "exonCount",
                "exonStarts",
                "exonEnds",
                "score",
                "name2",
                "cdsStartStat",
                "cdsEndStat",
                "exonFrames",
            ]
            column_count = data.shape[1]
            if column_count <= len(gene_pred_columns):
                data.columns = gene_pred_columns[:column_count]
            else:
                extra_columns = ["extra_{idx}".format(idx=idx) for idx in range(column_count - len(gene_pred_columns))]
                data.columns = gene_pred_columns + extra_columns

        rows = []
        for row in data.to_dict(orient="records"):
            rows.append({str(key): self._clean_annotation_value(value) for key, value in row.items()})

        ANNOTATION_CACHE[annotation_file] = rows
        return rows

    @classmethod
    def _get_any_annotation_value(cls, row: dict[str, Any], keys: Iterable[str], default: Any = None) -> Any:
        """Return a row value using flexible column aliases."""
        lookup = {cls._canonical_key(key): value for key, value in row.items()}
        for key in keys:
            value = lookup.get(cls._canonical_key(key))
            if value is not None:
                return value
        return default

    @classmethod
    def _annotation_identifiers(cls, row: dict[str, Any]) -> dict[str, str]:
        """Extract gene and transcript identifiers from one annotation row."""
        transcript_id = cls._get_any_annotation_value(
            row,
            ["transcript_id", "transcript", "transcript_name", "name", "id"],
            None,
        )
        gene_id = cls._get_any_annotation_value(
            row,
            ["gene_id", "gene", "gene_name", "name2", "geneid"],
            None,
        )
        transcript_id = str(transcript_id) if transcript_id is not None else "NA"
        gene_id = str(gene_id) if gene_id is not None else transcript_id
        return {"transcript_id": transcript_id, "gene_id": gene_id, "name": transcript_id}

    @classmethod
    def _annotation_exons(cls, row: dict[str, Any]) -> list[tuple[int, int]]:
        """Extract 0-based half-open exon intervals from one annotation row."""
        mapping: dict[str, Any] = {}
        for source_key, target_key in (
            ("exonStarts", "exon_starts"),
            ("exon_starts", "exon_starts"),
            ("exonStart", "exon_starts"),
            ("blockStarts", "block_starts"),
            ("block_starts", "block_starts"),
            ("exonEnds", "exon_ends"),
            ("exon_ends", "exon_ends"),
            ("exonEnd", "exon_ends"),
            ("blockEnds", "block_ends"),
            ("block_ends", "block_ends"),
            ("blockSizes", "block_sizes"),
            ("block_sizes", "block_sizes"),
        ):
            value = cls._get_any_annotation_value(row, [source_key], None)
            if value is not None:
                mapping[target_key] = value

        exons = cls._extract_exons(mapping)
        if exons:
            return exons

        starts = cls._parse_int_list(mapping.get("block_starts"))
        sizes = cls._parse_int_list(mapping.get("block_sizes"))
        chrom_start = cls._get_any_annotation_value(row, ["chromStart", "txStart", "tx_start"], None)
        if starts and sizes and chrom_start is not None:
            chrom_start_i = cls._to_int(chrom_start, 0)
            exons = [(chrom_start_i + start, chrom_start_i + start + size) for start, size in zip(starts, sizes)]
            return [(start, end) for start, end in exons if end > start]

        tx_start = cls._get_any_annotation_value(row, ["txStart", "tx_start", "start"], None)
        tx_end = cls._get_any_annotation_value(row, ["txEnd", "tx_end", "end"], None)
        if tx_start is not None and tx_end is not None:
            start_i = cls._to_int(tx_start, -1)
            end_i = cls._to_int(tx_end, -1)
            if start_i >= 0 and end_i > start_i:
                return [(start_i, end_i)]

        return []

    @classmethod
    def _genome_map_from_annotation_row(cls, row: dict[str, Any]) -> tuple[GenomeMap | None, dict[str, str]]:
        """Build a GenomeMap and metadata from one external annotation row."""
        ids = cls._annotation_identifiers(row)
        chrom = cls._get_any_annotation_value(row, ["chrom", "chromosome", "seqname", "chr"], None)
        strand = cls._get_any_annotation_value(row, ["strand", "gene_strand"], "+")
        exons = cls._annotation_exons(row)

        if chrom is None or not exons:
            return None, ids

        strand = str(strand)
        if strand not in {"+", "-"}:
            strand = "+"

        tx_to_genome: list[int] = []
        if strand == "+":
            for start, end in sorted(exons, key=lambda item: item[0]):
                tx_to_genome.extend(range(start + 1, end + 1))
        else:
            for start, end in sorted(exons, key=lambda item: item[0], reverse=True):
                tx_to_genome.extend(range(end, start, -1))

        if not tx_to_genome:
            return None, ids

        return GenomeMap(chrom=str(chrom), strand=strand, exons=sorted(exons, key=lambda item: item[0]), tx_to_genome=tx_to_genome), ids

    def _annotation_row_matches(self, ids: dict[str, str], target_meta: dict[str, str]) -> bool:
        """Return whether one external annotation row matches the selected target."""
        transcript_ids = {ids["transcript_id"], ids.get("name", ids["transcript_id"])}
        gene_ids = {ids["gene_id"]}

        target_transcripts = {
            str(target_meta.get("transcript_id", "")),
            str(target_meta.get("name", "")),
            str(target_meta.get("target_id", "")),
        }
        target_genes = {str(target_meta.get("gene_id", "")), str(target_meta.get("target_id", ""))}

        if self.id_type == "transcript":
            return bool(transcript_ids & target_transcripts)
        if self.id_type == "gene":
            return bool(gene_ids & target_genes)
        return bool((transcript_ids & target_transcripts) or (gene_ids & target_genes))

    @staticmethod
    def _annotation_length_score(item: tuple[GenomeMap, dict[str, str]]) -> int:
        """Return transcript length score for external annotation selection."""
        genome_map, _ = item
        return len(genome_map.tx_to_genome)

    def _external_genome_map_for_target(self, target_meta: dict[str, str]) -> tuple[GenomeMap | None, dict[str, str] | None]:
        """Return external genome mapping for the selected target when available."""
        if not self.annotation:
            return None, None

        matched: list[tuple[GenomeMap, dict[str, str]]] = []
        for row in self._read_annotation_rows():
            genome_map, ids = self._genome_map_from_annotation_row(row)
            if genome_map is None:
                continue
            if self._annotation_row_matches(ids, target_meta):
                matched.append((genome_map, ids))

        if not matched:
            return None, None

        if len(matched) == 1 or self.select_transcript == "first":
            selected = matched[0]
        else:
            selected = max(matched, key=self._annotation_length_score)
            if self.select_transcript == "highest":
                self._warn("--select-transcript highest is not available for external annotation; selected the longest annotated transcript.")

        if len(matched) > 1:
            self._warn(
                "External annotation matched {count} transcript records; selected {transcript}.".format(
                    count=len(matched),
                    transcript=selected[1]["transcript_id"],
                )
            )

        return selected

    def _resolve_coordinate_mode(self, requested: str, file_format: str, genome_map: GenomeMap | None) -> str:
        """Resolve the final coordinate mode from user request and available data."""
        if requested == "transcript":
            return "transcript"

        if genome_map is not None:
            return "genome"

        if requested in {"genome", "both"}:
            if file_format == "txt":
                self._warn(
                    "TXT input has no embedded genome_mapping and no usable external annotation was found; "
                    "switched to transcript-coordinate mode."
                )
            else:
                self._warn(
                    "No usable genome mapping was found for this target; switched to transcript-coordinate mode."
                )
        return "transcript"

    # ------------------------------------------------------------------
    # Import JSON and TXT density
    # ------------------------------------------------------------------

    def import_rpf(self) -> None:
        """Import the target RPF density profile."""
        self.file_format = RPFs._detect_rpf_format(self.rpf)
        if self.file_format == "json":
            self._import_json_profile()
        elif self.file_format == "txt":
            self._import_txt_profile()
        else:
            raise ValueError("Unsupported RPF density format: {fmt}".format(fmt=self.file_format))

        if self.profile is None or self.profile.empty:
            raise ValueError("Target RPF profile is empty after import.")

        print(
            "Imported target {target}: format={fmt}, coordinate={coord}, samples={samples}, rows={rows:,}.".format(
                target=self.target,
                fmt=self.file_format,
                coord=self.resolved_coordinate,
                samples=len(self.sample_name),
                rows=len(self.profile),
            ),
            flush=True,
        )

    def _import_json_profile(self) -> None:
        """Import a target profile from compact JSONL by streaming records."""
        requested_samples = self._parse_sample_names(self.sample_arg)
        all_sample_names: list[str] | None = None
        selected_samples: list[str] | None = None
        totals: dict[str, float] = {}
        candidates: list[dict[str, Any]] = []
        record_count = 0

        for record in RPFs.iter_json_records(self.rpf):
            current_samples = [str(sample) for sample in record.get("samples", {}).keys()]
            if all_sample_names is None:
                all_sample_names = current_samples
                if requested_samples is None:
                    selected_samples = list(all_sample_names)
                else:
                    missing = sorted(set(requested_samples) - set(all_sample_names))
                    if missing:
                        raise ValueError("Requested sample(s) not found: " + ", ".join(missing))
                    selected_samples = list(requested_samples)
                totals = {sample: 0.0 for sample in selected_samples}
            elif current_samples != all_sample_names:
                raise ValueError("Inconsistent sample order in JSON density records.")

            assert selected_samples is not None
            for sample in selected_samples:
                totals[sample] += self._density_sum(record, sample)

            if self._record_matches_target(record):
                candidates.append(record)

            record_count += 1

        if all_sample_names is None or selected_samples is None:
            raise ValueError("No transcript records were found in JSON density file: {file}".format(file=self.rpf))

        selected_record = self._select_json_record(candidates, selected_samples)
        meta = self._record_meta(selected_record)
        genome_map = self._extract_genome_map(selected_record)
        external_ids = None
        if genome_map is None and self.annotation:
            genome_map, external_ids = self._external_genome_map_for_target(
                {
                    "target_id": str(self.target),
                    "gene_id": meta["gene_id"],
                    "transcript_id": meta["transcript_id"],
                    "name": meta["name"],
                }
            )

        resolved_coordinate = self._resolve_coordinate_mode(self.coordinate, "json", genome_map)

        if external_ids is not None:
            meta = {
                "gene_id": external_ids.get("gene_id", meta["gene_id"]),
                "transcript_id": external_ids.get("transcript_id", meta["transcript_id"]),
                "name": meta["name"],
            }

        self.sample_name = selected_samples
        self.total_rpf_num = totals
        self.genome_map = genome_map if resolved_coordinate == "genome" else None
        self.resolved_coordinate = resolved_coordinate
        self.target_meta = TargetMeta(
            target_id=str(self.target),
            gene_id=meta["gene_id"],
            transcript_id=meta["transcript_id"],
            name=meta["name"],
            chrom=genome_map.chrom if genome_map else None,
            strand=genome_map.strand if genome_map else None,
            candidate_count=len(candidates),
            selected_by=self.select_transcript,
        )
        self.profile = self._build_profile_from_json_record(selected_record)
        self.feature_segments = self._build_feature_segments_from_profile(self.profile)

        print("Scanned JSON transcript records: {count:,}.".format(count=record_count), flush=True)

    def _import_txt_profile(self) -> None:
        """Import a target profile from legacy TXT in chunks."""
        requested_samples = self._parse_sample_names(self.sample_arg)
        header = pd.read_csv(self.rpf, sep="\t", nrows=0).columns.tolist()
        missing_base = [column for column in BASE_COLUMNS if column not in header]
        if missing_base:
            raise ValueError("TXT density file is missing required column(s): " + ", ".join(missing_base))

        all_sample_names = []
        seen = set()
        for column in header:
            if column.endswith(("_f0", "_f1", "_f2")):
                sample = column[:-3]
                if sample not in seen:
                    all_sample_names.append(sample)
                    seen.add(sample)

        if requested_samples is None:
            selected_samples = all_sample_names
        else:
            missing = sorted(set(requested_samples) - set(all_sample_names))
            if missing:
                raise ValueError("Requested sample(s) not found: " + ", ".join(missing))
            selected_samples = requested_samples

        frame_columns = [f"{sample}_f{frame}" for sample in selected_samples for frame in range(3)]
        totals = {sample: 0.0 for sample in selected_samples}
        target_chunks: list[pd.DataFrame] = []
        chunk_count = 0
        row_count = 0

        for chunk in pd.read_csv(self.rpf, sep="\t", chunksize=CHUNK_SIZE):
            for sample in selected_samples:
                totals[sample] += float(chunk.loc[:, [f"{sample}_f0", f"{sample}_f1", f"{sample}_f2"]].sum().sum())

            target_chunk = chunk.loc[chunk["name"].astype(str) == str(self.target), BASE_COLUMNS + frame_columns].copy()
            if not target_chunk.empty:
                target_chunks.append(target_chunk)

            chunk_count += 1
            row_count += len(chunk)

        if not target_chunks:
            raise ValueError("Target transcript was not found in the TXT density file: {target}".format(target=self.target))

        target_table = pd.concat(target_chunks, axis=0, ignore_index=True)
        target_table["codon"] = target_table["codon"].astype(str).str.upper()
        target_table = target_table.loc[~target_table["codon"].str.contains("N", regex=False), :].copy()
        target_table.sort_values("now_nt", inplace=True)
        target_table.reset_index(drop=True, inplace=True)

        external_target = {
            "target_id": str(self.target),
            "gene_id": str(self.target),
            "transcript_id": str(self.target),
            "name": str(self.target),
        }
        genome_map, external_ids = self._external_genome_map_for_target(external_target)
        resolved_coordinate = self._resolve_coordinate_mode(self.coordinate, "txt", genome_map)

        gene_id = external_ids.get("gene_id", str(self.target)) if external_ids else str(self.target)
        transcript_id = external_ids.get("transcript_id", str(self.target)) if external_ids else str(self.target)

        self.sample_name = selected_samples
        self.total_rpf_num = totals
        self.resolved_coordinate = resolved_coordinate
        self.genome_map = genome_map if resolved_coordinate == "genome" else None
        self.target_meta = TargetMeta(
            target_id=str(self.target),
            gene_id=gene_id,
            transcript_id=transcript_id,
            name=str(self.target),
            chrom=genome_map.chrom if genome_map else None,
            strand=genome_map.strand if genome_map else None,
            candidate_count=1,
            selected_by="txt_name",
        )
        self.profile = self._build_profile_from_txt_table(target_table)
        self.feature_segments = self._build_feature_segments_from_profile(self.profile)

        print(
            "Scanned TXT chunks={chunks:,}, rows={rows:,}.".format(chunks=chunk_count, rows=row_count),
            flush=True,
        )

    # ------------------------------------------------------------------
    # Profile construction
    # ------------------------------------------------------------------

    @staticmethod
    def _region_from_codon_index(codon_index: int, utr5_codons: int, cds_codons: int) -> str:
        """Return transcript region from codon index."""
        if codon_index < utr5_codons:
            return "5utr"
        if codon_index < utr5_codons + cds_codons:
            return "cds"
        return "3utr"

    def _coordinate_mapper(self) -> Callable[[int, int | None], tuple[float, int | None]]:
        """Return a function mapping transcript nucleotide index to plotting coordinate."""
        if self.resolved_coordinate != "genome" or self.genome_map is None:
            return lambda tx_nt0, genome_pos=None: (float(tx_nt0 + 1), None)

        genome_map = self.genome_map
        genome_to_display = self._build_genome_display_mapper(genome_map)

        def mapper(tx_nt0: int, genome_pos: int | None = None) -> tuple[float, int | None]:
            if tx_nt0 < 0 or tx_nt0 >= len(genome_map.tx_to_genome):
                return (float("nan"), None)
            pos = genome_map.tx_to_genome[tx_nt0]
            return (genome_to_display(pos), pos)

        return mapper

    def _build_genome_display_mapper(self, genome_map: GenomeMap) -> Callable[[int], float]:
        """Build a genomic-coordinate display mapper with optional intron compression."""
        exons_1based = [(start + 1, end) for start, end in genome_map.exons]
        exons_1based = sorted(exons_1based, key=lambda item: item[0])

        if self.intron_scale >= 0.999:
            return lambda pos: float(pos)

        display_ranges: list[tuple[int, int, float]] = []
        cursor = float(exons_1based[0][0])
        previous_end: int | None = None
        for start, end in exons_1based:
            if previous_end is None:
                display_start = cursor
            else:
                intron_length = max(0, start - previous_end - 1)
                cursor += intron_length * self.intron_scale
                display_start = cursor
            display_ranges.append((start, end, display_start))
            cursor = display_start + (end - start + 1)
            previous_end = end

        def mapper(pos: int) -> float:
            for start, end, display_start in display_ranges:
                if start <= pos <= end:
                    return display_start + (pos - start)
            return float(pos)

        return mapper

    def _build_profile_from_json_record(self, record: dict[str, Any]) -> pd.DataFrame:
        """Build a long profile table from one JSON transcript record."""
        if self.target_meta is None:
            raise ValueError("Target metadata has not been initialized.")

        trim = record.get("trim") if isinstance(record.get("trim"), dict) else {}
        codon_count = self._to_int(trim.get("codon_count"), 0)
        if codon_count <= 0:
            raise ValueError("Selected JSON transcript has no valid codon_count.")

        start_nt0 = self._to_int(trim.get("start_nt0"), 0)
        utr5_codons = self._to_int(trim.get("utr5_codons"), 0)
        cds_codons = self._to_int(trim.get("cds_codons"), 0)
        utr3_nt = self._to_int(trim.get("utr3_nt"), 0)
        trim_length_nt = self._to_int(trim.get("trim_length_nt"), codon_count * 3)
        from_tts_start = (utr3_nt - trim_length_nt) // 3 + 1
        sequence = str(record.get("sequence", "")).upper()
        coordinate_mapper = self._coordinate_mapper()

        sample_arrays = {
            sample: RPFs.density_to_frame_arrays(record, sample, codon_count) for sample in self.sample_name
        }

        rows = []
        for codon_index in range(codon_count):
            codon = sequence[codon_index * 3 : codon_index * 3 + 3].upper()
            if "N" in codon:
                continue
            region = self._region_from_codon_index(codon_index, utr5_codons, cds_codons)
            from_tis = codon_index - utr5_codons
            from_tts = from_tts_start + codon_index

            for frame_index in range(3):
                tx_nt0 = start_nt0 + codon_index * 3 + frame_index
                x_value, genome_pos = coordinate_mapper(tx_nt0)
                if not np.isfinite(x_value):
                    continue
                transcript_nt = tx_nt0 + 1

                for sample, frame_arrays in sample_arrays.items():
                    raw_count = int(frame_arrays[frame_index][codon_index])
                    density = self._normalize_count(sample, raw_count)
                    rows.append(
                        {
                            "Target": self.target_meta.target_id,
                            "GeneID": self.target_meta.gene_id,
                            "TranscriptID": self.target_meta.transcript_id,
                            "Name": self.target_meta.name,
                            "Sample": sample,
                            "CodonIndex": codon_index,
                            "Frame": str(frame_index),
                            "TranscriptNt": transcript_nt,
                            "GenomePos": genome_pos,
                            "X": x_value,
                            "Region": region,
                            "Codon": codon,
                            "FromTIS": from_tis,
                            "FromTTS": from_tts,
                            "RawCount": raw_count,
                            "Density": density,
                        }
                    )

        return pd.DataFrame.from_records(rows)

    def _build_profile_from_txt_table(self, target_table: pd.DataFrame) -> pd.DataFrame:
        """Build a long profile table from a legacy TXT target table."""
        if self.target_meta is None:
            raise ValueError("Target metadata has not been initialized.")

        coordinate_mapper = self._coordinate_mapper()
        rows = []
        for codon_index, row in enumerate(target_table.itertuples(index=False)):
            row_data = row._asdict()
            codon = str(row_data["codon"]).upper()
            if "N" in codon:
                continue
            now_nt = self._to_int(row_data["now_nt"], codon_index * 3 + 1)
            region = str(row_data["region"]).lower()

            for frame_index in range(3):
                transcript_nt = now_nt + frame_index
                tx_nt0 = transcript_nt - 1
                x_value, genome_pos = coordinate_mapper(tx_nt0)
                if not np.isfinite(x_value):
                    continue
                for sample in self.sample_name:
                    column = f"{sample}_f{frame_index}"
                    raw_count = self._to_int(row_data[column], 0)
                    density = self._normalize_count(sample, raw_count)
                    rows.append(
                        {
                            "Target": self.target_meta.target_id,
                            "GeneID": self.target_meta.gene_id,
                            "TranscriptID": self.target_meta.transcript_id,
                            "Name": self.target_meta.name,
                            "Sample": sample,
                            "CodonIndex": codon_index,
                            "Frame": str(frame_index),
                            "TranscriptNt": transcript_nt,
                            "GenomePos": genome_pos,
                            "X": x_value,
                            "Region": region,
                            "Codon": codon,
                            "FromTIS": self._to_int(row_data["from_tis"], 0),
                            "FromTTS": self._to_int(row_data["from_tts"], 0),
                            "RawCount": raw_count,
                            "Density": density,
                        }
                    )

        return pd.DataFrame.from_records(rows)

    def _normalize_count(self, sample: str, raw_count: float) -> float:
        """Normalize one raw count to RPM if requested."""
        if not self.norm:
            return float(raw_count)
        total = float(self.total_rpf_num.get(sample, 0.0))
        if total <= 0:
            return 0.0
        return float(raw_count) * RPM_SCALE / total

    def _build_feature_segments_from_profile(self, profile: pd.DataFrame) -> pd.DataFrame:
        """Build collapsed feature segments for the gene-structure track."""
        if profile.empty:
            return pd.DataFrame(columns=["Region", "XStart", "XEnd", "RawStart", "RawEnd"])

        base = (
            profile.loc[:, ["TranscriptNt", "GenomePos", "X", "Region"]]
            .drop_duplicates()
            .sort_values("X")
            .reset_index(drop=True)
        )
        if base.empty:
            return pd.DataFrame(columns=["Region", "XStart", "XEnd", "RawStart", "RawEnd"])

        raw_column = "GenomePos" if self.resolved_coordinate == "genome" and base["GenomePos"].notna().any() else "TranscriptNt"
        segments = []
        start_idx = 0
        for idx in range(1, len(base)):
            same_region = base.loc[idx, "Region"] == base.loc[idx - 1, "Region"]
            close_x = abs(float(base.loc[idx, "X"]) - float(base.loc[idx - 1, "X"])) <= 1.5
            if same_region and close_x:
                continue
            segments.append(self._segment_from_base(base.iloc[start_idx:idx], raw_column))
            start_idx = idx
        segments.append(self._segment_from_base(base.iloc[start_idx:], raw_column))

        return pd.DataFrame.from_records(segments)

    @staticmethod
    def _segment_from_base(base: pd.DataFrame, raw_column: str) -> dict[str, Any]:
        """Create one feature segment from base-resolution records."""
        raw_values = pd.to_numeric(base[raw_column], errors="coerce").dropna()
        return {
            "Region": str(base["Region"].iloc[0]),
            "XStart": float(base["X"].min()),
            "XEnd": float(base["X"].max()),
            "RawStart": int(raw_values.min()) if not raw_values.empty else None,
            "RawEnd": int(raw_values.max()) if not raw_values.empty else None,
        }

    # ------------------------------------------------------------------
    # Output tables
    # ------------------------------------------------------------------

    def output_profile(self) -> None:
        """Write the raw gene-level density profile and optional track table."""
        if self.profile is None:
            raise ValueError("Geneplot profile has not been imported yet.")

        profile_out = self.output + "_geneplot.profile.txt"
        self.profile.to_csv(profile_out, sep="	", index=False)
        self.output_files["profile_table"] = profile_out

        if self.export_track:
            self.output_track_table()

    def output_track_table(self) -> None:
        """Write a multi-sample bedGraph-like track table.

        Notes
        -----
        The output keeps the four standard bedGraph columns first and appends
        sample/frame/annotation columns so that one file can preserve all tracks.
        For strict bigWig conversion, split this table by ``Track`` and keep
        only ``chrom``, ``start``, ``end``, and ``value`` after sorting.
        """
        if self.profile is None or self.target_meta is None:
            raise ValueError("Geneplot profile has not been imported yet.")

        data = self.profile.copy()
        if self.frame != "all":
            data = data.loc[data["Frame"] == self.frame, :].copy()
        if data.empty:
            return

        if self.resolved_coordinate == "genome" and self.target_meta.chrom:
            data = data.loc[data["GenomePos"].notna(), :].copy()
            data["chrom"] = self.target_meta.chrom
            data["start"] = data["GenomePos"].astype(int) - 1
            data["end"] = data["GenomePos"].astype(int)
        else:
            data["chrom"] = self.target_meta.transcript_id
            data["start"] = data["TranscriptNt"].astype(int) - 1
            data["end"] = data["TranscriptNt"].astype(int)

        data["value"] = data["Density"].astype(float)
        if self.data_type == "rna":
            data["Track"] = data["Sample"].astype(str)
        else:
            data["Track"] = data["Sample"].astype(str) + "_frame" + data["Frame"].astype(str)
        data["DataType"] = self.data_type.upper()
        data["Coordinate"] = str(self.resolved_coordinate)
        keep_columns = [
            "chrom",
            "start",
            "end",
            "value",
            "Track",
            "Sample",
            "Frame",
            "Region",
            "RawCount",
            "Density",
            "DataType",
            "Coordinate",
            "GeneID",
            "TranscriptID",
        ]
        data = data.loc[:, keep_columns].sort_values(["Track", "chrom", "start", "end"])

        track_out = self.output + "_geneplot.bedgraph_like.txt"
        data.to_csv(track_out, sep="	", index=False)
        self.output_files["bedgraph_like_table"] = track_out

    # ------------------------------------------------------------------
    # Plot helpers
    # ------------------------------------------------------------------

    def _transform_density_values(self, values) -> pd.Series:
        """Transform density values for plotting only."""
        transformed = pd.Series(values, dtype="float64").replace([np.inf, -np.inf], np.nan).fillna(0.0)
        transformed = transformed.clip(lower=0.0)

        if self.plot_transform == "none":
            return transformed
        if self.plot_transform == "sqrt":
            return np.sqrt(transformed)
        if self.plot_transform in {"log", "log1p"}:
            return np.log1p(transformed)
        if self.plot_transform == "log2":
            return np.log2(transformed + 1.0)
        if self.plot_transform == "log10":
            return np.log10(transformed + 1.0)

        raise ValueError("Unsupported plot transform: {transform}".format(transform=self.plot_transform))

    def _density_axis_label(self) -> str:
        """Return y-axis label for plotted density values."""
        density_name = "RNA density" if self.data_type == "rna" else "RPF density"
        base = "{name} (RPM)".format(name=density_name) if self.norm else density_name
        transform_label = PLOT_TRANSFORM_LABELS.get(self.plot_transform, "")
        if transform_label:
            base = "{base}, {label} transformed".format(base=base, label=transform_label)
        if self.spike_clip:
            base += " (spike clipped)"
        return base

    def _plot_density_values(self, values) -> pd.Series:
        """Return transformed density values used for plotting.

        Spike clipping is applied only to plotted values. The exported
        ``profile.txt`` table always keeps the original raw and normalized
        density values.
        """
        transformed = self._transform_density_values(values)
        if not self.spike_clip:
            return transformed

        finite = transformed.replace([np.inf, -np.inf], np.nan).dropna()
        finite = finite.loc[finite > 0]
        if finite.empty:
            return transformed

        if self.clip_value is not None:
            cutoff = float(self.clip_value)
        else:
            cutoff = float(finite.quantile(self.clip_quantile))
        if not np.isfinite(cutoff) or cutoff <= 0:
            return transformed
        return transformed.clip(upper=cutoff)

    @staticmethod
    def _set_axis_style(ax) -> None:
        """Apply common axis style without grid."""
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    @staticmethod
    def _nice_tick_step(span: float, max_ticks: int = 8) -> float:
        """Return a regular tick interval for transcript/genome coordinates."""
        span = abs(float(span))
        if span <= 0:
            return 1.0
        raw_step = span / max(1, max_ticks - 1)
        preferred_steps = [1, 3, 10, 30, 50, 100, 300, 500, 1_000, 3_000, 5_000, 10_000, 30_000, 50_000, 100_000, 300_000, 500_000, 1_000_000]
        for step in preferred_steps:
            if step >= raw_step:
                return float(step)
        exponent = int(np.ceil(np.log10(raw_step)))
        return float(10 ** exponent)

    @classmethod
    def _regular_ticks(cls, x_min: float, x_max: float, max_ticks: int = 8) -> list[float]:
        """Return regular ticks within one x-axis range."""
        if x_min > x_max:
            x_min, x_max = x_max, x_min
        step = cls._nice_tick_step(x_max - x_min, max_ticks=max_ticks)
        start = np.ceil(x_min / step) * step
        end = np.floor(x_max / step) * step
        ticks = list(np.arange(start, end + step * 0.5, step)) if start <= end else []
        if not ticks:
            ticks = [x_min, x_max]
        while len(ticks) > max_ticks:
            ticks = ticks[::2]
        return [float(tick) for tick in ticks]

    @staticmethod
    def _format_tick_label(value: float) -> str:
        """Format x-axis tick labels."""
        if abs(value) >= 1_000_000:
            return "{:.2f}M".format(value / 1_000_000).rstrip("0").rstrip(".")
        if abs(value) >= 10_000:
            return "{:.0f}K".format(value / 1_000)
        if abs(value - round(value)) < 1e-6:
            return str(int(round(value)))
        return "{:.1f}".format(value)

    def _get_plot_data(self) -> pd.DataFrame:
        """Return profile rows used for the requested plotting frame."""
        if self.profile is None:
            raise ValueError("Geneplot profile has not been imported yet.")

        data = self.profile.copy()
        if self.frame != "all":
            data = data.loc[data["Frame"] == self.frame, :].copy()
        if data.empty:
            raise ValueError("No density rows remain after frame filtering.")
        return data

    def _plot_bar_sample(self, ax, sample_data: pd.DataFrame) -> None:
        """Draw one sample as bars.

        Ribo mode uses frame-colored bars when all frames are plotted. RNA mode
        uses one density color because codon-frame separation is not biologically
        meaningful for ordinary RNA coverage tracks.
        """
        sample_data = sample_data.sort_values(["X", "Frame"]).copy()
        sample_data["PlotDensity"] = self._plot_density_values(sample_data["Density"]).to_numpy(dtype=float)

        if self.data_type == "rna":
            colors = []
            for region in sample_data["Region"].astype(str).str.lower():
                if self.utr_gray and region != "cds":
                    colors.append(UTR_COLOR)
                else:
                    colors.append(RNA_COLOR)
            ax.bar(
                sample_data["X"],
                sample_data["PlotDensity"],
                width=self.bar_width,
                color=colors,
                edgecolor="none",
                linewidth=0,
            )
            return

        for frame in FRAME_ORDER:
            if self.frame != "all" and frame != self.frame:
                continue
            frame_df = sample_data.loc[sample_data["Frame"] == frame, :].copy()
            if frame_df.empty:
                continue
            colors = []
            for region in frame_df["Region"].astype(str).str.lower():
                if self.utr_gray and region != "cds":
                    colors.append(UTR_COLOR)
                else:
                    colors.append(FRAME_COLORS[frame])
            ax.bar(
                frame_df["X"],
                frame_df["PlotDensity"],
                width=self.bar_width,
                color=colors,
                edgecolor="none",
                linewidth=0,
            )

    def _plot_line_sample(self, ax, sample_data: pd.DataFrame) -> None:
        """Draw one sample as a line profile."""
        sample_data = sample_data.sort_values(["X", "Frame"]).copy()
        sample_data["PlotDensity"] = self._plot_density_values(sample_data["Density"]).to_numpy(dtype=float)

        if self.utr_gray:
            group_id = (
                (sample_data["Region"].astype(str).str.lower() != sample_data["Region"].astype(str).str.lower().shift())
                | (sample_data["X"].diff().abs() > 1.5)
            ).cumsum()
            for _, part in sample_data.groupby(group_id, sort=False):
                region = str(part["Region"].iloc[0]).lower()
                color = UTR_COLOR if region != "cds" else RNA_COLOR if self.data_type == "rna" else LINE_COLOR
                ax.plot(part["X"], part["PlotDensity"], color=color, linewidth=self.line_width)
        else:
            line_color = RNA_COLOR if self.data_type == "rna" else LINE_COLOR
            ax.plot(sample_data["X"], sample_data["PlotDensity"], color=line_color, linewidth=self.line_width)

    def _sync_yaxis(self, axes, plot_data: pd.DataFrame) -> None:
        """Use the same y-axis range for all sample panels."""
        values = self._plot_density_values(plot_data["Density"]).to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if self.y_max is not None:
            y_top = float(self.y_max)
        elif values.size > 0 and float(np.max(values)) > 0:
            y_top = float(np.max(values)) * 1.08
        else:
            y_top = 1.0

        for ax in axes:
            ax.set_ylim(0, y_top)

    @staticmethod
    def _format_track_sum(value: float) -> str:
        """Format a gene-level density sum for compact panel annotation."""
        if not np.isfinite(value):
            return "0"
        abs_value = abs(float(value))
        if abs_value >= 1_000_000:
            return "{:.2f}M".format(value / 1_000_000).rstrip("0").rstrip(".")
        if abs_value >= 10_000:
            return "{:.1f}K".format(value / 1_000).rstrip("0").rstrip(".")
        if abs_value >= 100:
            return "{:.0f}".format(value)
        if abs_value >= 10:
            return "{:.1f}".format(value).rstrip("0").rstrip(".")
        return "{:.2f}".format(value).rstrip("0").rstrip(".")

    def _sample_total_label(self, sample_data: pd.DataFrame) -> str:
        """Return the gene-level count/RPM label for one sample panel."""
        density_name = "RNA" if self.data_type == "rna" else "RPF"
        if self.norm:
            total_value = float(sample_data["Density"].sum()) if not sample_data.empty else 0.0
            return "{name} RPM sum = {value}".format(
                name=density_name,
                value=self._format_track_sum(total_value),
            )

        total_value = float(sample_data["RawCount"].sum()) if not sample_data.empty else 0.0
        return "{name} count = {value}".format(
            name=density_name,
            value=self._format_track_sum(total_value),
        )

    def _annotate_sample_total(self, ax, sample_data: pd.DataFrame) -> None:
        """Print the gene-level count/RPM sum in the upper-right corner of one sample panel."""
        label = self._sample_total_label(sample_data)
        ax.text(
            0.985,
            0.92,
            label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=max(self.tick_size, 1.0),
            color="#222222",
        )

    def _draw_strand_arrows(self, ax, x_min: float, x_max: float, strand: str | None) -> None:
        """Draw strand direction arrows on the gene structure line."""
        if strand not in {"+", "-"}:
            return
        span = max(float(x_max - x_min), 1.0)
        arrow_count = max(1, min(8, int(span / max(span / 6.0, 1.0))))
        arrow_length = span * 0.025
        positions = np.linspace(x_min + span * 0.12, x_max - span * 0.12, arrow_count)
        for pos in positions:
            if strand == "+":
                start, end = pos - arrow_length / 2.0, pos + arrow_length / 2.0
            else:
                start, end = pos + arrow_length / 2.0, pos - arrow_length / 2.0
            ax.annotate(
                "",
                xy=(end, 0.5),
                xytext=(start, 0.5),
                arrowprops={"arrowstyle": "-|>", "color": INTRON_COLOR, "linewidth": 0.7, "shrinkA": 0, "shrinkB": 0},
            )

    def _draw_gene_structure(self, ax, x_min: float, x_max: float) -> None:
        """Draw the 5UTR/CDS/3UTR structure track."""
        if self.feature_segments is None or self.feature_segments.empty:
            ax.set_axis_off()
            return

        self._set_axis_style(ax)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel("Structure", rotation=0, ha="right", va="center", labelpad=35, fontsize=self.label_size)
        ax.hlines(0.5, x_min, x_max, color=INTRON_COLOR, linewidth=0.7, zorder=1)

        for row in self.feature_segments.itertuples(index=False):
            region = str(getattr(row, "Region")).lower()
            x_start = float(getattr(row, "XStart"))
            x_end = float(getattr(row, "XEnd"))
            width = max(x_end - x_start + 1.0, 0.8)
            rect_start = x_start - 0.5
            if region == "cds":
                y0, height, color = 0.34, 0.32, CDS_COLOR
            else:
                y0, height, color = 0.42, 0.16, UTR_COLOR
            ax.add_patch(
                plt.Rectangle(
                    (rect_start, y0),
                    width,
                    height,
                    facecolor=color,
                    edgecolor="none",
                    zorder=2,
                )
            )

        strand = self.target_meta.strand if self.target_meta is not None else None
        self._draw_strand_arrows(ax, x_min, x_max, strand)

        if self.resolved_coordinate == "genome" and self.target_meta is not None and self.target_meta.chrom:
            xlabel = "Genomic coordinate on {chrom} (bp)".format(chrom=self.target_meta.chrom)
            if self.intron_scale < 0.999:
                xlabel += "; introns compressed by {scale:g}".format(scale=self.intron_scale)
        else:
            xlabel = "Transcript coordinate (nt)"
        ax.set_xlabel(xlabel, fontsize=self.label_size)

    def draw_geneplot(self) -> None:
        """Draw the IGV-like gene-level density plot."""
        if self.profile is None or self.target_meta is None:
            raise ValueError("Geneplot profile has not been imported yet.")

        plot_data = self._get_plot_data()
        x_min = float(plot_data["X"].min())
        x_max = float(plot_data["X"].max())
        x_span = max(x_max - x_min, 1.0)
        x_left = x_min - x_span * self.x_margin
        x_right = x_max + x_span * self.x_margin

        sample_count = len(self.sample_name)
        figure_height = min(
            max(3.2, sample_count * self.per_sample_height + self.structure_height + 1.2),
            self.max_height,
        )
        height_ratios = [1.0] * sample_count + [self.structure_height]

        fig, axes = plt.subplots(
            nrows=sample_count + 1,
            ncols=1,
            figsize=(self.figure_width, figure_height),
            sharex=True,
            gridspec_kw={"height_ratios": height_ratios, "hspace": 0.05},
        )
        axes = np.atleast_1d(axes)
        sample_axes = axes[:-1]
        structure_ax = axes[-1]

        for ax, sample in zip(sample_axes, self.sample_name):
            sample_data = plot_data.loc[plot_data["Sample"] == sample, :]
            if self.mode == "bar":
                self._plot_bar_sample(ax, sample_data)
            elif self.mode == "line":
                self._plot_line_sample(ax, sample_data)
            else:
                raise ValueError("mode must be 'line' or 'bar'.")

            self._set_axis_style(ax)
            self._annotate_sample_total(ax, sample_data)
            ax.set_ylabel(sample, rotation=0, ha="left", va="center", labelpad=12, fontsize=self.sample_label_size)
            ax.yaxis.set_label_position("right")
            ax.tick_params(axis="x", labelbottom=False, labelsize=self.tick_size)
            ax.tick_params(axis="y", labelsize=self.tick_size)

        self._sync_yaxis(sample_axes, plot_data)
        sample_axes[0].set_title(self._plot_title(), fontsize=self.title_size, pad=8)
        fig.text(0.035, 0.56, self._density_axis_label(), rotation=90, ha="center", va="center", fontsize=self.label_size)

        self._draw_gene_structure(structure_ax, x_left, x_right)
        ticks = self._regular_ticks(x_left, x_right, max_ticks=9)
        structure_ax.set_xticks(ticks)
        structure_ax.set_xticklabels([self._format_tick_label(tick) for tick in ticks], rotation=0, ha="center", fontsize=self.tick_size)

        for ax in axes:
            ax.set_xlim(x_left, x_right)

        self._add_legend(fig)
        fig.subplots_adjust(left=0.11, right=0.86, top=0.90, bottom=0.13)

        if self.output_format in {"pdf", "both"}:
            out_pdf = self.output + "_geneplot.pdf"
            fig.savefig(out_pdf, bbox_inches="tight")
            self.output_files["figure_pdf"] = out_pdf
        if self.output_format in {"png", "both"}:
            out_png = self.output + "_geneplot.png"
            fig.savefig(out_png, dpi=self.dpi, bbox_inches="tight")
            self.output_files["figure_png"] = out_png
        plt.close(fig)

    def _plot_title(self) -> str:
        """Return plot title."""
        if self.title:
            return self.title
        if self.target_meta is None:
            return str(self.target)
        title = "{gene} / {transcript}".format(
            gene=self.target_meta.gene_id,
            transcript=self.target_meta.transcript_id,
        )
        if self.resolved_coordinate == "genome" and self.target_meta.chrom:
            title += " ({chrom}, {strand})".format(
                chrom=self.target_meta.chrom,
                strand=self.target_meta.strand or "+",
            )
        return title

    def _add_legend(self, fig) -> None:
        """Add a compact figure legend."""
        handles = []
        labels = []
        density_label = "RNA density" if self.data_type == "rna" else "RPF density"
        if self.data_type == "rna":
            if self.mode == "bar":
                handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=RNA_COLOR, edgecolor="none"))
            else:
                handles.append(plt.Line2D([0], [0], color=RNA_COLOR, linewidth=self.line_width))
            labels.append(density_label)
        elif self.mode == "bar":
            for frame in FRAME_ORDER:
                if self.frame != "all" and frame != self.frame:
                    continue
                handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=FRAME_COLORS[frame], edgecolor="none"))
                labels.append("Frame " + frame)
        else:
            handles.append(plt.Line2D([0], [0], color=LINE_COLOR, linewidth=self.line_width))
            labels.append(density_label)

        if self.utr_gray:
            handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=UTR_COLOR, edgecolor="none"))
            labels.append("UTR")
        handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=CDS_COLOR, edgecolor="none"))
        labels.append("CDS")

        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=min(len(labels), 6),
            frameon=False,
            fontsize=self.legend_size,
            bbox_to_anchor=(0.5, 0.005),
        )

    def run(self) -> None:
        """Run the full geneplot workflow."""
        self.import_rpf()
        self.output_profile()
        self.draw_geneplot()
