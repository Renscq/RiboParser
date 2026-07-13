#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-13
# Version: 0.2.8-dev.001
# Function: Shuffle RPF density profiles within each transcript to generate randomized control data.
# Input: RPF density file in JSONL or legacy TXT format and an optional transcript filter.
# Output: Shuffled density file in JSONL and/or TXT format and a summary JSON file.

"""Fast and reproducible transcript-level shuffling of RPF density data.

The shuffle preserves the total RPF count of every transcript and sample. In
shared mode, one common position permutation is applied to all samples so that
cross-sample site-level covariance is retained. In independent mode, each
sample receives an independent permutation.

Random number streams are derived from the global seed and transcript name.
Consequently, results are reproducible and independent of worker count or task
completion order.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import time
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
from typing import Any, Iterable, Iterator, Sequence

import numpy as np
import pandas as pd

from . import RPFs

try:
    import orjson  # type: ignore
except ImportError:  # pragma: no cover
    orjson = None


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
SPARSE_ENCODING = "sparse_codon_frame"
DENSE_ENCODING = "dense_codon_frame"
PROGRESS_EVERY = 1000
JSON_COMPRESS_LEVEL = 1


class Shuffle(object):
    """Shuffle RPF density within each transcript.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Shuffle``.
    """

    def __init__(self, args):
        self.rpf_file = args.rpf
        self.gene_file = args.list
        self.output_prefix = args.output
        self.seed = int(args.seed)
        self.individual = bool(args.individual)
        self.thread = max(1, int(args.thread))
        self.input_format = self._detect_format(self.rpf_file)
        self.output_format = self._resolve_output_format(args.output_format)
        self.density_encoding = args.density_encoding

        self.gene_filter = self._read_gene_filter(self.gene_file)
        self.sample_name: list[str] = []
        self.output_files: OrderedDict[str, str] = OrderedDict()
        self.transcript_count = 0
        self.row_count = 0
        self.started_at = time.time()

    # ------------------------------------------------------------------
    # Input helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _detect_format(path: str) -> str:
        """Detect density file format."""
        lowered = path.lower()
        if lowered.endswith((".jsonl.gz", ".jsonl", ".json.gz", ".json")):
            return "json"
        if lowered.endswith((".txt", ".txt.gz", ".tsv", ".tab")):
            return "txt"
        return RPFs._detect_rpf_format(path)

    def _resolve_output_format(self, requested: str) -> str:
        """Resolve automatic output format."""
        if requested != "auto":
            return requested
        return self.input_format

    @staticmethod
    def _read_gene_filter(path: str | None) -> set[str] | None:
        """Read transcript identifiers from an optional filter table."""
        if not path:
            return None
        table = pd.read_csv(path, sep="\t", header=0, dtype=str)
        if table.empty:
            return set()
        column = "transcript_id" if "transcript_id" in table.columns else table.columns[0]
        return set(table[column].dropna().astype(str))

    @staticmethod
    def _record_name(record: dict[str, Any]) -> str:
        """Return transcript name from one JSON record."""
        value = record.get("name") or record.get("transcript_id")
        if value is None:
            raise ValueError("JSON transcript record lacks name or transcript_id.")
        return str(value)

    def import_rpf(self) -> None:
        """Inspect input density metadata before shuffling."""
        if self.input_format == "json":
            self.sample_name = RPFs.get_json_sample_names(self.rpf_file)
        else:
            raw = RPFs.read_txt_rpf_file(
                rpf_file=self.rpf_file,
                gene=self.gene_file,
                tis=None,
                tts=None,
                sample_name=None,
            )
            self.sample_name = RPFs.retrieve_sample_name(raw, None)

        if not self.sample_name:
            raise ValueError("No sample was found in the RPF density file.")

        print(
            "Detected input format={fmt}, samples={samples}, shuffle_mode={mode}.".format(
                fmt=self.input_format,
                samples=len(self.sample_name),
                mode="independent" if self.individual else "shared",
            ),
            flush=True,
        )

    # ------------------------------------------------------------------
    # Randomization helpers
    # ------------------------------------------------------------------

    def _stable_seed(self, transcript: str, sample: str | None = None) -> int:
        """Return a deterministic seed for one transcript/sample stream."""
        token = f"{self.seed}\0{transcript}\0{sample or 'shared'}".encode("utf-8")
        digest = hashlib.blake2b(token, digest_size=8).digest()
        return int.from_bytes(digest, byteorder="little", signed=False)

    def _shuffle_sample_arrays(
        self,
        transcript: str,
        sample_arrays: dict[str, tuple[list[int], list[int], list[int]]],
    ) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Shuffle codon-frame arrays while preserving sample totals."""
        if not sample_arrays:
            return {}

        first = next(iter(sample_arrays.values()))
        site_count = len(first[0]) * 3
        if site_count <= 1:
            return {
                sample: tuple(np.asarray(frame, dtype=np.int64) for frame in arrays)  # type: ignore[return-value]
                for sample, arrays in sample_arrays.items()
            }

        common_perm = None
        if not self.individual:
            common_rng = np.random.default_rng(self._stable_seed(transcript))
            common_perm = common_rng.permutation(site_count)

        shuffled: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for sample, arrays in sample_arrays.items():
            matrix = np.column_stack(arrays).astype(np.int64, copy=False)
            flat = matrix.reshape(-1)

            if self.individual:
                rng = np.random.default_rng(self._stable_seed(transcript, sample))
                perm = rng.permutation(site_count)
            else:
                perm = common_perm

            shuffled_flat = flat[perm]
            shuffled_matrix = shuffled_flat.reshape(matrix.shape)
            shuffled[sample] = (
                shuffled_matrix[:, 0].copy(),
                shuffled_matrix[:, 1].copy(),
                shuffled_matrix[:, 2].copy(),
            )

        return shuffled

    @staticmethod
    def _arrays_to_density(
        arrays: tuple[np.ndarray, np.ndarray, np.ndarray],
        encoding: str,
    ) -> dict[str, Any]:
        """Encode shuffled frame arrays as sparse or dense JSON density."""
        f0, f1, f2 = arrays
        if encoding == "dense":
            return {
                "encoding": DENSE_ENCODING,
                "index_base": 0,
                "f0": f0.astype(int).tolist(),
                "f1": f1.astype(int).tolist(),
                "f2": f2.astype(int).tolist(),
            }

        matrix = np.column_stack((f0, f1, f2))
        codon_index, frame = np.nonzero(matrix)
        count = matrix[codon_index, frame]
        return {
            "encoding": SPARSE_ENCODING,
            "index_base": 0,
            "codon_index": codon_index.astype(int).tolist(),
            "frame": frame.astype(int).tolist(),
            "count": count.astype(int).tolist(),
        }

    def _output_encoding(self, source_density: dict[str, Any]) -> str:
        """Resolve JSON density encoding for one output record."""
        if self.density_encoding in {"sparse", "dense"}:
            return self.density_encoding
        source = str(source_density.get("encoding", SPARSE_ENCODING))
        return "dense" if source in {DENSE_ENCODING, "dense"} else "sparse"

    # ------------------------------------------------------------------
    # JSON workflow
    # ------------------------------------------------------------------

    def _shuffle_json_record(self, record: dict[str, Any]) -> dict[str, Any]:
        """Shuffle one JSON transcript record."""
        transcript = self._record_name(record)
        codon_count = int(record.get("trim", {}).get("codon_count", 0))
        samples = record.get("samples")
        if not isinstance(samples, dict) or not samples:
            raise ValueError(f"Transcript {transcript} has no valid samples dictionary.")

        sample_arrays = {
            sample: RPFs.density_to_frame_arrays(record, sample, codon_count)
            for sample in self.sample_name
        }
        shuffled = self._shuffle_sample_arrays(transcript, sample_arrays)

        output_record = deepcopy(record)
        output_samples = output_record["samples"]
        for sample in self.sample_name:
            sample_entry = output_samples[sample]
            source_density = sample_entry.get("density", {})
            encoding = self._output_encoding(source_density)
            sample_entry["density"] = self._arrays_to_density(shuffled[sample], encoding)

        output_record["shuffle"] = {
            "seed": self.seed,
            "mode": "independent" if self.individual else "shared",
            "unit": "transcript_codon_frame_position",
        }
        return output_record

    @staticmethod
    def _json_dumps(record: dict[str, Any]) -> bytes:
        """Serialize one compact JSON record."""
        if orjson is not None:
            return orjson.dumps(record)
        return json.dumps(record, ensure_ascii=False, separators=(",", ":")).encode("utf-8")

    @staticmethod
    def _legacy_header(sample_names: Sequence[str]) -> list[str]:
        """Return legacy TXT output header."""
        header = BASE_COLUMNS.copy()
        for sample in sample_names:
            header.extend([f"{sample}_f0", f"{sample}_f1", f"{sample}_f2"])
        return header

    @staticmethod
    def _json_record_to_txt_rows(record: dict[str, Any], sample_names: Sequence[str]) -> Iterator[list[Any]]:
        """Convert one shuffled JSON record to legacy codon rows."""
        trim = record.get("trim", {})
        codon_count = int(trim.get("codon_count", 0))
        start_nt0 = int(trim.get("start_nt0", 0))
        utr5_codons = int(trim.get("utr5_codons", 0))
        cds_codons = int(trim.get("cds_codons", 0))
        utr3_nt = int(trim.get("utr3_nt", 0))
        trim_length_nt = int(trim.get("trim_length_nt", codon_count * 3))
        from_tts_start = (utr3_nt - trim_length_nt) // 3 + 1
        sequence = str(record.get("sequence", "")).upper()
        name = str(record.get("name") or record.get("transcript_id") or "NA")
        arrays = [RPFs.density_to_frame_arrays(record, sample, codon_count) for sample in sample_names]

        for idx in range(codon_count):
            codon = sequence[idx * 3 : idx * 3 + 3]
            if "N" in codon:
                continue
            if idx < utr5_codons:
                region = "5utr"
            elif idx < utr5_codons + cds_codons:
                region = "cds"
            else:
                region = "3utr"
            row: list[Any] = [
                name,
                start_nt0 + 1 + idx * 3,
                idx - utr5_codons,
                from_tts_start + idx,
                region,
                codon,
            ]
            for f0, f1, f2 in arrays:
                row.extend([f0[idx], f1[idx], f2[idx]])
            yield row

    def _process_json_batch(self, records: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """Shuffle a bounded JSON record batch in parallel."""
        if self.thread == 1 or len(records) <= 1:
            return [self._shuffle_json_record(record) for record in records]
        with ThreadPoolExecutor(max_workers=min(self.thread, len(records))) as executor:
            return list(executor.map(self._shuffle_json_record, records))

    def _shuffle_json(self) -> None:
        """Run streaming JSON shuffle with bounded parallel batches."""
        output_json = self.output_prefix + "_shuffle.jsonl.gz"
        output_txt = self.output_prefix + "_shuffle.txt"
        batch_size = max(16, self.thread * 4)
        batch: list[dict[str, Any]] = []

        json_handle = None
        txt_handle = None
        txt_writer = None
        try:
            if self.output_format in {"json", "both"}:
                json_handle = gzip.open(output_json, "wb", compresslevel=JSON_COMPRESS_LEVEL)
                self.output_files["json"] = output_json
            if self.output_format in {"txt", "both"}:
                txt_handle = open(output_txt, "w", encoding="utf-8", newline="")
                txt_writer = csv.writer(txt_handle, delimiter="\t", lineterminator="\n")
                txt_writer.writerow(self._legacy_header(self.sample_name))
                self.output_files["txt"] = output_txt

            def write_records(shuffled_records: Iterable[dict[str, Any]]) -> None:
                for shuffled_record in shuffled_records:
                    if json_handle is not None:
                        json_handle.write(self._json_dumps(shuffled_record))
                        json_handle.write(b"\n")
                    if txt_writer is not None:
                        for row in self._json_record_to_txt_rows(shuffled_record, self.sample_name):
                            txt_writer.writerow(row)
                            self.row_count += 1
                    self.transcript_count += 1
                    if self.transcript_count % PROGRESS_EVERY == 0:
                        print(f"transcripts={self.transcript_count:,}", flush=True)

            for record in RPFs.iter_json_records(self.rpf_file):
                name = self._record_name(record)
                if self.gene_filter is not None and name not in self.gene_filter:
                    continue
                batch.append(record)
                if len(batch) >= batch_size:
                    write_records(self._process_json_batch(batch))
                    batch = []

            if batch:
                write_records(self._process_json_batch(batch))
        finally:
            if json_handle is not None:
                json_handle.close()
            if txt_handle is not None:
                txt_handle.close()

    # ------------------------------------------------------------------
    # TXT workflow
    # ------------------------------------------------------------------

    def _shuffle_txt_group(
        self,
        transcript: str,
        positions: np.ndarray,
        values: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Shuffle one TXT transcript density matrix."""
        sample_arrays: dict[str, tuple[list[int], list[int], list[int]]] = {}
        for sample_idx, sample in enumerate(self.sample_name):
            offset = sample_idx * 3
            sample_arrays[sample] = (
                values[:, offset].astype(int).tolist(),
                values[:, offset + 1].astype(int).tolist(),
                values[:, offset + 2].astype(int).tolist(),
            )
        shuffled = self._shuffle_sample_arrays(transcript, sample_arrays)
        output = np.empty_like(values)
        for sample_idx, sample in enumerate(self.sample_name):
            offset = sample_idx * 3
            output[:, offset] = shuffled[sample][0]
            output[:, offset + 1] = shuffled[sample][1]
            output[:, offset + 2] = shuffled[sample][2]
        return positions, output

    @staticmethod
    def _txt_to_json_record(group: pd.DataFrame, sample_names: Sequence[str]) -> dict[str, Any]:
        """Convert one legacy TXT transcript group to compact JSON."""
        name = str(group["name"].iloc[0])
        codons = group["codon"].astype(str).str.upper().tolist()
        regions = group["region"].astype(str).tolist()
        codon_count = len(group)
        utr5_codons = sum(region == "5utr" for region in regions)
        cds_codons = sum(region == "cds" for region in regions)
        utr3_codons = sum(region == "3utr" for region in regions)
        start_nt0 = int(group["now_nt"].iloc[0]) - 1

        samples: OrderedDict[str, dict[str, Any]] = OrderedDict()
        for sample in sample_names:
            arrays = (
                group[f"{sample}_f0"].to_numpy(dtype=np.int64),
                group[f"{sample}_f1"].to_numpy(dtype=np.int64),
                group[f"{sample}_f2"].to_numpy(dtype=np.int64),
            )
            samples[sample] = {
                "profile": "unknown",
                "density": Shuffle._arrays_to_density(arrays, "sparse"),
            }

        return {
            "record_type": "transcript",
            "transcript_id": name,
            "gene_id": name,
            "name": name,
            "annotation": {"source": "legacy_txt"},
            "genome_mapping": None,
            "trim": {
                "shift5_nt": 0,
                "shift3_nt": 0,
                "start_nt0": start_nt0,
                "end_nt0": start_nt0 + codon_count * 3,
                "trim_length_nt": codon_count * 3,
                "codon_count": codon_count,
                "utr5_nt": utr5_codons * 3,
                "cds_nt": cds_codons * 3,
                "utr3_nt": utr3_codons * 3,
                "utr5_codons": utr5_codons,
                "cds_codons": cds_codons,
                "utr3_codons": utr3_codons,
            },
            "sequence": "".join(codons),
            "samples": samples,
        }

    def _shuffle_txt(self) -> None:
        """Shuffle legacy TXT input using transcript-level NumPy arrays."""
        raw = RPFs.read_txt_rpf_file(
            rpf_file=self.rpf_file,
            gene=self.gene_file,
            tis=None,
            tts=None,
            sample_name=None,
        ).to_pandas()
        if raw.empty:
            raise ValueError("RPF density table is empty after transcript filtering.")

        frame_columns = [
            f"{sample}_f{frame}"
            for sample in self.sample_name
            for frame in range(3)
        ]
        groups = [
            (str(name), indices.to_numpy(dtype=np.int64), raw.loc[indices, frame_columns].to_numpy(dtype=np.int64))
            for name, indices in raw.groupby("name", sort=False).groups.items()
        ]

        def task(item):
            return self._shuffle_txt_group(*item)

        if self.thread > 1 and len(groups) > 1:
            with ThreadPoolExecutor(max_workers=min(self.thread, len(groups))) as executor:
                results = list(executor.map(task, groups))
        else:
            results = [task(item) for item in groups]

        for positions, values in results:
            raw.loc[positions, frame_columns] = values

        self.transcript_count = len(groups)
        self.row_count = len(raw)

        if self.output_format in {"txt", "both"}:
            output_txt = self.output_prefix + "_shuffle.txt"
            raw.to_csv(output_txt, sep="\t", index=False)
            self.output_files["txt"] = output_txt

        if self.output_format in {"json", "both"}:
            output_json = self.output_prefix + "_shuffle.jsonl.gz"
            with gzip.open(output_json, "wb", compresslevel=JSON_COMPRESS_LEVEL) as out:
                for _, group in raw.groupby("name", sort=False):
                    record = self._txt_to_json_record(group, self.sample_name)
                    record["shuffle"] = {
                        "seed": self.seed,
                        "mode": "independent" if self.individual else "shared",
                        "unit": "transcript_codon_frame_position",
                    }
                    out.write(self._json_dumps(record))
                    out.write(b"\n")
            self.output_files["json"] = output_json

    # ------------------------------------------------------------------
    # Public workflow and summary
    # ------------------------------------------------------------------

    def shuffle_rpfs(self) -> None:
        """Shuffle the input RPF density file."""
        if self.input_format == "json":
            self._shuffle_json()
        else:
            self._shuffle_txt()

        if self.transcript_count == 0:
            raise ValueError("No transcript was retained for shuffling.")

    def output_rpfs(self) -> None:
        """Write a machine-readable shuffle summary."""
        elapsed = max(time.time() - self.started_at, 1e-6)
        summary_path = self.output_prefix + "_shuffle.summary.json"
        summary = OrderedDict(
            [
                ("tool", "rpf_Shuffle"),
                ("version", "0.2.8-dev.001"),
                ("input_rpf", os.path.abspath(self.rpf_file)),
                ("input_format", self.input_format),
                ("output_format", self.output_format),
                ("transcript_filter", os.path.abspath(self.gene_file) if self.gene_file else None),
                ("sample_count", len(self.sample_name)),
                ("samples", self.sample_name),
                ("transcript_count", self.transcript_count),
                ("row_count", self.row_count if self.row_count else None),
                ("parameters", OrderedDict([
                    ("seed", self.seed),
                    ("shuffle_mode", "independent" if self.individual else "shared"),
                    ("thread", self.thread),
                    ("density_encoding", self.density_encoding),
                ])),
                ("preserved_properties", [
                    "transcript_sample_total_count",
                    "sample_global_total_count",
                    "metadata_and_transcript_order",
                ]),
                ("output_files", self.output_files),
                ("elapsed_seconds", round(elapsed, 3)),
            ]
        )
        with open(summary_path, "w", encoding="utf-8") as out:
            json.dump(summary, out, ensure_ascii=False, indent=2)
            out.write("\n")
        self.output_files["summary_json"] = summary_path

        print(
            "Shuffled transcripts={transcripts:,}, elapsed={elapsed:.2f}s.".format(
                transcripts=self.transcript_count,
                elapsed=elapsed,
            ),
            flush=True,
        )
