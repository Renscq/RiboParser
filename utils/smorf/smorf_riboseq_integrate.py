#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Integrate sample-specific smORF evidence with early NoEvidence removal.
# Input: One evidence table per sample and optional ORF genePred annotation.
# Output: Integrated evidence, metric matrices, and filtered ORF genePred.

"""Simplified smORF evidence integration.

The input contract is intentionally narrow:

* every evidence file contains exactly one sample;
* every file contains the same unique ``orf_id`` set;
* row order may differ between files;
* every file was generated with all evidence levels retained.

The workflow has two passes. The first pass reads only ``sample``, ``orf_id``,
and ``translation_evidence``. ORFs labelled ``NoEvidence`` in every sample are
removed before full metric tables are read. The second pass reads only columns
required by the integrated table and selected matrices.

Sample files are read concurrently with threads. No dedicated temporary
working directory, position-map file, cache, or memory-mapped metric file is
created. Atomic output uses only ``<final>.tmp`` beside each final output.
"""

from __future__ import annotations

import glob
import gzip
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TextIO

import numpy as np
import pandas as pd

from utils.ribo.ArgsParser import message_print, progress_print


DEFAULT_THREADS = 4
OUTPUT_CHUNK_ROWS = 50_000
MISSING_TOKENS = frozenset({"", ".", "NA", "N/A", "NONE", "NAN"})
TRANSLATION_LABELS = (
    "NoEvidence",
    "LowConfidence",
    "MediumConfidence",
    "HighConfidence",
)
TRANSLATION_RANK = {
    label: rank
    for rank, label in enumerate(TRANSLATION_LABELS)
}
TRANSLATION_DECODE = np.asarray(TRANSLATION_LABELS, dtype=object)
PERIODICITY_LABELS = ("NA", "Weak", "Moderate", "Strong")
RELEASE_LABELS = ("NA", "Weak", "Moderate", "Strong")
SHAPE_LABELS = ("NA", "Skewed", "Disperse", "Intermediate", "Uniform")
PERIODICITY_DECODE = np.asarray(PERIODICITY_LABELS, dtype=object)
RELEASE_DECODE = np.asarray(RELEASE_LABELS, dtype=object)
SHAPE_DECODE = np.asarray(SHAPE_LABELS, dtype=object)

STATIC_COLUMNS = (
    "orf_id",
    "gene_id",
    "transcript_id",
    "chrom",
    "strand",
    "category",
    "genomic_start",
    "genomic_end",
    "nt_length",
    "coding_nt_length",
    "coding_codon_count",
)
STATIC_INTEGER_COLUMNS = frozenset(
    {
        "genomic_start",
        "genomic_end",
        "nt_length",
        "coding_nt_length",
        "coding_codon_count",
    }
)

INTEGRATION_NUMERIC_COLUMNS = (
    "rpf_sum",
    "covered_codon",
    "covered_codon_ratio",
    "frame0_density",
    "frame1_density",
    "frame2_density",
    "frame0_ratio",
    "frame1_ratio",
    "frame2_ratio",
    "periodicity_score",
    "release_ratio",
)
INTEGRATION_TEXT_COLUMNS = (
    "profile_status",
    "periodicity_label",
    "release_label",
    "coverage_shape",
    "translation_evidence",
)

MATRIX_GROUPS: dict[str, tuple[str, ...]] = {
    "rpf_sum": ("rpf_sum",),
    "density": (
        "frame0_density",
        "frame1_density",
        "frame2_density",
    ),
    "ratio": (
        "frame0_ratio",
        "frame1_ratio",
        "frame2_ratio",
    ),
    "covered_codon_ratio": ("covered_codon_ratio",),
    "translation_evidence": ("translation_evidence",),
}
MATRIX_SUFFIXES = {
    "rpf_sum": "rpf_sum_matrix",
    "density": "frame_density_matrix",
    "ratio": "frame_ratio_matrix",
    "covered_codon_ratio": "covered_codon_ratio_matrix",
    "translation_evidence": "translation_evidence_matrix",
}
KNOWN_MATRIX_SUFFIXES = tuple(MATRIX_SUFFIXES.values())


@dataclass(frozen=True, slots=True)
class SampleSource:
    """Describe one sample-specific evidence table.

    Attributes:
        path: Evidence table path.
        sample: Unique sample name stored in the table.
        row_to_master: Input-row to master-ORF position mapping.
    """

    path: str
    sample: str
    row_to_master: np.ndarray


@dataclass(slots=True)
class IntegrationConfig:
    """Store smORF integration settings.

    Attributes:
        input_files: Evidence table paths or glob patterns.
        output_integrated: Optional integrated table path.
        output_matrix: Optional matrix output template.
        output_density_matrix: Export frame-density matrix.
        output_ratio_matrix: Export frame-ratio matrix.
        output_covered_codon_ratio_matrix: Export covered-codon-ratio matrix.
        output_translation_evidence_matrix: Export evidence-label matrix.
        orf_genepred: Optional ORF genePred input.
        capture_labels: Evidence labels counted as captured.
        pass_labels: Evidence labels counted as reliable.
        excellent_min_samples: Reliable sample count defining Excellent.
        reliable_only: Keep only ORFs supported by at least one pass sample.
        threads: Concurrent sample-file readers.
    """

    input_files: tuple[str, ...]
    output_integrated: str | None = None
    output_matrix: str | None = None
    output_density_matrix: bool = False
    output_ratio_matrix: bool = False
    output_covered_codon_ratio_matrix: bool = False
    output_translation_evidence_matrix: bool = False
    orf_genepred: str | None = None
    capture_labels: tuple[str, ...] = field(
        default_factory=lambda: (
            "LowConfidence",
            "MediumConfidence",
            "HighConfidence",
        )
    )
    pass_labels: tuple[str, ...] = field(
        default_factory=lambda: (
            "MediumConfidence",
            "HighConfidence",
        )
    )
    excellent_min_samples: int = 2
    reliable_only: bool = False
    threads: int = DEFAULT_THREADS


@dataclass(frozen=True, slots=True)
class IntegrationResult:
    """Store final integration outputs.

    Attributes:
        sample_count: Number of input samples.
        input_orf_count: ORFs before all-NoEvidence removal.
        evidence_orf_count: ORFs after all-NoEvidence removal.
        retained_orf_count: ORFs after optional reliable-only filtering.
        integrated_output: Integrated output path.
        matrix_outputs: Generated matrix paths.
        genepred_output: Filtered genePred path.
    """

    sample_count: int
    input_orf_count: int
    evidence_orf_count: int
    retained_orf_count: int
    integrated_output: str | None
    matrix_outputs: tuple[str, ...]
    genepred_output: str | None


@dataclass(slots=True)
class SampleData:
    """Store one sample aligned to the reduced master ORF table.

    Attributes:
        sample: Sample name.
        numeric: Numeric metric arrays.
        translation_code: Encoded translation-evidence labels.
        profile_invalid: Invalid-profile mask.
        periodicity_code: Encoded periodicity labels.
        release_code: Encoded release labels.
        shape_code: Encoded coverage-shape labels.
    """

    sample: str
    numeric: dict[str, np.ndarray]
    translation_code: np.ndarray
    profile_invalid: np.ndarray
    periodicity_code: np.ndarray
    release_code: np.ndarray
    shape_code: np.ndarray


def eprint(message: str) -> None:
    """Print one plain integration message.

    Args:
        message: Message text.
    """
    message_print(message)


def smart_open(path: str | Path, mode: str = "rt") -> TextIO:
    """Open a plain or gzip-compressed file.

    Args:
        path: Input or output path.
        mode: Open mode.

    Returns:
        File handle.
    """
    file_path = Path(path)
    if file_path.name.lower().endswith(".gz"):
        return gzip.open(file_path, mode, encoding="utf-8", newline="")
    return file_path.open(mode, encoding="utf-8", newline="")


def split_items(value: str | None) -> tuple[str, ...]:
    """Split a comma-separated option.

    Args:
        value: Comma-separated text.

    Returns:
        Non-empty values.
    """
    if value is None:
        return ()
    return tuple(
        item.strip()
        for item in str(value).split(",")
        if item.strip()
    )


def resolve_input_files(values: tuple[str, ...]) -> tuple[str, ...]:
    """Resolve input paths and glob patterns.

    Args:
        values: Paths or glob patterns.

    Returns:
        Unique existing files in deterministic order.

    Raises:
        FileNotFoundError: If one path or pattern matches nothing.
        ValueError: If no file remains.
    """
    resolved: list[str] = []
    seen: set[str] = set()
    for value in values:
        matches = sorted(glob.glob(value))
        if not matches and Path(value).is_file():
            matches = [value]
        if not matches:
            raise FileNotFoundError(f"No evidence file matched: {value}")
        for match in matches:
            path = Path(match)
            if not path.is_file():
                continue
            canonical = str(path.resolve())
            if canonical in seen:
                continue
            seen.add(canonical)
            resolved.append(str(path))
    if not resolved:
        raise ValueError("No evidence input file was resolved.")
    return tuple(resolved)


def _read_header(path: str | Path) -> list[str]:
    """Read and validate a TSV header."""
    with smart_open(path, "rt") as handle:
        line = handle.readline().rstrip("\r\n")
    if not line:
        raise ValueError(f"Empty evidence table: {path}")
    header = line.split("\t")
    duplicated = sorted(
        name
        for name, count in pd.Series(header).value_counts().items()
        if count > 1
    )
    if duplicated:
        raise ValueError(
            f"Duplicate columns in {path}: {', '.join(duplicated)}"
        )
    return header


def _normalize_text(series: pd.Series) -> pd.Series:
    """Normalize text and missing tokens."""
    text = series.astype(str).str.strip()
    return text.mask(text.str.upper().isin(MISSING_TOKENS), ".")


def _strict_numeric(
    series: pd.Series,
    column: str,
    path: str,
) -> np.ndarray:
    """Convert one numeric column and reject malformed values."""
    text = series.astype(str).str.strip()
    missing = text.str.upper().isin(MISSING_TOKENS)
    numeric = pd.to_numeric(text.mask(missing), errors="coerce")
    invalid = numeric.isna() & ~missing
    if invalid.any():
        examples = text.loc[invalid].head(5).tolist()
        raise ValueError(
            f"Invalid numeric values in {column} from {path}: {examples}"
        )
    return numeric.to_numpy(dtype=np.float64)


def _encode_labels(
    series: pd.Series,
    labels: tuple[str, ...],
    column: str,
    path: str,
) -> np.ndarray:
    """Encode one categorical evidence column."""
    normalized = _normalize_text(series).replace({".": "NA"})
    mapping = {label: index for index, label in enumerate(labels)}
    unknown = sorted(set(normalized.unique()).difference(mapping))
    if unknown:
        raise ValueError(
            f"Unknown {column} labels in {path}: {', '.join(unknown)}"
        )
    return normalized.map(mapping).to_numpy(dtype=np.int8)


def _table_suffix(path: str | Path) -> tuple[str, str]:
    """Split a table path into stem and compound suffix."""
    text = str(path)
    for suffix in (".tsv.gz", ".txt.gz", ".tab.gz", ".tsv", ".txt", ".tab"):
        if text.lower().endswith(suffix):
            return text[: -len(suffix)], text[-len(suffix):]
    return text, ".txt"


def matrix_output_paths(
    template: str | Path,
    groups: tuple[str, ...],
) -> dict[str, str]:
    """Derive metric-specific matrix output paths.

    Args:
        template: User matrix output template.
        groups: Enabled matrix groups.

    Returns:
        Matrix-group to output-path mapping.
    """
    stem, suffix = _table_suffix(template)
    for known_suffix in KNOWN_MATRIX_SUFFIXES:
        token = f".{known_suffix}"
        if stem.endswith(token):
            stem = stem[: -len(token)]
            break
    return {
        group: f"{stem}.{MATRIX_SUFFIXES[group]}{suffix}"
        for group in groups
    }


def filtered_genepred_output_path(
    integrated_output: str | None,
    matrix_template: str | None,
) -> str:
    """Derive the filtered ORF genePred output path."""
    source = integrated_output or matrix_template
    if source is None:
        raise ValueError(
            "Cannot derive filtered genePred output without another output."
        )
    stem, _ = _table_suffix(source)
    for known_suffix in KNOWN_MATRIX_SUFFIXES:
        token = f".{known_suffix}"
        if stem.endswith(token):
            stem = stem[: -len(token)]
            break
    return f"{stem}.filtered.genePred"


class _AtomicTableWriter:
    """Write one tabular output atomically beside its final path."""

    def __init__(self, path: str | Path) -> None:
        """Initialize an atomic writer.

        Args:
            path: Final output path.
        """
        self.final_path = Path(path)
        self.temp_path = self.final_path.with_name(
            self.final_path.name + ".tmp"
        )
        self.handle: TextIO | None = None
        self.header_written = False

    def __enter__(self) -> "_AtomicTableWriter":
        """Open the temporary output."""
        self.final_path.parent.mkdir(parents=True, exist_ok=True)
        if self.final_path.name.lower().endswith(".gz"):
            self.handle = gzip.open(
                self.temp_path,
                "wt",
                encoding="utf-8",
                newline="",
            )
        else:
            self.handle = self.temp_path.open(
                "w",
                encoding="utf-8",
                newline="",
            )
        return self

    def write(self, table: pd.DataFrame) -> None:
        """Append one DataFrame block.

        Args:
            table: Output block.
        """
        if self.handle is None:
            raise RuntimeError("Atomic writer is not active.")
        table.to_csv(
            self.handle,
            sep="\t",
            index=False,
            header=not self.header_written,
            lineterminator="\n",
        )
        self.header_written = True

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> bool:
        """Commit a successful output or remove the temporary file."""
        if self.handle is not None:
            self.handle.close()
            self.handle = None
        if exc_type is None:
            os.replace(self.temp_path, self.final_path)
        else:
            self.temp_path.unlink(missing_ok=True)
        return False


class SmorfEvidenceIntegrator:
    """Integrate complete sample-specific smORF evidence tables."""

    def __init__(self, config: IntegrationConfig) -> None:
        """Initialize the integrator.

        Args:
            config: Integration settings.
        """
        self.config = config
        self.files = resolve_input_files(config.input_files)
        self.sources: list[SampleSource] = []
        self.master_ids: pd.Index | None = None
        self.master_to_reduced: np.ndarray | None = None
        self.static_table: pd.DataFrame | None = None
        self.sample_names: list[str] = []
        self.input_orf_count = 0
        self.evidence_orf_count = 0

    def run(self) -> IntegrationResult:
        """Run early filtering, integration, matrices, and genePred export.

        Returns:
            Integration result.
        """
        self._validate_config()
        self._scan_all_no_evidence()
        self._load_static_table()

        if self.static_table is None:
            raise RuntimeError("Reduced static ORF table was not created.")

        accumulators = self._initialize_accumulators(
            len(self.static_table)
        )
        matrix_store = self._initialize_matrix_store(
            len(self.static_table),
            len(self.sources),
        )
        self._read_and_integrate_samples(accumulators, matrix_store)

        final_mask = self._final_mask(accumulators)
        retained_count = int(final_mask.sum())
        integrated_output = None
        matrix_outputs: tuple[str, ...] = ()
        genepred_output = None

        if self.config.output_integrated:
            integrated = self._build_integrated_table(accumulators)
            self._write_dataframe(
                integrated.loc[final_mask].reset_index(drop=True),
                self.config.output_integrated,
            )
            integrated_output = self.config.output_integrated

        if self.config.output_matrix:
            paths = matrix_output_paths(
                self.config.output_matrix,
                tuple(matrix_store),
            )
            for group, path in paths.items():
                self._write_matrix(
                    group=group,
                    values=matrix_store[group],
                    final_mask=final_mask,
                    path=path,
                )
            matrix_outputs = tuple(paths.values())

        if self.config.orf_genepred:
            genepred_output = filtered_genepred_output_path(
                self.config.output_integrated,
                self.config.output_matrix,
            )
            final_ids = set(
                self.static_table.loc[final_mask, "orf_id"].astype(str)
            )
            self._write_filtered_genepred(
                input_path=self.config.orf_genepred,
                output_path=genepred_output,
                retained_ids=final_ids,
            )

        return IntegrationResult(
            sample_count=len(self.sources),
            input_orf_count=self.input_orf_count,
            evidence_orf_count=self.evidence_orf_count,
            retained_orf_count=retained_count,
            integrated_output=integrated_output,
            matrix_outputs=matrix_outputs,
            genepred_output=genepred_output,
        )

    def _validate_config(self) -> None:
        """Validate configuration values."""
        if not self.config.output_integrated and not self.config.output_matrix:
            raise ValueError(
                "At least one integrated or matrix output is required."
            )
        if self.config.threads < 1:
            raise ValueError("threads must be >= 1.")
        if self.config.excellent_min_samples < 1:
            raise ValueError("excellent_min_samples must be >= 1.")

        known = set(TRANSLATION_LABELS)
        unknown_capture = set(self.config.capture_labels).difference(known)
        unknown_pass = set(self.config.pass_labels).difference(known)
        if unknown_capture:
            raise ValueError(
                "Unknown capture labels: "
                + ", ".join(sorted(unknown_capture))
            )
        if unknown_pass:
            raise ValueError(
                "Unknown pass labels: "
                + ", ".join(sorted(unknown_pass))
            )
        if not set(self.config.pass_labels).issubset(
            set(self.config.capture_labels)
        ):
            raise ValueError(
                "pass_labels must be a subset of capture_labels."
            )

    def _scan_first_file(self) -> tuple[str, pd.Index, np.ndarray]:
        """Read the first file and define the master ORF index."""
        path = self.files[0]
        header = _read_header(path)
        required = {"sample", "orf_id", "translation_evidence"}
        missing = required.difference(header)
        if missing:
            raise ValueError(
                f"Evidence table {path} is missing: {', '.join(sorted(missing))}"
            )
        table = pd.read_csv(
            path,
            sep="\t",
            usecols=["sample", "orf_id", "translation_evidence"],
            dtype=str,
            keep_default_na=False,
            low_memory=False,
        )
        sample_values = _normalize_text(table["sample"]).unique()
        if len(sample_values) != 1 or sample_values[0] == ".":
            raise ValueError(
                f"Each evidence file must contain one sample: {path}"
            )
        ids = _normalize_text(table["orf_id"])
        if ids.eq(".").any() or ids.duplicated().any():
            raise ValueError(f"Invalid or duplicated orf_id values in {path}.")
        evidence = _encode_labels(
            table["translation_evidence"],
            TRANSLATION_LABELS,
            "translation_evidence",
            path,
        )
        return str(sample_values[0]), pd.Index(ids.astype(str)), evidence

    def _scan_other_file(
        self,
        path: str,
        master_index: pd.Index,
    ) -> tuple[str, np.ndarray, np.ndarray]:
        """Read one file's keys and map rows to the master ORF index."""
        header = _read_header(path)
        required = {"sample", "orf_id", "translation_evidence"}
        missing = required.difference(header)
        if missing:
            raise ValueError(
                f"Evidence table {path} is missing: {', '.join(sorted(missing))}"
            )
        table = pd.read_csv(
            path,
            sep="\t",
            usecols=["sample", "orf_id", "translation_evidence"],
            dtype=str,
            keep_default_na=False,
            low_memory=False,
        )
        sample_values = _normalize_text(table["sample"]).unique()
        if len(sample_values) != 1 or sample_values[0] == ".":
            raise ValueError(
                f"Each evidence file must contain one sample: {path}"
            )
        ids = _normalize_text(table["orf_id"])
        if ids.eq(".").any() or ids.duplicated().any():
            raise ValueError(f"Invalid or duplicated orf_id values in {path}.")
        if len(ids) != len(master_index):
            raise ValueError(
                f"ORF count differs in {path}: {len(ids):,} versus "
                f"{len(master_index):,}. Rerun smorf_evidence with all ORFs."
            )
        row_to_master = master_index.get_indexer(ids.astype(str))
        if np.any(row_to_master < 0) or np.unique(row_to_master).size != len(master_index):
            raise ValueError(
                f"ORF set differs in {path}. All evidence files must contain "
                "the same complete ORF set."
            )
        evidence = _encode_labels(
            table["translation_evidence"],
            TRANSLATION_LABELS,
            "translation_evidence",
            path,
        )
        return (
            str(sample_values[0]),
            row_to_master.astype(np.int32, copy=False),
            evidence,
        )

    def _scan_all_no_evidence(self) -> None:
        """Remove ORFs labelled NoEvidence in every sample before full reads."""
        eprint(
            "Scan minimal evidence columns: "
            f"files={len(self.files)}, threads={min(self.config.threads, len(self.files))}."
        )
        first_sample, master_index, first_evidence = self._scan_first_file()
        self.master_ids = master_index
        self.input_orf_count = len(master_index)
        any_evidence = first_evidence != TRANSLATION_RANK["NoEvidence"]
        sources_by_path: dict[str, SampleSource] = {
            self.files[0]: SampleSource(
                path=self.files[0],
                sample=first_sample,
                row_to_master=np.arange(
                    len(master_index),
                    dtype=np.int32,
                ),
            )
        }
        seen_samples = {first_sample}

        worker_count = min(
            self.config.threads,
            max(1, len(self.files) - 1),
        )
        if len(self.files) > 1:
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                future_to_path = {
                    executor.submit(
                        self._scan_other_file,
                        path,
                        master_index,
                    ): path
                    for path in self.files[1:]
                }
                for future in as_completed(future_to_path):
                    path = future_to_path[future]
                    sample, row_to_master, evidence = future.result()
                    if sample in seen_samples:
                        raise ValueError(
                            f"Sample appears in more than one file: {sample}"
                        )
                    seen_samples.add(sample)
                    any_evidence[
                        row_to_master[
                            evidence != TRANSLATION_RANK["NoEvidence"]
                        ]
                    ] = True
                    sources_by_path[path] = SampleSource(
                        path=path,
                        sample=sample,
                        row_to_master=row_to_master,
                    )

        retained_master = np.flatnonzero(any_evidence)
        master_to_reduced = np.full(
            len(master_index),
            -1,
            dtype=np.int32,
        )
        master_to_reduced[retained_master] = np.arange(
            len(retained_master),
            dtype=np.int32,
        )
        self.master_to_reduced = master_to_reduced
        self.evidence_orf_count = len(retained_master)
        self.sources = [sources_by_path[path] for path in self.files]
        self.sample_names = [source.sample for source in self.sources]

        removed = self.input_orf_count - self.evidence_orf_count
        removed_ratio = removed / max(self.input_orf_count, 1)
        eprint(
            "Remove all-sample NoEvidence ORFs: "
            f"input={self.input_orf_count:,}, removed={removed:,} "
            f"({removed_ratio:.1%}), retained={self.evidence_orf_count:,}."
        )

    def _load_static_table(self) -> None:
        """Read static annotation from the first evidence file."""
        if self.master_to_reduced is None:
            raise RuntimeError("NoEvidence filter was not initialized.")
        path = self.sources[0].path
        header = _read_header(path)
        columns = [column for column in STATIC_COLUMNS if column in header]
        if "orf_id" not in columns:
            raise ValueError(f"Evidence table lacks orf_id: {path}")
        table = pd.read_csv(
            path,
            sep="\t",
            usecols=columns,
            dtype=str,
            keep_default_na=False,
            low_memory=False,
        )
        if len(table) != self.input_orf_count:
            raise ValueError(f"Evidence file changed during integration: {path}")
        reduced = self.master_to_reduced >= 0
        table = table.loc[reduced].reset_index(drop=True)
        for column in table.columns:
            if column in STATIC_INTEGER_COLUMNS:
                numeric = _strict_numeric(
                    table[column],
                    column,
                    path,
                )
                if np.isnan(numeric).any():
                    raise ValueError(
                        f"Missing static integer values in {column} from {path}."
                    )
                table[column] = numeric.astype(np.int64)
            else:
                table[column] = _normalize_text(table[column])
        self.static_table = table

    def _required_full_columns(self, path: str) -> list[str]:
        """Return columns needed from one full sample evidence file."""
        header = set(_read_header(path))
        required = set(INTEGRATION_NUMERIC_COLUMNS)
        required.update(INTEGRATION_TEXT_COLUMNS)
        for group in self._enabled_matrix_groups():
            required.update(MATRIX_GROUPS[group])
        required.discard("sample")
        required.discard("orf_id")
        columns = [column for column in required if column in header]
        missing_core = {
            "rpf_sum",
            "covered_codon_ratio",
            "frame0_ratio",
            "periodicity_score",
            "translation_evidence",
        }.difference(columns)
        if missing_core:
            raise ValueError(
                f"Evidence table {path} lacks required integration columns: "
                + ", ".join(sorted(missing_core))
            )
        return columns

    def _read_sample_data(self, source: SampleSource) -> SampleData:
        """Read one full sample and align retained ORFs to reduced positions."""
        if self.master_to_reduced is None:
            raise RuntimeError("Reduced ORF mapping is unavailable.")
        columns = self._required_full_columns(source.path)
        numeric_columns = set(INTEGRATION_NUMERIC_COLUMNS)
        dtype_map = {
            column: (np.float32 if column in numeric_columns else str)
            for column in columns
        }
        try:
            table = pd.read_csv(
                source.path,
                sep="\t",
                usecols=columns,
                dtype=dtype_map,
                keep_default_na=False,
                low_memory=False,
            )
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"Failed to parse evidence metrics from {source.path}: {error}"
            ) from error
        if len(table) != len(source.row_to_master):
            raise ValueError(
                f"Evidence file changed during integration: {source.path}"
            )
        target = self.master_to_reduced[source.row_to_master]
        keep = target >= 0
        target = target[keep]
        size = self.evidence_orf_count

        numeric: dict[str, np.ndarray] = {}
        for column in INTEGRATION_NUMERIC_COLUMNS:
            output = np.full(size, np.nan, dtype=np.float32)
            if column in table.columns:
                values = table[column].to_numpy(dtype=np.float32, copy=False)
                output[target] = values[keep]
            numeric[column] = output

        translation_values = _encode_labels(
            table["translation_evidence"],
            TRANSLATION_LABELS,
            "translation_evidence",
            source.path,
        )
        translation_code = np.full(size, -1, dtype=np.int8)
        translation_code[target] = translation_values[keep]

        profile_invalid = np.zeros(size, dtype=bool)
        if "profile_status" in table.columns:
            invalid = (
                _normalize_text(table["profile_status"])
                .eq("InvalidProfile")
                .to_numpy()
            )
            profile_invalid[target] = invalid[keep]

        periodicity_code = np.zeros(size, dtype=np.int8)
        if "periodicity_label" in table.columns:
            values = _encode_labels(
                table["periodicity_label"],
                PERIODICITY_LABELS,
                "periodicity_label",
                source.path,
            )
            periodicity_code[target] = values[keep]

        release_code = np.zeros(size, dtype=np.int8)
        if "release_label" in table.columns:
            values = _encode_labels(
                table["release_label"],
                RELEASE_LABELS,
                "release_label",
                source.path,
            )
            release_code[target] = values[keep]

        shape_code = np.zeros(size, dtype=np.int8)
        if "coverage_shape" in table.columns:
            values = _encode_labels(
                table["coverage_shape"],
                SHAPE_LABELS,
                "coverage_shape",
                source.path,
            )
            shape_code[target] = values[keep]

        return SampleData(
            sample=source.sample,
            numeric=numeric,
            translation_code=translation_code,
            profile_invalid=profile_invalid,
            periodicity_code=periodicity_code,
            release_code=release_code,
            shape_code=shape_code,
        )

    def _enabled_matrix_groups(self) -> tuple[str, ...]:
        """Return enabled matrix groups in output order."""
        if not self.config.output_matrix:
            return ()
        groups = ["rpf_sum"]
        if self.config.output_density_matrix:
            groups.append("density")
        if self.config.output_ratio_matrix:
            groups.append("ratio")
        if self.config.output_covered_codon_ratio_matrix:
            groups.append("covered_codon_ratio")
        if self.config.output_translation_evidence_matrix:
            groups.append("translation_evidence")
        return tuple(groups)

    def _initialize_accumulators(self, size: int) -> dict[str, np.ndarray]:
        """Create compact ORF-level accumulators."""
        return {
            "capture_count": np.zeros(size, dtype=np.int16),
            "pass_count": np.zeros(size, dtype=np.int16),
            "no_count": np.zeros(size, dtype=np.int16),
            "low_count": np.zeros(size, dtype=np.int16),
            "medium_count": np.zeros(size, dtype=np.int16),
            "high_count": np.zeros(size, dtype=np.int16),
            "invalid_count": np.zeros(size, dtype=np.int16),
            "rpf_sum": np.zeros(size, dtype=np.float64),
            "rpf_max": np.zeros(size, dtype=np.float64),
            "covered_codon_ratio_max": np.zeros(size, dtype=np.float64),
            "frame0_ratio_max": np.zeros(size, dtype=np.float64),
            "periodicity_score_max": np.zeros(size, dtype=np.float64),
            "release_ratio_max": np.zeros(size, dtype=np.float64),
            "best_rank": np.full(size, -1, dtype=np.int8),
            "best_rpf": np.full(size, -np.inf, dtype=np.float64),
            "best_periodicity": np.full(size, -np.inf, dtype=np.float64),
            "best_coverage": np.full(size, -np.inf, dtype=np.float64),
            "best_frame0": np.full(size, -np.inf, dtype=np.float64),
            "best_sample_index": np.full(size, -1, dtype=np.int16),
            "best_translation_code": np.full(size, -1, dtype=np.int8),
            "best_rpf_sum": np.full(size, np.nan, dtype=np.float64),
            "best_covered_codon_ratio": np.full(size, np.nan, dtype=np.float64),
            "best_frame0_density": np.full(size, np.nan, dtype=np.float64),
            "best_frame1_density": np.full(size, np.nan, dtype=np.float64),
            "best_frame2_density": np.full(size, np.nan, dtype=np.float64),
            "best_frame0_ratio": np.full(size, np.nan, dtype=np.float64),
            "best_frame1_ratio": np.full(size, np.nan, dtype=np.float64),
            "best_frame2_ratio": np.full(size, np.nan, dtype=np.float64),
            "best_periodicity_score": np.full(size, np.nan, dtype=np.float64),
            "best_release_ratio": np.full(size, np.nan, dtype=np.float64),
            "best_periodicity_code": np.zeros(size, dtype=np.int8),
            "best_release_code": np.zeros(size, dtype=np.int8),
            "best_shape_code": np.zeros(size, dtype=np.int8),
        }

    def _initialize_matrix_store(
        self,
        orf_count: int,
        sample_count: int,
    ) -> dict[str, dict[str, np.ndarray]]:
        """Allocate in-memory matrices only for requested groups."""
        store: dict[str, dict[str, np.ndarray]] = {}
        for group in self._enabled_matrix_groups():
            metrics: dict[str, np.ndarray] = {}
            for metric in MATRIX_GROUPS[group]:
                if metric == "translation_evidence":
                    metrics[metric] = np.full(
                        (orf_count, sample_count),
                        -1,
                        dtype=np.int8,
                    )
                else:
                    metrics[metric] = np.zeros(
                        (orf_count, sample_count),
                        dtype=np.float32,
                    )
            store[group] = metrics
        return store

    def _read_and_integrate_samples(
        self,
        accum: dict[str, np.ndarray],
        matrix_store: dict[str, dict[str, np.ndarray]],
    ) -> None:
        """Read sample files concurrently and update outputs immediately."""
        worker_count = min(self.config.threads, len(self.sources))
        eprint(
            "Read full evidence tables: "
            f"samples={len(self.sources)}, threads={worker_count}, "
            f"ORFs={self.evidence_orf_count:,}."
        )
        source_index = {
            source.path: index
            for index, source in enumerate(self.sources)
        }
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            future_to_source = {
                executor.submit(self._read_sample_data, source): source
                for source in self.sources
            }
            completed = 0
            for future in as_completed(future_to_source):
                source = future_to_source[future]
                data = future.result()
                sample_index = source_index[source.path]
                self._accumulate_sample(
                    data=data,
                    sample_index=sample_index,
                    accum=accum,
                    matrix_store=matrix_store,
                )
                completed += 1
                progress_print(
                    "integrated sample="
                    f"{data.sample} ({completed}/{len(self.sources)})."
                )

    def _accumulate_sample(
        self,
        data: SampleData,
        sample_index: int,
        accum: dict[str, np.ndarray],
        matrix_store: dict[str, dict[str, np.ndarray]],
    ) -> None:
        """Update ORF-level accumulators from one aligned sample."""
        evidence = data.translation_code
        rpf = np.nan_to_num(data.numeric["rpf_sum"], nan=0.0)
        valid_profile = ~data.profile_invalid

        capture_codes = np.asarray(
            [TRANSLATION_RANK[label] for label in self.config.capture_labels]
        )
        pass_codes = np.asarray(
            [TRANSLATION_RANK[label] for label in self.config.pass_labels]
        )
        captured = np.isin(evidence, capture_codes) & (rpf > 0) & valid_profile
        passed = np.isin(evidence, pass_codes) & (rpf > 0) & valid_profile
        accum["capture_count"] += captured
        accum["pass_count"] += passed
        accum["invalid_count"] += data.profile_invalid
        accum["no_count"] += evidence == TRANSLATION_RANK["NoEvidence"]
        accum["low_count"] += evidence == TRANSLATION_RANK["LowConfidence"]
        accum["medium_count"] += evidence == TRANSLATION_RANK["MediumConfidence"]
        accum["high_count"] += evidence == TRANSLATION_RANK["HighConfidence"]

        accum["rpf_sum"] += rpf
        accum["rpf_max"] = np.maximum(accum["rpf_max"], rpf)
        for key, metric in (
            ("covered_codon_ratio_max", "covered_codon_ratio"),
            ("frame0_ratio_max", "frame0_ratio"),
            ("periodicity_score_max", "periodicity_score"),
            ("release_ratio_max", "release_ratio"),
        ):
            values = np.nan_to_num(data.numeric[metric], nan=0.0)
            accum[key] = np.maximum(accum[key], values)

        periodicity = np.nan_to_num(
            data.numeric["periodicity_score"],
            nan=-np.inf,
        )
        coverage = np.nan_to_num(
            data.numeric["covered_codon_ratio"],
            nan=-np.inf,
        )
        frame0 = np.nan_to_num(
            data.numeric["frame0_ratio"],
            nan=-np.inf,
        )
        better = (
            (evidence > accum["best_rank"])
            | ((evidence == accum["best_rank"]) & (rpf > accum["best_rpf"]))
            | (
                (evidence == accum["best_rank"])
                & (rpf == accum["best_rpf"])
                & (periodicity > accum["best_periodicity"])
            )
            | (
                (evidence == accum["best_rank"])
                & (rpf == accum["best_rpf"])
                & (periodicity == accum["best_periodicity"])
                & (coverage > accum["best_coverage"])
            )
            | (
                (evidence == accum["best_rank"])
                & (rpf == accum["best_rpf"])
                & (periodicity == accum["best_periodicity"])
                & (coverage == accum["best_coverage"])
                & (frame0 > accum["best_frame0"])
            )
        )
        accum["best_rank"][better] = evidence[better]
        accum["best_rpf"][better] = rpf[better]
        accum["best_periodicity"][better] = periodicity[better]
        accum["best_coverage"][better] = coverage[better]
        accum["best_frame0"][better] = frame0[better]
        accum["best_sample_index"][better] = sample_index
        accum["best_translation_code"][better] = evidence[better]
        for output_key, metric in (
            ("best_rpf_sum", "rpf_sum"),
            ("best_covered_codon_ratio", "covered_codon_ratio"),
            ("best_frame0_density", "frame0_density"),
            ("best_frame1_density", "frame1_density"),
            ("best_frame2_density", "frame2_density"),
            ("best_frame0_ratio", "frame0_ratio"),
            ("best_frame1_ratio", "frame1_ratio"),
            ("best_frame2_ratio", "frame2_ratio"),
            ("best_periodicity_score", "periodicity_score"),
            ("best_release_ratio", "release_ratio"),
        ):
            accum[output_key][better] = data.numeric[metric][better]
        accum["best_periodicity_code"][better] = data.periodicity_code[better]
        accum["best_release_code"][better] = data.release_code[better]
        accum["best_shape_code"][better] = data.shape_code[better]

        for group, metrics in matrix_store.items():
            for metric, matrix in metrics.items():
                if metric == "translation_evidence":
                    matrix[:, sample_index] = evidence
                else:
                    matrix[:, sample_index] = np.nan_to_num(
                        data.numeric[metric],
                        nan=0.0,
                    ).astype(np.float32, copy=False)

    def _final_mask(self, accum: dict[str, np.ndarray]) -> np.ndarray:
        """Return the final output mask after optional reliable filtering."""
        if self.config.reliable_only:
            return accum["pass_count"] > 0
        return np.ones(self.evidence_orf_count, dtype=bool)

    def _consensus(self, accum: dict[str, np.ndarray]) -> np.ndarray:
        """Calculate cross-sample consensus evidence labels."""
        high = accum["high_count"]
        medium = accum["medium_count"]
        low = accum["low_count"]
        consensus = np.full(
            self.evidence_orf_count,
            "NoEvidence",
            dtype=object,
        )
        high_mask = (high >= 2) | ((high >= 1) & (medium >= 1))
        medium_mask = ~high_mask & ((high >= 1) | (medium >= 1))
        low_mask = ~high_mask & ~medium_mask & (low >= 1)
        consensus[high_mask] = "HighConfidence"
        consensus[medium_mask] = "MediumConfidence"
        consensus[low_mask] = "LowConfidence"
        return consensus

    def _multi_sample_status(self, accum: dict[str, np.ndarray]) -> np.ndarray:
        """Calculate the compatibility multi-sample status."""
        captured = accum["capture_count"]
        passed = accum["pass_count"]
        status = np.full(
            self.evidence_orf_count,
            "NoEvidence",
            dtype=object,
        )
        status[captured > 0] = "CapturedButWeak"
        status[passed == 1] = "SingleSample"
        status[passed >= self.config.excellent_min_samples] = "Excellent"
        return status

    def _build_integrated_table(
        self,
        accum: dict[str, np.ndarray],
    ) -> pd.DataFrame:
        """Build the compact integrated ORF table."""
        if self.static_table is None:
            raise RuntimeError("Static table is unavailable.")
        sample_count = len(self.sources)
        best_sample = np.asarray(
            [
                self.sample_names[index] if index >= 0 else "."
                for index in accum["best_sample_index"]
            ],
            dtype=object,
        )
        best_translation = np.asarray(
            [
                TRANSLATION_DECODE[code] if code >= 0 else "NoEvidence"
                for code in accum["best_translation_code"]
            ],
            dtype=object,
        )
        columns: dict[str, Any] = {
            "total_sample_count": np.full(
                self.evidence_orf_count,
                sample_count,
                dtype=np.int16,
            ),
            "captured_sample_count": accum["capture_count"],
            "pass_sample_count": accum["pass_count"],
            "captured_sample_ratio": accum["capture_count"] / sample_count,
            "pass_sample_ratio": accum["pass_count"] / sample_count,
            "no_evidence_sample_count": accum["no_count"],
            "low_confidence_sample_count": accum["low_count"],
            "medium_confidence_sample_count": accum["medium_count"],
            "high_confidence_sample_count": accum["high_count"],
            "invalid_profile_sample_count": accum["invalid_count"],
            "rpf_sum_sum": accum["rpf_sum"],
            "rpf_sum_mean": accum["rpf_sum"] / sample_count,
            "rpf_sum_max": accum["rpf_max"],
            "covered_codon_ratio_max": accum["covered_codon_ratio_max"],
            "frame0_ratio_max": accum["frame0_ratio_max"],
            "periodicity_score_max": accum["periodicity_score_max"],
            "release_ratio_max": accum["release_ratio_max"],
            "best_sample": best_sample,
            "best_translation_evidence": best_translation,
            "best_rpf_sum": accum["best_rpf_sum"],
            "best_covered_codon_ratio": accum[
                "best_covered_codon_ratio"
            ],
            "best_frame0_density": accum["best_frame0_density"],
            "best_frame1_density": accum["best_frame1_density"],
            "best_frame2_density": accum["best_frame2_density"],
            "best_frame0_ratio": accum["best_frame0_ratio"],
            "best_frame1_ratio": accum["best_frame1_ratio"],
            "best_frame2_ratio": accum["best_frame2_ratio"],
            "best_periodicity_score": accum["best_periodicity_score"],
            "best_periodicity_label": PERIODICITY_DECODE[
                accum["best_periodicity_code"]
            ],
            "best_release_ratio": accum["best_release_ratio"],
            "best_release_label": RELEASE_DECODE[
                accum["best_release_code"]
            ],
            "best_coverage_shape": SHAPE_DECODE[
                accum["best_shape_code"]
            ],
            "consensus_translation_evidence": self._consensus(accum),
            "multi_sample_status": self._multi_sample_status(accum),
        }
        return pd.concat(
            [
                self.static_table.reset_index(drop=True),
                pd.DataFrame(columns),
            ],
            axis=1,
            copy=False,
        )

    def _write_dataframe(self, table: pd.DataFrame, path: str | Path) -> None:
        """Write one DataFrame atomically."""
        with _AtomicTableWriter(path) as writer:
            for start in range(0, len(table), OUTPUT_CHUNK_ROWS):
                writer.write(table.iloc[start:start + OUTPUT_CHUNK_ROWS])

    def _write_matrix(
        self,
        group: str,
        values: dict[str, np.ndarray],
        final_mask: np.ndarray,
        path: str | Path,
    ) -> None:
        """Write one metric matrix in bounded DataFrame blocks."""
        if self.static_table is None:
            raise RuntimeError("Static table is unavailable.")
        selected = np.flatnonzero(final_mask)
        static = self.static_table.loc[final_mask].reset_index(drop=True)
        static_columns = [
            column
            for column in STATIC_COLUMNS
            if column in static.columns
        ]
        with _AtomicTableWriter(path) as writer:
            for start in range(0, len(selected), OUTPUT_CHUNK_ROWS):
                stop = min(start + OUTPUT_CHUNK_ROWS, len(selected))
                output = static.loc[start:stop - 1, static_columns].copy()
                row_positions = selected[start:stop]
                for metric, matrix in values.items():
                    if group in {
                        "rpf_sum",
                        "covered_codon_ratio",
                        "translation_evidence",
                    }:
                        names = self.sample_names
                    else:
                        names = [
                            f"{sample}__{metric}"
                            for sample in self.sample_names
                        ]
                    if group in {
                        "rpf_sum",
                        "covered_codon_ratio",
                        "translation_evidence",
                    }:
                        for sample_index, name in enumerate(names):
                            column_values = matrix[row_positions, sample_index]
                            if metric == "translation_evidence":
                                output[name] = [
                                    TRANSLATION_DECODE[value]
                                    if value >= 0
                                    else "Missing"
                                    for value in column_values
                                ]
                            else:
                                output[name] = column_values
                    else:
                        for sample_index, name in enumerate(names):
                            output[name] = matrix[row_positions, sample_index]
                writer.write(output)

    def _write_filtered_genepred(
        self,
        input_path: str | Path,
        output_path: str | Path,
        retained_ids: set[str],
    ) -> None:
        """Filter ORF genePred records by the final ORF set."""
        found: set[str] = set()
        duplicates: set[str] = set()
        final_path = Path(output_path)
        temp_path = final_path.with_name(final_path.name + ".tmp")
        final_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with smart_open(input_path, "rt") as source, smart_open(
                temp_path,
                "wt",
            ) as target:
                for raw_line in source:
                    if raw_line.startswith("#"):
                        target.write(raw_line)
                        continue
                    text = raw_line.rstrip("\r\n")
                    if not text:
                        continue
                    orf_id = text.split("\t", 1)[0]
                    if orf_id not in retained_ids:
                        continue
                    if orf_id in found:
                        duplicates.add(orf_id)
                        continue
                    found.add(orf_id)
                    target.write(raw_line)
            missing = retained_ids.difference(found)
            if duplicates:
                raise ValueError(
                    "Duplicate retained ORFs in genePred: "
                    + ", ".join(sorted(duplicates)[:10])
                )
            if missing:
                raise ValueError(
                    f"Filtered genePred is missing {len(missing):,} retained "
                    "ORFs. Examples: "
                    + ", ".join(sorted(missing)[:10])
                )
            os.replace(temp_path, final_path)
        except BaseException:
            temp_path.unlink(missing_ok=True)
            raise


def run_integration(args: Any) -> IntegrationResult:
    """Run integration from parsed command-line arguments.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Integration result.
    """
    config = IntegrationConfig(
        input_files=tuple(args.input),
        output_integrated=args.output_integrated,
        output_matrix=args.output_matrix,
        output_density_matrix=args.output_density_matrix,
        output_ratio_matrix=args.output_ratio_matrix,
        output_covered_codon_ratio_matrix=(
            args.output_covered_codon_ratio_matrix
        ),
        output_translation_evidence_matrix=(
            args.output_translation_evidence_matrix
        ),
        orf_genepred=args.orf_genepred,
        capture_labels=split_items(args.capture_labels)
        or (
            "LowConfidence",
            "MediumConfidence",
            "HighConfidence",
        ),
        pass_labels=split_items(args.pass_labels)
        or ("MediumConfidence", "HighConfidence"),
        excellent_min_samples=args.excellent_min_samples,
        reliable_only=args.reliable_only,
        threads=args.threads,
    )
    return SmorfEvidenceIntegrator(config).run()
