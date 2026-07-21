#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-21
# Version: 0.2.8.18
# Function: Calculate and visualize transcript-aware meta-codon RPF density profiles.
# Input: RPF density file in JSONL or TXT format and optional codon/transcript lists.
# Output: Meta-codon density tables, sequence-context tables, and figures.

"""Core functions for RiboParser meta-codon analysis.

The module reads compact JSONL and legacy TXT RPF density files through
``RPFs.RPFData``. Target codon windows are extracted independently within each
transcript using consecutive ``from_tis`` coordinates, which prevents windows
from crossing transcript boundaries or spanning missing codon positions.

When ``--scale`` is enabled, each codon-level density is divided by the mean CDS
density of the corresponding transcript and sample. Therefore, a value of 1
represents the transcript-average CDS density, consistent with the codon
occupancy definition used by RiboParser.

Workflow
--------
1. Import and validate target codons.
2. Import JSONL or TXT RPF density through the shared RPF reader.
3. Select high-expression transcripts using CDS RPF counts.
4. Optionally convert density to RPM.
5. Optionally normalize each transcript by its mean CDS density.
6. Extract complete transcript-local windows around each target motif.
7. Optionally remove overlapping target windows.
8. Aggregate, smooth, output, and visualize meta-codon profiles.
"""

from __future__ import annotations

from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Any

import matplotlib

matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

from . import RPFs


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
RPM_SCALE = 1_000_000.0


@dataclass(slots=True)
class _GeneProfile:
    """Store one transcript-local codon profile."""

    name: str
    coordinates: np.ndarray
    codons: np.ndarray
    regions: np.ndarray
    density: np.ndarray


@dataclass(slots=True)
class _CodonResult:
    """Store one meta-codon calculation result."""

    codon: str
    raw_site_count: int
    retained_site_count: int
    density: pd.DataFrame
    sequence: pd.DataFrame


class MetaCodon:
    """Calculate transcript-aware meta-codon RPF density profiles.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_Meta_Codon``.
    """

    def __init__(self, args: Any):
        # Input and output.
        self.ribo = args.rpf
        self.gene = args.list
        self.output = args.output
        self.codon = args.codon

        # RPF import and filtering.
        self.frame = str(args.frame)
        self.rpf_num = float(args.min)
        self.tis = int(args.tis)
        self.tts = int(args.tts)
        self.norm = bool(args.normal)
        self.thread = max(1, int(getattr(args, "thread", 1)))

        # Meta-codon parameters.
        self.around = int(args.around)
        self.scale = bool(args.scale)
        self.unique = bool(args.unique)
        self.smooth = args.smooth
        self.fig = bool(args.fig)

        # Imported data.
        self.rpf_data: RPFs.RPFData | None = None
        self.file_format: str | None = None
        self.sample_name: list[str] = []
        self.sample_num = 0
        self.total_rpf_num: pd.Series | None = None
        self.gene_rpf_sum: pd.DataFrame | None = None
        self.high_gene: pd.Index = pd.Index([])
        self.high_rpf: pd.DataFrame | None = None
        self._gene_profiles: list[_GeneProfile] = []

        # Codon and output data.
        self.codon_list: pd.DataFrame | None = None
        self.density_df: OrderedDict[str, list[Any]] = OrderedDict()
        self.sequence_df: OrderedDict[str, list[Any]] = OrderedDict()

        # Smoothing configuration.
        self.smooth_window: int | None = None
        self.smooth_power: int | None = None

    # ------------------------------------------------------------------
    # Input and validation
    # ------------------------------------------------------------------
    @staticmethod
    def _normalize_codon(value: Any) -> str:
        """Return an uppercase DNA codon or codon motif."""
        return str(value).strip().upper().replace("U", "T").replace(" ", "")

    @staticmethod
    def _is_valid_codon_motif(codon: str) -> bool:
        """Return whether a codon motif is non-empty and 3-nt periodic."""
        return (
            bool(codon)
            and len(codon) % 3 == 0
            and not set(codon).difference("ACGT")
        )

    def import_codon(self) -> None:
        """Import, normalize, validate, and deduplicate target codons."""
        if self.codon:
            codon_table = pd.read_csv(
                self.codon,
                sep="\t",
                header=None,
                comment="#",
                dtype=str,
                usecols=[0],
            )
            codon_table.columns = ["codon"]
        else:
            _, annotation = RPFs.codon_table()
            codon_table = annotation.reset_index().loc[:, ["codon"]]

        codon_table["codon"] = codon_table["codon"].map(self._normalize_codon)
        codon_table = codon_table.loc[
            codon_table["codon"].str.lower().ne("codon")
        ].copy()

        valid_mask = codon_table["codon"].map(self._is_valid_codon_motif)
        invalid_codons = codon_table.loc[~valid_mask, "codon"].tolist()
        if invalid_codons:
            print(
                "Skip invalid codon motif(s): {codons}.".format(
                    codons=", ".join(invalid_codons)
                ),
                flush=True,
            )

        codon_table = (
            codon_table.loc[valid_mask, ["codon"]]
            .drop_duplicates(subset=["codon"], keep="first")
            .reset_index(drop=True)
        )
        if codon_table.empty:
            raise ValueError("No valid codon motif remains after validation.")

        self.codon_list = codon_table
        print(
            "Imported {count:,} unique codon motif(s).".format(
                count=len(self.codon_list)
            ),
            flush=True,
        )

    def _validate_imported_table(self, table: pd.DataFrame) -> None:
        """Validate the imported codon-level density table."""
        if table.empty:
            raise ValueError("RPF density table is empty after import and filtering.")

        missing = [
            column
            for column in BASE_COLUMNS + self.sample_name
            if column not in table.columns
        ]
        if missing:
            raise ValueError(
                "RPF density table is missing required column(s): "
                + ", ".join(missing)
            )

        if not self.sample_name:
            raise ValueError("No sample density columns were detected.")

    def _select_high_expression_transcripts(
        self,
        merged_rpf: pd.DataFrame,
    ) -> pd.DataFrame:
        """Select high-expression transcripts using mean CDS counts."""
        cds_mask = merged_rpf["region"].astype(str).str.lower().eq("cds")
        cds_table = merged_rpf.loc[
            cds_mask,
            ["name"] + self.sample_name,
        ].copy()
        if cds_table.empty:
            raise ValueError("No CDS codons remain after TIS/TTS filtering.")

        self.gene_rpf_sum = (
            cds_table.groupby("name", sort=False)[self.sample_name].sum()
        )
        mean_cds_count = self.gene_rpf_sum.mean(axis=1)
        self.high_gene = self.gene_rpf_sum.index[
            mean_cds_count >= self.rpf_num
        ]

        if len(self.high_gene) == 0:
            raise ValueError(
                "No transcript passed the minimum mean CDS RPF threshold: "
                f"{self.rpf_num:g}."
            )

        high_rpf = merged_rpf.loc[
            merged_rpf["name"].isin(self.high_gene),
            BASE_COLUMNS + self.sample_name,
        ].copy()
        high_rpf = high_rpf.astype(
            {sample: "float64" for sample in self.sample_name}
        )

        print(
            "Retained {retained:,} of {total:,} transcript(s) with mean CDS "
            "RPF count >= {minimum:g}.".format(
                retained=len(self.high_gene),
                total=len(self.gene_rpf_sum),
                minimum=self.rpf_num,
            ),
            flush=True,
        )
        return high_rpf

    def _normalize_to_rpm(self, table: pd.DataFrame) -> None:
        """Convert sample density columns to RPM in place."""
        if not self.norm:
            return
        if self.total_rpf_num is None:
            raise ValueError("Total RPF counts are unavailable for RPM normalization.")

        for sample in self.sample_name:
            total = float(self.total_rpf_num.get(sample, 0.0))
            if total > 0:
                table.loc[:, sample] *= RPM_SCALE / total
            else:
                table.loc[:, sample] = 0.0
                print(
                    f"Warning: sample {sample} has zero total RPF count; "
                    "its RPM density was set to zero.",
                    flush=True,
                )

    def _normalize_by_gene_mean(self, table: pd.DataFrame) -> None:
        """Normalize each position by transcript-specific mean CDS density."""
        if not self.scale:
            return

        cds_mask = table["region"].astype(str).str.lower().eq("cds")
        gene_mean = (
            table.loc[cds_mask, ["name"] + self.sample_name]
            .groupby("name", sort=False)[self.sample_name]
            .mean()
        )

        gene_names = table["name"].astype(str).to_numpy()
        denominator_values = gene_mean.reindex(gene_names).to_numpy(dtype=float)
        values = table.loc[:, self.sample_name].to_numpy(dtype=float)

        normalized = np.zeros_like(values, dtype=float)
        valid = np.isfinite(denominator_values) & (denominator_values > 0)
        np.divide(
            values,
            denominator_values,
            out=normalized,
            where=valid,
        )
        table.loc[:, self.sample_name] = normalized

    def _build_gene_profiles(self, table: pd.DataFrame) -> None:
        """Convert the retained table into compact transcript-local arrays."""
        profiles: list[_GeneProfile] = []

        for gene_name, gene_table in table.groupby("name", sort=False):
            gene_table = (
                gene_table.sort_values(
                    ["from_tis", "now_nt"],
                    kind="stable",
                )
                .drop_duplicates(subset=["from_tis"], keep="first")
                .reset_index(drop=True)
            )

            coordinates = pd.to_numeric(
                gene_table["from_tis"],
                errors="coerce",
            ).to_numpy(dtype=float)
            valid_coordinate = np.isfinite(coordinates)
            if not valid_coordinate.all():
                gene_table = gene_table.loc[valid_coordinate, :].reset_index(
                    drop=True
                )
                coordinates = coordinates[valid_coordinate]

            if gene_table.empty:
                continue

            density = (
                gene_table.loc[:, self.sample_name]
                .apply(pd.to_numeric, errors="coerce")
                .fillna(0.0)
                .to_numpy(dtype=float)
            )
            profiles.append(
                _GeneProfile(
                    name=str(gene_name),
                    coordinates=coordinates.astype(np.int64, copy=False),
                    codons=(
                        gene_table["codon"]
                        .astype(str)
                        .str.upper()
                        .str.replace("U", "T", regex=False)
                        .to_numpy(dtype="U3")
                    ),
                    regions=(
                        gene_table["region"]
                        .astype(str)
                        .str.lower()
                        .to_numpy(dtype="U4")
                    ),
                    density=density,
                )
            )

        if not profiles:
            raise ValueError("No transcript-local density profile could be built.")

        self._gene_profiles = profiles

    def import_rpf(self) -> None:
        """Import JSONL or TXT RPF density through the shared RPF reader."""
        self.rpf_data = RPFs.RPFData.from_file(
            rpf_file=self.ribo,
            sample_name=None,
            gene=self.gene,
            tis=self.tis,
            tts=self.tts,
            json_thread=self.thread,
        )
        self.file_format = str(self.rpf_data.file_format)
        self.sample_name = list(self.rpf_data.sample_name)
        self.sample_num = int(self.rpf_data.sample_num)
        self.total_rpf_num = self.rpf_data.total_rpf_num.astype(float)

        merged_rpf = self.rpf_data.get_frame(frame=self.frame)
        merged_rpf = RPFs.shift_site(
            merged_rpf=merged_rpf,
            sample_name=self.sample_name,
            shift_num=RPFs.set_codon_shift("A"),
        )
        self._validate_imported_table(merged_rpf)

        for sample in self.sample_name:
            merged_rpf.loc[:, sample] = pd.to_numeric(
                merged_rpf[sample],
                errors="coerce",
            ).fillna(0.0)

        high_rpf = self._select_high_expression_transcripts(merged_rpf)
        self._normalize_to_rpm(high_rpf)
        self._normalize_by_gene_mean(high_rpf)
        self._build_gene_profiles(high_rpf)

        # Keep a compact public reference without duplicating the large table.
        self.high_rpf = None

        print(
            "Imported RPF density: format={fmt}, samples={samples:,}, "
            "retained transcripts={genes:,}, retained rows={rows:,}.".format(
                fmt=self.file_format,
                samples=self.sample_num,
                genes=len(self._gene_profiles),
                rows=sum(
                    profile.coordinates.size
                    for profile in self._gene_profiles
                ),
            ),
            flush=True,
        )

    # ------------------------------------------------------------------
    # Smoothing
    # ------------------------------------------------------------------
    def smooth_rpf_density(self) -> None:
        """Parse and validate optional Savitzky-Golay smoothing parameters.

        Smoothing is applied to the final aggregated meta-codon profile rather
        than to every transcript. This avoids repeated whole-table copies and
        prevents target peaks from being mixed between unrelated transcript
        groups.
        """
        self.smooth_window = None
        self.smooth_power = None

        if self.smooth is None or str(self.smooth).strip() == "":
            return

        fields = [item.strip() for item in str(self.smooth).split(",")]
        if len(fields) != 2:
            raise ValueError(
                "--smooth must contain two integers: window,polynomial_order."
            )

        try:
            window = int(fields[0])
            power = int(fields[1])
        except ValueError as error:
            raise ValueError(
                "--smooth must contain two integers: window,polynomial_order."
            ) from error

        if window < 3:
            raise ValueError("The smoothing window must be >= 3.")
        if window % 2 == 0:
            raise ValueError("The smoothing window must be an odd integer.")
        if power < 0:
            raise ValueError("The smoothing polynomial order must be >= 0.")
        if power >= window:
            raise ValueError(
                "The smoothing polynomial order must be smaller than the window."
            )

        self.smooth_window = window
        self.smooth_power = power
        print(
            "Configured Savitzky-Golay smoothing: window={window}, power={power}.".format(
                window=window,
                power=power,
            ),
            flush=True,
        )

    def _smooth_profile(self, values: np.ndarray) -> np.ndarray:
        """Smooth one aggregated meta-codon matrix."""
        if self.smooth_window is None or self.smooth_power is None:
            return values

        profile_length = int(values.shape[0])
        window = min(self.smooth_window, profile_length)
        if window % 2 == 0:
            window -= 1

        minimum_window = self.smooth_power + 2
        if minimum_window % 2 == 0:
            minimum_window += 1

        if window < minimum_window:
            print(
                "Warning: skip smoothing because the meta-codon profile is "
                "shorter than the requested polynomial configuration.",
                flush=True,
            )
            return values

        smoothed = savgol_filter(
            values,
            window_length=window,
            polyorder=self.smooth_power,
            axis=0,
            mode="interp",
        )
        # RPF density cannot be negative; remove Savitzky-Golay edge artifacts.
        return np.clip(smoothed, 0.0, None)

    # ------------------------------------------------------------------
    # Transcript-local motif matching
    # ------------------------------------------------------------------
    @staticmethod
    def _motif_codons(codon: str) -> np.ndarray:
        """Split one DNA motif into codon tokens."""
        return np.asarray(
            [codon[index:index + 3] for index in range(0, len(codon), 3)],
            dtype="U3",
        )

    @staticmethod
    def _complete_window_indices(
        profile: _GeneProfile,
        start_coordinate: int,
        motif_length: int,
        around: int,
    ) -> np.ndarray | None:
        """Return indices for one complete consecutive transcript-local window."""
        expected_coordinates = np.arange(
            start_coordinate - around,
            start_coordinate + motif_length + around,
            dtype=np.int64,
        )
        indices = np.searchsorted(profile.coordinates, expected_coordinates)

        if np.any(indices >= profile.coordinates.size):
            return None
        if not np.array_equal(
            profile.coordinates[indices],
            expected_coordinates,
        ):
            return None
        if not np.all(profile.regions[indices] == "cds"):
            return None
        return indices

    @staticmethod
    def _isolated_site_mask(
        starts: np.ndarray,
        extraction_span: int,
    ) -> np.ndarray:
        """Return sites whose extraction windows do not overlap neighboring sites."""
        starts = np.asarray(starts, dtype=np.int64)
        if starts.size <= 1:
            return np.ones(starts.size, dtype=bool)

        delta = np.diff(starts)
        left_distance = np.concatenate(([np.inf], delta))
        right_distance = np.concatenate((delta, [np.inf]))
        return (
            (left_distance > extraction_span)
            & (right_distance > extraction_span)
        )

    def _find_gene_occurrences(
        self,
        profile: _GeneProfile,
        motif_codons: np.ndarray,
        apply_unique: bool,
    ) -> list[tuple[int, np.ndarray]]:
        """Find complete motif windows in one transcript."""
        motif_length = int(motif_codons.size)
        candidate_indices = np.flatnonzero(profile.codons == motif_codons[0])
        occurrences: list[tuple[int, np.ndarray]] = []

        for candidate_index in candidate_indices:
            start_coordinate = int(profile.coordinates[candidate_index])
            motif_coordinates = np.arange(
                start_coordinate,
                start_coordinate + motif_length,
                dtype=np.int64,
            )
            motif_indices = np.searchsorted(
                profile.coordinates,
                motif_coordinates,
            )

            if np.any(motif_indices >= profile.coordinates.size):
                continue
            if not np.array_equal(
                profile.coordinates[motif_indices],
                motif_coordinates,
            ):
                continue
            if not np.array_equal(
                profile.codons[motif_indices],
                motif_codons,
            ):
                continue
            if not np.all(profile.regions[motif_indices] == "cds"):
                continue

            window_indices = self._complete_window_indices(
                profile=profile,
                start_coordinate=start_coordinate,
                motif_length=motif_length,
                around=self.around,
            )
            if window_indices is None:
                continue

            occurrences.append((start_coordinate, window_indices))

        if not apply_unique or len(occurrences) <= 1:
            return occurrences

        starts = np.asarray(
            [start for start, _ in occurrences],
            dtype=np.int64,
        )
        extraction_span = 2 * self.around + motif_length - 1
        keep = self._isolated_site_mask(starts, extraction_span)
        return [
            occurrence
            for occurrence, retained in zip(occurrences, keep)
            if retained
        ]

    @staticmethod
    def _sequence_row(
        profile: _GeneProfile,
        start_coordinate: int,
        window_indices: np.ndarray,
        relative_positions: np.ndarray,
    ) -> dict[Any, Any]:
        """Build one sequence-context output row."""
        row: dict[Any, Any] = {
            "name": profile.name,
            "site_from_tis": int(start_coordinate),
        }
        for relative_position, codon in zip(
            relative_positions,
            profile.codons[window_indices],
        ):
            row[int(relative_position)] = str(codon)
        return row

    def _calculate_codon(self, codon: str) -> _CodonResult | None:
        """Calculate one single- or multi-codon meta profile."""
        motif_codons = self._motif_codons(codon)
        motif_length = int(motif_codons.size)
        relative_positions = np.arange(
            -self.around,
            self.around + motif_length,
            dtype=int,
        )

        density_sum = np.zeros(
            (relative_positions.size, self.sample_num),
            dtype=float,
        )
        raw_site_count = 0
        retained_site_count = 0
        sequence_rows: list[dict[Any, Any]] = []

        for profile in self._gene_profiles:
            raw_occurrences = self._find_gene_occurrences(
                profile=profile,
                motif_codons=motif_codons,
                apply_unique=False,
            )
            if self.unique and len(raw_occurrences) > 1:
                starts = np.asarray(
                    [start for start, _ in raw_occurrences],
                    dtype=np.int64,
                )
                extraction_span = 2 * self.around + motif_length - 1
                keep = self._isolated_site_mask(starts, extraction_span)
                retained_occurrences = [
                    occurrence
                    for occurrence, retained in zip(raw_occurrences, keep)
                    if retained
                ]
            else:
                retained_occurrences = raw_occurrences

            raw_site_count += len(raw_occurrences)
            retained_site_count += len(retained_occurrences)

            for start_coordinate, window_indices in retained_occurrences:
                density_sum += profile.density[window_indices, :]
                sequence_rows.append(
                    self._sequence_row(
                        profile=profile,
                        start_coordinate=start_coordinate,
                        window_indices=window_indices,
                        relative_positions=relative_positions,
                    )
                )

        if retained_site_count == 0:
            return None

        density_values = density_sum / float(retained_site_count)
        density_values = self._smooth_profile(density_values)

        density_table = pd.DataFrame(
            density_values,
            index=relative_positions,
            columns=self.sample_name,
        )
        density_table.index.name = "Codon"

        if self.frame == "all":
            nucleotide_position = relative_positions * 3
        else:
            nucleotide_position = relative_positions * 3 + int(self.frame)

        density_table.insert(0, "Frame", self.frame)
        density_table.insert(0, "Nucleotide", nucleotide_position)

        sequence_table = pd.DataFrame.from_records(sequence_rows)
        sequence_columns: list[Any] = [
            "name",
            "site_from_tis",
            *relative_positions.tolist(),
        ]
        sequence_table = sequence_table.reindex(columns=sequence_columns)
        sequence_table.index = [
            f"{row['name']}:{int(row['site_from_tis'])}"
            for _, row in sequence_table.iterrows()
        ]
        sequence_table.index.name = "Site"

        return _CodonResult(
            codon=codon,
            raw_site_count=int(raw_site_count),
            retained_site_count=int(retained_site_count),
            density=density_table,
            sequence=sequence_table,
        )

    def _store_result(self, result: _CodonResult | None) -> None:
        """Store one codon result in backward-compatible dictionaries."""
        if result is None:
            return

        self.density_df[result.codon] = [
            result.raw_site_count,
            result.retained_site_count,
            result.density,
        ]
        self.sequence_df[result.codon] = [
            result.raw_site_count,
            result.retained_site_count,
            result.sequence,
        ]

    def single_codon(self, codon: str) -> None:
        """Calculate and store one single-codon profile."""
        codon = self._normalize_codon(codon)
        if len(codon) != 3:
            raise ValueError("single_codon requires one 3-nt codon.")
        self._store_result(self._calculate_codon(codon))

    def multiple_codon(self, codon: str, codon_num: int | None = None) -> None:
        """Calculate and store one multi-codon motif profile."""
        codon = self._normalize_codon(codon)
        detected_num = len(codon) // 3
        if detected_num <= 1:
            raise ValueError("multiple_codon requires a motif longer than 3 nt.")
        if codon_num is not None and int(codon_num) != detected_num:
            raise ValueError("codon_num does not match the motif length.")
        self._store_result(self._calculate_codon(codon))

    def retrieve_codon_density(self) -> None:
        """Calculate all requested codon or codon-motif profiles."""
        if self.codon_list is None:
            raise ValueError("Codon motifs have not been imported.")
        if not self._gene_profiles:
            raise ValueError("RPF density has not been imported.")

        codons = self.codon_list["codon"].tolist()
        self.density_df.clear()
        self.sequence_df.clear()

        worker_count = min(self.thread, max(1, len(codons)))
        if worker_count == 1:
            results = [self._calculate_codon(codon) for codon in codons]
        else:
            print(
                f"Calculate meta-codon profiles with {worker_count} codon workers.",
                flush=True,
            )
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                results = list(executor.map(self._calculate_codon, codons))

        for codon, result in zip(codons, results):
            if result is None:
                print(
                    f"Codon {codon}: no complete retained window.",
                    flush=True,
                )
                continue
            self._store_result(result)
            print(
                "Codon {codon}: raw sites={raw:,}, retained sites={retained:,}.".format(
                    codon=codon,
                    raw=result.raw_site_count,
                    retained=result.retained_site_count,
                ),
                flush=True,
            )

        if not self.density_df:
            raise ValueError(
                "No complete meta-codon window was retained for any target motif."
            )

    # Backward-compatible misspelled method used by the current CLI entry.
    reterieve_codon_density = retrieve_codon_density

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    def output_meta_codon_density(self) -> None:
        """Write one meta-codon density table per target motif."""
        for codon, message in self.density_df.items():
            raw_site_count, retained_site_count, density_table = message
            output_file = (
                f"{self.output}_{codon}_{raw_site_count}_{retained_site_count}"
                "_meta_density.txt"
            )
            density_table.round(6).to_csv(
                output_file,
                sep="\t",
                index=True,
            )

    def output_meta_codon_seq(self) -> None:
        """Write one sequence-context table per target motif."""
        for codon, message in self.sequence_df.items():
            raw_site_count, retained_site_count, sequence_table = message
            output_file = (
                f"{self.output}_{codon}_{raw_site_count}_{retained_site_count}"
                "_meta_sequence.txt"
            )
            sequence_table.to_csv(
                output_file,
                sep="\t",
                index=True,
            )

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------
    @staticmethod
    def _plot_ticks(
        minimum: int,
        maximum: int,
        motif_length: int,
    ) -> list[int]:
        """Return readable integer ticks including the target motif."""
        span = maximum - minimum
        if span <= 12:
            ticks = list(range(minimum, maximum + 1))
        else:
            step = max(1, int(np.ceil(span / 8)))
            ticks = list(range(minimum, maximum + 1, step))

        ticks.extend([0, motif_length - 1, minimum, maximum])
        return sorted(
            {
                int(value)
                for value in ticks
                if minimum <= int(value) <= maximum
            }
        )

    def _density_ylabel(self) -> str:
        """Return the density-axis label."""
        if self.scale:
            return "Mean gene-normalized density"
        if self.norm:
            return "Mean RPM density"
        return "Mean RPF density"

    def draw_meta_codon(self) -> None:
        """Draw polished line profiles for all retained target motifs."""
        for codon, message in self.density_df.items():
            raw_site_count, retained_site_count, density_table = message
            motif_length = len(codon) // 3

            x_values = density_table.index.to_numpy(dtype=int)
            sample_table = density_table.loc[:, self.sample_name]

            figure_width = max(7.0, min(12.0, 6.5 + self.sample_num * 0.15))
            fig, ax = plt.subplots(
                figsize=(figure_width, 4.8),
                dpi=150,
            )

            for sample in self.sample_name:
                ax.plot(
                    x_values,
                    sample_table[sample].to_numpy(dtype=float),
                    linewidth=1.4,
                    label=sample,
                )

            ax.axvspan(
                -0.5,
                motif_length - 0.5,
                alpha=0.08,
                linewidth=0,
            )
            ax.axvline(
                0,
                linewidth=0.8,
                linestyle="--",
                alpha=0.65,
            )
            if motif_length > 1:
                ax.axvline(
                    motif_length - 1,
                    linewidth=0.8,
                    linestyle=":",
                    alpha=0.65,
                )

            ax.set_xlim(int(x_values.min()), int(x_values.max()))
            ax.set_xticks(
                self._plot_ticks(
                    minimum=int(x_values.min()),
                    maximum=int(x_values.max()),
                    motif_length=motif_length,
                )
            )
            ax.tick_params(axis="x", labelrotation=0)
            ax.set_xlabel("Relative codon position")
            ax.set_ylabel(self._density_ylabel())
            ax.set_title(
                "{codon} | raw sites={raw:,}, retained sites={retained:,}".format(
                    codon=codon,
                    raw=raw_site_count,
                    retained=retained_site_count,
                )
            )
            ax.grid(axis="y", linewidth=0.4, alpha=0.25)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            values = sample_table.to_numpy(dtype=float)
            finite_values = values[np.isfinite(values)]
            if finite_values.size and np.nanmin(finite_values) >= 0:
                ax.set_ylim(bottom=0)

            if self.sample_num > 0:
                ax.legend(
                    title="Sample",
                    loc="center left",
                    bbox_to_anchor=(1.02, 0.5),
                    frameon=False,
                    fontsize=8,
                    title_fontsize=9,
                )

            fig.tight_layout()
            fig.savefig(
                f"{self.output}_{codon}.pdf",
                bbox_inches="tight",
            )
            fig.savefig(
                f"{self.output}_{codon}.png",
                dpi=300,
                bbox_inches="tight",
            )
            plt.close(fig)
