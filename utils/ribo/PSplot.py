#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-06
# Version: dev004
# Function: Draw all pausing sites from selected genes on sample-resolved RPF profiles.
# Input: rpf_Odd_Ratio site table and compact RPF density JSONL records.
# Output: Gene-level pausing-site figures and optional sample metrics.

"""Ribosome pausing-site visualization utilities."""

from __future__ import annotations

import re
from collections import OrderedDict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("AGG")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import RPFs


FRAME_COLORS = {
    "0": "#D55E00",
    "1": "#0072B2",
    "2": "#009E73",
}
PAUSE_COLOR = "#CC3311"
PAUSE_REGION_COLOR = "#FFF9C4"
NON_PAUSE_COLOR = "#777777"
CONTROL_COLOR = "#4477AA"
TREAT_COLOR = "#EE7733"


class PSplot(object):
    """Draw sample-resolved RPF figures for selected genes.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``rpf_PSplot``.
    """

    def __init__(self, args):
        self.input_file = args.input
        self.rpf_file = args.rpf
        self.output_prefix = args.output
        self.target = args.target
        self.target_list = args.target_list
        self.local_window = 101
        self.pause_score = 20.0
        self.min_site_rpf = 3
        self.min_local_coverage = 0.10
        self.normal = bool(args.normal)
        self.view = args.view
        self.flank = int(args.flank)
        self.plot_transform = args.plot_transform
        self.y_scale = args.y_scale
        self.y_max = args.y_max
        self.pause_region = bool(args.pause_region)
        self.output_format = args.output_format
        self.export_metrics = bool(args.export_metrics)
        self.dpi = int(args.dpi)
        self.figure_width = 12.0
        self.per_sample_height = 1.15
        self.max_height = 20.0
        self.font_size = float(args.font_size)

        self.control: list[str] = []
        self.treat: list[str] = []
        self.samples: list[str] = []
        self.site = "P"
        self.frame = "all"
        self.events: pd.DataFrame | None = None
        self.records: dict[str, dict[str, Any]] = {}
        self.total_rpf: dict[str, float] = {}
        self.metrics: list[dict[str, Any]] = []
        self.output_files: OrderedDict[str, str] = OrderedDict()
        self.target_ids = self._resolve_target_ids()

    @staticmethod
    def _parse_csv(value: Any) -> list[str]:
        """Parse a comma-separated value into unique non-empty strings."""
        if value is None or str(value).strip() == "":
            return []
        return list(dict.fromkeys(item.strip() for item in str(value).split(",") if item.strip()))

    @staticmethod
    def _safe_name(value: Any) -> str:
        """Return a filesystem-safe identifier."""
        return re.sub(r"[^0-9A-Za-z._-]+", "_", str(value)).strip("_") or "event"

    @staticmethod
    def _first_value(table: pd.DataFrame, column: str) -> str | None:
        """Return one consistent non-empty value from a table column."""
        if column not in table.columns:
            return None
        values = table[column].dropna().astype(str).str.strip()
        values = values.loc[values != ""]
        unique_values = values.drop_duplicates().tolist()
        if len(unique_values) > 1:
            raise ValueError(
                "Input table contains inconsistent {column} values: {values}".format(
                    column=column,
                    values=", ".join(unique_values[:5]),
                )
            )
        return unique_values[0] if unique_values else None

    @staticmethod
    def _read_target_file(path: str) -> set[str]:
        """Read target identifiers from the first non-empty column."""
        targets: set[str] = set()
        header_tokens = {
            "target",
            "target_id",
            "name",
            "id",
            "gene",
            "gene_id",
            "transcript",
            "transcript_id",
        }
        first_data_line = True
        with open(path, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                token = line.split()[0]
                if first_data_line and token.lower() in header_tokens:
                    first_data_line = False
                    continue
                targets.add(token)
                first_data_line = False
        if not targets:
            raise ValueError("No target identifiers were found in: {path}".format(path=path))
        return targets

    def _resolve_target_ids(self) -> set[str] | None:
        """Return requested gene/transcript identifiers, if any."""
        if self.target:
            return {str(self.target)}
        if self.target_list:
            return self._read_target_file(self.target_list)
        return None

    @staticmethod
    def _record_meta(record: dict[str, Any]) -> dict[str, str]:
        """Extract common gene and transcript identifiers from one record."""
        annotation = record.get("annotation") if isinstance(record.get("annotation"), dict) else {}
        name = str(record.get("name") or record.get("transcript_id") or annotation.get("transcript_id") or "NA")
        transcript_id = str(record.get("transcript_id") or annotation.get("transcript_id") or name)
        gene_id = str(record.get("gene_id") or annotation.get("gene_id") or name)
        return {"name": name, "transcript_id": transcript_id, "gene_id": gene_id}

    def _matched_target_ids(self, record: dict[str, Any]) -> set[str]:
        """Return requested identifiers matched by one JSON transcript record."""
        if not self.target_ids:
            return set()
        meta = self._record_meta(record)
        record_ids = {meta["name"], meta["transcript_id"], meta["gene_id"]}
        return self.target_ids.intersection(record_ids)

    def _resolve_groups(self, table: pd.DataFrame) -> None:
        """Resolve control and treatment samples from input metadata."""
        control_value = self._first_value(table, "control_group")
        treat_value = self._first_value(table, "treat_group")

        self.control = self._parse_csv(control_value)
        self.treat = self._parse_csv(treat_value)
        if not self.control or not self.treat:
            raise ValueError(
                "Sample groups are unavailable. Use an rpf_Odd_Ratio dev002 or "
                "newer output containing control_group and treat_group metadata."
            )
        overlap = sorted(set(self.control).intersection(self.treat))
        if overlap:
            raise ValueError("Samples occur in both groups: {samples}".format(samples=", ".join(overlap)))
        self.samples = self.control + self.treat

    def _resolve_analysis_options(self, table: pd.DataFrame) -> None:
        """Resolve analysis settings from table metadata or defaults."""
        inferred_site = self._first_value(table, "analysis_site")
        inferred_frame = self._first_value(table, "analysis_frame")
        self.site = str(inferred_site or "P").upper()
        self.frame = str(inferred_frame or "all").lower()
        if self.site not in {"E", "P", "A"}:
            raise ValueError("Unsupported analysis site in input table: {site}".format(site=self.site))
        if self.frame not in {"0", "1", "2", "all"}:
            raise ValueError("Unsupported analysis frame in input table: {frame}".format(frame=self.frame))

        inferred_window = self._first_value(table, "local_window_codons")
        inferred_score = self._first_value(table, "pause_score_threshold")
        inferred_site_rpf = self._first_value(table, "min_site_rpf_threshold")
        inferred_coverage = self._first_value(table, "min_local_coverage_threshold")
        self.local_window = int(inferred_window or 101)
        self.pause_score = float(inferred_score or 20.0)
        self.min_site_rpf = int(inferred_site_rpf or 3)
        self.min_local_coverage = float(inferred_coverage or 0.10)
        if self.local_window < 3 or self.local_window % 2 == 0:
            raise ValueError("Resolved local_window_codons must be an odd integer >= 3.")
        if self.pause_score <= 1:
            raise ValueError("Resolved pause_score_threshold must be > 1.")
        if self.min_site_rpf < 1:
            raise ValueError("Resolved min_site_rpf_threshold must be >= 1.")
        if not 0 <= self.min_local_coverage <= 1:
            raise ValueError("Resolved min_local_coverage_threshold must be in [0, 1].")

    @staticmethod
    def _prepare_events(table: pd.DataFrame) -> pd.DataFrame:
        """Preserve every input pause-site row for target-gene plotting."""
        result = table.copy()
        result.reset_index(drop=True, inplace=True)
        result["event_number"] = np.arange(1, len(result) + 1, dtype=int)
        return result

    def read_events(self) -> None:
        """Read and filter the rpf_Odd_Ratio site table."""
        table = pd.read_csv(self.input_file, sep="\t", low_memory=False)
        required = {"name", "from_tis"}
        missing = sorted(required - set(table.columns))
        if missing:
            raise ValueError("Pausing-site table is missing required columns: {columns}".format(columns=", ".join(missing)))
        if table.empty:
            raise ValueError("Pausing-site table contains no rows: {path}".format(path=self.input_file))

        self._resolve_groups(table)
        self._resolve_analysis_options(table)
        self.events = self._prepare_events(table)
        if self.events.empty:
            raise ValueError("No pausing events remain after filtering.")

        print(
            "Candidate events={events:,}, transcripts={transcripts:,}, samples={samples:,}; "
            "site={site}, frame={frame}, local_window={window}.".format(
                events=len(self.events),
                transcripts=self.events["name"].astype(str).nunique(),
                samples=len(self.samples),
                site=self.site,
                frame=self.frame,
                window=self.local_window,
            ),
            flush=True,
        )

    @staticmethod
    def _density_sum(record: dict[str, Any], sample: str) -> float:
        """Return total density encoded for one sample in one record."""
        sample_entry = record.get("samples", {}).get(sample, {})
        density = sample_entry.get("density", {}) if isinstance(sample_entry, dict) else {}
        encoding = str(density.get("encoding", RPFs.SPARSE_ENCODING))
        if encoding in {RPFs.DENSE_ENCODING, "dense"}:
            return float(sum(density.get("f0", [])) + sum(density.get("f1", [])) + sum(density.get("f2", [])))
        return float(sum(density.get("count", [])))

    @staticmethod
    def _record_keys(record: dict[str, Any]) -> set[str]:
        """Return identifiers that can map one JSON record to an event name."""
        keys = set()
        for value in (record.get("name"), record.get("transcript_id")):
            if value is not None and str(value) != "":
                keys.add(str(value))
        annotation = record.get("annotation")
        if isinstance(annotation, dict) and annotation.get("transcript_id"):
            keys.add(str(annotation["transcript_id"]))
        return keys

    def read_rpf_records(self) -> None:
        """Stream the JSONL file once and retain only event transcripts."""
        if self.events is None:
            raise ValueError("Pausing events have not been imported.")
        target_names = set(self.events["name"].astype(str))
        self.total_rpf = {sample: 0.0 for sample in self.samples}
        all_samples: list[str] | None = None
        record_count = 0
        matched_requested_ids: set[str] = set()

        for record in RPFs.iter_json_records(self.rpf_file):
            current_samples = [str(sample) for sample in record.get("samples", {}).keys()]
            if all_samples is None:
                all_samples = current_samples
                missing_samples = sorted(set(self.samples) - set(all_samples))
                if missing_samples:
                    raise ValueError("Samples not found in JSON density: {samples}".format(samples=", ".join(missing_samples)))
            elif current_samples != all_samples:
                raise ValueError("Inconsistent sample order in JSON density records.")

            if self.normal:
                for sample in self.samples:
                    self.total_rpf[sample] += self._density_sum(record, sample)

            matched_ids = self._matched_target_ids(record)
            matched_requested_ids.update(matched_ids)
            matched = target_names.intersection(self._record_keys(record))
            if self.target_ids and not matched_ids:
                matched = set()
            for target_name in matched:
                if target_name in self.records:
                    raise ValueError("Multiple JSON records match event transcript: {target}".format(target=target_name))
                self.records[target_name] = record

            record_count += 1
            if not self.normal and not self.target_ids and len(self.records) == len(target_names):
                break

        if self.target_ids:
            missing_requested = sorted(self.target_ids - matched_requested_ids)
            if missing_requested:
                raise ValueError(
                    "Requested target ID(s) were not found in JSON density: {targets}".format(
                        targets=", ".join(missing_requested)
                    )
                )
            self.events = self.events.loc[
                self.events["name"].astype(str).isin(self.records)
            ].copy()
            if self.events.empty:
                raise ValueError("No selected pausing events were found for the requested target ID(s).")
            self.events.reset_index(drop=True, inplace=True)
            self.events["event_number"] = np.arange(1, len(self.events) + 1, dtype=int)
        else:
            missing_targets = sorted(target_names - set(self.records))
            if missing_targets:
                preview = ", ".join(missing_targets[:10])
                suffix = " ..." if len(missing_targets) > 10 else ""
                raise ValueError(
                    "Event transcript(s) not found in JSON density: {targets}{suffix}".format(
                        targets=preview,
                        suffix=suffix,
                    )
                )
        if self.normal:
            zero_totals = [sample for sample, value in self.total_rpf.items() if value <= 0]
            if zero_totals:
                raise ValueError("Cannot normalize samples with zero total RPF: {samples}".format(samples=", ".join(zero_totals)))

        print(
            "Scanned JSON records={records:,}; retained target records={targets:,}; "
            "selected events={events:,}.".format(
                records=record_count,
                targets=len(self.records),
                events=len(self.events),
            ),
            flush=True,
        )

    @staticmethod
    def _shift_matrix(matrix: np.ndarray, site: str) -> np.ndarray:
        """Shift codon-frame density using the same E/P/A convention as RPFs."""
        shift = RPFs.set_codon_shift(site)
        if shift == 0:
            return matrix.copy()
        shifted = np.zeros_like(matrix)
        if shift > 0:
            shifted[shift:, :] = matrix[:-shift, :]
        else:
            shifted[:shift, :] = matrix[-shift:, :]
        return shifted

    def _analysis_vector(self, matrix: np.ndarray) -> np.ndarray:
        """Collapse the selected frame(s) to the vector used for pause scoring."""
        if self.frame == "all":
            return matrix.sum(axis=1)
        return matrix[:, int(self.frame)]

    @staticmethod
    def _event_codon_index(event: pd.Series, record: dict[str, Any]) -> int:
        """Map an Odd Ratio event row to the JSON codon index."""
        trim = record.get("trim", {}) if isinstance(record.get("trim"), dict) else {}
        codon_count = int(trim.get("codon_count", 0))
        utr5_codons = int(trim.get("utr5_codons", 0))
        from_tis = int(float(event["from_tis"]))
        candidate = from_tis + utr5_codons
        if 0 <= candidate < codon_count:
            return candidate

        if "now_nt" in event.index and pd.notna(event["now_nt"]):
            start_nt0 = int(trim.get("start_nt0", 0))
            now_nt = int(float(event["now_nt"]))
            candidate = (now_nt - (start_nt0 + 1)) // 3
            if 0 <= candidate < codon_count:
                return candidate
        raise ValueError(
            "Cannot map event {name}:{position} to JSON codon coordinates.".format(
                name=event["name"],
                position=event["from_tis"],
            )
        )

    def _sample_metrics(self, vector: np.ndarray, event_index: int) -> dict[str, Any]:
        """Calculate sample-level pause evidence over the full local window."""
        half_window = self.local_window // 2
        start = max(0, event_index - half_window)
        end = min(vector.size, event_index + half_window + 1)
        local = vector[start:end]
        local_mean = float(local.mean()) if local.size else 0.0
        local_coverage = float(np.count_nonzero(local) / local.size) if local.size else 0.0
        site_rpf = int(vector[event_index])
        if local_mean > 0:
            score = float(site_rpf / local_mean)
        elif site_rpf > 0:
            score = float("inf")
        else:
            score = 0.0
        is_pause = (
            site_rpf >= self.min_site_rpf
            and score >= self.pause_score
            and local_coverage >= self.min_local_coverage
        )
        return {
            "site_rpf": site_rpf,
            "local_mean": local_mean,
            "local_coverage": local_coverage,
            "pause_score": score,
            "sample_pause": bool(is_pause),
        }

    def _transform(self, values: np.ndarray) -> np.ndarray:
        """Apply a non-negative display-only transformation."""
        values = np.nan_to_num(values.astype(float), nan=0.0, posinf=0.0, neginf=0.0)
        values = np.clip(values, 0.0, None)
        if self.plot_transform == "none":
            return values
        if self.plot_transform == "sqrt":
            return np.sqrt(values)
        if self.plot_transform == "log1p":
            return np.log1p(values)
        if self.plot_transform == "log2":
            return np.log2(values + 1.0)
        if self.plot_transform == "log10":
            return np.log10(values + 1.0)
        raise ValueError("Unsupported plot transformation: {value}".format(value=self.plot_transform))

    def _display_matrix(self, raw_matrix: np.ndarray, sample: str) -> np.ndarray:
        """Return raw or RPM-normalized frame density for plotting."""
        values = raw_matrix.astype(float)
        if self.normal:
            values = values * 1_000_000.0 / self.total_rpf[sample]
        return self._transform(values)

    @staticmethod
    def _format_number(value: float) -> str:
        """Format a compact numeric label."""
        if not np.isfinite(value):
            return "Inf"
        if abs(value) >= 1000:
            return "{:.1f}K".format(value / 1000.0).rstrip("0").rstrip(".")
        if abs(value) >= 10:
            return "{:.1f}".format(value).rstrip("0").rstrip(".")
        return "{:.2f}".format(value).rstrip("0").rstrip(".")

    def _gene_output_stem(self, record: dict[str, Any]) -> str:
        """Return a unique output stem for one gene/transcript profile."""
        meta = self._record_meta(record)
        return "{prefix}_{gene}_{transcript}_PSplot".format(
            prefix=self.output_prefix,
            gene=self._safe_name(meta["gene_id"]),
            transcript=self._safe_name(meta["transcript_id"]),
        )

    def _metric_row(
        self,
        event: pd.Series,
        sample: str,
        group: str,
        metrics: dict[str, Any],
        record: dict[str, Any],
        site_number: int | None = None,
    ) -> dict[str, Any]:
        """Build one exported sample-by-site metric row."""
        meta = self._record_meta(record)
        metric_row = {
            "event_number": int(event["event_number"]),
            "site_number": site_number,
            "gene_id": meta["gene_id"],
            "transcript_id": meta["transcript_id"],
            "name": str(event["name"]),
            "from_tis": int(float(event["from_tis"])),
            "codon": str(event.get("codon", "NA")),
            "pause_class": str(event.get("pause_class", "NA")),
            "sample": sample,
            "group": group.lower(),
            "analysis_site": self.site,
            "analysis_frame": self.frame,
        }
        metric_row.update(metrics)
        return metric_row

    def _plot_frame_bars(self, ax, x: np.ndarray, matrix: np.ndarray) -> None:
        """Draw selected frame-resolved RPF bars."""
        if self.frame == "all":
            offsets = {"0": -0.26, "1": 0.0, "2": 0.26}
            for frame in ("0", "1", "2"):
                ax.bar(
                    x + offsets[frame],
                    matrix[:, int(frame)],
                    width=0.25,
                    color=FRAME_COLORS[frame],
                    edgecolor="none",
                    linewidth=0,
                    label="Frame " + frame,
                    zorder=2,
                )
        else:
            frame = self.frame
            ax.bar(
                x,
                matrix[:, int(frame)],
                width=0.78,
                color=FRAME_COLORS[frame],
                edgecolor="none",
                linewidth=0,
                label="Frame " + frame,
                zorder=2,
            )

    def _draw_gene_events(
        self,
        target_events: pd.DataFrame,
        record: dict[str, Any],
        matrices: dict[str, np.ndarray],
    ) -> None:
        """Draw all selected pause sites from one transcript in one figure."""
        unique_events = (
            target_events.sort_values(["from_tis", "event_number"], kind="mergesort")
            .drop_duplicates(subset=["from_tis"], keep="first")
            .reset_index(drop=True)
        )
        event_details: list[tuple[int, pd.Series, int]] = []
        for row_index, event in unique_events.iterrows():
            event_details.append((row_index + 1, event, self._event_codon_index(event, record)))

        trim = record.get("trim", {}) if isinstance(record.get("trim"), dict) else {}
        codon_count = int(trim.get("codon_count", 0))
        utr5_codons = int(trim.get("utr5_codons", 0))
        event_indices = [item[2] for item in event_details]
        if self.view == "region":
            start = max(0, min(event_indices) - self.flank)
            end = min(codon_count, max(event_indices) + self.flank + 1)
        else:
            start = 0
            end = codon_count
        indices = np.arange(start, end, dtype=int)
        x = indices - utr5_codons

        plot_matrices = {
            sample: self._display_matrix(matrices[sample][indices, :], sample)
            for sample in self.samples
        }
        metrics_by_sample: dict[str, list[dict[str, Any]]] = {}
        for sample in self.samples:
            vector = self._analysis_vector(matrices[sample])
            metrics_by_sample[sample] = [
                self._sample_metrics(vector, event_index)
                for _, _, event_index in event_details
            ]

        sample_count = len(self.samples)
        figure_height = min(max(3.2, sample_count * self.per_sample_height + 1.7), self.max_height)
        fig, axes = plt.subplots(
            nrows=sample_count,
            ncols=1,
            figsize=(self.figure_width, figure_height),
            sharex=True,
            squeeze=False,
            gridspec_kw={"hspace": 0.06},
        )
        axes_list = axes[:, 0]

        if self.y_max is not None:
            shared_top = float(self.y_max)
        else:
            maxima = [float(np.max(values)) if values.size else 0.0 for values in plot_matrices.values()]
            shared_top = max(max(maxima, default=0.0) * 1.12, 1.0)

        for sample_index, (ax, sample) in enumerate(zip(axes_list, self.samples)):
            plot_matrix = plot_matrices[sample]
            sample_metrics = metrics_by_sample[sample]
            group = "Control" if sample in self.control else "Treatment"
            group_color = CONTROL_COLOR if sample in self.control else TREAT_COLOR

            if self.pause_region:
                for _, _, event_index in event_details:
                    event_x = float(event_index - utr5_codons)
                    ax.axvspan(
                        event_x - 1.5,
                        event_x + 1.5,
                        facecolor=PAUSE_REGION_COLOR,
                        edgecolor="none",
                        linewidth=0,
                        alpha=0.3,
                        zorder=0,
                    )
            self._plot_frame_bars(ax, x.astype(float), plot_matrix)
            ax.grid(False)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_color(group_color)
            ax.spines["left"].set_linewidth(2.2)
            ax.tick_params(axis="both", labelsize=max(self.font_size - 1.0, 1.0))

            if self.y_scale == "shared" or self.y_max is not None:
                y_top = shared_top
            else:
                y_top = max(float(np.max(plot_matrix)) * 1.12 if plot_matrix.size else 0.0, 1.0)
            ax.set_ylim(0, y_top)

            pause_count = 0
            for detail_index, (site_number, event, event_index) in enumerate(event_details):
                event_x = float(event_index - utr5_codons)
                metrics = sample_metrics[detail_index]
                is_pause = bool(metrics["sample_pause"])
                pause_count += int(is_pause)
                marker_y = y_top * 0.94
                ax.scatter(
                    [event_x],
                    [marker_y],
                    marker="v",
                    s=34,
                    facecolor=PAUSE_COLOR if is_pause else "white",
                    edgecolor=PAUSE_COLOR if is_pause else NON_PAUSE_COLOR,
                    linewidth=0.8,
                    zorder=7,
                    clip_on=False,
                )
                if sample_index == 0:
                    ax.text(
                        event_x,
                        y_top * (0.78 if site_number % 2 == 0 else 0.98),
                        "P{number}\n{position:+d} {codon}".format(
                            number=site_number,
                            position=int(float(event["from_tis"])),
                            codon=str(event.get("codon", "NA")),
                        ),
                        ha="center",
                        va="top",
                        rotation=90,
                        fontsize=max(self.font_size - 2.0, 1.0),
                        color=PAUSE_COLOR,
                        zorder=8,
                    )
                self.metrics.append(
                    self._metric_row(
                        event,
                        sample,
                        group,
                        metrics,
                        record,
                        site_number=site_number,
                    )
                )

            ax.text(
                1.01,
                0.58,
                "{sample}\n{group} | pauses={pause}/{total}".format(
                    sample=sample,
                    group=group,
                    pause=pause_count,
                    total=len(event_details),
                ),
                transform=ax.transAxes,
                ha="left",
                va="center",
                fontsize=self.font_size,
                color="#222222",
                fontweight="bold" if pause_count else "normal",
            )
            if sample_index < sample_count - 1:
                ax.tick_params(axis="x", labelbottom=False)

        meta = self._record_meta(record)
        view_label = "whole transcript" if self.view == "gene" else "pause span with flanks"
        fig.suptitle(
            "{gene} / {transcript} | {count} selected pause sites\n"
            "{view}; site={site}, frame={frame}".format(
                gene=meta["gene_id"],
                transcript=meta["transcript_id"],
                count=len(event_details),
                view=view_label,
                site=self.site,
                frame=self.frame,
            ),
            fontsize=self.font_size + 2.0,
            y=0.995,
        )
        axes_list[-1].set_xlabel("Codon position relative to TIS", fontsize=self.font_size)
        density_label = "RPF density (RPM)" if self.normal else "Raw RPF count"
        if self.plot_transform != "none":
            density_label += ", {value} transformed".format(value=self.plot_transform)
        fig.text(0.035, 0.50, density_label, rotation=90, ha="center", va="center", fontsize=self.font_size)

        handles = []
        labels = []
        frames = ("0", "1", "2") if self.frame == "all" else (self.frame,)
        for frame in frames:
            handles.append(plt.Rectangle((0, 0), 1, 1, color=FRAME_COLORS[frame], linewidth=0))
            labels.append("Frame " + frame)
        handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="v",
                linestyle="none",
                markerfacecolor=PAUSE_COLOR,
                markeredgecolor=PAUSE_COLOR,
            )
        )
        labels.append("Sample-level pause")
        handles.append(
            plt.Line2D([0], [0], marker="v", linestyle="none", markerfacecolor="white", markeredgecolor=NON_PAUSE_COLOR)
        )
        labels.append("Sample-level non-pause")
        if self.pause_region:
            handles.append(
                plt.Rectangle(
                    (0, 0),
                    1,
                    1,
                    facecolor=PAUSE_REGION_COLOR,
                    edgecolor="none",
                    linewidth=0,
                    alpha=0.3,
                )
            )
            labels.append("Pause-site region")
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=min(len(labels), 6),
            frameon=False,
            fontsize=max(self.font_size - 1.0, 1.0),
            bbox_to_anchor=(0.5, 0.005),
        )
        fig.subplots_adjust(left=0.10, right=0.78, top=0.91, bottom=0.12)

        output_stem = self._gene_output_stem(record)
        output_path = Path(output_stem)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_key = self._safe_name(meta["transcript_id"])
        if self.output_format in {"pdf", "both"}:
            pdf_file = output_stem + ".pdf"
            fig.savefig(pdf_file, bbox_inches="tight")
            self.output_files[output_key + "_pdf"] = pdf_file
        if self.output_format in {"png", "both"}:
            png_file = output_stem + ".png"
            fig.savefig(png_file, dpi=self.dpi, bbox_inches="tight")
            self.output_files[output_key + "_png"] = png_file
        plt.close(fig)

    def draw_events(self) -> None:
        """Decode each target transcript once and draw all its pause sites."""
        if self.events is None:
            raise ValueError("Pausing events have not been imported.")
        for target_name, target_events in self.events.groupby("name", sort=False):
            record = self.records[str(target_name)]
            codon_count = int(record.get("trim", {}).get("codon_count", 0))
            if codon_count <= 0:
                raise ValueError("Target JSON record has no valid codon_count: {target}".format(target=target_name))
            matrices = {
                sample: self._shift_matrix(
                    RPFs.density_to_frame_matrix(record, sample, codon_count),
                    self.site,
                )
                for sample in self.samples
            }
            print(
                "Plot transcript {name}: pause sites={sites:,}.".format(
                    name=target_name,
                    sites=target_events["from_tis"].nunique(),
                ),
                flush=True,
            )
            self._draw_gene_events(target_events, record, matrices)

        transcript_count = self.events["name"].astype(str).nunique()
        print("Generated transcript figures={count:,}.".format(count=transcript_count), flush=True)

    def output_metrics(self) -> None:
        """Write per-event, per-sample metrics used for track labels."""
        if not self.metrics:
            raise ValueError("No sample-level pause metrics are available.")
        metrics_file = self.output_prefix + "_PSplot.metrics.txt"
        metrics_path = Path(metrics_file)
        metrics_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame.from_records(self.metrics).to_csv(metrics_file, sep="\t", index=False)
        self.output_files["metrics"] = metrics_file
        print("Sample-level pause metrics written: {path}".format(path=metrics_file), flush=True)
