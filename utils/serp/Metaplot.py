#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-10
# Version: dev002
# Function: Calculate sample-first SeRP enrichment metaplots around TIS and TTS.
# Input: TXT/JSON RPF density data and matched control/IP sample groups.
# Output: Sample metaprofiles, paired enrichment profiles, summary tables, and figures.

"""Sample-first metaplot workflow for selective ribosome profiling.

The workflow first builds one normalized TIS/TTS metagene profile for every
sample from a common set of transcripts. Matched IP/control enrichment is then
calculated from those sample-level metaprofiles, followed by aggregation across
biological replicate pairs.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from utils.serp.Enrichment import calculate_matched_local_enrichment


BASE_COLUMNS = ["name", "now_nt", "from_tis", "from_tts", "region", "codon"]
ENRICH_COLOR = "#0072B2"
REPLICATE_COLOR = "#8A8A8A"
NEUTRAL_COLOR = "#7A7A7A"


class SeRPMetaplot:
    """Calculate sample-first SeRP enrichment around TIS and TTS.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from ``serp_metaplot``.

    Notes
    -----
    All samples use the same retained transcript set. For each sample, RPF
    density is first normalized by the sample library size and then averaged
    across retained transcripts at each relative codon position. Enrichment is
    subsequently calculated for matched IP/control sample metaprofiles. This
    differs from a transcript-first ratio workflow in which IP/control ratios
    are calculated independently for every transcript before aggregation.
    """

    def __init__(self, args) -> None:
        """Initialize metaplot parameters and result containers."""
        self.rpf_file = str(args.rpf)
        self.output = str(args.output)
        self.total_rpf_file = args.norm
        self.control_samples = self._parse_samples(args.control, "--ck")
        self.ip_samples = self._parse_samples(args.ip, "--ip")
        self.sample_names = self.control_samples + self.ip_samples
        self.label = (
            str(args.label)
            if args.label
            else "{0} / {1}".format(
                ",".join(self.ip_samples),
                ",".join(self.control_samples),
            )
        )

        self.min_rpf = int(args.min)
        self.scale = float(args.scale)
        self.window = int(args.window)
        self.pseudocount = float(args.pseudocount)
        self.tis_length = int(args.tis)
        self.tts_length = int(args.tts)
        self.aggregate = str(args.aggregate)
        self.detail = bool(args.detail)
        self.show_replicates = bool(args.show_replicates)

        self.output_format = str(args.output_format)
        self.dpi = int(args.dpi)
        self.font_size = float(args.font_size)

        self.rpf_data = None
        self.all_rpf = None
        self.total_rpf_num = None
        self.file_format = None
        self.input_gene_count = 0
        self.coverage_gene_count = 0
        self.retained_genes: list[str] = []

        self.sample_meta_table = None
        self.pair_enrichment_table = None
        self.metaplot_table = None
        self.pair_names: list[str] = []

    @staticmethod
    def _parse_samples(value: str, argument_name: str) -> list[str]:
        """Parse one comma-separated sample argument and reject duplicates."""
        samples = [item.strip() for item in str(value).split(",") if item.strip()]
        if not samples:
            raise ValueError("{0} requires at least one sample.".format(argument_name))
        if len(samples) != len(set(samples)):
            raise ValueError(
                "{0} contains duplicated sample names: {1}".format(
                    argument_name,
                    ", ".join(samples),
                )
            )
        return samples

    def import_rpf(self) -> None:
        """Import TXT/JSON RPF density data through the common RPF reader."""
        from utils.ribo.RPFs import RPFData

        if len(self.control_samples) != len(self.ip_samples):
            raise ValueError(
                "serp_metaplot requires equal numbers of control and IP samples so "
                "biological replicates can be paired by sample order."
            )

        rpf_data = RPFData.from_file(
            rpf_file=self.rpf_file,
            sample_name=self.sample_names,
        )
        merged_rpf = rpf_data.get_frame(frame="all")

        missing_samples = [
            sample for sample in self.sample_names if sample not in merged_rpf.columns
        ]
        if missing_samples:
            raise ValueError(
                "RPF density data are missing requested sample(s): {0}".format(
                    ", ".join(missing_samples)
                )
            )

        self.rpf_data = rpf_data
        self.file_format = str(rpf_data.file_format)
        self.all_rpf = merged_rpf.loc[:, BASE_COLUMNS + self.sample_names].copy()
        self.input_gene_count = int(self.all_rpf["name"].nunique())

        if self.total_rpf_file:
            total_table = pd.read_csv(
                self.total_rpf_file,
                sep="\t",
                header=None,
                usecols=[0, 1],
                names=["sample", "total_rpf"],
            )
            total_series = pd.Series(
                pd.to_numeric(total_table["total_rpf"], errors="raise").values,
                index=total_table["sample"].astype(str),
                dtype="float64",
            )
            missing_totals = [
                sample for sample in self.sample_names if sample not in total_series.index
            ]
            if missing_totals:
                raise ValueError(
                    "Normalization file is missing sample(s): {0}".format(
                        ", ".join(missing_totals)
                    )
                )
            self.total_rpf_num = total_series.reindex(self.sample_names)
        else:
            self.total_rpf_num = rpf_data.total_rpf_num.reindex(self.sample_names)

        if self.total_rpf_num.isna().any():
            missing_totals = self.total_rpf_num[self.total_rpf_num.isna()].index.tolist()
            raise ValueError(
                "Cannot determine total RPF counts for sample(s): {0}".format(
                    ", ".join(missing_totals)
                )
            )
        if (self.total_rpf_num <= 0).any():
            invalid_samples = self.total_rpf_num[
                self.total_rpf_num <= 0
            ].index.tolist()
            raise ValueError(
                "Total RPF counts must be positive for sample(s): {0}".format(
                    ", ".join(invalid_samples)
                )
            )

        print("Input density format: {0}".format(self.file_format), flush=True)
        print("Transcripts: {0:,}".format(self.input_gene_count), flush=True)
        print(
            "Replicate pairs: {0}".format(
                ", ".join(
                    "{0}/{1}".format(ip, ck)
                    for ck, ip in zip(self.control_samples, self.ip_samples)
                )
            ),
            flush=True,
        )

    def select_genes(self) -> None:
        """Select one common transcript set for all sample metaprofiles."""
        cds = self.all_rpf.loc[self.all_rpf["region"].eq("cds")].copy()
        if cds.empty:
            raise ValueError("No CDS records are available in the RPF density file.")

        gene_counts = cds.groupby("name", sort=False)[self.sample_names].sum()
        count_pass = gene_counts.index[(gene_counts >= self.min_rpf).all(axis=1)]
        self.coverage_gene_count = int(len(count_pass))

        candidate = cds.loc[cds["name"].isin(count_pass)]
        tis_positions = candidate.loc[
            candidate["from_tis"].between(0, self.tis_length - 1),
            ["name", "from_tis"],
        ].drop_duplicates()
        tts_positions = candidate.loc[
            candidate["from_tts"].between(-self.tts_length + 1, 0),
            ["name", "from_tts"],
        ].drop_duplicates()

        tis_complete = tis_positions.groupby("name", sort=False)["from_tis"].size()
        tts_complete = tts_positions.groupby("name", sort=False)["from_tts"].size()
        tis_complete = tis_complete.index[tis_complete.eq(self.tis_length)]
        tts_complete = tts_complete.index[tts_complete.eq(self.tts_length)]

        retained = count_pass.intersection(tis_complete).intersection(tts_complete)
        self.retained_genes = [str(item) for item in retained.tolist()]
        if not self.retained_genes:
            raise ValueError(
                "No transcript passed both the sample-level RPF filter and the complete "
                "TIS/TTS window requirements. Reduce --min, --tis, or --tts."
            )

        print(
            "Transcripts passing all-sample RPF filter: {0:,}".format(
                self.coverage_gene_count
            ),
            flush=True,
        )
        print(
            "Transcripts with complete metaplot windows: {0:,}".format(
                len(self.retained_genes)
            ),
            flush=True,
        )

    def _normalized_window_data(
        self,
        meta_name: str,
        coord_column: str,
        positions: np.ndarray,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return per-transcript normalized densities and sample metaprofiles."""
        retained_set = set(self.retained_genes)
        window = self.all_rpf.loc[
            self.all_rpf["name"].isin(retained_set)
            & self.all_rpf["region"].eq("cds")
            & self.all_rpf[coord_column].isin(positions),
            ["name", coord_column] + self.sample_names,
        ].copy()

        if window.empty:
            raise ValueError("No {0} records are available after filtering.".format(meta_name))
        duplicated = window.duplicated(["name", coord_column], keep=False)
        if bool(duplicated.any()):
            examples = (
                window.loc[duplicated, ["name", coord_column]]
                .drop_duplicates()
                .head(10)
            )
            raise ValueError(
                "Duplicated transcript/position rows were detected in {0}: {1}".format(
                    meta_name,
                    examples.to_dict(orient="records"),
                )
            )

        normalized = window.loc[:, self.sample_names].astype(float).div(
            self.total_rpf_num.loc[self.sample_names].values,
            axis=1,
        ) * self.scale
        normalized.insert(0, "position", window[coord_column].astype(int).to_numpy())
        normalized.insert(0, "transcript", window["name"].astype(str).to_numpy())

        expected_rows = len(self.retained_genes) * len(positions)
        if len(normalized) != expected_rows:
            raise ValueError(
                "{0} window is incomplete after transcript filtering: expected {1:,} rows, "
                "observed {2:,}.".format(meta_name, expected_rows, len(normalized))
            )

        sample_meta = normalized.groupby("position", sort=False)[self.sample_names].mean()
        sample_meta = sample_meta.reindex(positions)
        if sample_meta.isna().any().any():
            raise ValueError(
                "Missing sample metaprofile values were detected in {0}.".format(meta_name)
            )
        sample_meta.index = sample_meta.index.astype(int)
        return normalized, sample_meta

    @staticmethod
    def _summarize_pair_enrichment(
        pair_ratio: pd.DataFrame,
        positions: np.ndarray,
        meta_name: str,
        aggregate: str,
    ) -> pd.DataFrame:
        """Summarize enrichment across matched biological replicate pairs."""
        matrix = pair_ratio.reindex(positions).to_numpy(dtype=float)
        mean_values = np.nanmean(matrix, axis=1)
        median_values = np.nanmedian(matrix, axis=1)
        min_values = np.nanmin(matrix, axis=1)
        max_values = np.nanmax(matrix, axis=1)
        if matrix.shape[1] >= 3:
            q25_values = np.nanpercentile(matrix, 25, axis=1)
            q75_values = np.nanpercentile(matrix, 75, axis=1)
        else:
            q25_values = np.full(len(positions), np.nan, dtype=float)
            q75_values = np.full(len(positions), np.nan, dtype=float)
        selected = median_values if aggregate == "median" else mean_values

        return pd.DataFrame(
            {
                "meta": meta_name,
                "position": positions.astype(int),
                "gene_number": int(0),
                "pair_number": int(pair_ratio.shape[1]),
                "enrichment": selected,
                "mean_enrichment": mean_values,
                "median_enrichment": median_values,
                "min_enrichment": min_values,
                "max_enrichment": max_values,
                "q25_enrichment": q25_values,
                "q75_enrichment": q75_values,
            }
        )

    def _write_detail_rows(
        self,
        handle,
        meta_name: str,
        normalized: pd.DataFrame,
    ) -> None:
        """Write optional per-transcript normalized sample densities."""
        if handle is None:
            return
        sample_index = {sample: normalized.columns.get_loc(sample) for sample in self.sample_names}
        transcript_index = normalized.columns.get_loc("transcript")
        position_index = normalized.columns.get_loc("position")
        for row in normalized.itertuples(index=False, name=None):
            transcript = str(row[transcript_index])
            position = int(row[position_index])
            for sample in self.sample_names:
                value = float(row[sample_index[sample]])
                handle.write(
                    "{0}\t{1}\t{2}\t{3}\t{4:.10g}\n".format(
                        meta_name,
                        transcript,
                        position,
                        sample,
                        value,
                    )
                )

    def calculate_metaplot(self) -> None:
        """Build sample metaprofiles first, then calculate paired enrichment."""
        tis_positions = np.arange(0, self.tis_length, dtype=int)
        tts_positions = np.arange(-self.tts_length + 1, 1, dtype=int)

        sample_records = []
        pair_records = []
        summary_tables = []

        detail_handle = None
        if self.detail:
            detail_path = Path(self.output + ".metaplot_gene.txt.gz")
            detail_handle = gzip.open(detail_path, "wt", encoding="utf-8", newline="")
            detail_handle.write(
                "meta\ttranscript\tposition\tsample\tnormalized_density\n"
            )

        try:
            for meta_name, coord_column, positions in (
                ("TIS", "from_tis", tis_positions),
                ("TTS", "from_tts", tts_positions),
            ):
                normalized, sample_meta = self._normalized_window_data(
                    meta_name=meta_name,
                    coord_column=coord_column,
                    positions=positions,
                )
                self._write_detail_rows(detail_handle, meta_name, normalized)

                for position, values in sample_meta.iterrows():
                    for sample in self.sample_names:
                        sample_records.append(
                            {
                                "meta": meta_name,
                                "position": int(position),
                                "sample": sample,
                                "gene_number": len(self.retained_genes),
                                "meta_density": float(values[sample]),
                            }
                        )

                pair_ratio, _ = calculate_matched_local_enrichment(
                    raw_gene_rpm=sample_meta,
                    control_samples=self.control_samples,
                    ip_samples=self.ip_samples,
                    window=self.window,
                    pseudocount=self.pseudocount,
                )
                self.pair_names = [str(column) for column in pair_ratio.columns]

                for pair_name, control_sample, ip_sample in zip(
                    pair_ratio.columns,
                    self.control_samples,
                    self.ip_samples,
                ):
                    pair_series = pair_ratio[pair_name].reindex(positions)
                    for position, value in pair_series.items():
                        pair_records.append(
                            {
                                "meta": meta_name,
                                "position": int(position),
                                "pair": str(pair_name),
                                "control_sample": control_sample,
                                "ip_sample": ip_sample,
                                "enrichment": float(value),
                            }
                        )

                summary = self._summarize_pair_enrichment(
                    pair_ratio=pair_ratio,
                    positions=positions,
                    meta_name=meta_name,
                    aggregate=self.aggregate,
                )
                summary["gene_number"] = len(self.retained_genes)
                summary_tables.append(summary)
        finally:
            if detail_handle is not None:
                detail_handle.close()

        self.sample_meta_table = pd.DataFrame.from_records(sample_records)
        self.pair_enrichment_table = pd.DataFrame.from_records(pair_records)
        self.metaplot_table = pd.concat(summary_tables, axis=0, ignore_index=True)

    def write_tables(self) -> None:
        """Write sample, paired-enrichment, aggregate, and run-summary tables."""
        self.sample_meta_table.to_csv(
            self.output + ".metaplot_sample.txt",
            sep="\t",
            index=False,
            float_format="%.8g",
        )
        self.pair_enrichment_table.to_csv(
            self.output + ".metaplot_pair.txt",
            sep="\t",
            index=False,
            float_format="%.8g",
        )
        self.metaplot_table.to_csv(
            self.output + ".metaplot.txt",
            sep="\t",
            index=False,
            float_format="%.8g",
        )

        summary = pd.DataFrame(
            [
                ("input_rpf", self.rpf_file),
                ("input_format", self.file_format),
                ("control_samples", ",".join(self.control_samples)),
                ("ip_samples", ",".join(self.ip_samples)),
                ("replicate_pairs", len(self.control_samples)),
                ("label", self.label),
                ("input_transcripts", self.input_gene_count),
                ("rpf_filter_transcripts", self.coverage_gene_count),
                ("retained_transcripts", len(self.retained_genes)),
                ("minimum_gene_rpf", self.min_rpf),
                ("normalization_scale", self.scale),
                ("gene_aggregation", "mean_normalized_density"),
                ("enrichment_order", "sample_metaplot_then_matched_ratio"),
                ("local_enrichment_window_codons", self.window),
                ("pseudocount", self.pseudocount),
                ("tis_window_codons", self.tis_length),
                ("tts_window_codons", self.tts_length),
                ("replicate_pair_aggregation", self.aggregate),
            ],
            columns=["metric", "value"],
        )
        summary.to_csv(
            self.output + ".metaplot_summary.txt",
            sep="\t",
            index=False,
        )

    def _save_figure(self, fig) -> None:
        """Save the combined TIS/TTS figure in requested formats."""
        if self.output_format in {"png", "both"}:
            fig.savefig(
                self.output + ".metaplot.png",
                dpi=self.dpi,
                bbox_inches="tight",
            )
        if self.output_format in {"pdf", "both"}:
            fig.savefig(
                self.output + ".metaplot.pdf",
                bbox_inches="tight",
            )
        plt.close(fig)

    def draw_metaplot(self) -> None:
        """Draw paired-replicate and aggregate TIS/TTS enrichment profiles."""
        if self.metaplot_table is None or self.metaplot_table.empty:
            raise RuntimeError("Metaplot table has not been calculated.")

        fig, axes = plt.subplots(
            1,
            2,
            figsize=(7.0, 3.1),
            sharey=True,
            dpi=self.dpi,
            constrained_layout=True,
        )

        for ax, meta_name in zip(axes, ["TIS", "TTS"]):
            summary = self.metaplot_table.loc[
                self.metaplot_table["meta"].eq(meta_name)
            ]
            pairs = self.pair_enrichment_table.loc[
                self.pair_enrichment_table["meta"].eq(meta_name)
            ]

            if self.show_replicates:
                for pair_name, pair_df in pairs.groupby("pair", sort=False):
                    ax.plot(
                        pair_df["position"].to_numpy(dtype=float),
                        pair_df["enrichment"].to_numpy(dtype=float),
                        color=REPLICATE_COLOR,
                        linewidth=0.8,
                        alpha=0.45,
                        label="Replicate pair" if pair_name == self.pair_names[0] else None,
                    )

            ax.plot(
                summary["position"].to_numpy(dtype=float),
                summary["enrichment"].to_numpy(dtype=float),
                color=ENRICH_COLOR,
                linewidth=1.6,
                label="{0} enrichment".format(self.aggregate.title()),
            )
            ax.axhline(1.0, color=NEUTRAL_COLOR, linewidth=0.7, linestyle="--")
            ax.axvline(0.0, color=NEUTRAL_COLOR, linewidth=0.6)
            ax.set_title(
                "{0} (n={1:,})".format(meta_name, len(self.retained_genes)),
                fontsize=self.font_size,
            )
            ax.set_xlabel("Position relative to {0} (codons)".format(meta_name))
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=max(self.font_size - 1.0, 1.0), length=3.0)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=6))

        axes[0].set_ylabel("IP / control enrichment")
        axes[0].yaxis.label.set_size(self.font_size)
        for ax in axes:
            ax.xaxis.label.set_size(self.font_size)

        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            unique = []
            seen = set()
            for handle, label in zip(handles, labels):
                if not label or label in seen:
                    continue
                unique.append((handle, label))
                seen.add(label)
            fig.legend(
                [item[0] for item in unique],
                [item[1] for item in unique],
                loc="upper center",
                bbox_to_anchor=(0.5, 1.04),
                ncol=max(1, len(unique)),
                frameon=False,
                fontsize=max(self.font_size - 1.0, 1.0),
            )
        fig.suptitle(self.label, fontsize=self.font_size + 0.5, y=1.09)
        self._save_figure(fig)

    def result_counts(self) -> dict[str, int | float | str]:
        """Return compact workflow result counts for terminal reporting."""
        return {
            "input_transcripts": self.input_gene_count,
            "rpf_filter_transcripts": self.coverage_gene_count,
            "retained_transcripts": len(self.retained_genes),
            "replicate_pairs": len(self.control_samples),
            "enrichment_order": "sample_metaplot_then_matched_ratio",
            "replicate_pair_aggregation": self.aggregate,
        }
