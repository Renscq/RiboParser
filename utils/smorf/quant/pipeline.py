#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-08-01
# Version: dev003
# Function: Orchestrate reliable smORF density quantification.
# Input: Validated CLI arguments.
# Output: Density matrix and quantification summary.

"""Orchestrate reliable smORF density quantification."""

from __future__ import annotations

import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from utils.ribo.ArgsParser import progress_print

from .counter import (
    _initialize_sample_worker,
    _quantify_sample,
    _records_by_chromosome,
    _sample_worker,
)
from .input import (
    read_density_design,
    read_smorf_genepred,
)
from .models import (
    VALID_FRAMES,
    QuantConfig,
    QuantResult,
)
from .output import (
    _write_matrix,
    output_path_from_prefix,
)


class SmORFQuantifier:
    """Quantify reliable smORFs across genomic P-site density samples."""

    def __init__(self, config: QuantConfig) -> None:
        """Initialize the quantifier and validate scalar configuration."""
        frame = str(config.frame).lower()
        if frame not in VALID_FRAMES:
            raise ValueError(f"Unsupported frame: {config.frame}. Use all, 0, 1, or 2.")
        if int(config.threads) < 1:
            raise ValueError("threads must be >= 1.")
        self.config = QuantConfig(
            genepred=str(config.genepred),
            density_list=str(config.density_list),
            output_prefix=str(config.output_prefix),
            frame=frame,
            include_stop=bool(config.include_stop),
            threads=int(config.threads),
        )

    def run(self) -> QuantResult:
        """Run annotation parsing, sample quantification, and atomic output."""
        records = read_smorf_genepred(
            self.config.genepred,
            include_stop=self.config.include_stop,
        )
        samples = read_density_design(self.config.density_list)
        records_by_chrom = _records_by_chromosome(records)
        sample_order = tuple(sample.sample for sample in samples)
        workers = max(1, min(self.config.threads, len(samples)))
        counts_by_sample: dict[str, np.ndarray] = {}

        with tempfile.TemporaryDirectory(prefix="smorf_quant.") as work_root:
            if workers == 1:
                for number, sample in enumerate(samples, start=1):
                    result = _quantify_sample(
                        sample,
                        records_by_chrom,
                        len(records),
                        self.config.frame,
                        str(Path(work_root) / f"sample_{number}"),
                    )
                    if result.observed_target_chromosomes == 0:
                        raise ValueError(
                            f"No density chromosome matched the smORF genePred "
                            f"for sample {result.sample}."
                        )
                    counts_by_sample[result.sample] = result.counts
                    progress_print(
                        f"quantified sample {number:,}/{len(samples):,}: {result.sample}"
                    )
            else:
                with ProcessPoolExecutor(
                    max_workers=workers,
                    initializer=_initialize_sample_worker,
                    initargs=(records_by_chrom, len(records)),
                ) as executor:
                    futures = {
                        executor.submit(
                            _sample_worker,
                            sample,
                            self.config.frame,
                            str(Path(work_root) / f"sample_{number}"),
                        ): sample.sample
                        for number, sample in enumerate(samples, start=1)
                    }
                    completed = 0
                    for future in as_completed(futures):
                        result = future.result()
                        if result.observed_target_chromosomes == 0:
                            raise ValueError(
                                f"No density chromosome matched the smORF "
                                f"genePred for sample {result.sample}."
                            )
                        counts_by_sample[result.sample] = result.counts
                        completed += 1
                        progress_print(
                            f"quantified sample {completed:,}/{len(samples):,}: {result.sample}"
                        )

        output_path = output_path_from_prefix(self.config.output_prefix)
        _write_matrix(
            output_path,
            records,
            sample_order,
            counts_by_sample,
        )
        return QuantResult(
            output=str(output_path),
            orf_count=len(records),
            sample_count=len(samples),
            effective_workers=workers,
            frame=self.config.frame,
            include_stop=self.config.include_stop,
        )


def run_smorf_quant(args: object) -> QuantResult:
    """Build a quantification config from argparse arguments and run it."""
    quantifier = SmORFQuantifier(
        QuantConfig(
            genepred=str(args.genepred),
            density_list=str(args.density_list),
            output_prefix=str(args.output),
            frame=str(getattr(args, "frame", "all")),
            include_stop=bool(getattr(args, "include_stop", False)),
            threads=int(getattr(args, "threads", 1)),
        )
    )
    return quantifier.run()
