#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Rensc
# Date: 2026-07-22
# Version: 0.2.8.18
# Function: Define data models and evidence thresholds for smORF Ribo-seq analysis.
# Input: Density-track metadata and evidence-mode settings.
# Output: Validated threshold and track objects.

"""Data models and presets for smORF Ribo-seq evidence analysis."""

from __future__ import annotations

from dataclasses import dataclass, replace


SUPPORTED_DENSITY_FORMATS = frozenset({"auto", "wig", "bedgraph"})
VALID_STRANDS = frozenset({"+", "-", "."})
STOP_CODONS = frozenset({"TAA", "TAG", "TGA", "UAA", "UAG", "UGA"})
EVIDENCE_LEVELS = {
    "NoEvidence": 0,
    "LowConfidence": 1,
    "MediumConfidence": 2,
    "HighConfidence": 3,
}


@dataclass(frozen=True, slots=True)
class DensityTrack:
    """Describe one P-site density track.

    Attributes:
        sample: Biological sample identifier.
        strand: ``+``, ``-``, or ``.`` for unstranded density.
        path: WIG/bedGraph path.
        file_format: ``auto``, ``wig``, or ``bedgraph``.
    """

    sample: str
    strand: str
    path: str
    file_format: str = "auto"


@dataclass(frozen=True, slots=True)
class EvidenceThresholds:
    """Store sample-level evidence thresholds.

    Attributes:
        min_rpf_sum: Minimum coding-region P-site sum.
        min_covered_codon: Minimum covered coding codons.
        min_codon_coverage: Minimum covered-codon fraction.
        strong_periodicity: Frame-0 ratio for strong periodicity.
        moderate_periodicity: Frame-0 ratio for moderate periodicity.
        strong_start_pause: Strong start-pause ratio.
        moderate_start_pause: Moderate start-pause ratio.
        strong_stop_pause: Strong stop-pause ratio.
        moderate_stop_pause: Moderate stop-pause ratio.
        strong_release: Strong pre-stop/post-stop ratio.
        moderate_release: Moderate pre-stop/post-stop ratio.
        uniform_coverage_ratio: Minimum codon coverage for uniform shape.
        uniform_gini: Maximum Gini index for uniform shape.
        uniform_max_to_mean: Maximum peak-to-mean ratio for uniform shape.
        skewed_max_to_mean: Peak-to-mean ratio defining skewed coverage.
        skewed_top_fraction: Top-10% density fraction defining skew.
        min_shape_codons: Minimum ORF codons for shape classification.
        min_boundary_codons: Minimum ORF codons for pause classification.
        pseudocount: Ratio pseudocount.
    """

    min_rpf_sum: float
    min_covered_codon: int
    min_codon_coverage: float
    strong_periodicity: float
    moderate_periodicity: float
    strong_start_pause: float = 1.50
    moderate_start_pause: float = 1.20
    strong_stop_pause: float = 1.50
    moderate_stop_pause: float = 1.20
    strong_release: float = 3.00
    moderate_release: float = 1.50
    uniform_coverage_ratio: float = 0.40
    uniform_gini: float = 0.50
    uniform_max_to_mean: float = 5.00
    skewed_max_to_mean: float = 10.00
    skewed_top_fraction: float = 0.70
    min_shape_codons: int = 10
    min_boundary_codons: int = 9
    pseudocount: float = 0.10

    def validate(self) -> None:
        """Validate threshold ordering and ranges.

        Raises:
            ValueError: If a threshold is invalid.
        """
        if self.min_rpf_sum < 0:
            raise ValueError("min_rpf_sum must be >= 0.")
        if self.min_covered_codon < 1:
            raise ValueError("min_covered_codon must be >= 1.")
        if not 0 <= self.min_codon_coverage <= 1:
            raise ValueError("min_codon_coverage must be in [0, 1].")
        if not (
            0 <= self.moderate_periodicity
            <= self.strong_periodicity
            <= 1
        ):
            raise ValueError(
                "Require 0 <= moderate_periodicity <= "
                "strong_periodicity <= 1."
            )
        if self.moderate_start_pause > self.strong_start_pause:
            raise ValueError(
                "moderate_start_pause must be <= strong_start_pause."
            )
        if self.moderate_stop_pause > self.strong_stop_pause:
            raise ValueError(
                "moderate_stop_pause must be <= strong_stop_pause."
            )
        if self.moderate_release > self.strong_release:
            raise ValueError(
                "moderate_release must be <= strong_release."
            )
        if self.min_shape_codons < 3:
            raise ValueError("min_shape_codons must be >= 3.")
        if self.min_boundary_codons < 9:
            raise ValueError("min_boundary_codons must be >= 9.")
        if self.pseudocount < 0:
            raise ValueError("pseudocount must be >= 0.")

    @classmethod
    def from_mode(cls, mode: str) -> "EvidenceThresholds":
        """Build one evidence-threshold preset.

        Args:
            mode: ``sensitive``, ``balanced``, or ``strict``.

        Returns:
            Validated threshold object.

        Raises:
            ValueError: If the mode is unknown.
        """
        presets = {
            "sensitive": cls(
                min_rpf_sum=3.0,
                min_covered_codon=2,
                min_codon_coverage=0.10,
                moderate_periodicity=0.50,
                strong_periodicity=0.65,
            ),
            "balanced": cls(
                min_rpf_sum=5.0,
                min_covered_codon=3,
                min_codon_coverage=0.15,
                moderate_periodicity=0.55,
                strong_periodicity=0.70,
            ),
            "strict": cls(
                min_rpf_sum=10.0,
                min_covered_codon=5,
                min_codon_coverage=0.25,
                moderate_periodicity=0.60,
                strong_periodicity=0.75,
            ),
        }
        key = str(mode).lower()
        if key not in presets:
            raise ValueError(
                f"Unknown evidence mode: {mode}. "
                f"Available: {', '.join(sorted(presets))}"
            )
        thresholds = presets[key]
        thresholds.validate()
        return thresholds

    def with_overrides(
        self,
        *,
        min_rpf_sum: float | None = None,
        min_covered_codon: int | None = None,
        min_codon_coverage: float | None = None,
        moderate_periodicity: float | None = None,
        strong_periodicity: float | None = None,
    ) -> "EvidenceThresholds":
        """Return a copy with selected user overrides.

        Args:
            min_rpf_sum: Optional abundance override.
            min_covered_codon: Optional covered-codon override.
            min_codon_coverage: Optional codon-coverage override.
            moderate_periodicity: Optional moderate periodicity override.
            strong_periodicity: Optional strong periodicity override.

        Returns:
            Validated threshold copy.
        """
        updated = replace(
            self,
            min_rpf_sum=(
                self.min_rpf_sum
                if min_rpf_sum is None
                else float(min_rpf_sum)
            ),
            min_covered_codon=(
                self.min_covered_codon
                if min_covered_codon is None
                else int(min_covered_codon)
            ),
            min_codon_coverage=(
                self.min_codon_coverage
                if min_codon_coverage is None
                else float(min_codon_coverage)
            ),
            moderate_periodicity=(
                self.moderate_periodicity
                if moderate_periodicity is None
                else float(moderate_periodicity)
            ),
            strong_periodicity=(
                self.strong_periodicity
                if strong_periodicity is None
                else float(strong_periodicity)
            ),
        )
        updated.validate()
        return updated
