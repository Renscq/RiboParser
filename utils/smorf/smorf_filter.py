# Author: Rensc
# date: 2026-05-22

"""
Basic ORF filtering module.

This module filters ORFs based on rule-based criteria from ORF.message.txt.

Supported filtering evidence:
1. ORF source strand.
2. ORF priority.
3. ORF completeness.
4. ORF category.
5. Start codon type.
6. ORF peptide length.
7. Kozak PWM similarity score.
"""
from typing import Dict, List, Tuple, Optional
from .smorf_kozak import KozakPWM


class ORFFilter:
    """
    Rule-based ORF filter with Kozak PWM scoring.
    """

    def __init__(
        self,
        keep_start_codons: str = "ATG,CTG,GTG,TTG",
        min_aa: int = 8,
        max_aa: int = 10000,
        min_kozak_pwm_score: float = 0.0,
        kozak_pwm: Optional[KozakPWM] = None,
        keep_categories: str = "uORF,dORF,lncORF,iORF,emORF,overlap_uORF,overlap_dORF,other_ORF,annotated_ORF,annotated_mORF",
        remove_categories: str = "same_frame_iORF,antisense_ORF",
        require_sense: bool = True,
        require_primary: bool = True,
        require_complete: bool = True,
    ):
        """
        Initialize ORF filter.
        """

        self.keep_start_codons = self._parse_set(keep_start_codons)
        self.min_aa = min_aa
        self.max_aa = max_aa
        self.min_kozak_pwm_score = min_kozak_pwm_score
        self.kozak_pwm = kozak_pwm
        self.keep_categories = self._parse_set(keep_categories)
        self.remove_categories = self._parse_set(remove_categories)
        self.require_sense = require_sense
        self.require_primary = require_primary
        self.require_complete = require_complete

    @staticmethod
    def _parse_set(value: str) -> set:
        """
        Convert comma-separated string to set.
        """

        if value is None or value == "":
            return set()

        return {x.strip() for x in value.split(",") if x.strip()}

    @staticmethod
    def read_table(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
        """
        Read tab-delimited message table.
        """

        records = []

        with open(path, "r") as handle:
            header = handle.readline().rstrip("\n").split("\t")

            for line in handle:
                line = line.rstrip("\n")

                if not line:
                    continue

                fields = line.split("\t")
                record = dict(zip(header, fields))
                records.append(record)

        return header, records

    @staticmethod
    def write_table(
        path: str,
        header: List[str],
        records: List[Dict[str, str]],
    ) -> None:
        """
        Write tab-delimited message table.
        """

        with open(path, "w") as out:
            out.write("\t".join(header) + "\n")

            for record in records:
                out.write("\t".join(record.get(col, "") for col in header) + "\n")

    def add_kozak_pwm_fields(
        self,
        header: List[str],
        records: List[Dict[str, str]],
    ) -> List[str]:
        """
        Add Kozak PWM score fields.
        """

        new_header = list(header)

        for col in [
            "kozak_pwm_name",
            "kozak_pwm_source",
            "kozak_pwm_consensus",
            "kozak_pwm_score",
            "kozak_pwm_level",
            "kozak_consensus_identity",
            "filter_status",
            "filter_reason",
        ]:
            if col not in new_header:
                new_header.append(col)

        for record in records:
            if self.kozak_pwm is None:
                record["kozak_pwm_name"] = "none"
                record["kozak_pwm_source"] = "none"
                record["kozak_pwm_consensus"] = "NA"
                record["kozak_pwm_score"] = "NA"
                record["kozak_pwm_level"] = "NA"
                record["kozak_consensus_identity"] = "NA"
                continue

            kozak_seq = record.get("kozak_seq", "")
            score = self.kozak_pwm.score(kozak_seq)
            identity = self.kozak_pwm.identity_to_consensus(kozak_seq)

            record["kozak_pwm_name"] = self.kozak_pwm.name
            record["kozak_pwm_source"] = self.kozak_pwm.source
            record["kozak_pwm_consensus"] = self.kozak_pwm.consensus
            record["kozak_pwm_score"] = "{:.6f}".format(score)
            record["kozak_pwm_level"] = self.kozak_pwm.level(score)
            record["kozak_consensus_identity"] = "{:.6f}".format(identity)

        return new_header

    def filter_records(
        self,
        records: List[Dict[str, str]],
        report_every: int = 1000,
    ) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
        """
        Filter ORF records.
        """

        passed = []
        removed = []
        total_orfs = len(records)

        for idx, record in enumerate(records, start=1):
            # Print filtering progress every N ORFs.
            if report_every > 0 and idx % report_every == 0:
                print(
                    "[smORFFilter] [{}/{}] Processing ORF: {}".format(
                        idx,
                        total_orfs,
                        record.get("orf_id", "NA"),
                    ),
                    flush=True,
                )

            keep, reason = self.check_record(record)
            record["filter_status"] = "PASS" if keep else "FAIL"
            record["filter_reason"] = reason

            if keep:
                passed.append(record)
            else:
                removed.append(record)

        return passed, removed

    def check_record(self, record: Dict[str, str]) -> Tuple[bool, str]:
        """
        Check one ORF record.
        """

        reasons = []

        source_strand = record.get("source_strand", "")
        priority = record.get("priority", "")
        completeness = record.get("completeness", "")
        category = record.get("category", "")
        start_codon = record.get("start_codon", "")
        aa_length = self._safe_int(record.get("aa_length", "0"))

        if self.require_sense and source_strand != "sense":
            reasons.append("non_sense_strand")

        if self.require_primary and priority != "primary":
            reasons.append("non_primary")

        if self.require_complete and completeness != "complete":
            reasons.append("incomplete_orf")

        if self.remove_categories and category in self.remove_categories:
            reasons.append("removed_category:{}".format(category))

        if self.keep_categories and category not in self.keep_categories:
            reasons.append("category_not_allowed:{}".format(category))

        if self.keep_start_codons and start_codon not in self.keep_start_codons:
            reasons.append("start_codon_not_allowed:{}".format(start_codon))

        if aa_length < self.min_aa:
            reasons.append("too_short")

        if aa_length > self.max_aa:
            reasons.append("too_long")

        if self.kozak_pwm is not None:
            kozak_pwm_score = self._safe_float(record.get("kozak_pwm_score", "0"))

            if kozak_pwm_score < self.min_kozak_pwm_score:
                reasons.append("weak_kozak_pwm")

        if reasons:
            return False, ";".join(reasons)

        return True, "PASS"

    @staticmethod
    def _safe_int(value: str) -> int:
        """
        Convert string to int safely.
        """

        try:
            return int(value)
        except ValueError:
            return 0

    @staticmethod
    def _safe_float(value: str) -> float:
        """
        Convert string to float safely.
        """

        try:
            return float(value)
        except ValueError:
            return 0.0
