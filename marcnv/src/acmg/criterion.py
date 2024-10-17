import csv
import os
from dataclasses import dataclass

from marcnv.src.acmg import core

_ACMG_HEADER = [
    "Evidence Type",
    "Pretext",
    "Index",
    "Evidence",
    "Suggested points",
    "Min Score",
    "Max Score",
    "Tooltip",
]


@dataclass
class ACMGCriterion:
    evidence_type: str
    pretext: str
    evidence: str
    suggested_points: float | None
    min_score: float | None
    max_score: float | None
    tooltip: str


def _load_acmg_file(filepath: str) -> dict[str, ACMGCriterion]:
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"ACMG file not found: {filepath}")

    criteria: dict[str, ACMGCriterion] = {}
    with open(filepath, "r") as f:
        reader = csv.DictReader(
            f,
            delimiter="\t",
        )

        if reader.fieldnames != _ACMG_HEADER:
            raise ValueError(f"ACMG file header mismatch: {reader.fieldnames}. Expected: {_ACMG_HEADER}")

        for row in reader:
            criterion = ACMGCriterion(
                evidence_type=row["Evidence Type"] if row["Evidence Type"] else "",
                pretext=row["Pretext"] if row["Pretext"] else "",
                evidence=row["Evidence"] if row["Evidence"] else "",
                suggested_points=float(row["Suggested points"]) if row["Suggested points"] else None,
                min_score=float(row["Min Score"]) if row["Min Score"] else None,
                max_score=float(row["Max Score"]) if row["Max Score"] else None,
                tooltip=row["Tooltip"] if row["Tooltip"] else "",
            )
            if row["Index"] in criteria:
                raise ValueError(f'Duplicate index in ACMG criteria: {row['Index']}')
            criteria[row["Index"]] = criterion
    return criteria


def get_acmg_criteria(duplication: bool) -> dict[str, ACMGCriterion]:
    return _load_acmg_file(core.ACMG_GAIN_TSV_FILEPATH if duplication else core.ACMG_LOSS_TSV_FILEPATH)
