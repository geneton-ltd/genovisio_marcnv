import bisect
import enum
from dataclasses import dataclass, field


class Severity(enum.StrEnum):
    BENIGN = "Benign"
    LBENIGN = "Likely benign"
    VOUS = "Uncertain"
    LPATHOGENIC = "Likely pathogenic"
    PATHOGENIC = "Pathogenic"

    @classmethod
    def from_score(cls, score: float) -> "Severity":
        """Returns the enum value from score."""
        epsilon = 0.00000001
        severity_thresholds = [-1.0 + epsilon, -0.9 + epsilon, 0.9, 1.0]
        severity_index = bisect.bisect(severity_thresholds, score)
        vals = list(cls)
        return cls(vals[severity_index])


@dataclass
class SectionResult:
    section: str
    option: str
    reason: str
    score: float
    evidence: str


@dataclass
class Prediction:
    score: float
    criteria: list[SectionResult]
    severity: Severity = field(init=False)

    def __post_init__(self) -> None:
        self.severity = Severity.from_score(self.score)
