import annotation

from marcnv.src.acmg import classification
from marcnv.src.acmg.acmg_classify import MarCNVClassifier


def test_classify():
    annot = annotation.Annotation.load_from_json("tests/annotation_test.json.gz")

    classifier = MarCNVClassifier(annot)

    result = classifier.classify()
    assert result == classification.Prediction(
        score=0.15,
        criteria=[
            classification.SectionResult(
                section="1",
                option="1A",
                reason="The number of overlapping protein-coding genes (7) or enhancers (68) is more than zero.",
                score=0.0,
                evidence="1A. Contains protein-coding or other known functionally important elements",
            ),
            classification.SectionResult(
                section="2",
                option="2H",
                reason="Overlaps a gene with a predicted high risk of haplo-insufficiency (Gene name: MFSD14CP, Predictors: GHIS, ExAC)",
                score=0.15,
                evidence="2H. Multiple HI predictors suggest that AT LEAST ONE gene in the interval is haploinsufficient (HI)",
            ),
            classification.SectionResult(
                section="3",
                option="3A",
                reason="Overlaps 7 protein-coding genes.",
                score=0.0,
                evidence="3A. 0-24 genes",
            ),
            classification.SectionResult(
                section="4",
                option="4Skip",
                reason="Manual decision needed.",
                score=0.0,
                evidence="Skip if either your CNV overlapped with an established HI gene/region in Section 2, OR there have been no reports associating either the CNV or any genes within the CNV with human phenotypes caused by loss of function (LOF) or copy number loss",
            ),
            classification.SectionResult(
                section="5",
                option="5F",
                reason="No family history is available.",
                score=0.0,
                evidence="5F. Inheritance information is unavailable or uninformative.",
            ),
        ],
    )

    assert result.severity == classification.Severity.VOUS
