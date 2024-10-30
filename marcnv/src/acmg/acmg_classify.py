import sys
from dataclasses import dataclass, field
from typing import Any

import annotation

from marcnv.src.acmg import classification, core, criterion


def evaluate_transcript(transcript: annotation.TranscriptRegion) -> tuple[str | None, str | None]:
    # If transcript's 5' is in CNV
    if transcript["flag_five_inside"]:
        if len(transcript["five_prime_utr_overlaps"]) > 0:
            text_5p = f'Overlaps {len(transcript["five_prime_utr_overlaps"])} start (5\') regions with a total length of {sum(transcript["five_prime_utr_overlaps"])}bp.'
        else:
            text_5p = "Overlaps no start (5') regions."

        # Count the number overlapping CDS
        if len(transcript["cds_overlaps"]) == 0:
            return "2C-2", f"{text_5p}\nOverlaps no CDS regions."
        else:
            return (
                "2C-1",
                f'{text_5p}\nOverlaps {len(transcript["cds_overlaps"])} CDS regions with a total length of {sum(transcript["cds_overlaps"])}bp.',
            )

    # If transcript's 3' is in CNV
    if transcript["flag_three_inside"]:
        if len(transcript["three_prime_utr_overlaps"]) > 0:
            text_3p = f'Overlaps {len(transcript["three_prime_utr_overlaps"])} end (3\') regions with total length of {sum(transcript["three_prime_utr_overlaps"])}bp.'
        else:
            text_3p = "Overlaps no end (3') regions."

        # Count the number overlapping CDS
        if len(transcript["cds_overlaps"]) == 0:
            return "2D-1", f"{text_3p}\nOverlaps no CDS regions."
        elif len(transcript["cds_overlaps"]) == 1:
            return "2D-3", f'{text_3p}\nOnly last CDS is involved with overlap of {transcript["cds_overlaps"][0]}bp.'
        else:
            return (
                "2D-4",
                f'{text_3p}\nOverlaps {len(transcript["cds_overlaps"])} CDS regions with a total length of {sum(transcript["cds_overlaps"])}bp.',
            )

    # If CNV is in the gene (not very likely)
    if transcript["flag_contained"]:
        return (
            "2E",
            f'A transcript completely contains the target CNV.\n'
            f'Overlaps {len(transcript["cds_overlaps"])} CDS regions with a total length of {sum(transcript["cds_overlaps"])}bp.',
        )

    # TODO ?
    # print(f'WARNING: transcript {transcript["ID"]} is not evaluated - it is wholly inside the CNV, although the gene is not.')
    return None, None


def evaluate_gene(transcript_regions: list[annotation.TranscriptRegion]) -> tuple[str, str, str]:
    # Evaluate all transcripts and pick the most severe (evaluate from the longest)
    most_severe_transcript_name = "UNKNOWN"
    most_severe_transcript_reason = ""
    most_severe_transcript_option = ""
    severity_seq = ["2C-1", "2D-4", "2D-3", "2C-2", "2E", "2D-2", "2D-1", ""]
    for transcript in sorted(transcript_regions, key=lambda x: x["length"], reverse=True):
        option, reason = evaluate_transcript(transcript)
        if option and reason:
            rev_severity = severity_seq.index(option)
            if rev_severity < severity_seq.index(most_severe_transcript_option):
                most_severe_transcript_name = transcript["identifier"]
                most_severe_transcript_reason = reason
                most_severe_transcript_option = option

    return most_severe_transcript_option, most_severe_transcript_reason, most_severe_transcript_name


def get_reason_from_benigncnvs(cnv: annotation.CNVRegionAnnotation, benign_cnvs: list[Any]) -> str:
    reasons: list[str] = []
    for benign_cnv in benign_cnvs:
        overlap = cnv.get_overlap_with_region(benign_cnv["start"], benign_cnv["end"])
        reasons.append(
            f'Accession number: {benign_cnv["variantaccession"]}, pubmedid: {benign_cnv.get("pubmedid", "UNKNOWN")}. '
            f'Overlap: {overlap}bp ({overlap / float(cnv.length) * 100.0:.1f}%).'
        )
    reason = f"Overlaps {len(benign_cnvs)} established benign CNV, but contains additional protein coding genes."
    reason += " Benign CNVs: (displaying only first 5)\n" if len(reasons) > 5 else " Benign CNVs:\n"
    reason += "\n".join(reasons[:5])
    return reason


def evaluate_section2_duplication(annot: annotation.Annotation) -> tuple[str, str]:
    cont_eval: tuple[str, str] | None = None

    inside_only_regions = annot.get_triplosensitivity_regions(
        annotation.enums.Overlap.CONTAINED_INSIDE, core.HI_TS_SCORES
    )
    print(f"Evaluating section 2: {inside_only_regions=}", file=sys.stderr)
    if (count := len(inside_only_regions)) > 0:
        scores = [int(r["Triplosensitivity Score"]) for r in inside_only_regions]
        names = [r["ISCA Region Name"] for r in inside_only_regions]
        if len(set(scores)) == 1:
            detail = f'Found {count} such regions: {', '.join(names)}; all with TS score {scores[0]}.'
        else:
            detail = f'Found {count} such regions: {', '.join(names)} with TS scores: {scores}, respectively.'
        return "2A", "Completely contains at least one established TS region. " + detail

    inside_only_genes = annot.get_triplosensitivity_genes(annotation.enums.Overlap.CONTAINED_INSIDE, core.HI_TS_SCORES)
    print(f"Evaluating section 2: {inside_only_genes=}", file=sys.stderr)
    if (count := len(inside_only_genes)) > 0:
        scores = [int(g["Triplosensitivity Score"]) for g in inside_only_genes]
        names = [g["Gene Symbol"] for g in inside_only_genes]
        if len(set(scores)) == 1:
            detail = f'Found {count} such genes: {', '.join(names)}; all with TS score {scores[0]}.'
        else:
            detail = f'Found {count} such genes: {', '.join(names)} with TS scores: {scores}, respectively.'
        return "2A", "Completely contains at least one established TS gene. " + detail

    all_ts_regions = annot.get_triplosensitivity_regions(annotation.enums.Overlap.ANY, core.HI_TS_SCORES)
    print(f"Evaluating section 2: {all_ts_regions=}", file=sys.stderr)
    if len(all_ts_regions) > 0:
        max_overlap = max([annot.cnv.get_overlap_with_region(r["start"], r["end"]) for r in all_ts_regions])
        names = [r["ISCA Region Name"] for r in all_ts_regions]
        names_str = (", ".join(names[:10]) + ", ...") if len(names) > 10 else ", ".join(names)
        reason = (
            f"Partially overlaps established TS region(s) ({len(names)} - {names_str}, max. overlap {max_overlap}bp)."
        )
        return "2B", reason

    # Load (protein) coding genes and genes on breakpoints
    gene_names = sorted([gene["gene_name"] for gene in annot.get_genes(overlap=annotation.enums.Overlap.ANY)])
    protein_genes_on_breakpoints = annot.get_genes(
        gene_type="protein_coding", overlap=annotation.enums.Overlap.START_OR_END
    )
    print(f"Evaluating section 2: {gene_names=}, {protein_genes_on_breakpoints=}", file=sys.stderr)

    # Compare protein coding genes in benign CNV - searching for identical gene content
    benign_cnvs = annot.get_benign_cnvs_gs_outer(frequency_threshold=core.MIN_FREQUENCY_BENIGN)
    print(f"Evaluating section 2: {benign_cnvs=}", file=sys.stderr)
    for benign_cnv in benign_cnvs:
        genes_cnv = sorted([gene["gene_name"] for gene in benign_cnv["genes"]])
        # here we assume that benign CNV start and end in between genes
        if genes_cnv == gene_names and len(protein_genes_on_breakpoints) == 0:
            reason = f'Identical in gene content ({len(gene_names)} genes) to a benign CNV gain (variant_accesion={benign_cnv["variantaccession"]}).'
            return "2C", reason

    # Smaller than established benign CNV, breakpoints are OK
    for benign_cnv in benign_cnvs:
        if (
            annot.cnv.is_overlapping(benign_cnv["start"], benign_cnv["end"], annotation.enums.Overlap.SPAN_ENTIRE)
            and len(protein_genes_on_breakpoints) == 0
        ):
            reason = (
                f'Smaller than an established benign CNV gain (variant_accesion={benign_cnv["variantaccession"]}), breakpoints do not '
                f'interrupt protein-coding genes.'
            )
            return "2D", reason

    # Smaller than established benign CNV, breakpoints are NOT OK
    for benign_cnv in benign_cnvs:
        if (
            annot.cnv.is_overlapping(benign_cnv["start"], benign_cnv["end"], annotation.enums.Overlap.SPAN_ENTIRE)
            and len(protein_genes_on_breakpoints) > 0
        ):
            gene_names = [g["gene_name"] for g in protein_genes_on_breakpoints]
            reason = (
                f'Smaller than an established benign CNV gain (variant_accesion={benign_cnv["variantaccession"]}), but breakpoints '
                f'potentially interrupt protein-coding gene(s) ({", ".join(gene_names)}).'
            )
            if cont_eval is None:
                cont_eval = "2E", reason

    # Larger than established benign CNV, identical protein coding genes
    protein_genes = annot.get_genes(gene_type="protein_coding")
    print(f"Evaluating section 2: {protein_genes=}", file=sys.stderr)
    for benign_cnv in benign_cnvs:
        if annot.cnv.is_overlapping(benign_cnv["start"], benign_cnv["end"], annotation.enums.Overlap.CONTAINED_INSIDE):
            other_genes = [
                g for g in protein_genes if g["start"] >= benign_cnv["start"] or g["end"] <= benign_cnv["end"]
            ]
            if len(other_genes) == 0:
                reason = (
                    f'Larger than an established benign CNV gain (variant_accesion={benign_cnv["variantaccession"]}), does not include '
                    f'additional protein-coding genes.'
                )
                return "2F", reason

    # Overlapping a benign CNV
    if len(benign_cnvs) > 0 and cont_eval is None:
        cont_eval = "2G", get_reason_from_benigncnvs(annot.cnv, benign_cnvs)

    # Complete containment of an HI gene
    for hi_gene in annot.get_haploinsufficient_genes(annotation.enums.Overlap.CONTAINED_INSIDE, core.HI_TS_SCORES):
        print(f"Evaluating section 2: inside_only{hi_gene=}", file=sys.stderr)
        reason = f'Completely contains an established HI gene {hi_gene["Gene Symbol"]} with HI score {hi_gene["Haploinsufficiency Score"]}.'
        if cont_eval is None:
            cont_eval = "2H", reason

    # Breakpoints with a hi_gene
    for hi_gene in annot.get_haploinsufficient_genes(annotation.enums.Overlap.SPAN_ENTIRE, core.HI_TS_SCORES):
        print(f"Evaluating section 2: span_whole_only{hi_gene=}", file=sys.stderr)
        reason = (
            f'Both breakpoints are within the same HI gene {hi_gene["Gene Symbol"]} - gene-level sequence variant, possibly resulting '
            f'in loss of function (LOF).'
        )
        return "2I", reason

    for hi_gene in annot.get_haploinsufficient_genes(annotation.enums.Overlap.START_OR_END, core.HI_TS_SCORES):
        print(f"Evaluating section 2: partial_both{hi_gene=}", file=sys.stderr)
        reason = f'One breakpoint is within an established HI gene {hi_gene["Gene Symbol"]}, the patient’s phenotype is unknown.'
        if cont_eval is None:
            cont_eval = "2J", reason

    # Breakpoint within other gene
    genes_on_breakpoints = annot.get_genes(overlap=annotation.enums.Overlap.START_OR_END)
    if len(genes_on_breakpoints) > 0:
        gene_names = [g["gene_name"] for g in genes_on_breakpoints]
        reason = f"One or both breakpoints are within gene(s) of no established clinical significance ({gene_names})."
        if cont_eval is None:
            cont_eval = "2L", reason

    # Skip section if there are no supporting data
    if cont_eval:
        return cont_eval
    else:
        reason = "The section is skipped due to lack of supporting data (no TS/HI regions/genes and benign CNVs)."
        return "2Skip", reason


def evaluate_section2(annot: annotation.Annotation) -> tuple[str, str]:
    # cont_evaluation clauses:
    cont_eval: tuple[str, str] | None = None

    if annot.cnv.is_duplication:
        return evaluate_section2_duplication(annot)

    inside_only_regions = annot.get_haploinsufficient_regions(
        annotation.enums.Overlap.CONTAINED_INSIDE, core.HI_TS_SCORES
    )
    print(f"Evaluating section 2: {inside_only_regions=}", file=sys.stderr)
    if (count := len(inside_only_regions)) > 0:
        scores = [int(r["Haploinsufficiency Score"]) for r in inside_only_regions]
        names = [r["ISCA Region Name"] for r in inside_only_regions]
        if len(set(scores)) == 1:
            detail = f'Found {count} such regions: {', '.join(names)}; all with HI score {scores[0]}.'
        else:
            detail = f'Found {count} such regions: {', '.join(names)} with HI scores: {scores}, respectively.'
        return "2A", "Completely contains at least one established HI region. " + detail

    inside_only_genes = annot.get_haploinsufficient_genes(annotation.enums.Overlap.CONTAINED_INSIDE, core.HI_TS_SCORES)
    print(f"Evaluating section 2: {inside_only_genes=}", file=sys.stderr)
    if (count := len(inside_only_genes)) > 0:
        scores = [int(g["Haploinsufficiency Score"]) for g in inside_only_genes]
        names = [g["Gene Symbol"] for g in inside_only_genes]
        if len(set(scores)) == 1:
            detail = f'Found {count} such genes: {', '.join(names)}; all with HI score {scores[0]}.'
        else:
            detail = f'Found {count} such genes: {', '.join(names)} with HI scores: {scores}, respectively.'
        return "2A", "Completely contains at least one established HI gene. " + detail

    # Evaluation of every single HI gene:
    for hi_gene in annot.get_haploinsufficient_genes(annotation.enums.Overlap.ANY, core.HI_TS_SCORES):
        gene_info = annot.get_gene_by_name(hi_gene["Gene Symbol"])
        if gene_info is None:
            print(
                f'WARNING: gene {hi_gene["Gene Symbol"]} NOT FOUND in GenCode, probably mismatch in coordinates GenCode/Clingen.'
            )
            continue
        if gene_info["gene_type"] != "protein_coding":
            print(f'WARNING: evaluated GENE TYPE is {gene_info["gene_type"]}')

        transcript_regions = annot.get_gene_transcript_regions(hi_gene["Gene Symbol"])
        print(f"Evaluating section 2: {transcript_regions=} for {hi_gene}", file=sys.stderr)

        # Evaluate the gene and all its transcripts, return the worst value
        option, partial_reason, transcript_name = evaluate_gene(transcript_regions)
        reason = f'{partial_reason} (Gene name: {hi_gene['Gene Symbol']}, Transcript ID: {transcript_name})'

        if option != "" and option != "2D-1":  # 2D-1 is "Continue evaluation"
            return option, reason
        if option == "2D-1":
            cont_eval = option, reason

    # Partial overlap with hi_range
    hi_regions = annot.get_haploinsufficient_regions(annotation.enums.Overlap.ANY, core.HI_TS_SCORES)
    print(f"Evaluating section 2: {hi_regions=}", file=sys.stderr)
    if len(hi_regions) > 0:
        names = [r["ISCA Region Name"] for r in hi_regions]
        max_overlap = max([annot.cnv.get_overlap_with_region(r["start"], r["end"]) for r in hi_regions])
        names_str = (", ".join(names[:10]) + ", ...") if len(names) > 10 else ", ".join(names)
        reason = (
            f"Partially overlaps established HI region(s) ({len(names)} - {names_str}, max. overlap {max_overlap}bp)."
        )
        return "2B", reason

    # Completely contained within a benign CNV
    benign_cnvs = annot.get_benign_cnvs_gs_outer(core.MIN_FREQUENCY_BENIGN)
    print(f"Evaluating section 2: {benign_cnvs=}", file=sys.stderr)
    for benign_cnv in benign_cnvs:
        if annot.cnv.is_overlapping(benign_cnv["start"], benign_cnv["end"], annotation.enums.Overlap.SPAN_ENTIRE):
            reason = (
                f'An established benign CNV {benign_cnv["variantaccession"]} with pubmedid {benign_cnv.get("pubmedid", "UNKNOWN")} '
                f'completely contains the target CNV.'
            )
            return "2F", reason

    # Overlapping a benign CNV
    if len(benign_cnvs) > 0 and cont_eval is None:
        reason = get_reason_from_benigncnvs(annot.cnv, benign_cnvs)
        cont_eval = "2G", reason

    # check HI predictors
    high_risk_genes = annot.get_high_risk_loss_genes()
    print(f"Evaluating section 2: {high_risk_genes=}", file=sys.stderr)
    if len(high_risk_genes) > 0:
        gene = high_risk_genes[0]
        reason = (
            f'Overlaps a gene with a predicted high risk of haplo-insufficiency (Gene name: {gene["gene_name"]}, Predictors: '
            f'{", ".join(gene["risk_predictors"])})'
        )
        return "2H", reason

    # Skip section if no data available
    if cont_eval is not None:
        return cont_eval
    else:
        reason = "The section is skipped due to lack of supporting data (no HI regions/genes and benign CNVs)."
        return "2Skip", reason


@dataclass
class MarCNVClassifier:
    annot: annotation.Annotation
    acmg_criteria: dict[str, criterion.ACMGCriterion] = field(init=False)

    def __post_init__(self) -> None:
        self.acmg_criteria = criterion.get_acmg_criteria(self.annot.cnv.is_duplication)

    def _build_section_result(self, section: str, option: str, reason: str) -> classification.SectionResult:
        criterion = self.acmg_criteria[option]
        return classification.SectionResult(
            section=section,
            option=option,
            reason=reason,
            evidence=criterion.evidence,
            score=criterion.suggested_points if criterion.suggested_points is not None else 0,
        )

    def _evaluate_section1(self) -> classification.SectionResult:
        # Find number of genes and regulatory elements
        gene_count = len(self.annot.get_genes(gene_type="protein_coding"))
        enhancers_count = self.annot.count_regulatory_types()["enhancer"]

        print(f"Evaluating section 1: {gene_count=}, {enhancers_count=}", file=sys.stderr)
        # Pick final option and assign reason
        if gene_count + enhancers_count == 0:
            option = "1B"
            reason = "The number of overlapping protein-coding genes and regulatory elements is zero."
        else:
            option = "1A"
            reason = f"The number of overlapping protein-coding genes ({gene_count}) or enhancers ({enhancers_count}) is more than zero."
        return self._build_section_result(
            section="1",
            option=option,
            reason=reason,
        )

    def _evaluate_section2(self) -> classification.SectionResult:
        section = evaluate_section2(self.annot)
        return self._build_section_result(
            section="2",
            option=section[0],
            reason=section[1],
        )

    def _evaluate_section3(self) -> classification.SectionResult:
        protein_genes = self.annot.get_genes(gene_type="protein_coding")

        if self.annot.cnv.is_duplication:
            thresholds = core.DUPLICATION_GENES_THRESHOLDS
        else:
            thresholds = core.DELETION_GENES_THRESHOLDS

        print(f"Evaluating section 3: {len(protein_genes)=}", file=sys.stderr)
        if len(protein_genes) < thresholds[0]:
            option = "3A"
        elif len(protein_genes) < thresholds[1]:
            option = "3B"
        else:
            option = "3C"

        return self._build_section_result(
            section="3",
            option=option,
            reason=f"Overlaps {len(protein_genes)} protein-coding genes.",
        )

    def _evaluate_section4(self) -> classification.SectionResult:
        reason = "Manual decision needed."
        option = "4Skip"

        common_variability_regions = self.annot.get_common_variability_regions()
        print(f"Evaluating section 4: {common_variability_regions=}", file=sys.stderr)
        if len(common_variability_regions) >= 1:
            region = [r for r in common_variability_regions if r["population"] == "nfe"][0]
            reason = f'Common population variation {self.annot.cnv.genomic_coord} for population {region["population"]} has frequency of {region["frequency"] * 100.0}%.'
            option = "4O"

        return self._build_section_result(
            section="4",
            option=option,
            reason=reason,
        )

    def _evaluate_section5(self) -> classification.SectionResult:
        return self._build_section_result(
            section="5",
            option="5F",
            reason="No family history is available.",
        )

    def classify(self) -> classification.Prediction:
        sections = [
            self._evaluate_section1(),
            self._evaluate_section2(),
            self._evaluate_section3(),
            self._evaluate_section4(),
            self._evaluate_section5(),
        ]

        final_score = round(sum([section.score for section in sections]), 8)

        return classification.Prediction(
            score=final_score,
            criteria=sections,
        )
