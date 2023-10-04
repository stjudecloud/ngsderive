import csv
import logging
import sys
from collections import defaultdict

from ..utils import GFF, NGSFile, NGSFileType, get_reads_rg, validate_read_group_info

logger = logging.getLogger("strandedness")


def get_filtered_reads_from_region(samfile, gene, min_quality, apply_filters=True):
    for read in samfile.fetch(gene["seqname"], gene["start"], gene["end"]):
        if apply_filters and (
            read.is_qcfail
            or read.is_duplicate
            or read.is_secondary
            or read.is_unmapped
            or read.mapq < min_quality
        ):
            continue
        yield read


def disqualify_gene(gene, gff, samfile):
    if gene["seqname"] not in samfile.references:
        return True

    # if there are overlapping features on the positive and negative strand
    # ignore this gene.
    hits = gff.query(gene["seqname"], gene["start"], gene["end"])

    has_positive_gene = None
    has_negative_gene = None

    for hit in hits:
        if hit["feature"] == "gene":
            continue
        if hit["strand"] == "+":
            has_positive_gene = True
        elif hit["strand"] == "-":
            has_negative_gene = True

        if has_positive_gene and has_negative_gene:
            break

    if has_positive_gene and has_negative_gene:
        return True

    return False


def get_predicted_strandedness(forward_evidence_pct, reverse_evidence_pct):
    if (
        forward_evidence_pct == 0 and reverse_evidence_pct == 0
    ):  # no evidence for either
        return "Unknown"

    if 40 <= forward_evidence_pct <= 60:
        return "Unstranded"
    # This second Unstranded check is redundant with the first check,
    # but more explicit
    if 40 <= reverse_evidence_pct <= 60:
        return "Unstranded"
    if 80 <= forward_evidence_pct:
        return "Stranded-Forward"
    if 80 <= reverse_evidence_pct:
        return "Stranded-Reverse"

    return "Inconclusive"


def determine_strandedness(
    ngsfilepath,
    gff,
    n_genes,
    min_mapq,
    minimum_reads_per_gene,
    split_by_rg,
    max_iterations_per_try,
    checked_genes,
    overall_evidence,
):
    try:
        ngsfile = NGSFile(ngsfilepath)
    except FileNotFoundError:
        result = {
            "File": ngsfilepath,
            "TotalReads": "N/A",
            "ForwardPct": "N/A",
            "ReversePct": "N/A",
            "Predicted": "Error opening file.",
        }

        if split_by_rg:
            result["ReadGroup"] = "N/A"

        return [result]

    if ngsfile.filetype != NGSFileType.BAM:
        raise RuntimeError(
            f"Invalid file: {ngsfilepath}. `strandedness` only supports aligned BAM files!"
        )
    samfile = ngsfile.handle

    n_tested_genes = 0
    n_reads_observed = 0
    checked_genes = set()

    logger.debug("Starting sampling...")
    total_iterations = 0
    while True:
        total_iterations += 1
        if total_iterations > max_iterations_per_try:
            logger.warning("Max iterations reached! Moving forward with prediction.")
            break

        if n_tested_genes >= n_genes:
            break

        gene = gff.sample()

        if gene["gene_id"] in checked_genes:
            continue

        if disqualify_gene(gene, gff, samfile):
            checked_genes.add(gene["gene_id"])
            continue

        logger.debug("== Candidate Gene ==")
        logger.debug(f"  [*] ID: {gene['gene_id']}")
        logger.debug(
            f"  [*] Location: {gene['seqname']}:{gene['start']}-{gene['end']} ({gene['strand']})"
        )

        checked_genes.add(gene["gene_id"])
        relevant_reads = get_filtered_reads_from_region(
            samfile, gene, min_quality=min_mapq
        )

        reads_in_gene = 0
        this_genes_evidence = defaultdict(lambda: defaultdict(int))

        for read in relevant_reads:
            reads_in_gene += 1
            if read.is_paired:
                if read.is_read1:
                    read_id = "1"
                elif read.is_read2:
                    read_id = "2"
                else:
                    raise RuntimeError("Read is not read 1 or read 2?")
            else:
                # SE reads are equivalent to just assuming the read is read 1.
                read_id = "1"

            if not read.is_reverse:
                read_strand_id = "+"
            else:
                read_strand_id = "-"

            gene_strand_id = gene["strand"]
            reads_observed_state = read_id + read_strand_id + gene_strand_id
            this_reads_rg = get_reads_rg(read)
            this_genes_evidence["overall"][reads_observed_state] += 1
            this_genes_evidence[this_reads_rg][reads_observed_state] += 1

        if reads_in_gene >= minimum_reads_per_gene:
            logger.debug(
                f"    - Sufficient read count ({reads_in_gene} >= {minimum_reads_per_gene})"
            )
            rg_log = ""
            if not this_genes_evidence["unknown_read_group"]:
                this_genes_evidence.pop("unknown_read_group")
            for rg, states in this_genes_evidence.items():
                state_logs = ";".join([f"{state}={n}" for state, n in states.items()])
                rg_log += f"{rg}:{state_logs} "

            logger.debug(f"    - {rg_log}")

            for rg, evidence in this_genes_evidence.items():
                for state, n in evidence.items():
                    overall_evidence[rg][state] += n

            n_tested_genes += 1
            n_reads_observed += reads_in_gene
        else:
            logger.debug(
                f"    - Read count too low ({reads_in_gene} < {minimum_reads_per_gene})"
            )

    rgs_in_header_not_in_seq = validate_read_group_info(
        set(overall_evidence.keys()),
        samfile.header,
    )
    for rg in rgs_in_header_not_in_seq:
        overall_evidence[rg] = defaultdict(int)  # init rg to all zeroes

    if split_by_rg:
        results = []
        for rg, rg_evidence in overall_evidence.items():
            evidence_stranded_forward = (
                rg_evidence["1++"]
                + rg_evidence["1--"]
                + rg_evidence["2+-"]
                + rg_evidence["2-+"]
            )
            evidence_stranded_reverse = (
                rg_evidence["1+-"]
                + rg_evidence["1-+"]
                + rg_evidence["2++"]
                + rg_evidence["2--"]
            )
            total_reads = evidence_stranded_forward + evidence_stranded_reverse
            if total_reads == 0 and rg == "unknown_read_group":
                continue

            forward_pct = (
                0
                if total_reads <= 0
                else round(evidence_stranded_forward / total_reads * 100, 2)
            )
            reverse_pct = (
                0
                if total_reads <= 0
                else round(evidence_stranded_reverse / total_reads * 100, 2)
            )
            predicted = get_predicted_strandedness(forward_pct, reverse_pct)

            results.append(
                {
                    "File": ngsfilepath,
                    "ReadGroup": rg,
                    "TotalReads": total_reads,
                    "ForwardPct": str(forward_pct) + "%",
                    "ReversePct": str(reverse_pct) + "%",
                    "Predicted": predicted,
                }
            )
        return (results, checked_genes, overall_evidence)

    evidence_stranded_forward = (
        overall_evidence["overall"]["1++"]
        + overall_evidence["overall"]["1--"]
        + overall_evidence["overall"]["2+-"]
        + overall_evidence["overall"]["2-+"]
    )
    evidence_stranded_reverse = (
        overall_evidence["overall"]["1+-"]
        + overall_evidence["overall"]["1-+"]
        + overall_evidence["overall"]["2++"]
        + overall_evidence["overall"]["2--"]
    )
    total_reads = evidence_stranded_forward + evidence_stranded_reverse
    forward_pct = (
        0
        if total_reads <= 0
        else round(evidence_stranded_forward / total_reads * 100, 2)
    )
    reverse_pct = (
        0
        if total_reads <= 0
        else round(evidence_stranded_reverse / total_reads * 100, 2)
    )
    predicted = get_predicted_strandedness(forward_pct, reverse_pct)

    return (
        [
            {
                "File": ngsfilepath,
                "TotalReads": total_reads,
                "ForwardPct": str(forward_pct) + "%",
                "ReversePct": str(reverse_pct) + "%",
                "Predicted": predicted,
            }
        ],
        checked_genes,
        overall_evidence,
    )


def main(
    ngsfiles,
    gene_model_file,
    outfile,
    n_genes,
    minimum_reads_per_gene,
    only_protein_coding_genes,
    min_mapq,
    split_by_rg,
    max_tries,
    max_iterations_per_try,
):
    logger.info("Arguments:")
    logger.info(f"  - Gene model file: {gene_model_file}")
    logger.info(f"  - Number of genes: {n_genes}")
    logger.info(f"  - Minimum reads per gene: {minimum_reads_per_gene}")
    logger.info(f"  - Only consider protein coding genes: {only_protein_coding_genes}")
    logger.info(f"  - Minimum MAPQ: {min_mapq}")
    logger.info(f"  - Split by RG: {split_by_rg}")

    if max_iterations_per_try < n_genes:
        logger.error(
            "Max iteration per try cannot be less than number of genes to search!"
        )
        sys.exit(1)

    logger.info("Reading gene model...")
    gff = GFF(
        gene_model_file,
        feature_type="gene",
        store_results=True,
        need_tabix=True,
        only_protein_coding_genes=only_protein_coding_genes,
    )

    logger.info(f"  - {len(gff.entries)} features processed.")

    fieldnames = ["TotalReads", "ForwardPct", "ReversePct", "Predicted"]
    if split_by_rg:
        fieldnames = ["ReadGroup"] + fieldnames
    fieldnames = ["File"] + fieldnames

    writer = None

    for ngsfilepath in ngsfiles:
        tries_for_file = 0
        checked_genes = set()
        overall_evidence = defaultdict(lambda: defaultdict(int))

        while True:
            tries_for_file += 1
            logger.info(f"Processing {ngsfilepath}, try #{tries_for_file}...")

            entries, checked_genes, overall_evidence = determine_strandedness(
                ngsfilepath,
                gff,
                n_genes=n_genes,
                min_mapq=min_mapq,
                minimum_reads_per_gene=minimum_reads_per_gene,
                split_by_rg=split_by_rg,
                max_iterations_per_try=max_iterations_per_try,
                checked_genes=checked_genes,
                overall_evidence=overall_evidence,
            )

            entries_contains_inconclusive = False
            for entry in entries:
                if entry["Predicted"] == "Inconclusive":
                    entries_contains_inconclusive = True
                    break

            if entries_contains_inconclusive and tries_for_file < max_tries:
                continue

            if not writer:
                writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
            for entry in entries:
                writer.writerow(entry)
                outfile.flush()

            break
