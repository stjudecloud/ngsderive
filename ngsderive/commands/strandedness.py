import csv
import itertools
import logging
import sys

from collections import defaultdict

from ..utils import NGSFile, NGSFileType, GFF

logger = logging.getLogger("strandedness")


def get_filtered_reads_from_region(samfile, gene, min_quality=30, apply_filters=True):
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


def disqualify_gene(gene, gff):
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


def get_reads_rg(read, default="unknown_read_group"):
    for (k, v) in read.tags:
        if k == "RG":
            return v

    return default


def get_predicted_strandedness(forward_evidence_pct, reverse_evidence_pct):
    predicted = "Inconclusive"
    if 0.4 <= forward_evidence_pct and forward_evidence_pct <= 0.6:
        predicted = "Unstranded"
    if 0.8 <= forward_evidence_pct:
        predicted = "Stranded-Forward"
    elif 0.8 <= reverse_evidence_pct:
        predicted = "Stranded-Reverse"

    return predicted


def determine_strandedness(
    ngsfilepath,
    gff,
    n_genes=100,
    min_mapq=30,
    minimum_reads_per_gene=10,
    split_by_rg=False,
    only_protein_coding_genes=True,
    max_iterations_per_try=1000,
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
            "Invalid file: {}. `strandedness` only supports aligned BAM files!".format(
                ngsfilepath
            )
        )
    samfile = ngsfile.handle

    n_tested_genes = 0
    n_reads_observed = 0
    checked_genes = set()

    read_groups = ["unknown_read_group"]
    if "RG" in samfile.header:
        read_groups += [rg["ID"] for rg in samfile.header["RG"]]

    overall_evidence = defaultdict(dict)

    for rg in read_groups:
        overall_evidence[rg] = defaultdict(int)

    logger.debug("Starting sampling...")
    total_iterations = 0
    while True:
        total_iterations += 1
        if total_iterations > max_iterations_per_try:
            logger.warn("Max iterations reached! Moving forward with prediction.")
            break

        if n_tested_genes >= n_genes:
            break

        gene = gff.sample()

        if gene["gene_id"] in checked_genes:
            continue

        if disqualify_gene(gene, gff):
            continue

        logging.debug("== Candidate Gene ==")
        logging.debug("  [*] ID: {}".format(gene["gene_id"]))
        logging.debug(
            "  [*] Location: {}:{}-{} ({})".format(
                gene["seqname"], gene["start"], gene["end"], gene["strand"]
            )
        )
        logging.debug("  [*] Actions:")

        checked_genes.add(gene["gene_id"])
        relevant_reads = get_filtered_reads_from_region(
            samfile, gene, min_quality=min_mapq
        )

        reads_in_gene = 0
        this_genes_evidence = defaultdict(dict)
        for rg in read_groups:
            this_genes_evidence[rg] = defaultdict(int)

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
            this_genes_evidence[this_reads_rg][reads_observed_state] += 1

        if reads_in_gene >= minimum_reads_per_gene:
            logging.debug(
                "    - Sufficient read count ({} >= {})".format(
                    reads_in_gene, minimum_reads_per_gene
                )
            )
            rg_log = ""
            if not this_genes_evidence["unknown_read_group"]:
                this_genes_evidence.pop("unknown_read_group")
            for rg, states in this_genes_evidence.items():
                state_logs = ";".join([f"{state}={n}" for state, n in states.items()])
                rg_log += f"{rg}:{state_logs} "

            logging.debug("    - {}".format(rg_log))

            for rg in this_genes_evidence.keys():
                for state in this_genes_evidence[rg].keys():
                    overall_evidence[rg][state] += this_genes_evidence[rg][state]

            n_tested_genes += 1
            n_reads_observed += reads_in_gene
        else:
            logging.debug(
                "    - Read count too low ({} < {})".format(
                    reads_in_gene, minimum_reads_per_gene
                )
            )

    if split_by_rg:
        results = []
        for rg in read_groups:
            evidence_stranded_forward = (
                overall_evidence[rg]["1++"]
                + overall_evidence[rg]["1--"]
                + overall_evidence[rg]["2+-"]
                + overall_evidence[rg]["2-+"]
            )
            evidence_stranded_reverse = (
                overall_evidence[rg]["1+-"]
                + overall_evidence[rg]["1-+"]
                + overall_evidence[rg]["2++"]
                + overall_evidence[rg]["2--"]
            )
            total_reads = evidence_stranded_forward + evidence_stranded_reverse
            if total_reads == 0 and rg == "unknown_read_group":
                continue

            forward_pct = (
                0
                if total_reads <= 0
                else round(evidence_stranded_forward / total_reads, 4)
            )
            reverse_pct = (
                0
                if total_reads <= 0
                else round(evidence_stranded_reverse / total_reads, 4)
            )
            predicted = get_predicted_strandedness(forward_pct, reverse_pct)

            results.append(
                {
                    "File": ngsfilepath,
                    "ReadGroup": rg,
                    "TotalReads": total_reads,
                    "ForwardPct": forward_pct,
                    "ReversePct": reverse_pct,
                    "Predicted": predicted,
                }
            )
        return results
    else:
        evidence_stranded_forward = 0
        evidence_stranded_reverse = 0

        for rg in read_groups:
            evidence_stranded_forward += (
                overall_evidence[rg]["1++"]
                + overall_evidence[rg]["1--"]
                + overall_evidence[rg]["2+-"]
                + overall_evidence[rg]["2-+"]
            )
            evidence_stranded_reverse += (
                overall_evidence[rg]["1+-"]
                + overall_evidence[rg]["1-+"]
                + overall_evidence[rg]["2++"]
                + overall_evidence[rg]["2--"]
            )

        total_reads = evidence_stranded_forward + evidence_stranded_reverse
        forward_pct = (
            0 if total_reads <= 0 else round(evidence_stranded_forward / total_reads, 4)
        )
        reverse_pct = (
            0 if total_reads <= 0 else round(evidence_stranded_reverse / total_reads, 4)
        )
        predicted = get_predicted_strandedness(forward_pct, reverse_pct)

        return [
            {
                "File": ngsfilepath,
                "TotalReads": total_reads,
                "ForwardPct": forward_pct,
                "ReversePct": reverse_pct,
                "Predicted": predicted,
            }
        ]


def main(
    ngsfiles,
    gene_model_file,
    outfile=sys.stdout,
    delimiter="\t",
    n_genes=100,
    minimum_reads_per_gene=10,
    only_protein_coding_genes=True,
    min_mapq=30,
    split_by_rg=False,
    max_tries=3,
    max_iterations_per_try=1000,
):
    logger.info("Arguments:")
    logger.info("  - Gene model file: {}".format(gene_model_file))
    logger.info("  - Number of genes: {}".format(n_genes))
    logger.info("  - Minimum reads per gene: {}".format(minimum_reads_per_gene))
    logger.info(
        "  - Only consider protein coding genes: {}".format(only_protein_coding_genes)
    )
    logger.info("  - Minimum MAPQ: {}".format(min_mapq))
    logger.info("  - Split by RG: {}".format(split_by_rg))

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

    logger.info("  - {} features processed.".format(len(gff.entries)))

    fieldnames = ["TotalReads", "ForwardPct", "ReversePct", "Predicted"]
    if split_by_rg:
        fieldnames = ["ReadGroup"] + fieldnames
    fieldnames = ["File"] + fieldnames

    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=delimiter)
    writer.writeheader()
    outfile.flush()

    for ngsfilepath in ngsfiles:
        tries_for_file = 0

        while True:
            tries_for_file += 1
            logger.info("Processing {}, try #{}...".format(ngsfilepath, tries_for_file))

            entries = determine_strandedness(
                ngsfilepath,
                gff,
                n_genes=n_genes,
                min_mapq=min_mapq,
                minimum_reads_per_gene=minimum_reads_per_gene,
                split_by_rg=split_by_rg,
                only_protein_coding_genes=only_protein_coding_genes,
                max_iterations_per_try=max_iterations_per_try,
            )

            entries_contains_inconclusive = False
            for entry in entries:
                if entry["Predicted"] == "Inconclusive":
                    entries_contains_inconclusive = True
                    break

            if entries_contains_inconclusive and tries_for_file < max_tries:
                continue

            for entry in entries:
                writer.writerow(entry)
                outfile.flush()

            break
