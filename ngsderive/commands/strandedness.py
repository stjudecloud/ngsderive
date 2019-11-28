import csv
import itertools
import logging
import pysam
import random
import tabix

from collections import defaultdict

from .. import utils

logger = logging.getLogger('strandedness')

def get_filtered_reads_from_region(samfile, gene, min_quality=30, apply_filters=True):
  for read in samfile.fetch(gene['seqname'], gene["start"], gene["end"]):
    if apply_filters and (read.is_qcfail or read.is_duplicate or read.is_secondary or read.is_unmapped or read.mapq < min_quality):
      continue
    yield read


def disqualify_gene(gene, gtf_tabix, only_consider_protein_genes=True):
  # potentially only consider protein coding genes.
  if only_consider_protein_genes:
    if not "attr_gene_type" in gene: 
      return True
    if not "protein" in gene['attr_gene_type']:
      return True 

  # if there are overlapping features on the positive and negative strand
  # ignore this gene.
  hits = gtf_tabix.query(gene['seqname'], gene['start'], gene['end'])
  has_positive_gene = None
  has_negative_gene = None

  for hit in hits:
    # must be a gene
    if not hit[2] == "gene":
      continue

    if hit[6] == "+":
      has_positive_gene = hit
    elif hit[6] == "-":
      has_negative_gene = hit
    
    if has_positive_gene and has_negative_gene:
      break
  
  if has_positive_gene and has_negative_gene:
    return True 

  return False


def main(ngsfiles, gtf, gene_model_file, outfile, n_genes=100, minimum_reads_per_gene=50):
  writer = csv.DictWriter(outfile, fieldnames=["File", "TotalReads", "ForwardPct", "ReversePct", "Predicted"], delimiter="\t")
  writer.writeheader()

  for ngsfile in ngsfiles:
    logger.debug("Loading SAM file.")
    samfile = pysam.AlignmentFile(ngsfile, "rb")
    logger.debug("Creating opening tabix version of gene model.")
    gtf_tabix = tabix.open(gene_model_file)

    n_tested_genes = 0
    n_reads_observed = 0
    overall_evidence = defaultdict(int)
    gene_blacklist = set()

    logger.debug("Starting sampling (n={}).".format(n_genes))
    while True:
      if n_tested_genes >= n_genes:
        break

      gene = random.choice(gtf.entries)

      if gene["attr_gene_id"] in gene_blacklist:
        continue

      if disqualify_gene(gene, gtf_tabix):
        continue

      logging.debug("== Candidate Gene ==")
      logging.debug("  [*] Name: {}".format(gene['attr_gene_name']))
      logging.debug("  [*] Location: {}:{}-{} ({})".format(gene['seqname'], gene['start'], gene['end'], gene['strand']))
      logging.debug("  [*] Actions:")

      gene_blacklist.add(gene["attr_gene_id"])
      relevant_reads = get_filtered_reads_from_region(samfile, gene)

      reads_in_gene = 0
      this_genes_evidence = defaultdict(int)

      for read in relevant_reads:
        reads_in_gene += 1
        if not read.is_paired:
          raise RuntimeError("This tool currently only works for paired-end data! Please contact the author if you'd like SE data to be supported.")

        if read.is_read1:
          read_id = "1"
        elif read.is_read2:
          read_id = "2"
        else:
          raise RuntimeError("Read is not read 1 or read 2?")
      
        if not read.is_reverse:
          read_strand_id = "+"
        else:
          read_strand_id = "-"

        gene_strand_id = gene['strand']
        this_genes_evidence[read_id + read_strand_id + gene_strand_id] += 1

      if reads_in_gene >= minimum_reads_per_gene:
        logging.debug("    - Sufficient read count ({} >= {})".format(reads_in_gene, minimum_reads_per_gene))
        logging.debug("    - {}".format(" ".join(["{}:{}".format(k, v) for k, v in this_genes_evidence.items()])))
        for key in this_genes_evidence.keys():
          overall_evidence[key] += this_genes_evidence[key]

        n_tested_genes += 1
        n_reads_observed += reads_in_gene
      else:
        logging.debug("    - Read count too low ({} < {})".format(reads_in_gene, minimum_reads_per_gene))


    evidence_stranded_forward = overall_evidence["1++"] + overall_evidence["1--"] + overall_evidence["2+-"] + overall_evidence["2-+"]
    evidence_stranded_reverse = overall_evidence["1+-"] + overall_evidence["1-+"] + overall_evidence["2++"] + overall_evidence["2--"]
    total_reads = evidence_stranded_forward + evidence_stranded_reverse

    forward_pct = round(evidence_stranded_forward / total_reads, 4)
    reverse_pct = round(evidence_stranded_reverse / total_reads, 4)

    predicted = "Inconclusive"
    if 0.4 <= forward_pct and forward_pct <= 0.6:
      predicted = "Unstranded"
    if 0.8 <= forward_pct:
      predicted = "Stranded-Forward"
    elif 0.8 <= reverse_pct:
      predicted = "Stranded-Reverse"

    result = {
      "File": ngsfile, 
      "TotalReads": total_reads, 
      "ForwardPct": forward_pct, 
      "ReversePct": reverse_pct,
      "Predicted": predicted
    }

    writer.writerow(result)
