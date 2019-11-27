import itertools
import pysam
import random
import tabix
from collections import defaultdict

from .. import utils

def get_filtered_reads_from_region(samfile, gene, min_quality=30, apply_filters=True):
  for read in samfile.fetch(gene['seqname'], gene["start"], gene["end"]):
    if apply_filters and (read.is_qcfail or read.is_duplicate or read.is_secondary or read.is_unmapped or read.mapq < min_quality):
      continue
    yield read


def filter_gene(gene, gtf_tabix, only_consider_protein_genes=True):
  # potentially only consider protein coding genes.
  if only_consider_protein_genes and not "protein" in gene['attribute']:
    return False

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
    return False

  return True


def main(ngsfile, gene_model_file, n_genes=500, minimum_reads_per_gene=10):
  samfile = pysam.AlignmentFile(ngsfile, "rb")
  gtf_tabix = tabix.open(gene_model_file)
  gtf = utils.GFF(gene_model_file, feature_filter=["gene"])

  n_genes_we_tried = 0
  n_tested_genes = 0
  n_reads_observed = 0
  overall_evidence = defaultdict(int)
  gene_blacklist = set()

  while True:
    if n_tested_genes >= n_genes:
      break

    gene = random.choice(gtf.entries)

    if gene["seqname"] in gene_blacklist:
      continue

    if filter_gene(gene, gtf_tabix):
      continue

    gene_blacklist.add(gene["seqname"])
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
      for key in this_genes_evidence.keys():
        overall_evidence[key] += this_genes_evidence[key]

      n_tested_genes += 1
      n_reads_observed += reads_in_gene

    n_genes_we_tried += 1

  print("Number of genes we tried: {}".format(n_genes_we_tried))
  evidence_stranded_forward = overall_evidence["1++"] + overall_evidence["1--"] + overall_evidence["2+-"] + overall_evidence["2-+"]
  evidence_stranded_reverse = overall_evidence["1+-"] + overall_evidence["1-+"] + overall_evidence["2++"] + overall_evidence["2--"]
  total = evidence_stranded_forward + evidence_stranded_reverse
  return evidence_stranded_forward, evidence_stranded_reverse, total
