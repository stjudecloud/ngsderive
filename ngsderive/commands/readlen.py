import csv
import itertools
import pysam
from collections import defaultdict

def main(ngsfiles, outfile, n_samples=10000, cutoff=0.7):
  writer = csv.DictWriter(outfile, fieldnames=["File", "ReadLength"], delimiter="\t")
  writer.writeheader()

  for ngsfile in ngsfiles:
    readlengths = defaultdict(int)
    samfile = pysam.AlignmentFile(ngsfile, "rb")

    # accumulate read lengths
    for read in itertools.islice(samfile, n_samples):
      readlengths[len(read.query)] += 1

    # normalize values
    for readlen in readlengths:
      readlengths[readlen] /= n_samples

    max_readlen = sorted(readlengths.keys(), reverse=True)[0]

    # ensure % of max readlength > [cutoff]
    # if not, cannot determine
    result = {
      "File": ngsfile,
      "ReadLength": max_readlen if readlengths[max_readlen] > cutoff else -1
    }

    writer.writerow(result)
