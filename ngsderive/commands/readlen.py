import itertools
import pysam
from collections import defaultdict

def main(ngsfile, n_samples=10000, cutoff=0.8):
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
  if readlengths[max_readlen] > cutoff:
    return max_readlen
  else:
    return -1
