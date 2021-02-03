import csv
import itertools
import pysam
import sys

import logging
from collections import defaultdict

from ..utils import NGSFile, NGSFileType

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# sets are created by doing PHRED+33 decoding of each encoding's ASCII ranges
SANGER_SET = set([i for i in range(0, 93)])
ILLUMINA_1_0_SET = set([i for i in range(26, 93)])
ILLUMINA_1_3_SET = set([i for i in range(31, 93)])


def main(ngsfiles, outfile=sys.stdout, delimiter="\t", n_samples=1000000):

    writer = csv.DictWriter(
        outfile,
        fieldnames=["File", "Evidence", "ProbableEncoding"],
        delimiter=delimiter,
    )
    writer.writeheader()
    outfile.flush()

    if n_samples < 1:
        n_samples = None

    for ngsfilepath in ngsfiles:
        try:
            ngsfile = NGSFile(ngsfilepath, store_qualities=True)
        except FileNotFoundError:
            result = {
                "File": ngsfilepath,
                "Evidence": "File not found.",
                "ProbableEncoding": "N/A",
            }
            writer.writerow(result)
            outfile.flush()
            continue

        score_set = set()
        for read in itertools.islice(ngsfile, n_samples):
            score_set.update(read["quality"])

        max_phred_score = str(max(score_set) + 33)
        min_phred_score = str(min(score_set) + 33)
        if score_set <= ILLUMINA_1_3_SET:
            result = {
                "File": ngsfilepath,
                "Evidence": "ASCII range: {}-{}".format(
                    min_phred_score, max_phred_score
                ),
                "ProbableEncoding": "Illumina 1.3",
            }
        elif score_set <= ILLUMINA_1_0_SET:
            result = {
                "File": ngsfilepath,
                "Evidence": "ASCII range: {}-{}".format(
                    min_phred_score, max_phred_score
                ),
                "ProbableEncoding": "Solexa/Illumina 1.0",
            }
        elif score_set <= SANGER_SET:
            result = {
                "File": ngsfilepath,
                "Evidence": "ASCII range: {}-{}".format(
                    min_phred_score, max_phred_score
                ),
                "ProbableEncoding": "Sanger/Illumina 1.8",
            }
        else:
            result = {
                "File": ngsfilepath,
                "Evidence": "ASCII values outside known PHRED encoding ranges: {}-{}".format(
                    min_phred_score, max_phred_score
                ),
                "ProbableEncoding": "Unknown",
            }

        writer.writerow(result)
        outfile.flush()
