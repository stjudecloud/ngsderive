import csv
import itertools
import logging

from ..utils import NGSFile

logger = logging.getLogger("encoding")

# sets are created by doing PHRED+33 decoding of each encoding's ASCII ranges
SANGER_SET = set(i for i in range(0, 93))
ILLUMINA_1_0_SET = set(i for i in range(26, 93))
ILLUMINA_1_3_SET = set(i for i in range(31, 93))


def main(ngsfiles, outfile, n_reads):
    writer = csv.DictWriter(
        outfile,
        fieldnames=["File", "Evidence", "ProbableEncoding"],
        delimiter="\t",
    )
    writer.writeheader()
    outfile.flush()

    if n_reads < 1:
        n_reads = None

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
        for read in itertools.islice(ngsfile, n_reads):
            score_set.update(read["quality"])

        highest_ascii = str(max(score_set) + 33)
        lowest_ascii = str(min(score_set) + 33)
        result = {
            "File": ngsfilepath,
            "Evidence": f"ASCII range: {lowest_ascii}-{highest_ascii}",
        }
        if score_set <= ILLUMINA_1_3_SET:
            result["ProbableEncoding"] = "Illumina 1.3"
        elif score_set <= ILLUMINA_1_0_SET:
            result["ProbableEncoding"] = "Solexa/Illumina 1.0"
        elif score_set <= SANGER_SET:
            result["ProbableEncoding"] = "Sanger/Illumina 1.8"
        else:
            # overwrite result["Evidence"] with more info
            result[
                "Evidence"
            ] = f"ASCII values outside known PHRED encoding ranges: {lowest_ascii}-{highest_ascii}"
            result["ProbableEncoding"] = "Unknown"

        writer.writerow(result)
        outfile.flush()
