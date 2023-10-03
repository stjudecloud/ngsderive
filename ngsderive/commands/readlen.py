import csv
import itertools
import logging
from collections import defaultdict

from ..utils import NGSFile

logger = logging.getLogger("readlen")


def main(
    ngsfiles,
    outfile,
    n_reads,
    majority_vote_cutoff,
):
    writer = csv.DictWriter(
        outfile,
        fieldnames=["File", "Evidence", "MajorityPctDetected", "ConsensusReadLength"],
        delimiter="\t",
    )
    writer.writeheader()
    outfile.flush()

    if n_reads < 1:
        n_reads = None

    for ngsfilepath in ngsfiles:
        read_lengths = defaultdict(int)
        try:
            ngsfile = NGSFile(ngsfilepath)
        except FileNotFoundError:
            result = {
                "File": ngsfilepath,
                "Evidence": "Error opening file.",
                "MajorityPctDetected": "N/A",
                "ConsensusReadLength": "N/A",
            }

            writer.writerow(result)
            outfile.flush()
            continue

        # accumulate read lengths
        total_reads_sampled = 0
        for read in itertools.islice(ngsfile, n_reads):
            total_reads_sampled += 1
            read_lengths[len(read["query"])] += 1

        read_length_keys_sorted = sorted(
            [int(k) for k in read_lengths.keys()], reverse=True
        )
        putative_max_readlen = read_length_keys_sorted[0]

        # note that simply picking the read length with the highest amount of evidence
        # doesn't make sense things like adapter trimming might shorten the read length,
        # but the read length should never grow past the maximum value.

        # if not, cannot determine, return -1
        pct = round(read_lengths[putative_max_readlen] / total_reads_sampled * 100, 2)
        logger.info(f"Max read length percentage: {pct}")
        majority_readlen = putative_max_readlen if pct > majority_vote_cutoff else -1

        result = {
            "File": ngsfilepath,
            "Evidence": ";".join(
                [f"{k}={read_lengths[k]}" for k in read_length_keys_sorted]
            ),
            "MajorityPctDetected": str(pct) + "%",
            "ConsensusReadLength": majority_readlen,
        }

        writer.writerow(result)
        outfile.flush()
