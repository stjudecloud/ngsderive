import csv
import itertools

import logging
from collections import defaultdict

from ..utils import NGSFile, NGSFileType

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_reads_rg(read, default="unknown_read_group"):
    for k, v in read.tags:
        if k == "RG":
            return v

    return default


def main(ngsfiles, outfile, n_reads, lenient, split_by_rg):
    writer = csv.DictWriter(
        outfile,
        fieldnames=[
            "File",
            "Read1s",
            "Read2s",
            "Neither",
            "Both",
            "Mate state",
            "Endedness",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    outfile.flush()

    if n_reads < 1:
        n_reads = None

    for ngsfilepath in ngsfiles:
        try:
            ngsfile = NGSFile(ngsfilepath)
        except FileNotFoundError:
            result = {
                "File": ngsfilepath,
                "Read1s": "N/A",
                "Read2s": "N/A",
                "Neither": "N/A",
                "Both": "N/A",
                "Mate state": "N/A",
                "Endedness": "File not found.",
            }
            writer.writerow(result)
            outfile.flush()
            continue

        if ngsfile.filetype != NGSFileType.BAM and ngsfile.filetype != NGSFileType.SAM:
            raise RuntimeError(
                "Invalid file: {}. `endedness` only supports SAM/BAM files!".format(
                    ngsfilepath
                )
            )
        samfile = ngsfile.handle

        read_groups = ["unknown_read_group"]
        if "RG" in samfile.header:
            read_groups += [rg["ID"] for rg in samfile.header["RG"]]

        read1s = 0
        read2s = 0
        neither = 0
        both = 0
        for read in itertools.islice(samfile, n_reads):
            if read.is_read1 and not read.is_read2:
                read1s += 1
            elif not read.is_read1 and read.is_read2:
                read2s += 1
            elif not read.is_read1 and not read.is_read2:
                neither += 1
            elif read.is_read1 and read.is_read2:
                both += 1
            else:
                raise RuntimeError(
                    "This shouldn't be possible. Please contact the developers."
                )
        assert (read1s + read2s + neither + both) > 0
        # only read1s present
        if (read1s > 0) and (read2s == 0 and neither == 0 and both == 0):
            result = {
                "File": ngsfilepath,
                "Read1s": read1s,
                "Read2s": read2s,
                "Neither": neither,
                "Both": both,
                "Mate state": "Legal",
                "Endedness": "Single-End",
            }
        # only read2s present
        elif (read2s > 0) and (read1s == 0 and neither == 0 and both == 0):
            result = {
                "File": ngsfilepath,
                "Read1s": read1s,
                "Read2s": read2s,
                "Neither": neither,
                "Both": both,
                "Mate state": "Illegal",
                "Endedness": "Single-End",
            }
        # only neither or only both present
        elif (
            ((neither > 0) and (not both > 0)) or ((not neither > 0) and (both > 0))
        ) and (read1s == 0 and read2s == 0):
            result = {
                "File": ngsfilepath,
                "Read1s": read1s,
                "Read2s": read2s,
                "Neither": neither,
                "Both": both,
                "Mate state": "Illegal",
                "Endedness": "Single-End",
            }
        # legal reads mixed with illegal reads
        elif (read1s > 0 or read2s > 0) and (neither > 0 or both > 0):
            result = {
                "File": ngsfilepath,
                "Read1s": read1s,
                "Read2s": read2s,
                "Neither": neither,
                "Both": both,
                "Mate state": "Illegal",
                "Endedness": "Inconclusive",
            }
        # any mix of neither and both, regardless of read1/2s
        elif neither > 0 and both > 0:
            result = {
                "File": ngsfilepath,
                "Read1s": read1s,
                "Read2s": read2s,
                "Neither": neither,
                "Both": both,
                "Mate state": "Illegal",
                "Endedness": "Inconclusive",
            }
        else:
            assert neither == 0 and both == 0

            read1_frac = read1s / (read1s + read2s)
            if read1_frac > 0.4 and read1_frac < 0.6:
                result = {
                    "File": ngsfilepath,
                    "Read1s": read1s,
                    "Read2s": read2s,
                    "Neither": neither,
                    "Both": both,
                    "Mate state": "Legal",
                    "Endedness": "Paired-End",
                }
            else:
                result = {
                    "File": ngsfilepath,
                    "Read1s": read1s,
                    "Read2s": read2s,
                    "Neither": neither,
                    "Both": both,
                    "Mate state": "Legal",
                    "Endedness": "Inconclusive",
                }

        writer.writerow(result)
        outfile.flush()

        if result["Mate state"] == "Illegal":
            logger.warning("Illegal mate state detected!")
            raise SystemExit(2)
        if result["Endedness"] == "Inconclusive":
            logger.warning("Could not determine endedness!")
            raise SystemError(3)
