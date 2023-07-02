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


def resolve_flag_count(read1s, read2s, neither, both):
    # only read1s present
    if (read1s > 0) and (read2s == 0 and neither == 0 and both == 0):
        return {
            "Read1s": read1s,
            "Read2s": read2s,
            "Neither": neither,
            "Both": both,
            "Mate state": "Legal",
            "Endedness": "Single-End",
        }
    # only read2s present
    if (read2s > 0) and (read1s == 0 and neither == 0 and both == 0):
        return {
            "Read1s": read1s,
            "Read2s": read2s,
            "Neither": neither,
            "Both": both,
            "Mate state": "Illegal",
            "Endedness": "Single-End",
        }
    # only neither or only both present
    if (
        ((neither > 0) and (not both > 0)) or ((not neither > 0) and (both > 0))
    ) and (read1s == 0 and read2s == 0):
        return {
            "Read1s": read1s,
            "Read2s": read2s,
            "Neither": neither,
            "Both": both,
            "Mate state": "Illegal",
            "Endedness": "Single-End",
        }
    # legal reads mixed with illegal reads
    if (read1s > 0 or read2s > 0) and (neither > 0 or both > 0):
        return {
            "Read1s": read1s,
            "Read2s": read2s,
            "Neither": neither,
            "Both": both,
            "Mate state": "Illegal",
            "Endedness": "Inconclusive",
        }
    # any mix of neither and both, regardless of read1/2s
    if neither > 0 and both > 0:
        return {
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
        if read1_frac > 0.45 and read1_frac < 0.55:
            return {
                "Read1s": read1s,
                "Read2s": read2s,
                "Neither": neither,
                "Both": both,
                "Mate state": "Legal",
                "Endedness": "Paired-End",
            }
        return {
            "Read1s": read1s,
            "Read2s": read2s,
            "Neither": neither,
            "Both": both,
            "Mate state": "Legal",
            "Endedness": "Inconclusive",
        }




def main(ngsfiles, outfile, n_reads, lenient, split_by_rg):
    if not split_by_rg:
        fieldnames = [
            "File",
            "Read1s",
            "Read2s",
            "Neither",
            "Both",
            "Mate state",
            "Endedness",
        ]
    else:
        fieldnames = [
            "File",
            "Read group",
            "Read1s",
            "Read2s",
            "Neither",
            "Both",
            "Mate state",
            "Endedness",
        ]
    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t",
    )
    writer.writeheader()
    outfile.flush()

    if n_reads < 1:
        n_reads = None

    sysexit = 0
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
            if split_by_rg:
                result["Read group"] = "N/A"
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

        mate_flags = defaultdict(lambda: {"read1s": 0, "read2s": 0, "neither": 0, "both": 0})

        for read in itertools.islice(samfile, n_reads):
            # only count primary alignments
            if read.is_secondary:  # true if not primary alignment
                continue

            rg = get_reads_rg(read)

            if read.is_read1 and not read.is_read2:
                mate_flags["overall"]["read1s"] += 1
                mate_flags[rg]["read1s"] += 1
            elif not read.is_read1 and read.is_read2:
                mate_flags["overall"]["read2s"] += 1
                mate_flags[rg]["read2s"] += 1
            elif not read.is_read1 and not read.is_read2:
                mate_flags["overall"]["neither"] += 1
                mate_flags[rg]["neither"] += 1
            elif read.is_read1 and read.is_read2:
                mate_flags["overall"]["both"] += 1
                mate_flags[rg]["both"] += 1
            else:
                raise RuntimeError(
                    "This shouldn't be possible. Please contact the developers."
                )
        assert (mate_flags["overall"]["read1s"]
            + mate_flags["overall"]["read2s"]
            + mate_flags["overall"]["neither"]
            + mate_flags["overall"]["both"]
        ) > 0

        if not split_by_rg:
            result = resolve_flag_count(
                mate_flags["overall"]["read1s"],
                mate_flags["overall"]["read2s"],
                mate_flags["overall"]["neither"],
                mate_flags["overall"]["both"]
            )
            result["File"] = ngsfilepath
            writer.writerow(result)
            outfile.flush()

            if result["Mate state"] == "Illegal":
                logger.warning("Illegal mate state detected!")
                if not lenient:
                    sysexit = 2
            if result["Endedness"] == "Inconclusive":
                logger.warning("Could not determine endedness!")
                if not lenient and sysexit == 0:
                    sysexit = 3

        else:
            for rg in mate_flags:
                if rg == "unknown_read_group":
                    if (mate_flags[rg]["read1s"] + mate_flags[rg]["read2s"] + mate_flags[rg]["neither"] + mate_flags[rg]["both"]) == 0:
                        continue
                result = resolve_flag_count(
                    mate_flags[rg]["read1s"],
                    mate_flags[rg]["read2s"],
                    mate_flags[rg]["neither"],
                    mate_flags[rg]["both"]
                )
                result["File"] = ngsfilepath
                result["Read group"] = rg
                writer.writerow(result)
                outfile.flush()

                if result["Mate state"] == "Illegal":
                    logger.warning("Illegal mate state detected!")
                    if not lenient:
                        sysexit = 2
                if result["Endedness"] == "Inconclusive":
                    logger.warning("Could not determine endedness!")
                    if not lenient and sysexit == 0:
                        sysexit = 3
    if sysexit != 0:
        raise SystemExit(sysexit)

