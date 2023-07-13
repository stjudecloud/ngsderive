import csv
import itertools
import logging
from collections import defaultdict

from .strandedness import get_reads_rg

from ..utils import NGSFile, NGSFileType

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def resolve_endedness(
    firsts, lasts, neither, both, paired_deviance, reads_per_template=None
):
    result = {
        "f+l-": firsts,
        "f-l+": lasts,
        "f-l-": neither,
        "f+l+": both,
    }
    if reads_per_template is not None:
        result["Reads per template"] = reads_per_template

    # only firsts present
    if (firsts > 0) and (lasts == 0 and neither == 0 and both == 0):
        result["Endedness"] = "Unknown"
        return result
    # only lasts present
    if (lasts > 0) and (firsts == 0 and neither == 0 and both == 0):
        result["Endedness"] = "Unknown"
        return result
    # only neither present
    if (neither > 0) and (firsts == 0 and lasts == 0 and both == 0):
        result["Endedness"] = "Unknown"
        return result
    # only both present
    if (both > 0) and (firsts == 0 and lasts == 0 and neither == 0):
        if reads_per_template is None or round(reads_per_template) == 1:
            result["Endedness"] = "Single-End"
        else:
            result["Endedness"] = "Unknown"
        return result
    # first/lasts mixed with neither/both reads
    if (firsts > 0 or lasts > 0) and (neither > 0 or both > 0):
        result["Endedness"] = "Unknown"
        return result
    # any mix of neither and both, regardless of first/lasts
    if neither > 0 and both > 0:
        result["Endedness"] = "Unknown"
        return result
    else:
        assert neither == 0 and both == 0

        read1_frac = firsts / (firsts + lasts)
        if read1_frac > (0.5 - paired_deviance) and read1_frac < (
            0.5 + paired_deviance
        ):
            if reads_per_template is None or round(reads_per_template) == 2:
                result["Endedness"] = "Paired-End"
            else:
                result["Endedness"] = "Unknown"
            return result
        result["Endedness"] = "Unknown"
        return result


def find_reads_per_template(read_names):
    count = 0
    sum = 0
    for c in read_names.values():
        count += 1
        sum += c
    rpt = sum / count
    return rpt


def main(ngsfiles, outfile, n_reads, paired_deviance, lenient, no_rpt, split_by_rg):
    fieldnames = [
        "File",
        "f+l-",
        "f-l+",
        "f-l-",
        "f+l+",
        "Endedness",
    ]
    if split_by_rg:
        fieldnames.insert(1, "Read group")
    if not no_rpt:
        fieldnames.insert(-1, "Reads per template")

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
                "f+l-": "N/A",
                "f-l+": "N/A",
                "f-l-": "N/A",
                "f+l+": "N/A",
                "Endedness": "File not found.",
            }
            if split_by_rg:
                result["Read group"] = "N/A"
            if not no_rpt:
                result["Reads per template"] = "N/A"
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

        ordering_flags = defaultdict(
            lambda: {"firsts": 0, "lasts": 0, "neither": 0, "both": 0}
        )
        read_names = None
        if not no_rpt:
            read_names = defaultdict(lambda: defaultdict(lambda: 0))

        for read in itertools.islice(samfile, n_reads):
            # only count primary alignments and unmapped reads
            if read.is_secondary and not read.is_unmapped:
                continue

            rg = get_reads_rg(read)
            if read_names is not None:
                read_names["overall"][read.query_name] += 1
                read_names[rg][read.query_name] += 1

            if read.is_read1 and not read.is_read2:
                ordering_flags["overall"]["firsts"] += 1
                ordering_flags[rg]["firsts"] += 1
            elif not read.is_read1 and read.is_read2:
                ordering_flags["overall"]["lasts"] += 1
                ordering_flags[rg]["lasts"] += 1
            elif not read.is_read1 and not read.is_read2:
                ordering_flags["overall"]["neither"] += 1
                ordering_flags[rg]["neither"] += 1
            elif read.is_read1 and read.is_read2:
                ordering_flags["overall"]["both"] += 1
                ordering_flags[rg]["both"] += 1
            else:
                raise RuntimeError(
                    "This shouldn't be possible. Please contact the developers."
                )
        assert (
            ordering_flags["overall"]["firsts"]
            + ordering_flags["overall"]["lasts"]
            + ordering_flags["overall"]["neither"]
            + ordering_flags["overall"]["both"]
        ) > 0

        if not split_by_rg:
            if read_names is not None:
                reads_per_template = find_reads_per_template(read_names["overall"])
            else:
                reads_per_template = None
            result = resolve_endedness(
                ordering_flags["overall"]["firsts"],
                ordering_flags["overall"]["lasts"],
                ordering_flags["overall"]["neither"],
                ordering_flags["overall"]["both"],
                paired_deviance,
                reads_per_template=reads_per_template,
            )

            if result["Endedness"] == "Unknown":
                logger.warning("Could not determine endedness!")
                if not lenient:
                    sysexit = 2

            result["File"] = ngsfilepath
            writer.writerow(result)
            outfile.flush()

        else:
            for rg in ordering_flags:
                if (
                    rg == "unknown_read_group"
                    and (
                        ordering_flags[rg]["firsts"]
                        + ordering_flags[rg]["lasts"]
                        + ordering_flags[rg]["neither"]
                        + ordering_flags[rg]["both"]
                    )
                    == 0
                ):
                    continue

                if read_names is not None:
                    reads_per_template = find_reads_per_template(read_names[rg])
                else:
                    reads_per_template = None
                result = resolve_endedness(
                    ordering_flags[rg]["firsts"],
                    ordering_flags[rg]["lasts"],
                    ordering_flags[rg]["neither"],
                    ordering_flags[rg]["both"],
                    paired_deviance,
                    reads_per_template=reads_per_template,
                )

                if result["Endedness"] == "Unknown":
                    logger.warning("Could not determine endedness!")
                    if not lenient:
                        sysexit = 2

                result["File"] = ngsfilepath
                result["Read group"] = rg
                writer.writerow(result)
                outfile.flush()

    if sysexit != 0:
        raise SystemExit(sysexit)
