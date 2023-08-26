import csv
import itertools
import logging
from collections import defaultdict
from math import isclose
from sys import intern

from ..utils import NGSFile, NGSFileType
from .strandedness import get_reads_rg

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def resolve_endedness(
    firsts, lasts, neither, both, paired_deviance, round_rpt, reads_per_template=None
):
    result = {
        "f+l-": firsts,
        "f-l+": lasts,
        "f-l-": neither,
        "f+l+": both,
    }
    if reads_per_template is not None:
        result["Reads per template"] = reads_per_template
        if round_rpt:
            reads_per_template = round(reads_per_template)

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
        if reads_per_template is None or isclose(reads_per_template, 1.0):
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
        lower_limit = 0.5 - paired_deviance
        upper_limit = 0.5 + paired_deviance
        if isclose(read1_frac, 0.5) or (
            read1_frac < upper_limit and read1_frac > lower_limit
        ):
            if reads_per_template is None or isclose(reads_per_template, 2.0):
                result["Endedness"] = "Paired-End"
            else:
                result["Endedness"] = "Unknown"
            return result
        result["Endedness"] = "Unknown"
        return result


def find_reads_per_template(read_names):
    tot_reads = 0
    tot_templates = 0
    read_group_reads = defaultdict(lambda: 0)
    read_group_templates = defaultdict(lambda: 0)
    for read_name, rg_list in read_names.items():
        num_reads = len(rg_list)
        tot_reads += num_reads
        tot_templates += 1
        rg_set = set(rg_list)
        if len(rg_set) == 1:
            rg = rg_list[0]
            read_group_reads[rg] += num_reads
            read_group_templates[rg] += 1
        else:
            logger.warning(
                f"QNAME {read_name} in multiple read groups: {', '.join(rg_set)}"
            )
            for rg in rg_list:
                read_group_reads[rg] += 1
            for rg in rg_set:
                read_group_templates[rg] += 1

    read_group_rpt = {}
    read_group_rpt["overall"] = tot_reads / tot_templates
    for rg in read_group_reads:
        read_group_rpt[rg] = read_group_reads[rg] / read_group_templates[rg]

    return read_group_rpt


def main(
    ngsfiles,
    outfile,
    n_reads,
    paired_deviance,
    lenient,
    calc_rpt,
    round_rpt,
    split_by_rg,
):
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
    if calc_rpt:
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
            if calc_rpt:
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

        ordering_flags = defaultdict(
            lambda: {"firsts": 0, "lasts": 0, "neither": 0, "both": 0}
        )
        read_names = None
        if calc_rpt:
            read_names = defaultdict(list)

        for read in itertools.islice(samfile, n_reads):
            # only count primary alignments and unmapped reads
            if (read.is_secondary or read.is_supplementary) and not read.is_unmapped:
                continue

            rg = intern(get_reads_rg(read))
            if read_names is not None:
                read_names[read.query_name].append(rg)

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

        rg_rpt = None
        if read_names is not None:
            rg_rpt = find_reads_per_template(read_names)

        if not split_by_rg:
            if rg_rpt is not None:
                reads_per_template = rg_rpt["overall"]
            else:
                reads_per_template = None
            result = resolve_endedness(
                ordering_flags["overall"]["firsts"],
                ordering_flags["overall"]["lasts"],
                ordering_flags["overall"]["neither"],
                ordering_flags["overall"]["both"],
                paired_deviance,
                round_rpt,
                reads_per_template,
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

                if rg_rpt is not None:
                    reads_per_template = rg_rpt[rg]
                else:
                    reads_per_template = None
                result = resolve_endedness(
                    ordering_flags[rg]["firsts"],
                    ordering_flags[rg]["lasts"],
                    ordering_flags[rg]["neither"],
                    ordering_flags[rg]["both"],
                    paired_deviance,
                    round_rpt,
                    reads_per_template,
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
