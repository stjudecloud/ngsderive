import logging
import sys
import os
import csv

from collections import defaultdict
from pathlib import Path

from ..utils import NGSFile, NGSFileType, GFF, JunctionCache

logger = logging.getLogger("junction-annotation")


def annotate_event(
    found_pos,
    reference_positions,
    fuzzy_range,
):
    for pos in reference_positions.irange(
        found_pos - fuzzy_range, found_pos + fuzzy_range
    ):  # return first match even if multiple found
        return False, pos
    return True, None


def annotate_junctions(
    ngsfilepath,
    cache,
    min_intron,
    min_mapq,
    min_reads,
    fuzzy_range,
    consider_unannotated_references_novel,
    junction_dir,
    disable_junction_files,
):
    try:
        ngsfile = NGSFile(ngsfilepath)
    except FileNotFoundError:
        result = {
            "File": ngsfilepath,
            "total_junctions": "N/A",
            "total_splice_events": "N/A",
            "known_junctions": "N/A",
            "partial_novel_junctions": "N/A",
            "complete_novel_junctions": "N/A",
            "known_spliced_reads": "N/A",
            "complete_novel_spliced_reads": "N/A",
            "partial_novel_spliced_reads": "N/A",
        }
        return [result]

    if ngsfile.filetype != NGSFileType.BAM:
        raise RuntimeError(
            "Invalid file: {}. `junction-annotation` only supports aligned BAM files!".format(
                ngsfilepath
            )
        )
    samfile = ngsfile.handle

    junction_file = None
    if not disable_junction_files:
        junction_filename = os.path.join(junction_dir, ngsfile.basename)
        junction_file = open(f"{junction_filename}.junctions.tsv", "w")
        print(
            "\t".join(
                ["chrom", "intron_start", "intron_end", "read_count", "annotation"]
            ),
            file=junction_file,
        )

    num_known = 0
    num_novel = 0
    num_partial = 0
    num_known_spliced_reads = 0
    num_novel_spliced_reads = 0
    num_partial_spliced_reads = 0

    for contig in samfile.references:
        num_too_few_reads = 0
        logger.info(f"Searching {contig} for splice junctions...")
        found_introns = samfile.find_introns(
            [seg for seg in samfile.fetch(contig) if seg.mapping_quality >= min_mapq]
        )

        events = [
            (intron_start, intron_end, num_reads)
            for (intron_start, intron_end), num_reads in found_introns.items()
            if intron_end - intron_start >= min_intron
        ]
        if not events:
            logger.debug(
                f"No valid splice junctions in {contig}. {len(found_introns)} potential junctions too short."
            )
            continue
        logger.debug(
            f"Found {len(events)} potential splice junctions. {len(found_introns) - len(events)} potential junctions too short."
        )

        if contig not in cache.exon_starts:
            logger.info(f"{contig} not found in GFF. All events marked `unannotated_reference`.")
            annotation = "unannotated_reference"
            if consider_unannotated_references_novel:
                logger.info(f"Events being considered novel for summary report.")

            for intron_start, intron_end, num_reads in events:
                if num_reads < min_reads:
                    num_too_few_reads += 1
                    continue
                if consider_unannotated_references_novel:
                    num_novel += 1
                    num_novel_spliced_reads += num_reads

                if junction_file:
                    print(
                        "\t".join(
                            [
                                contig,
                                str(intron_start),
                                str(intron_end),
                                str(num_reads),
                                annotation,
                            ]
                        ),
                        file=junction_file,
                    )

            logger.debug(
                f"{num_too_few_reads} potential junctions didn't have enough read support."
            )
            logger.debug(f"{len(events) - num_too_few_reads} junctions annotated.")
            continue

        collapsed_junctions = defaultdict(int)

        for intron_start, intron_end, num_reads in events:
            start_novel, ref_start = annotate_event(
                intron_start, cache.exon_ends[contig], fuzzy_range
            )

            end_novel, ref_end = annotate_event(
                intron_end, cache.exon_starts[contig], fuzzy_range
            )

            if ref_start:
                start = ref_start
            else:
                start = intron_start
            if ref_end:
                end = ref_end
            else:
                end = intron_end

            # if fuzzy searching, collapse these reads into nearby events
            if fuzzy_range:
                collapsed_junctions[(start, end)] += num_reads
            # if not fuzzy searching, tally reads and write to outfile
            else:
                if num_reads < min_reads:
                    num_too_few_reads += 1
                    continue
                if start_novel and end_novel:
                    annotation = "complete_novel"
                    num_novel += 1
                    num_novel_spliced_reads += num_reads
                elif start_novel or end_novel:
                    annotation = "partial_novel"
                    num_partial += 1
                    num_partial_spliced_reads += num_reads
                else:
                    annotation = "annotated"
                    num_known += 1
                    num_known_spliced_reads += num_reads

                if junction_file:
                    print(
                        "\t".join(
                            [
                                contig,
                                str(start),
                                str(end),
                                str(num_reads),
                                annotation,
                            ]
                        ),
                        file=junction_file,
                    )

        # if not fuzzy searching, collapsed_junctions is empty and loop is skipped,
        # reads will have been tallied and written to outfile already
        for (intron_start, intron_end), num_reads in sorted(
            collapsed_junctions.items()
        ):
            if num_reads < min_reads:
                num_too_few_reads += 1
                continue
            start_novel, _ = annotate_event(intron_start, cache.exon_ends[contig], 0)

            end_novel, _ = annotate_event(intron_end, cache.exon_starts[contig], 0)

            if start_novel and end_novel:
                annotation = "complete_novel"
                num_novel += 1
                num_novel_spliced_reads += num_reads
            elif start_novel or end_novel:
                annotation = "partial_novel"
                num_partial += 1
                num_partial_spliced_reads += num_reads
            else:
                annotation = "annotated"
                num_known += 1
                num_known_spliced_reads += num_reads

            if not disable_junction_files:
                print(
                    "\t".join(
                        [
                            contig,
                            str(intron_start),
                            str(intron_end),
                            str(num_reads),
                            annotation,
                        ]
                    ),
                    file=junction_file,
                )

        logger.debug(
            f"{num_too_few_reads} potential junctions didn't have enough read support."
        )
        if not fuzzy_range:
            logger.debug(f"{len(events) - num_too_few_reads} junctions annotated.")
        else:
            logger.debug(
                f"{len(collapsed_junctions) - num_too_few_reads} junctions annotated."
            )
    if not disable_junction_files:
        junction_file.close()

    result = {
        "File": ngsfilepath,
        "total_junctions": num_known + num_novel + num_partial,
        "total_splice_events": num_known_spliced_reads
        + num_novel_spliced_reads
        + num_partial_spliced_reads,
        "known_junctions": num_known,
        "partial_novel_junctions": num_partial,
        "complete_novel_junctions": num_novel,
        "known_spliced_reads": num_known_spliced_reads,
        "partial_novel_spliced_reads": num_partial_spliced_reads,
        "complete_novel_spliced_reads": num_novel_spliced_reads,
    }
    return result


def main(
    ngsfiles,
    gene_model_file,
    outfile,
    min_intron,
    min_mapq,
    min_reads,
    fuzzy_range,
    consider_unannotated_references_novel,
    junction_dir,
    disable_junction_files,
):
    logger.info("Arguments:")
    logger.info("  - Gene model file: {}".format(gene_model_file))
    logger.info("  - Minimum intron length: {}".format(min_intron))
    logger.info("  - Minimum MAPQ: {}".format(min_mapq))
    logger.info("  - Minimum reads per junction: {}".format(min_reads))
    logger.info("  - Fuzzy junction range: +-{}".format(fuzzy_range))
    logger.info("  - Consider unannotated references novel: {}".format(consider_unannotated_references_novel))
    if not disable_junction_files:
        logger.info("  - Junction file directory: {}".format(junction_dir))
    else:
        logger.info("  - Junction file directory: <disabled>")

    logger.info("Processing gene model...")
    gff = GFF(
        gene_model_file,
        feature_type="exon",
        dataframe_mode=False,
    )
    cache = JunctionCache(gff)
    logger.info("Done")

    junction_dir = Path(junction_dir)
    if not disable_junction_files:
        junction_dir.mkdir(parents=True, exist_ok=True)

    writer = None
    for ngsfilepath in ngsfiles:
        entry = annotate_junctions(
            ngsfilepath,
            cache,
            min_intron=min_intron,
            min_mapq=min_mapq,
            min_reads=min_reads,
            fuzzy_range=fuzzy_range,
            consider_unannotated_references_novel=consider_unannotated_references_novel,
            junction_dir=junction_dir,
            disable_junction_files=disable_junction_files,
        )

        if not writer:
            fieldnames = [
                "File",
                "total_junctions",
                "total_splice_events",
                "known_junctions",
                "partial_novel_junctions",
                "complete_novel_junctions",
                "known_spliced_reads",
                "partial_novel_spliced_reads",
                "complete_novel_spliced_reads",
            ]

            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
        writer.writerow(entry)
        outfile.flush()
