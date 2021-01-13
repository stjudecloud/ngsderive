import logging
import sys
import csv

from collections import defaultdict

from ..utils import NGSFile, NGSFileType, GFF, JunctionCache, ContigEnd

logger = logging.getLogger("junction-annotation")


def annotate_event(
    found_pos,
    reference_positions,
    fuzzy_range,
):
    if fuzzy_range == 0:
        if found_pos in reference_positions:
            return False, found_pos
        else:
            return True, None
    for pos in reference_positions.irange(
        found_pos - fuzzy_range, found_pos + fuzzy_range
    ):  # return first match even if multiple found
        return False, pos
    return True, None


def annotate_junctions(
    ngsfilepath,
    gff,
    min_intron,
    min_mapq,
    min_reads,
    fuzzy_range,
    disable_junction_files,
):
    try:
        ngsfile = NGSFile(ngsfilepath)
    except FileNotFoundError:
        result = {
            "file": ngsfilepath,
            "total_junctions": "N/A",
            "total_splice_events": "N/A",
            "known_junctions": "N/A",
            "complete_novel_junctions": "N/A",
            "partial_novel_junctions": "N/A",
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

    if not disable_junction_files:
        junction_file = open(f"{ngsfile.basename}.junctions.tsv", "w")
        print(
            "\t".join(
                ["chrom", "intron_start", "intron_end", "read_count", "annotation"]
            ),
            file=junction_file,
        )

    cache = JunctionCache(gff)

    num_known = 0
    num_novel = 0
    num_partial = 0
    num_known_spliced_reads = 0
    num_novel_spliced_reads = 0
    num_partial_spliced_reads = 0

    for contig in samfile.references:

        logger.debug(f"Searching {contig} for splice junctions...")
        found_introns = samfile.find_introns(
            [seg for seg in samfile.fetch(contig) if seg.mapping_quality >= min_mapq]
        )

        events = [
            (intron_start, intron_end, num_reads)
            for (intron_start, intron_end), num_reads in found_introns.items()
            if num_reads >= min_reads and intron_end - intron_start >= min_intron
        ]
        if not events:
            logger.debug(
                f"No valid splice junctions on {contig}. {len(found_introns)} potential junctions discarded."
            )
            continue
        logger.debug(
            f"Found {len(events)} valid splice junctions. {len(found_introns) - len(events)} potential junctions filtered."
        )

        if cache.EOF:
            logger.info(f"{contig} not found in GFF. All events novel.")
            annotation = "complete_novel"
            for intron_start, intron_end, num_reads in events:
                num_novel += 1
                num_novel_spliced_reads += num_reads
            continue  # annotate rest of contigs as novel
        if contig != cache.cur_contig:
            try:
                cache.advance_contigs(contig)
            except StopIteration:
                logger.info(f"{contig} not found in GFF. All events novel.")
                annotation = "complete_novel"
                for intron_start, intron_end, num_reads in events:
                    num_novel += 1
                    num_novel_spliced_reads += num_reads
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
                continue  # annotate rest of contigs as novel

        collapsed_junctions = defaultdict(int)

        for intron_start, intron_end, num_reads in events:
            start_novel, ref_start = annotate_event(
                intron_start, cache.exon_ends, fuzzy_range
            )

            end_novel, ref_end = annotate_event(
                intron_end, cache.exon_starts, fuzzy_range
            )

            if ref_start:
                start = ref_start
            else:
                start = intron_start
            if ref_end:
                end = ref_end
            else:
                end = intron_end

            if fuzzy_range:
                collapsed_junctions[(start, end)] += num_reads
            else:  # only execute if not fuzzy searching
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
                                str(start),
                                str(end),
                                str(num_reads),
                                annotation,
                            ]
                        ),
                        file=junction_file,
                    )

        # if not fuzzy searching, collapsed_junctions is empty and loop is skipped
        for (intron_start, intron_end), num_reads in sorted(
            collapsed_junctions.items()
        ):
            start_novel, _ = annotate_event(intron_start, cache.exon_ends, 0)

            end_novel, _ = annotate_event(intron_end, cache.exon_starts, 0)

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

    junction_file.close()

    result = {
        "file": ngsfilepath,
        "total_junctions": num_known + num_novel + num_partial,
        "total_splice_events": num_known_spliced_reads
        + num_novel_spliced_reads
        + num_partial_spliced_reads,
        "known_junctions": num_known,
        "complete_novel_junctions": num_novel,
        "partial_novel_junctions": num_partial,
        "known_spliced_reads": num_known_spliced_reads,
        "complete_novel_spliced_reads": num_novel_spliced_reads,
        "partial_novel_spliced_reads": num_partial_spliced_reads,
    }
    return result


def main(
    ngsfiles,
    gene_model_file,
    outfile=sys.stdout,
    delimiter="\t",
    min_intron=50,
    min_mapq=30,
    min_reads=1,
    fuzzy_range=0,
):
    logger.info("Arguments:")
    logger.info("  - Gene model file: {}".format(gene_model_file))
    logger.info("  - Minimum intron length: {}".format(min_intron))
    logger.info("  - Minimum MAPQ: {}".format(min_mapq))
    logger.info("  - Minimum reads per junction: {}".format(min_reads))
    logger.info("  - Fuzzy junction range: +-{}".format(fuzzy_range))

    gff = GFF(
        gene_model_file,
        feature_type="exon",
        dataframe_mode=False,
    )
    logger.debug("Opened gene model")

    fieldnames = [
        "file",
        "total_junctions",
        "total_splice_events",
        "known_junctions",
        "complete_novel_junctions",
        "partial_novel_junctions",
        "known_spliced_reads",
        "complete_novel_spliced_reads",
        "partial_novel_spliced_reads",
    ]

    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=delimiter)
    writer.writeheader()
    outfile.flush()

    for ngsfilepath in ngsfiles:
        entry = annotate_junctions(
            ngsfilepath,
            gff,
            min_intron=min_intron,
            min_mapq=min_mapq,
            min_reads=min_reads,
            fuzzy_range=fuzzy_range,
        )
        writer.writerow(entry)
        outfile.flush()
