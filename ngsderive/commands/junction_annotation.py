import logging
import sys
import csv

from ..utils import NGSFile, NGSFileType, GFF, JunctionCache, ContigEnd

logger = logging.getLogger("junction-annotation")


def annotate_junctions(
    ngsfilepath,
    gff,
    min_intron=50,
    min_mapq=30,
    min_reads=1,
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

    cache = JunctionCache(gff)

    num_known = 0
    num_novel = 0
    num_partial = 0
    num_known_spliced_reads = 0
    num_novel_spliced_reads = 0
    num_partial_spliced_reads = 0
    num_start_novel = 0
    num_end_novel = 0

    for contig in samfile.references:
        if cache.EOF:
            logger.warning(
                f"Reached end of GTF! Did not find {contig} in GTF. Last contig was {cache.cur_contig}"
            )
            break
        if contig != cache.cur_contig:
            cache.advance_contigs(contig)

        logger.debug(f"Searching {contig} for splice junctions...")
        found_introns = samfile.find_introns(
            [seg for seg in samfile.fetch(contig) if seg.mapping_quality >= min_mapq]
        )

        events = [
            (intron_start - 1, intron_end, num_reads)
            for (intron_start, intron_end), num_reads in found_introns.items()
            if num_reads >= min_reads and intron_end - intron_start > min_intron
        ]
        if not events:
            logger.debug(f"No valid splice junctions on {contig}")
            continue
        logger.debug(
            f"Found {len(events)} valid splice junctions. {len(found_introns) - len(events)} potential junctions filtered."
        )

        for n, (intron_start, intron_end, num_reads) in enumerate(events):
            start_novel = None
            if intron_start in cache.exon_ends:
                start_novel = False
            else:
                start_novel = True
                num_start_novel += 1

            end_novel = None
            if intron_end in cache.exon_starts:
                end_novel = False
            else:
                end_novel = True
                num_end_novel += 1

            if start_novel and end_novel:
                num_novel += 1
                num_novel_spliced_reads += num_reads
            elif start_novel or end_novel:
                num_partial += 1
                num_partial_spliced_reads += num_reads
            else:
                num_known += 1
                num_known_spliced_reads += num_reads

    logger.info(f"starts novel={num_start_novel}")
    logger.info(f"ends novel={num_end_novel}")
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
    return [result]


def main(
    ngsfiles,
    gene_model_file,
    outfile=sys.stdout,
    delimiter="\t",
    min_intron=50,
    min_mapq=30,
    min_reads=1,
):
    logger.info("Arguments:")
    logger.info("  - Gene model file: {}".format(gene_model_file))
    logger.info("  - Minimum intron length: {}".format(min_intron))
    logger.info("  - Minimum MAPQ: {}".format(min_mapq))
    logger.info("  - Minimum reads per junction: {}".format(min_reads))

    logger.info("Opening gene model...")
    gff = GFF(
        gene_model_file,
        feature_type="exon",
        dataframe_mode=False,
    )

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
        entries = annotate_junctions(
            ngsfilepath,
            gff,
            min_intron=min_intron,
            min_mapq=min_mapq,
            min_reads=min_reads,
        )

        for entry in entries:
            writer.writerow(entry)
            outfile.flush()
