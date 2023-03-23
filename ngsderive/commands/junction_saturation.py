import logging
import sys
import os
import csv

from collections import defaultdict
from pathlib import Path

from .junction_annotation import annotate_junctions
from ..utils import NGSFile, NGSFileType, GFF, JunctionCache

logger = logging.getLogger("junction-saturation")


def main(
    ngsfiles,
    gene_model_file,
    outfile=sys.stdout,
    sample_start=5,
    sample_step=5,
    sample_end=100,
    min_intron=50,
    min_mapq=30,
    min_reads=2,
    fuzzy_range=0,
):
    logger.info("Arguments:")
    logger.info("  - Gene model file: {}".format(gene_model_file))
    logger.info("  - Sample start percentage: {}".format(sample_start))
    logger.info("  - Sample step percentage: {}".format(sample_step))
    logger.info("  - Sample end percentage : {}".format(sample_end))
    logger.info("  - Minimum intron length: {}".format(min_intron))
    logger.info("  - Minimum MAPQ: {}".format(min_mapq))
    logger.info("  - Minimum reads per junction: {}".format(min_reads))
    logger.info("  - Fuzzy junction range: +-{}".format(fuzzy_range))

    logger.debug("Processing gene model...")
    gff = GFF(
        gene_model_file,
        feature_type="exon",
        dataframe_mode=False,
    )
    cache = JunctionCache(gff)
    logger.debug("Done")

    sample_percents = list(range(sample_start, sample_end, sample_step)) + [sample_end]

    fieldnames = ["File"] + [
        fieldname
        for percent in sample_percents
        for fieldname in (
            str(percent) + "_percent_known_junctions",
            str(percent) + "_percent_partial_novel_junctions",
            str(percent) + "_percent_complete_novel_junctions",
        )
    ]

    print("\t".join(fieldnames), file=outfile, flush=True)

    for ngsfilepath in ngsfiles:
        results = dict()
        for sample_percent in sample_percents:
            logger.info(
                f"Beginning annotate_junctions() with a sample rate of {sample_percent / 100}"
            )
            result = annotate_junctions(
                ngsfilepath,
                cache,
                min_intron=min_intron,
                min_mapq=min_mapq,
                min_reads=min_reads,
                fuzzy_range=fuzzy_range,
                consider_unannotated_references_novel=False,
                junction_dir=None,
                disable_junction_files=True,
                sample_rate=sample_percent / 100,
            )
            results[sample_percent] = result

        entry = [ngsfilepath] + [
            value
            for percent in sample_percents
            for value in (
                str(results[percent]["known_junctions"]),
                str(results[percent]["partial_novel_junctions"]),
                str(results[percent]["complete_novel_junctions"]),
            )
        ]

        print("\t".join(entry), file=outfile, flush=True)
