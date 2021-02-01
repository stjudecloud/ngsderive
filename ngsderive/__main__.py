#!/usr/bin/env python3

import argparse
import csv
import logging
import sys

from ngsderive import utils
from ngsderive.commands import (
    readlen,
    instrument,
    strandedness,
    encoding,
    junction_annotation,
)

logger = logging.getLogger("ngsderive")


def get_args():
    class SaneFormatter(
        argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
    ):
        pass

    parser = argparse.ArgumentParser(
        description="Derive useful information (or best guess) from next-generation sequencing files.",
        formatter_class=SaneFormatter,
    )
    subparsers = parser.add_subparsers(dest="subcommand")

    common = argparse.ArgumentParser(add_help=False, formatter_class=SaneFormatter)
    common.add_argument(
        "ngsfiles",
        type=str,
        nargs="+",
        help="Next-generation sequencing files to process (BAM or FASTQ).",
    )
    common.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="Write to filename rather than standard out.",
        default="stdout",
    )
    common.add_argument(
        "--delimiter", default="<tab>", help="Delimiter for the outfile."
    )
    common.add_argument(
        "--debug", default=False, action="store_true", help="Enable DEBUG log level."
    )
    common.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="Enable INFO log level.",
    )
    common.add_argument("--version", action="version", version="%(prog)s 1.0.0")

    readlen = subparsers.add_parser(
        "readlen", parents=[common], formatter_class=SaneFormatter
    )
    readlen.add_argument(
        "-c",
        "--majority-vote-cutoff",
        type=float,
        help="To call a majority readlen, the maximum read length must have at least `majority-vote-cutoff`%% reads in support.",
        default=0.7,
    )
    readlen.add_argument(
        "-n",
        "--n-samples",
        type=int,
        help="How many reads to sample. Any n < 1 to parse whole file.",
        default=1000000,
    )

    instrument = subparsers.add_parser(
        "instrument", parents=[common], formatter_class=SaneFormatter
    )
    instrument.add_argument(
        "-n",
        "--n-samples",
        type=int,
        help="How many reads to sample. Any n < 1 to parse whole file.",
        default=10000,
    )

    strandedness = subparsers.add_parser(
        "strandedness", parents=[common], formatter_class=SaneFormatter
    )
    strandedness.add_argument(
        "-g", "--gene-model", help="Gene model as a GFF/GTF file.", required=True
    )
    strandedness.add_argument(
        "--max-tries",
        type=int,
        default=3,
        help="When inconclusive, the test will repeat until this many tries have been reached.",
    )
    strandedness.add_argument(
        "--max-iterations-per-try",
        type=int,
        default=1000,
        help="At most, search this many times for genes that satisfy our search criteria.",
    )
    strandedness.add_argument(
        "-m",
        "--minimum-reads-per-gene",
        type=int,
        help="Filter any genes that don't have at least `m` reads.",
        default=10,
    )
    strandedness.add_argument(
        "-n", "--n-genes", type=int, help="How many genes to sample.", default=100
    )
    strandedness.add_argument(
        "-q",
        "--min-mapq",
        type=int,
        help="Minimum MAPQ to consider for reads.",
        default=30,
    )
    protein_coding_parser = strandedness.add_mutually_exclusive_group(required=False)
    protein_coding_parser.add_argument(
        "--only-protein-coding-genes",
        dest="only_protein_coding_genes",
        action="store_true",
        help="Only consider protein coding genes",
    )
    protein_coding_parser.add_argument(
        "--no-only-protein-coding-genes",
        dest="only_protein_coding_genes",
        action="store_false",
    )
    split_by_rg_parser = strandedness.add_mutually_exclusive_group(required=False)
    split_by_rg_parser.add_argument(
        "--split-by-rg",
        dest="split_by_rg",
        action="store_true",
        help="Contain one entry per read group.",
    )
    split_by_rg_parser.add_argument(
        "--no-split-by-rg", dest="split_by_rg", action="store_false"
    )
    strandedness.set_defaults(only_protein_coding_genes=True, split_by_rg=False)

    encoding = subparsers.add_parser(
        "encoding", parents=[common], formatter_class=SaneFormatter
    )
    encoding.add_argument(
        "-n",
        "--n-samples",
        type=int,
        help="How many reads to sample. Any n < 1 to parse whole file.",
        default=1000000,
    )

    junction_annotation = subparsers.add_parser(
        "junction-annotation", parents=[common], formatter_class=SaneFormatter
    )
    junction_annotation.add_argument(
        "-g", "--gene-model", help="Gene model as a GFF/GTF file.", required=True
    )
    junction_annotation.add_argument(
        "-j",
        "--junction-files-dir",
        help="Directory to write annotated junction files to.",
        default="./",
    )
    junction_annotation.add_argument(
        "-d",
        "--disable-junction-files",
        help="Disable generating junction files in current working directory.",
        action="store_true",
    )
    junction_annotation.add_argument(
        "-i",
        "--min-intron",
        type=int,
        help="Minimum size of intron to be considered a splice.",
        default=50,
    )
    junction_annotation.add_argument(
        "-q",
        "--min-mapq",
        type=int,
        help="Minimum MAPQ to consider for reads.",
        default=30,
    )
    junction_annotation.add_argument(
        "-m",
        "--min-reads",
        type=int,
        help="Filter any junctions that don't have at least `m` reads.",
        default=2,
    )
    junction_annotation.add_argument(
        "-k",
        "--fuzzy-junction-match-range",
        type=int,
        help="Consider found splices within `+-k` bases of a known splice event annotated.",
        default=0,
    )

    args = parser.parse_args()
    if not args.subcommand:
        parser.print_help()
        sys.exit(1)

    return args


def setup_logging(log_level=logging.INFO):
    """Set up the logging.

    Forked from MIT code here: https://github.com/MisterWil/abodepy.
    """
    logging.basicConfig(level=log_level)
    fmt = "%(asctime)s %(levelname)s (%(threadName)s) " "[%(name)s] %(message)s"
    colorfmt = "%(log_color)s{}%(reset)s".format(fmt)
    datefmt = "%H:%M:%S"

    try:
        from colorlog import ColoredFormatter

        logging.getLogger().handlers[0].setFormatter(
            ColoredFormatter(
                colorfmt,
                datefmt=datefmt,
                reset=True,
                log_colors={
                    "DEBUG": "cyan",
                    "INFO": "green",
                    "WARNING": "yellow",
                    "ERROR": "red",
                    "CRITICAL": "red",
                },
            )
        )
    except ImportError:
        pass

    logger = logging.getLogger("")
    logger.setLevel(log_level)


def process_args(args):
    # setup logging
    log_level = logging.WARN

    if args.verbose:
        log_level = logging.INFO

    if args.debug:
        log_level = logging.DEBUG

    setup_logging(log_level)

    # set output file
    if args.outfile == "stdout":
        args.outfile = sys.stdout
    else:
        args.outfile = open(args.outfile, "w")

    # set delimiter
    if args.delimiter == "<tab>":
        args.delimiter = "\t"


def run():
    args = get_args()
    process_args(args)

    if args.subcommand == "readlen":
        readlen.main(
            args.ngsfiles,
            outfile=args.outfile,
            delimiter=args.delimiter,
            n_samples=args.n_samples,
            majority_vote_cutoff=args.majority_vote_cutoff,
        )
    if args.subcommand == "instrument":
        instrument.main(
            args.ngsfiles,
            outfile=args.outfile,
            delimiter=args.delimiter,
            n_samples=args.n_samples,
        )
    if args.subcommand == "strandedness":
        strandedness.main(
            args.ngsfiles,
            args.gene_model,
            outfile=args.outfile,
            delimiter=args.delimiter,
            n_genes=args.n_genes,
            minimum_reads_per_gene=args.minimum_reads_per_gene,
            only_protein_coding_genes=args.only_protein_coding_genes,
            min_mapq=args.min_mapq,
            split_by_rg=args.split_by_rg,
            max_iterations_per_try=args.max_iterations_per_try,
        )
    if args.subcommand == "encoding":
        encoding.main(
            args.ngsfiles,
            outfile=args.outfile,
            delimiter=args.delimiter,
            n_samples=args.n_samples,
        )
    if args.subcommand == "junction-annotation":
        junction_annotation.main(
            args.ngsfiles,
            args.gene_model,
            outfile=args.outfile,
            delimiter=args.delimiter,
            min_intron=args.min_intron,
            min_mapq=args.min_mapq,
            min_reads=args.min_reads,
            fuzzy_range=args.fuzzy_junction_match_range,
            junction_dir=args.junction_files_dir,
            disable_junction_files=args.disable_junction_files,
        )
