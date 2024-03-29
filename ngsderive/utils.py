import enum
import gzip
import logging
import os
import random
import re
import subprocess
from collections import defaultdict
from operator import itemgetter

import gtfparse
import pysam
import tabix
from sortedcontainers import SortedList

logger = logging.getLogger("utils")


class NGSFileType(enum.Enum):
    FASTQ = 1
    SAM = 2
    BAM = 3


class NGSFile:
    def __init__(self, filename, store_qualities=False):
        self.filename = filename
        self.store_qualities = store_qualities
        self.basename = os.path.basename(self.filename)
        self.ext = ".".join(self.basename.split(".")[1:])
        self.readmode = "r"
        self.gzipped = False
        self.read_num = 0

        if (
            self.ext.endswith(".gz")
            or self.ext.endswith(".bgz")
            or self.ext.endswith(".bam")
        ):
            self.readmode = "rb"
            self.gzipped = True

        if (
            self.ext.endswith("fastq")
            or self.ext.endswith("fq")
            or self.ext.endswith("fastq.gz")
            or self.ext.endswith("fq.gz")
        ):
            self.filetype = NGSFileType.FASTQ
            if self.gzipped:
                self.handle = gzip.open(self.filename, mode=self.readmode)
            else:
                self.handle = open(self.filename, mode=self.readmode, encoding="utf-8")
        elif self.ext.endswith("sam"):
            self.filetype = NGSFileType.SAM
            self.handle = pysam.AlignmentFile(self.filename, self.readmode)
        elif self.ext.endswith("bam"):
            self.filetype = NGSFileType.BAM
            self.handle = pysam.AlignmentFile(self.filename, self.readmode)
        else:
            raise RuntimeError(f"Could not determine NGS file type: {self.filename}")

        logger.debug(f"Opened NGS file '{self.filename}' as {self.filetype}")
        logger.debug(f"Gzipped: {self.gzipped}")
        logger.debug(f"Readmode: {self.readmode}")

    def __iter__(self):
        return self

    def __next__(self):
        query_name = None
        query = None
        read_group = None
        quality = None
        if self.filetype == NGSFileType.FASTQ:
            query_name = self.handle.readline().strip()
            if not query_name or query_name == "":
                raise StopIteration()

            query = self.handle.readline().strip()
            _plusline = self.handle.readline().strip()
            quality_string = self.handle.readline().strip()

            if self.gzipped:
                query_name = query_name.decode("utf-8")
                query = query.decode("utf-8")
                if self.store_qualities:
                    quality_string = quality_string.decode("utf-8")

            quality = []
            if self.store_qualities:
                for char in quality_string:
                    # PHRED+33 decoding
                    ascii_code = ord(char)
                    quality.append(ascii_code - 33)

            if query_name.startswith("@"):
                query_name = query_name[1:]

        elif self.filetype == NGSFileType.SAM or self.filetype == NGSFileType.BAM:
            read = next(self.handle)
            query_name = read.query_name
            query = read.query_alignment_sequence
            read_group = read.get_tag("RG")
            if self.store_qualities:
                quality = read.query_alignment_qualities

        self.read_num += 1

        result = {
            "query_name": query_name,
            "query": query,
        }
        if self.store_qualities:
            result["quality"] = quality
        if read_group:
            result["read_group"] = read_group

        return result


def sort_gff(filename):
    sorted_gff_name_tmp = filename
    ext = sorted_gff_name_tmp.split(".")[-1]
    gzipped = False
    if ext == "gz":
        gzipped = True
        sorted_gff_name_tmp = ".".join(sorted_gff_name_tmp.split(".")[:-1])
        ext = sorted_gff_name_tmp.rsplit(".", maxsplit=1)[-1]  # get ext behind `.gz`
        handle = gzip.open(filename, "r")
    else:
        handle = open(filename, "r", encoding="utf-8")
    sorted_gff_name_tmp = ".".join(
        sorted_gff_name_tmp.split(".")[:-1]
    )  # remove gtf/gff/gff3 ext
    sorted_gff_name = ".".join([sorted_gff_name_tmp, "sorted", ext])
    compressed_gff_name = ".".join([sorted_gff_name, "gz"])

    # Test if this was done in a previous run
    already_created = True
    try:
        # pytabix outputs unwanted text when file doesn't exist
        if os.path.isfile(compressed_gff_name):
            tabixed = tabix.open(compressed_gff_name)
            tabixed.queryi(0, 100, 200)
        else:
            raise FileNotFoundError
    except (FileNotFoundError, tabix.TabixError):
        already_created = False

    if already_created:
        logger.warning(f"Found and reusing {compressed_gff_name}.")
        return compressed_gff_name

    entries = []
    header_lines = []
    header_parsed = False  # only keep comments at start of file
    for line in handle:
        if gzipped:
            line = line.decode("utf-8")
        line = line.strip()

        if not line:
            continue
        if not header_parsed and line[0] == "#":
            header_lines.append(line)
            continue
        if line[0] == "#":
            continue
        header_parsed = True

        [
            seqname,
            source,
            feature,
            start,
            end,
            score,
            strand,
            frame,
            attributes,
        ] = line.split("\t")

        result = [
            seqname,
            source,
            feature,
            int(start),
            int(end),
            score,
            strand,
            frame,
            attributes.replace(';"', '"').replace(
                ";-", "-"
            ),  # correct error in ensemble 78 release
        ]
        entries.append(result)

    entries.sort(key=itemgetter(0, 3, 4))  # seqname, start, end
    new_gff = open(sorted_gff_name, "w", encoding="utf-8")
    for line in header_lines:
        print(line, file=new_gff)
    for entry in entries:
        print("\t".join([str(field) for field in entry]), file=new_gff)
    new_gff.close()
    subprocess.run(["bgzip", "-f", sorted_gff_name], check=True)
    subprocess.run(["tabix", "-p", "gff", compressed_gff_name], check=True)
    logger.warning(f"Created a sorted and tabixed GFF, `{compressed_gff_name}`.")
    return compressed_gff_name


class GFF:
    def __init__(
        self,
        filename,
        dataframe_mode=False,
        store_results=False,
        need_tabix=False,
        feature_type=None,
        gene_exclude_list=None,
        only_protein_coding_genes=False,
    ):
        if not os.path.isfile(filename):
            logger.error(f"Gene model {filename} does not exist!")
            raise SystemExit(1)
        if need_tabix:
            try:
                self.tabix = tabix.open(filename)
                self.tabix.queryi(0, 100, 200)
            except tabix.TabixError:
                logger.warning(f"{filename} not tabixed! Sorting and tabixing new GFF.")
                filename = sort_gff(filename)  # replace provided filename with new file
                self.tabix = tabix.open(filename)
        self.filename = filename
        self.basename = os.path.basename(self.filename)
        self.df = None
        self.entries = None
        self.gene_exclude_list = None
        if gene_exclude_list:
            self.gene_exclude_list = set(
                [
                    item.strip()
                    for item in open(
                        gene_exclude_list, "r", encoding="utf-8"
                    ).readlines()
                ]
            )
        self.feature_type = feature_type

        if dataframe_mode:
            if self.feature_type:
                self.df = gtfparse.read_gtf(filename, features=[self.feature_type])
            else:
                self.df = gtfparse.read_gtf(filename)
            if self.gene_exclude_list:
                if "gene_name" in self.df.columns:
                    self.df = self.df[
                        self.df["gene_name"] not in self.gene_exclude_list
                    ]
                else:
                    logger.warning(
                        "`gene_name` field missing from GFF; could not filter using provided gene_exclude_list."
                    )
            if only_protein_coding_genes:
                if "gene_type" in self.df:  # Gencode
                    self.df = self.df[self.df["gene_type"].str.contains("protein")]
                elif "gene_biotype" in self.df:  # ENSEMBL
                    self.df = self.df[self.df["gene_biotype"].str.contains("protein")]
                else:
                    logger.warning(
                        "Could not isolate protein coding genes. Using all genes."
                    )

        else:
            self.ext = ".".join(self.basename.split(".")[1:])
            self.gzipped = False
            if self.ext.endswith(".gz") or self.ext.endswith(".bgz"):
                self.gzipped = True
                self._handle = gzip.open(filename, "r")
            else:
                self._handle = open(filename, "r", encoding="utf-8")

            self._attr_regexes = [r"(\S+)=(\S+)", r"(\S+) \"(\S+)\""]

            if store_results:
                self.entries = [entry for entry in self]

    def __iter__(self):
        if self.df is not None:
            raise NotImplementedError("iteration in Dataframe mode not implemented")
        return self

    def __next__(self):
        if self.df is not None:
            raise NotImplementedError("iteration in Dataframe mode not implemented")
        while True:
            if self.gzipped:
                line = self._handle.readline().decode("utf-8").strip()
            else:
                line = self._handle.readline().strip()

            if not line:
                raise StopIteration

            if line.startswith("#"):
                continue
            line = line.strip()

            [
                seqname,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                attributes,
            ] = line.split("\t")

            if self.feature_type and feature != self.feature_type:
                continue

            if self.gene_exclude_list:
                selected_bad_gene = False
                for bad_gene in self.gene_exclude_list:
                    if bad_gene in attributes:
                        selected_bad_gene = True
                        break
                if selected_bad_gene:
                    continue

            result = {
                "seqname": seqname,
                "source": source,
                "feature": feature,
                "start": int(start),
                "end": int(end),
                "score": score,
                "strand": strand,
                "frame": frame,
            }

            attributes = attributes.replace(';"', '"').replace(
                ";-", "-"
            )  # correct error in ensemble 78 release
            for attr_raw in [s.strip() for s in attributes.split(";")]:
                if not attr_raw:
                    continue

                for regex in self._attr_regexes:
                    match = re.match(regex, attr_raw)
                    if match:
                        key, value = match.group(1), match.group(2)
                        result[key.strip()] = value.strip()
            return result

    def sample(self):
        if self.df is None and not self.entries:
            raise NotImplementedError("sample() not implemented in iterator mode")
        if self.df is not None:
            return self.df.sample().to_dict("records")[0]
        return random.choice(self.entries)

    def query(self, contig, start, end):
        if self.df is None and not self.entries:
            raise NotImplementedError("query() not implemented in iterator mode")
        if self.df is not None:
            return self.df.query(
                f"seqname=='{contig}' and start=={start} and end=={end}"
            ).to_dict("records")

        raw_hits = self.tabix.query(contig, start, end)
        hits = []
        for hit in raw_hits:
            if self.gene_exclude_list:
                selected_bad_gene = False
                for bad_gene in self.gene_exclude_list:
                    if bad_gene in hit[8]:
                        selected_bad_gene = True
                        break
                if selected_bad_gene:
                    continue

            result = {
                "seqname": hit[0],
                "source": hit[1],
                "feature": hit[2],
                "start": int(hit[3]),
                "end": int(hit[4]),
                "score": hit[5],
                "strand": hit[6],
                "frame": hit[7],
            }

            attributes = (
                hit[8].replace(';"', '"').replace(";-", "-")
            )  # correct error in ensemble 78 release
            for attr_raw in [s.strip() for s in attributes.split(";")]:
                if not attr_raw or attr_raw == "":
                    continue

                for regex in self._attr_regexes:
                    match = re.match(regex, attr_raw)
                    if match:
                        key, value = match.group(1), match.group(2)
                        result[key.strip()] = value.strip()

            hits.append(result)
        return hits


class JunctionCache:
    def __init__(self, gff):
        self.gff = gff
        self.exon_starts = defaultdict(SortedList)
        self.exon_ends = defaultdict(SortedList)
        while True:
            try:
                next(self)
            except StopIteration:
                contigs = ", ".join(self.exon_starts.keys())
                logger.debug(f"Cached {contigs}")
                break

    def __iter__(self):
        return self

    def __next__(self):
        exon = next(self.gff)

        # GFF is 1-based, end inclusive
        # PySam is 0-based, end exclusive
        # starts need to have 1 subtracted
        # ends are already equivelant
        start, end = exon["start"] - 1, exon["end"]
        self.exon_starts[exon["seqname"]].add(start)
        self.exon_ends[exon["seqname"]].add(end)


MOCK_READ_GROUPS = {"overall", "unknown_read_group"}


def get_reads_rg(read, default="unknown_read_group"):
    for k, v in read.tags:
        if k == "RG":
            return v

    return default


def validate_read_group_info(sequence_read_groups, header):
    sequence_read_groups -= MOCK_READ_GROUPS
    header_read_groups = set()
    if "RG" in header:
        header_read_groups = {rg["ID"] for rg in header["RG"]}
    rgs_in_seq_not_in_header = sequence_read_groups - header_read_groups
    for rg in rgs_in_seq_not_in_header:
        logger.warning(
            f"Read group '{rg}' was found in sequences but not the file header!"
        )
    rgs_in_header_not_in_seq = header_read_groups - sequence_read_groups
    for rg in rgs_in_header_not_in_seq:
        logger.warning(
            f"Read group '{rg}' is in the file header but was not found in sampled reads!"
        )

    return rgs_in_header_not_in_seq
