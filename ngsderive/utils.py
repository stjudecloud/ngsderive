import enum
import gzip
import logging
import os
import re
import pysam
from gtfparse import read_gtf
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
                self.handle = open(self.filename, mode=self.readmode)
        elif self.ext.endswith("sam"):
            self.filetype = NGSFileType.SAM
            self.handle = pysam.AlignmentFile(self.filename, self.readmode)
        elif self.ext.endswith("bam"):
            self.filetype = NGSFileType.BAM
            self.handle = pysam.AlignmentFile(self.filename, self.readmode)
        else:
            raise RuntimeError(
                "Could not determine NGS file type: {}".format(self.filename)
            )

        logger.debug("Opened NGS file '{}' as {}".format(self.filename, self.filetype))
        logger.debug("Gzipped: {}".format(self.gzipped))
        logger.debug("Readmode: {}".format(self.readmode))

    def __iter__(self):
        self.read_num = 0
        return self

    def __next__(self):
        if self.filetype == NGSFileType.FASTQ:
            query_name = self.handle.readline().strip()
            if not query_name or query_name == "":
                raise StopIteration()

            query = self.handle.readline().strip()
            plusline = self.handle.readline().strip()
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
            if self.store_qualities:
                quality = read.query_alignment_qualities

        self.read_num += 1

        if self.store_qualities:
            return {"query_name": query_name, "query": query, "quality": quality}
        else:
            return {"query_name": query_name, "query": query}


class GFF:
    def __init__(
        self,
        filename,
        dataframe_mode=False,
        feature_type=None,
        filters=None,
        gene_blacklist=None,
        only_protein_coding_genes=False,
    ):
        self.df = None
        self.gene_blacklist = None
        if gene_blacklist:
            self.gene_blacklist = set(
                [item.strip() for item in open(gene_blacklist, "r").readlines()]
            )
        self.feature_type = feature_type
        self.filters = filters  # TODO not sure what these are meant to be

        if dataframe_mode:
            self.df = read_gtf(filename)
            if self.feature_type:
                self.df = self.df[self.df["feature"] == self.feature_type]
            if self.gene_blacklist:
                self.df = self.df[self.df["gene_name"] not in self.gene_blacklist]
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
            self.filename = filename
            self.basename = os.path.basename(self.filename)
            self.ext = ".".join(self.basename.split(".")[1:])
            if self.ext.endswith(".gz") or self.ext.endswith(".bgz"):
                self._handle = gzip.open(filename, "r")
            else:
                self._handle = open(filename, "r")

            self._attr_regexes = [r"(\S+)=(\S+)", r"(\S+) \"(\S+)\""]

    def sample(self):
        if self.df is None:
            raise NotImplementedError("sample() not implemented in iterator mode")
        return self.df.sample().to_dict("records")[0]

    def query(self, chr, start, end):
        if self.df is None:
            raise NotImplementedError("query() not implemented in iterator mode")
        return self.df.query(
            f"seqname=='{chr}' and start=={start} and end=={end}"
        ).to_dict("records")

    def filter(self, func):
        # TODO unclear what this is for
        self.df = filter(func, self.df)

    def __iter__(self):
        if self.df:
            raise NotImplementedError("iteration in Dataframe mode not implemented")
        return self

    def __next__(self):
        if self.df:
            raise NotImplementedError("iteration in Dataframe mode not implemented")
        while True:
            line = self._handle.readline().decode("utf-8")
            if not line:
                raise StopIteration

            if line.startswith("#"):
                continue

            [
                seqname,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                attribute,
            ] = line.split("\t")

            if self.feature_type and feature != self.feature_type:
                continue

            if self.gene_blacklist:
                for bad_gene in self.gene_blacklist:
                    if bad_gene in attribute:
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

            attribute = attribute.replace(';"', '"').replace(
                ";-", "-"
            )  # correct error in ensemble 78 release
            for attr_raw in [s.strip() for s in attribute.split(";")]:
                if not attr_raw or attr_raw == "":
                    continue

                for regex in self._attr_regexes:
                    match = re.match(regex, attr_raw)
                    if match:
                        key, value = match.group(1), match.group(2)
                        result[key.strip()] = value.strip()

            # TODO not sure what this is supposed to be used for
            entry_passes_filters = True
            if self.filters:
                for (k, v) in self.filters.items():
                    if k not in result or result[k] != v:
                        entry_passes_filters = False
                        break

            if entry_passes_filters:
                return result

    def next(self):
        return self.__next__()


class ContigEnd(Exception):
    pass


class Junction_cache:
    # def __init__(self, gff, cache_size=10000000):
    def __init__(self, gff):
        self.gff = gff
        # self.max_junctions = cache_size

        exon = next(self.gff)
        self.cur_contig = exon["seqname"]
        logger.debug(f"Caching {self.cur_contig}...")

        tmp_junctions = []
        while True:
            tmp_junctions += (exon["start"] - 1, exon["end"] - 1)
            exon = gff.next()
            if self.cur_contig != exon["seqname"]:
                break
        self.junctions = SortedList(tmp_junctions)
        logger.debug("Done")
        self.last_exon = exon

    def __iter__(self):
        return self

    def __next__(self):
        self.junctions.pop(0)
        self.junctions.pop(1)
        exon = next(self.gff)
        if self.cur_contig != exon["seqname"]:
            self.last_exon = exon
            raise ContigEnd

        start, end = exon["start"] - 1, exon["end"] - 1
        self.junctions.update((start, end))
        return (start, end)

    def next(self):
        return self.__next__()

    def advance_contigs(self, contig=None):
        self.junctions.clear()
        exon = self.last_exon
        # if exon["seqname"] == self.cur_contig:
        #     logger.debug(f"Skipping rest of {self.cur_contig} in GTF")
        # while exon["seqname"] == self.cur_contig:  # pass through rest of cur_contig
        #     try:
        #         exon = next(self.gff)
        #     except ContigEnd:
        #         break
        if contig:
            while exon["seqname"] != contig:
                try:
                    exon = next(self.gff)
                except ContigEnd:
                    logger.warning(f"Skipped {self.cur_contig} searching for {contig}")
        self.cur_contig = exon["seqname"]
        logger.debug(f"Caching {self.cur_contig}...")

        # for _ in range(self.max_junctions):
        tmp_junctions = []
        while True:
            tmp_junctions += (exon["start"] - 1, exon["end"] - 1)
            try:
                exon = next(self.gff)
            except StopIteration or ContigEnd:
                if self.cur_contig != exon["seqname"]:
                    self.last_exon = exon
                break
        logger.debug("Done")
        self.junctions.update(tmp_junctions)
