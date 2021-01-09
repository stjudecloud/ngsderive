import enum
import gzip
import logging
import os
import re
import pysam

logger = logging.getLogger('utils')


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

        if self.ext.endswith(".gz") or self.ext.endswith(
                ".bgz") or self.ext.endswith(".bam"):
            self.readmode = "rb"
            self.gzipped = True

        if self.ext.endswith("fastq") or \
           self.ext.endswith("fq") or \
           self.ext.endswith("fastq.gz") or \
           self.ext.endswith("fq.gz"):
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
            raise RuntimeError("Could not determine NGS file type: {}".format(
                self.filename))

        logger.debug("Opened NGS file '{}' as {}".format(
            self.filename, self.filetype))
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
            quality = self.handle.readline().strip()

            if self.gzipped:
                query_name = query_name.decode("utf-8")
                query = query.decode("utf-8")
                if self.store_qualities:
                  quality = quality.decode("utf-8")

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
    def __init__(self,
                 filename,
                 feature_type=None,
                 filters=None,
                 gene_blacklist=None):
        self.gene_blacklist = None
        if gene_blacklist:
            self.gene_blacklist = set([
                item.strip() for item in open(gene_blacklist, 'r').readlines()
            ])
        self._handle = gzip.open(filename, "r")
        self.feature_type = feature_type
        self.filters = filters
        self.entries = []
        self.attr_regexes = [r"(\S+)=(\S+)", r"(\S+) \"(\S+)\""]

        while True:
            line = self._handle.readline().decode("utf-8")
            if not line:
                break

            if line.startswith("#"):
                continue

            [
                seqname, source, feature, start, end, score, strand, frame,
                attribute
            ] = line.split("\t")

            if self.feature_type and feature != self.feature_type:
                continue

            entry_passes_gene_blacklist = True
            if self.gene_blacklist:
                for bad_gene in self.gene_blacklist:
                    if bad_gene in attribute:
                        entry_passes_gene_blacklist = False
                        break

            if not entry_passes_gene_blacklist:
                continue

            result = {
                "seqname": seqname,
                "source": source,
                "feature": feature,
                "start": int(start),
                "end": int(end),
                "score": score,
                "strand": strand,
                "frame": frame
            }

            for attr_raw in [s.strip() for s in attribute.split(";")]:
                if not attr_raw or attr_raw == "":
                    continue

                for regex in self.attr_regexes:
                    match = re.match(regex, attr_raw)
                    if match:
                        key, value = match.group(1), match.group(2)
                        result["attr_" + key] = value.strip()

            entry_passes_filters = True
            if entry_passes_filters and self.filters:
                for (k, v) in self.filters.items():
                    if k not in result or result[k] != v:
                        entry_passes_filters = False
                        break

            if entry_passes_filters:
                self.entries.append(result)

    def filter(self, func):
        self.entries = filter(func, self.entries)

    def __iter__(self):
        self.i = 0
        return self

    def __next__(self):
        if self.i < len(self.entries):
            self.i += 1
            return self.entries[self.i - 1]
        else:
            raise StopIteration
