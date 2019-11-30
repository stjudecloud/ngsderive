import enum
import gzip
import logging
import os
import pysam

logger = logging.getLogger('utils')

class NGSFileType(enum.Enum):
  FASTQ = 1
  SAM = 2
  BAM = 3

class NGSFile:
  def __init__(self, filename):
    self.filename = filename
    self.basename = os.path.basename(self.filename)
    self.ext = ".".join(self.basename.split(".")[1:])
    self.readmode = "r"
    self.gzipped = False

    if self.ext.endswith(".gz") or self.ext.endswith(".bgz") or self.ext.endswith(".bam"):
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
      raise RuntimeError("Could not determine NGS file type: {}".format(self.filename))

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
      quality = self.handle.readline().strip()

      if self.gzipped:
        query_name = query_name.decode("utf-8")
        query = query.decode("utf-8")

      if query_name.startswith("@"):
        query_name = query_name[1:]

    elif self.filetype == NGSFileType.SAM or self.filetype == NGSFileType.BAM:
      read = next(self.handle)
      query_name = read.query_name
      query = read.query

    self.read_num += 1

    return {
      "query_name": query_name,
      "query": query
    }

class GFF:
  def __init__(self, filename, feature_filter=None):
    self._handle = gzip.open(filename, "r")
    self.feature_filter = feature_filter
    self.features = []

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

      if self.feature_filter and not feature in self.feature_filter:
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

        [key, value] = attr_raw.split("=")
        result["attr_" + key] = value.strip()

      self.features.append(result)

  def __iter__(self):
    self.i = 0
    return self

  def __next__(self):
    if self.i < len(self.features):
      self.i += 1
      return self.features[self.i - 1]
    else:
      raise StopIteration
