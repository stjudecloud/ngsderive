import gzip

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

      try:
        [
            seqname, source, feature, start, end, score, strand, frame,
            attribute
        ] = line.split("\t")
      except:
        print(line)

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

      for attr_raw in attribute.split(";"):
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
