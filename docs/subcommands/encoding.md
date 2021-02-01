# encoding

The `encoding` subcommand detects the likely PHRED score encoding scheme used for FASTQ and BAM files. PHRED scores are encoded as ASCII characters, but which ASCII characters encode for which relative quality scores depends on the build and generation of the sequencing machine used. The BAM specification calls for `PHRED+33` (AKA Illumina1.8+ or Sanger encoding), however that is also the most permissive encoding. It is possible for a stricter encoding to be mis-translated as `PHRED+33` and appear to be of higher quality than it is in truth. This subcommand can be used to check if base scores are suspiciously high throughout a BAM, and a mis-translation may have occured.

`ngsderive` obtained details of the different encoding schemes [here][https://en.wikipedia.org/wiki/FASTQ_format#Encoding].

## Limitations

* All 3 possible encodings have the same upper range, and only differ in terms of lower quality scores. Thus if the provided read data is all of high quality, it may be classified as a stricter encoding than was originally used to generate the data.

    If the downstream use of this value is needed to ensure all PHRED scores are within the derived encoding, use `--n-samples -1` to parse the entire file.
