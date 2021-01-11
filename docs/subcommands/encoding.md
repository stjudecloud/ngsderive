# encoding

The `encoding` subcommand detects the likely PHRED score encoding scheme used for FASTQ files. PHRED scores are encoded as ASCII characters, but which ASCII characters encode for which relative quality scores depends on the build and generation of the sequencing machine used. `ngsderive` obtained details of the different encoding schemes [here][https://en.wikipedia.org/wiki/FASTQ_format#Encoding].

## Limitations

* BAM files are not currently supported for encoding detection.

* All 3 possible encodings have the same upper range, and only differ in terms of lower quality scores. Thus if the provided read data is all of high quality, it may be classified as a stricter encoding than was originally used to generate the data.

    If the downstream use of this value is needed to ensure all PHRED scores are within the derived encoding, use `--n-samples -1` to parse the entire file.
