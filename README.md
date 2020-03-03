# ngsderive

[![PyPI Version](https://img.shields.io/pypi/v/ngsderive.svg)](https://pypi.python.org/pypi/ngsderive/)
[![PyPI Python Versions](https://img.shields.io/pypi/pyversions/ngsderive.svg)](https://pypi.python.org/pypi/ngsderive/)
[![PyPI Project Status](https://img.shields.io/pypi/status/ngsderive.svg)](https://pypi.python.org/pypi/ngsderive/)


[![CI
status](https://github.com/claymcleod/ngsderive/workflows/CI/badge.svg)](https://github.com/claymcleod/ngsderive/actions)

`ngsderive` is a forensic analysis tool useful in backwards computing information 
from next-generation sequencing data. Notably, results are provided as a 'best guess' — the tool does 
not claim 100% accuracy and results should be considered with that understanding.

Note that this utility only implements commands which were not available at the 
time of writing in common NGS utilities (e.g. [Picard](https://broadinstitute.github.io/picard/)).

## Getting Started

These instructions will get you a copy of the project up and running on your
local machine for development and testing purposes. See deployment for notes on
how to deploy the project on a live system.

### Installing

To get started with `ngsderive`, you can install it using pip:

```bash
pip install git+https://github.com/claymcleod/ngsderive.git
```

### Development

To get started on a development version of the code, just run the following:

```bash
conda create -n ngsderive-dev python=3.7 poetry -y
conda activate ngsderive-dev
git clone git@github.com:claymcleod/ngsderive.git
cd ngsderive
poetry install
```
## Usage

### Illumina machine type

The `instrument` subcommand will attempt to backward compute the machine that
generated a NGS file using (1) the instrument id(s) and (2) the flowcell id(s). 

### Limitations

* This command may not comprehensively detect the correct machines as there is
  no published catalog of Illumina serial numbers. As we encounter more serial
  numbers in practice, we update this code.

### Read length calculation

The `readlen` subcommand can be used to compute the read length used during
sequencing can be estimated (or reported as inconclusive). At the time of
writing, the algorithm used is roughly:

1. Compute distribution of read lengths for the first `--n-samples` reads in a
   file.
2. Assuming read length in the file can only decrease from the actual read
   length (from adapter trimming or similar), the putative maximum read length
   is considered to be the highest detected read length.
3. If the percentage of reads that are evidence for the putative maximum read
   length makes up at least `--majority-vote-cutoff`% of the reads, the putative
   read length is considered to be confirmed. If not, the consensus read length
   will be return as -1 (could not determine).
   1. For example, if 100bp is the maximum read length detected and 85% percent
      of the reads support that claim, then we considered 100bp as the consensus
      read length. If only 30% of the reads indicated 100bp, the tool cannot
      report a consensus).

### Strandedness inference

The `strandedness` subcommand can be used to determine strandedness of RNA-Seq
samples. Note that due to the need to (1) to examine the alignment status of a
read [not included in FASTQ] and (2) fetch regions using random access [not
available with SAM files], only BAM files are currently supported for this
subcommand.

Strandedness can estimated by observing the following characteristics of a
particular read:

* Whether the read is read 1 or read 2 ("read ordinal").
* Whether the read was aligned to the + or - strand ("read strand")
* Given a gene model, whether a feature of interest (usually a gene) falls on
  the + or - strand ("gene strand").

A shorthand notation for the state of a read can be achieved by simply
concatenating the three characteristics above (e.g., `1+-` means that a read 1
was aligned to the positive strand and a gene was observed at the same location
on the negative strand).

Given the notation above, the following lookup table can be used to see whether
a read is evidence for forward-strandedness or reverse-strandedness:

| Patterns                   | Evidence for strandedness type |
| -------------------------- | ------------------------------ |
| `1++`, `1--`, `2+-`, `2-+` | Forward                        |
| `2++`, `2--`, `1+-`, `1-+` | Reverse                        |

By default, the strandedness check is designed to work with the
[GENCODE][gencode-website] geneset. Either a GTF or GFF file can be used as the
gene model — you can use the following one liner to prepare the latest geneset
for hg38.

```bash

# hg38
curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gff3.gz | gunzip -c | sort -k1,1 -k4,4n -k5,5n | bgzip > gencode.v32.annotation.gff3.gz
tabix -p gff gencode.v32.annotation.gff3.gz
```

If you would like to use the script on hg19, it takes a little more finesse (given the different formats of the attribute column between versions):

```bash
# hg19
curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.annotation.gtf.gz | gunzip -c | sort -k1,1 -k4,4n -k5,5n | python <(cat 
<<END    
import re 
import sys

for line in [l.strip() for l in sys.stdin]:
  if line.startswith("#"):
    print(line)                     
  else:                       
    columns = line.split('\t')
    if len(columns) != 9:                                                    
      raise RuntimeError("Unexpected column number: {}".format(len(columns)))
    
    print('\t'.join(columns[0:8]), end="\t")
    
    attrs_post = []
    for attr in columns[8].split(";"):                        
      groups = re.match(r"\s?(\S+) (\S+)\s?", attr)
      if groups:             
        key = groups.group(1)                    
        value = groups.group(2).replace("\"", "").replace(" ", ",")
        attrs_post.append(key + "=" + value)

    print(";".join(attrs_post))
END
) | sed 's/^chrM/chrMT/g' | sed 's/^chr//g' | bgzip > gencode.v32lift37.annotation.gtf.gz
tabix -p gff gencode.v32lift37.annotation.gtf.gz
```

At the time of writing, the algorithm works roughly like this:

1. The gene model is read in and only `gene` features are retained.
2. For `--n-genes` times, a randomly sampled gene is selected from the gene
   model. The gene must pass a quality check. Of particular interest,
   1. The gene must not be an overlapping feature on the opposite strand which
      would present ambiguous results.
   2. *Optionally*, the gene must be a protein coding gene.
   3. *Optionally*, the gene must have at least `--minimum-reads-per-gene`
      minimum reads per gene.
3. All of the reads from that region of the genome are extracted and put through
   several quality filters including but not limited to:
   1. The read must not be marked as QC-failed.
   2. The read must not be marked as a duplicate.
   3. The read must not be marked as secondary.
   4. The read must not be unmapped.
   5. *Optionally*, the read have a minimum MAPQ score.
4. For all reads that pass the above filters, compute the evidence and tally
   results.

This lookup table is used for the classification of strandedness based on the evidence:

| Lookup                            | Value              |
| --------------------------------- | ------------------ |
| 40% <= `forward_reads_pct` <= 60% | `Unstranded`       |
| 80% <= `forward_reads_pct`        | `Stranded-Forward` |
| 80% <= `reverse_reads_pct`        | `Stranded-Reverse` |
| Else                              | `Inconclusive`     |

The tool will repeat the strandedness test at most `--max-tries` times to try to find a
non-`Inconclusive`prediction.

### Differences from popular methods

The most popular strandedness inference tool that the author is aware of is
RSeQC's
[infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py). The
main difference is that RSeQC starts at the beginning of the BAM file and takes
the first `n` reads that match its criteria. If the BAM is coordinate sorted,
this would mean its not uncommon to have all of the evidence at the beginning of
`chr1`. Anecdotally, this method differs in that it is slightly slower than
`infer_experiment.py` but is expected to be more robust to biases caused by
which reads are sampled.


### Limitations

* Does not yet work with single-end data (simply because the author doesn't have
  any on hand). The tool will throw an error if any unpaired reads are
  discovered (let us know in the issues if you need this supported).
* Though hundreds of Unstranded and Stranded-Reverse data has been tested and
  verified, Stranded-Forward data has not been tested to work with this tool
  (simply because the author doesn't have on hand). We do not anticipate any
  issues since Stranded-Reverse is working well.

## Running the tests

```bash
> py.test
```

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting pull requests to us.


## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

## License

This project is licensed as follows:
* All code related to the `instrument` subcommand is licensed under the [AGPL
  v2.0][agpl-v2]. This is not due any strict requirement, but out of deference
  to some [code][10x-inspiration] I drew inspiration from (and copied patterns
  from), the decision was made to license this code consistently.
* The rest of the project is licensed under the MIT License - see the
  [LICENSE.md](LICENSE.md) file for details.

[10x-inspiration]:
https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py
[agpl-v2]: http://www.affero.org/agpl2.html
[gencode-website]:
https://www.gencodegenes.org/
