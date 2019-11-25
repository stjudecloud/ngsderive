# ngsderive

[![CI status](https://github.com/claymcleod/ngsderive/workflows/CI/badge.svg)](https://github.com/claymcleod/ngsderive/actions)

`ngsderive` is a set of utilities developed to backwards compute various attributes from
next-generation sequencing data. This command line utility only implements
commands which were not available at the time of writing in common NGS utilities
(e.g., [Picard](https://broadinstitute.github.io/picard/)). All utilities are
provided as-is with no warranties: many tools just provide suggestions, and
though we've done the best we can + this toolkit evolves as we learn more, we 
don't claim 100% accuracy for all utilities provided.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Installing

To get started with `ngsderive`, you can install it using pip:

```bash
pip install git+https://github.com/claymcleod/ngsderive.git
```

Currently, the following commands are supported:

```bash
> ngsderive                                                                                                                                                      ✔  6699  12:12:51
Usage:
  ngsderive readlen <ngsfiles>... [--outfile=<outfile>]
  ngsderive machine <ngsfiles>... [--outfile=<outfile>]
  ngsderive (-h | --help)
  ngsderive --version
```

## Usage

### Illumina machine type

The `ngsderive instrument` subcommand will attempt to backward compute the
machine that generated a NGS file using (1) the instrument id(s) and (2) the
flowcell id(s). Note that this command may not comprehensively detect the
correct machines as there is no published catalog of Illumina serial numbers.
As we encounter more serial numbers in practice, we update this code.

### Read length calculation

Using the `ngsderive readlen` subcommand, one can backward compute the
readlength used during sequencing. Currently, the algorithm used is roughly:

1. Compute distribution of read lengths for the first `n` reads in a file
   (default: 10000).
2. If the maximum read length makes up `p`% of the reads, the read length is
   equal to that number (e.g., if 85% percent of the reads are 100bp, then the
   read length is considered 100bp.)
3. If #2 does not hold true, the read length cannot be computed confidently.

## Running the tests

```bash
> py.test
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## License

This project is licensed as follows:
* All code related to the `instrument` subcommand is licensed under the [AGPL
  v2.0][agpl-v2]. This is not due any strict requirement, but out of deference
  to some [code][10x-inspiration] I drew inspiration from (and copied patterns
  from), the decision was made to license this code consistently.
* The rest of the project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

[agpl-v2]: http://www.affero.org/agpl2.html
[10x-inspiration]:
https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py
