<p align="center">
  <h1 align="center">
    ngsderive
  </h1>

  <p align="center">
    <a href="https://actions-badge.atrox.dev/stjudecloud/ngsderive/goto" target="_blank">
      <img alt="Actions: CI Status"
          src="https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2Fstjudecloud%2Fngsderive%2Fbadge&style=flat" />
    </a>
    <a href="https://pypi.org/project/ngsderive/" target="_blank">
      <img alt="PyPI"
          src="https://img.shields.io/pypi/v/ngsderive?color=orange">
    </a>
    <a href="https://pypi.python.org/pypi/ngsderive/" target="_blank">
      <img alt="PyPI: Downloads"
          src="https://img.shields.io/pypi/dm/ngsderive?color=orange">
    </a>
    <a href="https://pypi.python.org/pypi/ngsderive/" target="_blank">
      <img alt="PyPI: Downloads"
          src="https://img.shields.io/pypi/pyversions/ngsderive?color=orange">
    </a>
    <a href="https://codecov.io/gh/stjudecloud/ngsderive" target="_blank">
      <img alt="Code Coverage"
          src="https://codecov.io/gh/stjudecloud/ngsderive/branch/master/graph/badge.svg" />
    </a>
    <a href="https://github.com/stjudecloud/ngsderive/blob/master/LICENSE.md" target="_blank">
    <img alt="License: MIT"
          src="https://img.shields.io/badge/License-MIT-blue.svg" />
    </a>
  </p>


  <p align="center">
    Forensic analysis tool useful in backwards computing information from next-generation sequencing data. 
    <br />
    <a href="https://stjudecloud.github.io/ngsderive/"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/stjudecloud/ngsderive/issues/new?assignees=&labels=&template=feature_request.md&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
    ·
    <a href="https://github.com/stjudecloud/ngsderive/issues/new?assignees=&labels=&template=bug_report.md&title=Descriptive%20Title&labels=bug">Report Bug</a>
    ·
    ⭐ Consider starring the repo! ⭐
    <br />
  </p>
</p>

Notice: `ngsderive` is a forensic analysis tool useful in backwards computing information 
from next-generation sequencing data. Notably, results are provided as a 'best guess' — the tool does not claim 100% accuracy and results should be considered with that understanding.

Note that this utility only implements commands which were not available at the 
time of writing in common NGS utilities (e.g. [Picard](https://broadinstitute.github.io/picard/)).

## 🎨 Features

The following attributes can be guessed using ngsderive:

* <b>Illumina Instrument.</b> Infer which Illumina instrument was used to generate the data by matching against known instrument and flowcell naming patterns. Each guess comes with a confidence score. 
* <b>RNA-Seq Strandedness.</b> Infer from the data whether RNA-Seq data was generated using a Stranded-Forward, Stranded-Reverse, or Unstranded protocol.
* <b>Pre-trimmed Read Length.</b> Compute the distribution of read lengths in the file and attempt to guess what the original read length of the experiment was.

## 📚 Getting Started

### Installation

#### Python Package Index

You can install ngsderive using the Python Package Index ([PyPI](https://pypi.org/)).

```bash
pip install ngsderive
```

## 🖥️ Development

If you are interested in contributing to the code, please first review
our [CONTRIBUTING.md][contributing-md] document. 

To bootstrap a development environment, please use the following commands.

```bash
# Clone the repository
git clone git@github.com:stjudecloud/ngsderive.git
cd ngsderive

# Install the project using poetry
poetry install
```

## 🚧️ Tests

ngsderive provides a (currently patchy) set of tests — both unit and end-to-end.

```bash
py.test
```

## 🤝 Contributing

Contributions, issues and feature requests are welcome!<br />Feel free to check [issues page](https://github.com/stjudecloud/ngsderive/issues). You can also take a look at the [contributing guide][contributing-md].

## 📝 License

Copyright © 2020 [St. Jude Cloud Team](https://github.com/stjudecloud).<br />
This project is [MIT][license-md] licensed.

[contributing-md]: https://github.com/stjudecloud/ngsderive/blob/master/CONTRIBUTING.md
[license-md]: https://github.com/stjudecloud/ngsderive/blob/master/LICENSE.md
