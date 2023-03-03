# junction-saturation

The `junction-saturation` command will assist in determining if sequencing depth for an RNA-Seq sample is sufficient for alternative splicing analysis. The command will run [`ngsderive junction-annotation`](./junction_annotation.md) with increasing sampling rates. If sequencing depth is sufficient, we'd expect the number of junctions discovered to plateau, and increasing sampling rate will only increase the number of reads supporting each junction. If the slope is still rising at 100% sampling rate, it's likely that there are junctions present in the sample which haven't been sequenced.

The following options are adjustable for junction detection, read filtering, and sampling behavior:

| Option                               | Help                                                                                                  |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------- |
| `-i`, `--min-intron`                 | Minimum intron length to be considered a junction (default: `50`)                                     |
| `-q`, `--min-mapq`                   | Minimum read quality to be considered a supporting read (default: `30`)                               |
| `-m`, `--min-reads`                  | Minimum number of reads supporting a junction for it to be reported (default: `2`)                    |
| `-k`, `--fuzzy-junction-match-range` | Consider exonic starts/ends within `+-k` bases of a known intron start/end to be known (default: `0`) |
| `-s`, `--sample-start`               | Percentage to begin sampling (default: `5`)                                                           |
| `-p`, `--sample-step`                | Percentage to step between `sample-start` and `sample-end` (default: `5`)                             |
| `-e`, `--sample-end`                 | Percentage to stop sampling (default: `100`)                                                          |
