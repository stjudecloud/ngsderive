# junction-annotation

The `junction-annotation` command will search a bam file for splice junctions, compare them to a gene model, and output whether the found junctions are novel, partially novel, or already annotated in a gene model. In addition to a summary of the junctions and spliced reads, per-junction annotations are output with their location and number of supporting reads.

The following criteria are adjustable for junction detection and read filtering:

| Option                               | Help                                                                                                  |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------- |
| `-i`, `--min-intron`                 | Minimum intron length to be considered a junction (default: `50`)                                     |
| `-q`, `--min-mapq`                   | Minimum read quality to be considered a supporting read (default: `30`)                               |
| `-m`, `--min-reads`                  | Minimum number of reads supporting a junction for it to be reported (default: `2`)                    |
| `-k`, `--fuzzy-junction-match-range` | Consider exonic starts/ends within `+-k` bases of a known intron start/end to be known (default: `0`) |
