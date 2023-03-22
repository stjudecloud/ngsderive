# junction-annotation

The `junction-annotation` command will search an RNA-Seq bam file for splice junctions, compare them to a gene model, and output whether the found junctions are novel, partially novel, or already annotated in a gene model. In addition to a summary of the junctions and spliced reads, per-junction annotations are output with their location and number of supporting reads.

The following criteria are adjustable for junction detection and read filtering:

| Option                               | Help                                                                                                  |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------- |
| `-i`, `--min-intron`                            | Minimum intron length to be considered a junction (default: `50`)                                     |
| `-q`, `--min-mapq`                              | Minimum read quality to be considered a supporting read (default: `30`)                               |
| `-m`, `--min-reads`                             | Minimum number of reads supporting a junction for it to be reported (default: `2`)                    |
| `-k`, `--fuzzy-junction-match-range`            | Consider exonic starts/ends within `+-k` bases of a known intron start/end to be known (default: `0`) |
| `-c`, `--consider-unannotated-references-novel` | For the summary report, consider all events on unannotated reference sequences `complete_novel`. Default is to exclude them from the summary. Either way, they will be annotated as `unannotated_reference` in the junctions file. (default: False) |

In the resulting `<basename>.junctions.tsv` files, note that the coordinates are 0-based, end-exclusive. Generation of these files can be disabled with `--disable-junction-files`. Instead only the summary of junction and splice counts will be generated.
