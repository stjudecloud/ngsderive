# strandedness

The `strandedness` subcommand can be used to determine strandedness of RNA-Seq samples. Note that due to the need to (1) examine the alignment status of a read [not included in FASTQ] and (2) fetch regions using random access [not available with SAM files], only BAM files are currently supported for this subcommand.

Strandedness can estimated by observing the following characteristics of a particular read:

* Whether the read is read 1 or read 2 ("read ordinal").
* Whether the read was aligned to the + or - strand ("read strand")
* Given a gene model, whether a feature of interest (usually a gene) falls on the + or - strand ("gene strand").

A shorthand notation for the state of a read can be achieved by simply concatenating the three characteristics above (e.g., `1+-` means that a read 1 was aligned to the positive strand and a gene was observed at the same location on the negative strand).

Given the notation above, the following lookup table can be used to see whether a read is evidence for forward-strandedness or reverse-strandedness:

| Patterns                   | Evidence for strandedness type |
| -------------------------- | ------------------------------ |
| `1++`, `1--`, `2+-`, `2-+` | Forward                        |
| `2++`, `2--`, `1+-`, `1-+` | Reverse                        |

Any valid GTF or GFF file should be compatible with `ngsderive`. The strandedness check requires random access of the gene model, which means the file must be sorted and tabixed. If you supply an unsorted or untabixed GTF/GFF, we will create one for you using the supplied file. If you encounter any errors related to your gene model choice, please let us know using [GitHub Issues][issues].

## Algorithm

At the time of writing, the algorithm works roughly like this:

1. The gene model is read in and only `gene` features are retained.
2. For `--n-genes` times, a randomly sampled gene is selected from the gene model. The gene must pass a quality check. Of particular interest,
   1. The gene must not be an overlapping feature on the opposite strand which would present ambiguous results.
   2. *Optionally*, the gene must be a protein coding gene. This defaults to `True`.
   3. *Optionally*, the gene must have at least `--minimum-reads-per-gene` minimum reads per gene. This defaults to 10 reads.
3. All of the reads from that region of the genome are extracted and put through several quality filters including but not limited to:
   1. The read must not be marked as QC-failed.
   2. The read must not be marked as a duplicate.
   3. The read must not be marked as secondary.
   4. The read must not be unmapped.
   5. *Optionally*, the read have a minimum MAPQ score. This defaults to a MAPQ score of 30.
4. For all reads that pass the above filters, compute the evidence and tally results.

This lookup table is used for the classification of strandedness based on the evidence:

| Lookup                            | Value              |
| --------------------------------- | ------------------ |
| 40% <= `forward_reads_pct` <= 60% | `Unstranded`       |
| 80% <= `forward_reads_pct`        | `Stranded-Forward` |
| 80% <= `reverse_reads_pct`        | `Stranded-Reverse` |
| Else                              | `Inconclusive`     |

The tool will repeat the strandedness test at most `--max-tries` times to try to find a non-`Inconclusive` prediction.

## Differences

The most popular strandedness inference tool that the author is aware of is RSeQC's [infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py). The main difference is that RSeQC starts at the beginning of the BAM file and takes the first `n` reads that match its criteria. If the BAM is coordinate sorted, this would mean its not uncommon to have all of the evidence at the beginning of `chr1`. Anecdotally, this method differs in that it is slightly slower than `infer_experiment.py` but is expected to be more robust to biases caused by which reads are sampled. Further, reads within RSeQC may fall into ambiguous regions for which a read is not decisively evidence for a particular strandedness. These cases are eliminated by the use of (a) filtering to only coding regions and (b) removal of any regions of the gene model for which genes overlap on each strand.

## Limitations

* Does not yet work with single-end data (simply because the author doesn't have any on hand). The tool will throw an error if any unpaired reads are discovered (let us know in [the issues][issues] if you need this supported).
* Though hundreds of Unstranded and Stranded-Reverse data has been tested and verified, Stranded-Forward data has not been tested to work with this tool (simply because the author doesn't have on hand). We do not anticipate any issues since Stranded-Reverse is working well.

[gencode-website]: https://www.gencodegenes.org
[issues]: https://github.com/stjudecloud/ngsderive/issues
