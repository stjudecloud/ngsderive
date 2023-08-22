# endedness

The `endedness` subcommand determines whether a BAM is made of Single-End or Paired-End reads **according to the BAM specification**. It will return `Unknown` if the BAM is in an ambiguous or unknown state.
If ngsderive is returning `Unknown` for a case you believe should be supported, please log an issue on our GitHub with details.

## Algorithm

By default, ngsderive will examine the bitwise FLAG field of each read primary or unmapped read in the input BAM or SAM. It keeps a tally of how the bits `0x40` (first segment in the template) and `0x80` (last segment in the template) are set.
For Single-End data, both bits should be set for every read. Each read is the only segment in the template, and is therefore both the first and the last segment.
For Paired-End data, _exactly_ half of all reads should be the first in the template (AKA read1, with `0x40` set) and _exactly_ half should be the last in the template (AKA read2, with `0x80` set).

With the `-r` or `--calc-rpt` option enabled, there is an additional check. ngsderive will count how many primary or unmapped reads belong to each unique query name (Reads-Per-Template, abbreviated RPT).
For Single-End data, RPT should be exactly 1. There should only be one primary or unmapped read mapping to each query name.
For Paired-End data, RPT should be exactly 2. Each query name should be used once for read1 and a second time for read2.
This check has been disabled by default due to its large memory requirements.

The default values for `--paired-deviance` (`0.0`) and `--round-rpt` (`False`) are suitable only if the default `--n-reads` (`-1`) is used. If only a subset of the BAM or SAM file is being processed, please set "paired deviance" to an appropriate 0\<x\<0.5 value and enable RPT rounding. An appropriate value for paired deviance depends on how much of the input file(s) is being processed.

## Limitations

Many "Single-End" BAMs do not adhere to the SAM file format specification. Specifically, bits `0x40` and `0x80` of the bitwise FLAG field are left unset by many aligners. The SAM file format specification details this case as representing a loss of information; therefore ngsderive can make no claims about the endedness of such files.
