# Changelog

<!--next-version-placeholder-->

## v2.2.1 (2023-02-21)
### Fix
* **strandedness** Make `--max-iterations-per-try` dynamic, based on `--n-genes` instead of a static default of `1000` ([#88](https://github.com/stjudecloud/ngsderive/pull/88))

## v2.2.0 (2021-06-24)
### Feature
* **instrument:** Attempt to recover when read names not in Illumina format ([#14](https://github.com/stjudecloud/ngsderive/issues/14)) ([`9ab0951`](https://github.com/stjudecloud/ngsderive/commit/9ab0951a04bdcffdd12f2291ec25c46d74d666d1))

## v2.1.0 (2021-04-05)
### Feature
* Instead of failing on an unknown contig, try new gene ([#11](https://github.com/stjudecloud/ngsderive/issues/11)) ([`6ff35a5`](https://github.com/stjudecloud/ngsderive/commit/6ff35a5c033c5d16145562cb85a2177c65ddd478))

## v2.0.0 (2021-03-11)
### Fix
* Remove the option to override the default tab output delimiter ([#12](https://github.com/stjudecloud/ngsderive/issues/12)) ([`fbe382d`](https://github.com/stjudecloud/ngsderive/commit/fbe382d67e22f3cc585a86da21b788d7e79f76a3))

### Breaking
* This removes the `--delimiter` commandline arg ([`fbe382d`](https://github.com/stjudecloud/ngsderive/commit/fbe382d67e22f3cc585a86da21b788d7e79f76a3))

## v1.4.1 (2021-02-05)
### Fix
* **junction-annotation:** Capitalize `file` in result outfile ([`fa01f8b`](https://github.com/stjudecloud/ngsderive/commit/fa01f8ba8312ad396bd6161e4d6ae81bbf1388b5))

## v1.4.0 (2021-02-03)
### Feature
* Junction annotation ([#7](https://github.com/stjudecloud/ngsderive/issues/7)) ([`f723731`](https://github.com/stjudecloud/ngsderive/commit/f723731343a3ab328643fa2d57391c91ef6efd43))

## v1.3.0 (2021-02-02)
### Feature
* Bam support for encoding ([#9](https://github.com/stjudecloud/ngsderive/issues/9)) ([`b1941a7`](https://github.com/stjudecloud/ngsderive/commit/b1941a74b36ab6b578e91273bc9ca41743dfaefa))
