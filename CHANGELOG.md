# Changelog

<!--next-version-placeholder-->

## v3.3.1 (2023-08-17)

### Fix

* **strandedness:** Pass in required max_tries arg to main() ([#136](https://github.com/stjudecloud/ngsderive/issues/136)) ([`544be1a`](https://github.com/stjudecloud/ngsderive/commit/544be1ab85bcbede256b3f6b74fb80002a3de459))

## v3.3.0 (2023-08-16)

### Feature

* Endedness subcommand ([#124](https://github.com/stjudecloud/ngsderive/issues/124)) ([`d61e82d`](https://github.com/stjudecloud/ngsderive/commit/d61e82dd53272929dca2e825d172aff2bb97bacc))

### Documentation

* **encoding:** `n-samples` -> `n-reads` ([`665d1d2`](https://github.com/stjudecloud/ngsderive/commit/665d1d2dac56a09bed4101a52ce2f3cddeecd738))

## v3.2.1 (2023-06-27)
### Fix
* -1 n-reads default for encoding. don't falsely imply random sampling ([#119](https://github.com/stjudecloud/ngsderive/issues/119)) ([`c83275a`](https://github.com/stjudecloud/ngsderive/commit/c83275a018703c57f0ce8c90cf4b94fbd658e300))

## v3.2.0 (2023-06-26)
### Feature
* Handle case of BAM without RG tag ([#100](https://github.com/stjudecloud/ngsderive/issues/100)) ([`4436abf`](https://github.com/stjudecloud/ngsderive/commit/4436abf70540a33b1825deba3c8b0b499b93bbc5))

## v3.1.1 (2023-05-08)
### Fix
* Revert "feat: junction-saturation (#16)" ([#110](https://github.com/stjudecloud/ngsderive/issues/110)) ([`9262af0`](https://github.com/stjudecloud/ngsderive/commit/9262af0d8290e3f1ae70ba47b58fbea9b9c15d35))

## v3.1.0 (2023-03-23)
### Feature
* Junction-saturation ([#16](https://github.com/stjudecloud/ngsderive/issues/16)) ([`e4811e7`](https://github.com/stjudecloud/ngsderive/commit/e4811e73e2164dd8f6150d6f3ae2f0f6a8af2f99))

## v3.0.0 (2023-03-22)
### Feature
* **junction-annotation:** By default discard unannotated contigs ([#96](https://github.com/stjudecloud/ngsderive/issues/96)) ([`d2f2e18`](https://github.com/stjudecloud/ngsderive/commit/d2f2e18338656da53e2d6d8c74f2a48850a1f125))
* Adds Dockerfile ([#97](https://github.com/stjudecloud/ngsderive/issues/97)) ([`26242b2`](https://github.com/stjudecloud/ngsderive/commit/26242b27fe2d73e2413f872a61496cffca3922f2))

### Breaking
* Always annotate as `unannotated_reference` in the junctions file. By default, exclude from the summary report. With flag, call them "novel" in the summary ([`d2f2e18`](https://github.com/stjudecloud/ngsderive/commit/d2f2e18338656da53e2d6d8c74f2a48850a1f125))

## v2.3.1 (2023-02-22)
### Fix
* Set default `n-genes` to 1000 for improved accuracy ([`914525d`](https://github.com/stjudecloud/ngsderive/commit/914525da2449ad6588c1518287a9b73cc91cf9d4))

## v2.3.0 (2023-02-21)
### Feature
* Make retries cumulative instead of a reset ([#91](https://github.com/stjudecloud/ngsderive/issues/91)) ([`72c7b0f`](https://github.com/stjudecloud/ngsderive/commit/72c7b0fca9c72b4caafb307451ff3b1c08c671bc))

### Documentation
* Clarify other important improvement regarding disambiguation of… ([#15](https://github.com/stjudecloud/ngsderive/issues/15)) ([`5d58156`](https://github.com/stjudecloud/ngsderive/commit/5d58156f79ddf57f8d0518bc45745e923d8f4c04))
* Refactor to exclude junction-annotation from "best guess" notice ([#13](https://github.com/stjudecloud/ngsderive/issues/13)) ([`c8e208e`](https://github.com/stjudecloud/ngsderive/commit/c8e208e2edf4c8d81fb5cddc322b38a5acc57a44))

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
