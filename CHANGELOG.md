# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.9.0] - 2021-09-02
### Added
* CHANGELOG.md will document notable changes
### Changed
* `--bam` based dehumanising checks for BAM index with pysam `AlignmentFile.has_index` to skip count pass and use `idxstats` if index is available for modest improvement
* `--bam` based dehumanising uses multiprocessing and breaks the input BAM into chunks
