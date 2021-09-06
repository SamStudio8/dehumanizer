# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [ZeroVer Versioning](https://0ver.org/).

## [0.9.0] - 2021-09-02
### Added
* CHANGELOG.md will document notable changes
* `dehumanise --version` will report version number
### Changed
* `--bam` based dehumanising checks for BAM index with pysam `AlignmentFile.has_index` to skip count pass and use `idxstats` if index is available for modest improvement
* Merged the default humref, decoy and HLA references into one FASTA available for download as most users do not care about which ref causes a read to be discarded
### Fixed
* Cleared up confusion surrounding pre-built mmi indexes by reintroducing manifest syntax to map each input file to a minimap2 preset
    * `dehumanise` will exit 65 to indicate that no references for a given preset are listed in the manifest
    * `dehumanise` will exit 78 to indicate that a reference in the manifest is not mapped to a minimap2 preset that matches `--preset`
* `dehumanise --help` now works
