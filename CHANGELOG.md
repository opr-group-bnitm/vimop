# Changelog

## [1.0.4] - 2025-
### Added

### Changed

### Fixed
- mainifest name (caused error when loading in EPI2ME)
- command line example run tutorial 


## [1.0.3] - 2025-11-17
### Added
- data base building module
- profiles for apptainer and conda
- tutorials for command line and EPI2ME installation and usage

### Changed
- installation and usage moved from README to separate tutorials

### Fixed
- reporting contig stats
- trimming 0 bases (before seqtk trim defaulted to quality based trimming designed for illumina reads)

## [1.0.2] - 2025-06-17
### Added

### Changed

### Fixed
- Fixed version number in manifest (required for correct version display in EPI2ME Desktop)


## [1.0.1] - 2025-06-17
### Added
- CHANGELOG.md

### Changed
- README.md edited.

### Fixed
- Fixed bug that led to crash when the sequence contains a non standard nucleotide (ACGT) that is also not an N.
