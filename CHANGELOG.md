# Changelog


## [1.0.4] - 2025-12-
### Added

### Changed

### Fixed
- manifest name (caused bug when loading in EPI2ME desktop version 5.3.0)
- tutorial error (--out_dir instead of -o)
- bug "Cannot invoke method optional() on null object" when using Nextflow 25.10.2
- conda profile works with "nextflow run" command

### Removed


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
