# Changelog
## [1.0.7]
### Added
- conda profile: 1 h timeout for creating the conda environments to avoid hanging on slow connections

### Changed
- code updated for the strict syntax parser of Nextflow 26 (default since 26.04); older Nextflow versions still work (checked with 24.10)
- ingress helpers that start processes are now sub-workflows
- container images are set in one place (`params.images`), with the same images and versions as before
- conda environments no longer list the `defaults` channel; channels are ordered conda-forge, bioconda, nanoporetech, and the centrifuge environment now includes conda-forge
- process directive syntax (`memory` in `classify_contigs`)

### Fixed
- boolean and numeric parameters given on the command line, which Nextflow 26 passes as text, are evaluated correctly (e.g. `--validate_params false`, `--assemble_notarget false`)
- crash with `--monochrome_logs` or `--min_*` resource parameters on the command line with Nextflow 26
- parameter validation rejected boolean and numeric parameters given on the command line with Nextflow 26 (e.g. `--min_cpus 8`, `--download_db_all`)
- `cleanup = true` had no effect because it was inside `process {}`; as a top-level setting it now deletes a run's files in the work directory when the run completes successfully (failed runs keep them for `-resume`)
- conda profile: the sample report failed with `No module named 'pkg_resources'` (removed in setuptools 82, still used by ezcharts 0.11.2); the report environment now pins `setuptools<82` and `pandas<3`, as when the Docker image was built

### Removed
- unsupported `maxParallel` setting in `download_parts` (downloads are still limited to one at a time by `maxForks`)
- unused `imageVersionDownload` setting
- `nextflow.enable.dsl = 2` statements (DSL2 is the default)

## [1.0.6] - 2026-06-29
### Added

### Changed

### Fixed
- remove unsupported symbolic alleles when using bcftools consensus (currently supported: <DEL>, <*>, <NON_REF>)

### Removed

## [1.0.5] - 2026-04-22
### Added
- paper citation
  
### Changed

### Fixed
- error when inversions denoted with <INV> are in the vcf-file before consensus calling with bcftools

### Removed

## [1.0.4] - 2025-12-10
### Added

### Changed
- tutorial revised

### Fixed
- manifest name (caused bug when loading in EPI2ME desktop version 5.3.0)
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
