
# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


[standard-readme]: https://github.com/RichardLitt/standard-readme


## Unreleased
### Changed
- Updated dependencies

### Fixed
- Compiler error and deprecation warnings


## [1.0.1] - 2021-02-22
### Added
- A wonderful logo :-)

### Changed
- Updated README and other docs
- Some jobs in the workflow are grouped to reduce the number of cluster jobs
- Workflow requires a minimum Snakemake version
- Ignoring unused parameter in `process-pile-ups`; will be removed in next
  major release
- Disentangled workflow configuration for better usability and less build time
  for Sakemake's DAG

### Removed
- Old documentation parts/details

### Fixed
- Sporadically lost masked regions in mask homogenization
- Handling of cyclic scaffolds
- Overly strict handling of types in DENTIST's config file
- Several minor bugs


## [1.0.0] - 2021-02-04
### Added

- A Docker container! This means you can just `--use-singularity` with
  Snakemake.
- Workflow rule to just produce all the repeat masks (this is used in the
  paper to calculate the repeat content of the assemblies)
- Automatic validation of the closed gaps with an alignment of the reads
  against a preliminary gap-closed assembly:
    - Added command `bed2mask`
    - Optionally write a BED file of closed gaps
    - Added command `validate-regions`
    - Added interface for reading/writing Dazzler track extras which is
      utilized to communicate the contig and read IDs between `output` and
      `validate-regions`
- Extensively documented the example workflow config `./snakemake/snakemake.yml`
- Local alignment chaining via command `chain-local-alignments` and internally
- Using chaining to filter/improve pile up alignments
- Added possibility to revert CLI options via `--revert`
- All multi-valued CLI options take their value from a comma-separated list
  and/or by giving the same option multiple times
- Added `full_validation` flag to workflow to keep the preliminary assembly
  and validation results
- Added `no_purge_output` flag to workflow to prevent the automatic skipping
  of invalid gaps; this also will not trigger the validation if not requested
  explicitly
- Possibility to lazily read local alignments from `.las` file
- Greatly improved performance of reading `.las` files by switching to binary
  interface
- Possibility to manually skip filling of gaps
- `DBdust` for improved sensitivity in alignments
- Homogenized masks implemented via new command `propagate-mask` which
  translates a given mask via an alignment from one DB/DAM to another. The
  masks are propagated from the assembly to the reads and back to gain
  sensitivity.

### Changed
- `dentist --dependencies` now reports the availability of the listed tools
  and exits non-zero if some dependency is missing
- Avoid loading the full alignment into memory when masking
- Refactored `getAlignments` and `getFlatLocalAlignments` such that the caller
  has full control over the buffering strategy
- Streamlined option passing for `daligner`, `damapper`, `datander`, `DBdust`
  and `daccord`. Also, removed `--reference-error` and `--reads-error` in favor of default value `-e.7` in all cases. Basically, one must not modify
  `-e` without adjusting the rest of the options accordingly. Setting it to the
  minimum value just ensures no alignments are discarded for no good reason.
- Masks can be created on a "block-level" and later merged with `merged-masks`
  (this is incompatible with Dazzler's `Catrack`)
- Avoid cryptic error message if alignment is not a valid `.las` file
- Distributed tandem repeat masking
- Renamed `workdir` → `tmpdir` to avoid confusion
- Raised default value of `--max-insertion-error` (experiments show a small
  drop in correctness but large gain in contiguity)
- Replaced obscure `damapper` argument from `block_alignments` by
  `block_a=FULL_DB` or `block_b=FULL_DB`
- Behavior of environment substitution in workflow config files:
    - The config file may contain `default_env` and/or `override_env`. This
      allows to create "template" config files and "instantiation" config files
      because snakemake allows the user to specify more than one config file.
    - If just the placeholder string is given (e.g. `$FOO`) then substitute
      the exact value and type in the current `env` dictionary. This means,
      the type of a value given via `default_env` or `override_env` is
      copied which in turn prevents type error in DENTIST's config file.
- Included `filter-mask` into the standard commands because it can be useful
  for adjusting repeat masks
- Removed some dead code

### Removed
- Dependency on `LAdump` and `dumpLA`

### Fixed
- Added `properAlignmentAllowance` to `completelyCovers`
- Fixed workflow for non-PacBio reads
- Snakemake can now run in single pass; before a separate call was required to
  create DENTIST's config file.
- Deactivated `daligner`'s bridging (`-B`) when self-alignments are requested
  (`-I`) to avoid a bug.
- Many more fixes to workflow and related files under `./snakemake/`
- Details in commands in README
- Remedied syntax highlighting errors
- Many technical bugs/errors
- Compilation with `ldc2`


## [1.0.0-beta.3] - 2020-07-23
### Added
- Always skip file locking with environment variable `SKIP_FILE_LOCKING=1`

## [1.0.0-beta.2] - 2020-07-23
### Added
- Allow use of environment variables in Snakemake workflow config
- Avoid appending to DBs by design
- Improved README:
    - Advice on how to choose parameters
    - Advice on how to run DENTIST with different read types
    - Version information to dependencies
- Log level information to log messages
- More logging on failed gap closing

### Changed
- Simplified usage of `--workdir`: no need to manually create
  the designated directory
- Improvements to close more gaps:
    - Custom pre-consensus alignment filtering
    - Add support sequence to cropped reads to ensure daligner finds alignments
    - Allow cropping in masked region if necessary
    - Selectively ignore repeat mask to allow post consensus alignments
    - Increased sensitivity in pileup alignments by adding the bridging option
      of `daligner`
- Select reference read for consensus by intrinsic QVs → better
  consensus quality
- Moved flag `--max-insertion-error` from `process` to `output` stage so
  trying different values becomes much faster
- Automatically deduce trace point spacing in all places
- Faster check if `.las` files are empty → faster CLI options checking
- Naming of temporary files for easier inspection
- Use `DBdust` for post consensus alignment
- Produce `.db` for cropped pileups (temporary files) to make `DAScover`
  and `DASqv` work
- Removed `-I` option from `daligner` calls (avoid useless alignment)

### Fixed
- Several bugs in Snakemake workflow
- Significantly improved number of closed gaps
- Coordinates in AGP output
- Bug in procedure that identifies a good cropping position
- Error that caused `--proper-alignment-allowance` to have no effect by default


## [1.0.0-beta.1] - 2020-03-17
### Added
- post-consensus alignment and validation with new parameter
  `--max-insertion-error`
- inserted sequences are highlighted by upper-case letters which can be
  turned off with `--no-highlight-insertions`
- batch ranges may end with a `$` indicating the end of the pileup DB
- some mechanisms for early error detection
- write duplicate contig IDs to contig alignment cache for easier debugging
- added support for complementary contig alignments in `check-results`
- allow `.db` databases as reference
- improved version reporting
- updated README with additional instructions

### Changed
- integrated Snakemake workflow into a single file and removed "testing"
  workflow
- cropping and splicing of insertions:
    - existing sequence is completely retained
    - moved from `process-pile-ups` to `output`
    - binary format of insertions DBs (breaking change) to gain more freedom
      in later steps
    - splice sites are chosen based on the post-consensus alignments
- ambiguities in the alignment of reads are now detected globally
- weakly anchored alignments are discarded early in the filtering pipeline
- the self- and read-alignment-based masks are now computed separately
- coverage values may now be fractional
- improved README by adhering to [Standard Readme][standard-readme]
- better (error) reporting
- temporary files have more informative names
- many minor refactorings and extensions

### Removed
- combined self- and read-alignment-based masking: old behvaiour can be copied
  by using the `--mask` parameter and supplying both masks to all commands

### Fixed
- trying all possible reference reads for consensus in order to find a
  non-failing reference
- corrected insertion splicing in case of reverse-complement alignment of the
  consensus
- bug that caused `check-results` to discard all alignments in certain loci
- added missing logic for cropped contigs in `getGapState` in `check-results`


## [0.0.1] - 2018-01-03
### Added
- work-around for `damapper` bug
- histograms generated by `check-results` include a column for
  `.999` sequence identity
- `check-results` optionally writes a detailed gap report

### Changed
- simplified the coverage bounds interface of `mask-repetitive-regions`: only max-values and/or the read coverage are required
- improve consensus quality by using `lasfilteralignments` to remove deteriorating local alignments
- reduced value of `--min-reads-per-pile-up` to `--min-spanning-reads` (default: 3) to better work with extremely low coverage
- reduced the default minimum anchor length to 500
- fix & simplify score function for read alignments
- suppress generation of two las files for reads alignment in
  `generate-dazzler-options`
- insertion DBs do not include information about existing contigs anymore
  which makes validation of results easier
- renamed `debug-graph` → `debug-scaffold` for clarity
- added logging for discarded pile ups
- `check-results` counts complementary alignments (inversions) as errors
- remove leading/trailing gaps from all checks by `check-results`
- improved quality of documentation

## [0.0.1-rc.1] - 2018-10-15
### Added
- implemented routine to generate exact alignments from `.las` alignment chains
- better conflict handling for scaffold graph (see `--min-spanning-reads` and
  `--best-pile-up-margin`)
- resolve conflicting cropping information attached to contigs
- accept both `.db` and `.dam` files
- new pipeline command `process-pile-ups` to compute consensus in parallel
- introduced build config `testing` which produces additional commands for
  evaluation of results:
    - `translocate-gaps`
    - `find-closable-gaps`
    - `translate-coords`
    - `check-results`
- improved logging in general
- output of assembly graph using `--debug-graph`
- sorted CLI options alphabetically
- inspect Dazzler masks using `show-mask` command
- separate command `mask-repetitive-regions` for masking
- `show-*` commands have `--json` switch
- CLI switch `--usage`
- runtime improvements
- improved type definitions for better debugging experience
- many small bug fixes
