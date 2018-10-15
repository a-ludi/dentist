# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed

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
