# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0.0] - 2026-07-25
Adjustments to allow flexibility in the way the workflow is run, as discssed in
- [GRD-1174](https://jira.oicr.on.ca/browse/GRD-1174)
- added various enables arguments to the runDragenGermline, that default to running snv calling only
- enableCnvSelfNormalization included, which run Cnv calls in a manner that generates a uniformity of coverage metric.  This is needed for QC

## [Unreleased] - 2025-01-16
Initial Implementation
- [GRD-855](https://jira.oicr.on.ca/browse/GRD-855)
