# `pynanalogue`

Nanalogue = *N*ucleic Acid *Analogue*
Nanalogue is a tool to parse or analyse BAM/Mod BAM files with a single-molecule focus.
We expose some of Nanalogue's functions through a python interface here.

Software is in the alpha stage.

## Requirements

- Python 3.9 or higher
- Rust toolchain (for building from source)

## Versioning

We use [Semantic Versioning](https://semver.org/) (SemVer) for version numbers.

**Current Status: Pre-1.0 (0.x.y)**

While in 0.x.y versions:
- The API may change without notice
- Breaking changes can occur in minor version updates
- This is a development phase with no stability guarantees

**After 1.0.0 Release:**

Once we reach version 1.0.0, we will guarantee:
- No breaking changes in minor (x.**Y**.z) or patch (x.y.**Z**) releases
- Clear migration guides for major version updates
- Deprecation warnings at least one minor version before removal of features

# Acknowledgments

This software was developed at the Earlham Institute in the UK.
This work was supported by the Biotechnology and Biological Sciences
Research Council (BBSRC), part of UK Research and Innovation,
through the Core Capability Grant BB/CCG2220/1 at the Earlham Institute
and the Earlham Institute Strategic Programme Grant Cellular Genomics
BBX011070/1 and its constituent work packages BBS/E/ER/230001B 
(CellGen WP2 Consequences of somatic genome variation on traits).
The work was also supported by the following response-mode project grants:
BB/W006014/1 (Single molecule detection of DNA replication errors) and
BB/Y00549X/1 (Single molecule analysis of Human DNA replication).
This research was supported in part by NBI Research Computing
through use of the High-Performance Computing system and Isilon storage.
