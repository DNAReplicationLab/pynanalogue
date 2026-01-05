# `pynanalogue`

PyNanalogue = *Py*thon *N*ucleic Acid *Analogue*
Nanalogue is a tool to parse or analyse BAM/Mod BAM files with a single-molecule focus.
We expose some of Nanalogue's functions through a python interface here.

[![Python Tests (3.9-3.14 Ubuntu & Mac), Benchmark, Linting](https://github.com/DNAReplicationLab/pynanalogue/actions/workflows/python-tests.yml/badge.svg)](https://github.com/DNAReplicationLab/pynanalogue/actions/workflows/python-tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A common pain point in genomics analyses is that BAM files are information-dense
which makes it difficult to gain insight from them. PyNanalogue hopes to make it easy
to extract and process this information, with a particular focus on single-molecule
aspects and DNA/RNA modifications. Despite this focus, some of pynanalogue's functions are
quite general and can be applied to almost any BAM file.

# Requirements

- Python 3.9 or higher
- Rust toolchain (for building from source)

# Installation

This section in development.

PyNanalogue should be available on PyPI.
You can run the following command or an equivalent to install it.

```bash
pip install nanalogue
```

Common wheels (`manylinux/mac`) are available there.
Please open an issue if you want more wheels!

## Musllinux note

Although we have a `musllinux` wheel, a package we depend on, `polars-runtime-32`, does not
have one (see this [issue](https://github.com/pola-rs/polars/issues/25568)), which means
you may have to build the wheel yourself during installation. Please open an issue if you
spot that this has changed! This issue affects only those that use an Alpine Linux distribution
or similar. If you are using Ubuntu/Debian etc., this issue should not affect you.

# Functions

Our package exposes the following python functions.
They usually have lots of optional arguments.
Among other operations, the options allow you to subsample the BAM file (`sample_fraction`),
restrict read and/or modification data to a specific genomic region (`region` or `mod_region`),
restrict by one or several read ids (`read_ids`),
a specific mapping type (`read_filter`), filter modification data suitably
(`min_mod_qual`, `reject_mod_qual_non_inclusive`) etc.

## Read info

Prints information about reads in JSON. A sample output snippet follows.

```python
import pynanalogue as pn
import json
result_bytes = pn.read_info("some_file.bam")
decoded_output = json.loads(result_bytes)
```

A record from the decoded output might look like
```json
[
{
   "read_id": "cd623d4a-510d-4c6c-9d88-10eb475ac59d",
   "sequence_length": 2104,
   "contig": "contig_0",
   "reference_start": 7369,
   "reference_end": 9473,
   "alignment_length": 2104,
   "alignment_type": "primary_reverse",
   "mod_count": "C-m:263;N+N:2104;(probabilities >= 0.5020, PHRED base qual >= 0)"
}
]
```

## Window reads

Documentation in development.

## Polars bam mods

Documentation in development.

## Simulate mod bam

Documentation in development.

# Further documentation

In addition to this repository, we are developing a
companion cookbook [here](https://www.nanalogue.com).

# Versioning

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
