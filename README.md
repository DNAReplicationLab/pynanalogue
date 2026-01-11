# `pynanalogue`

PyNanalogue = *Py*thon *N*ucleic Acid *Analogue*.

Nanalogue is a tool to parse or analyse BAM/Mod BAM files with a single-molecule focus.
We expose some of Nanalogue's functions through a python interface here.

[![Python Tests (3.9-3.14 Ubuntu & Mac), Benchmark, Linting](https://github.com/DNAReplicationLab/pynanalogue/actions/workflows/python-tests.yml/badge.svg)](https://github.com/DNAReplicationLab/pynanalogue/actions/workflows/python-tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A common pain point in genomics analyses is that BAM files are information-dense
which makes it difficult to gain insight from them. PyNanalogue hopes to make it easy
to extract and process this information, with a particular focus on single-molecule
aspects and DNA/RNA modifications. Despite this focus, some of pynanalogue's functions are
quite general and can be applied to almost any BAM file.

We can process any type of DNA/RNA modifications occuring in any pattern (single/multiple mods,
spatially-isolated/non-isolated etc.). All we require is that the data is stored
in a BAM file in the mod BAM format (i.e. using MM/ML tags as laid down in the
[specifications](https://samtools.github.io/hts-specs/SAMtags.pdf)).

## Table of Contents

Note: these links work in Github but may not work on PyPI.
But the table of contents is correct.

- [Requirements](#requirements)
- [Installation](#installation)
  - [More details](#more-details)
  - [Note on musllinux](#note-on-musllinux)
- [Functions](#functions)
  - [Read info](#read-info)
    - [Documentation](#documentation)
    - [Sample input and output](#sample-input-and-output)
  - [Window reads](#window-reads)
    - [Documentation](#documentation-1)
    - [Sample input and output](#sample-input-and-output-1)
  - [Polars bam mods](#polars-bam-mods)
    - [Documentation](#documentation-2)
    - [Sample input and output](#sample-input-and-output-2)
  - [Simulate mod bam](#simulate-mod-bam)
    - [Example](#example)
- [Further documentation](#further-documentation)
- [Versioning](#versioning)
- [Acknowledgments](#acknowledgments)

# Requirements

- Python 3.9 or higher
- Rust toolchain (for building from source)

# Installation

PyNanalogue should be available on PyPI.
You can run the following command or an equivalent to install it.

```bash
pip install pynanalogue
```

## More details

Common wheels (`manylinux/mac`) are available in PyPI.
Please open an issue if you want more wheels!

### Note on musllinux

Although we have a `musllinux` wheel, a package we depend on (`polars-runtime-32`) does not
have one (see this [issue](https://github.com/pola-rs/polars/issues/25568)), which means
you may have to build the wheel yourself during installation. Please open an issue if you
spot that this has changed! This issue affects only those that use an Alpine Linux distribution
or similar. If you are using Ubuntu/Debian/Fedora etc., this issue should not affect you.

# Functions

Our package exposes the following python functions.
They usually have lots of optional arguments.
Among other operations, the options allow you to subsample the BAM file (`sample_fraction`),
restrict read and/or modification data to a specific genomic region (`region` or `mod_region`),
restrict by one or several read ids (`read_ids`),
a specific mapping type (`read_filter`), filter modification data suitably
(`min_mod_qual`, `reject_mod_qual_non_inclusive`) etc.
Please see each section below for the list of options, which you can access by
reading the docstring of each function.

## Read info

Prints information about reads in JSON.
In this section, we show how to get documentation about the function,
a sample execution, and a sample output snippet.

### Documentation

The function has lots of helpful optional arguments.
Please run the command below to see what they are.

```python
import pynanalogue as pn
print(pn.read_info.__doc__)
```

### Sample input and output

A sample execution and output snippet follows.

```python
import pynanalogue as pn
import json
result_bytes = pn.read_info("tests/data/examples/example_1.bam")
decoded_output = json.loads(result_bytes)
```

A record from the decoded output might look like the following.
You will get one record per alignment.
```json
[
{
	"read_id": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
	"sequence_length": 48,
	"contig": "dummyIII",
	"reference_start": 23,
	"reference_end": 71,
	"alignment_length": 48,
	"alignment_type": "primary_forward",
	"mod_count": "T+T:3;(probabilities >= 0.5020, PHRED base qual >= 0)"
}
]
```

## Window reads

Output windowed modification densities of reads as a polars dataframe.
In this section, we show how to get documentation about the function,
a sample execution, and a sample output snippet.

### Documentation

The function has lots of helpful optional arguments and some required arguments
like window size and step size.
Please run the command below to see what they are.

```python
import pynanalogue as pn
print(pn.window_reads.__doc__)
```

### Sample input and output

A sample execution and output snippet follows.

```python
import pynanalogue as pn
import polars as pl
df = pn.window_reads("tests/data/examples/example_1.bam",win = 2,step = 1)
```

The output is a polars dataframe. If printed in tsv format, a few rows may look like this.
(This was generated from a file without basecalling quality information, which is why we show 255s here).
```text
#contig	ref_win_start	ref_win_end	read_id	win_val	strand	base	mod_strand	mod_type	win_start	win_end	basecall_qual
dummyI	9	13	5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	0	+	T	+	T	0	4	255
dummyI	12	14	5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	0	+	T	+	T	3	5	255
dummyI	13	17	5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	0	+	T	+	T	4	8	255
dummyIII	26	32	a4f36092-b4d5-47a9-813e-c22c3b477a0c	1	+	T	+	T	3	9	255
dummyIII	31	51	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	+	T	+	T	8	28	255
dummyIII	50	63	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0	+	T	+	T	27	40	255
dummyIII	62	71	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	+	T	+	T	39	48	255
dummyII	15	17	fffffff1-10d2-49cb-8ca3-e8d48979001b	0	-	T	+	T	12	14	255
dummyII	16	20	fffffff1-10d2-49cb-8ca3-e8d48979001b	0	-	T	+	T	13	17	255
dummyII	19	23	fffffff1-10d2-49cb-8ca3-e8d48979001b	0	-	T	+	T	16	20	255
dummyII	22	24	fffffff1-10d2-49cb-8ca3-e8d48979001b	0.5	-	T	+	T	19	21	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0	.	G	-	7200	28	30	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0	.	G	-	7200	29	31	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0	.	G	-	7200	30	33	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0	.	G	-	7200	32	44	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0	.	G	-	7200	43	45	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	1	.	T	+	T	3	9	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	.	T	+	T	8	28	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0	.	T	+	T	27	40	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	.	T	+	T	39	48	255
```

## Polars bam mods

Output raw modification data as a polars dataframe.
In this section, we show how to get documentation about the function,
a sample execution, and a sample output snippet.
Please note that as we report every modified position per molecule as a separate
row in a dataframe, the data size could get very big.
So, we recommend querying per region or subsampling the BAM file
in order to not run into memory issues -- there are options in this
function to do so.
We may develop an iterable version of this function in the future.
Please open an issue if you are interested in this, or a pull request if you can do this!

### Documentation

The function has lots of helpful optional arguments.
Please run the command below to see what they are.

```python
import pynanalogue as pn
print(pn.polars_bam_mods.__doc__)
```

### Sample input and output

A sample execution and output snippet follows.

```python
import pynanalogue as pn
import polars as pl
df = pn.polars_bam_mods("tests/data/examples/example_1.bam")
```

The output is a polars dataframe. If printed in tsv format, it might look like this.
Mod quality is a probability represented as a number between 0 and 255,
where 0 means not modified and 255 means modified with certainty.
This is how modification data is stored in the mod BAM format.

```text
read_id	seq_len	alignment_type	align_start	align_end	contig	contig_id	base	is_strand_plus	mod_code	position	ref_position	mod_quality
5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	8	primary_forward	9	17	dummyI	0	T	true	T	0	9	4
5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	8	primary_forward	9	17	dummyI	0	T	true	T	3	12	7
5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	8	primary_forward	9	17	dummyI	0	T	true	T	4	13	9
5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	8	primary_forward	9	17	dummyI	0	T	true	T	7	16	6
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	primary_forward	23	71	dummyIII	2	T	true	T	3	26	221
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	primary_forward	23	71	dummyIII	2	T	true	T	8	31	242
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	primary_forward	23	71	dummyIII	2	T	true	T	27	50	3
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	primary_forward	23	71	dummyIII	2	T	true	T	39	62	47
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	primary_forward	23	71	dummyIII	2	T	true	T	47	70	239
fffffff1-10d2-49cb-8ca3-e8d48979001b	33	primary_reverse	3	36	dummyII	1	T	true	T	12	15	3
fffffff1-10d2-49cb-8ca3-e8d48979001b	33	primary_reverse	3	36	dummyII	1	T	true	T	13	16	3
fffffff1-10d2-49cb-8ca3-e8d48979001b	33	primary_reverse	3	36	dummyII	1	T	true	T	16	19	4
fffffff1-10d2-49cb-8ca3-e8d48979001b	33	primary_reverse	3	36	dummyII	1	T	true	T	19	22	3
fffffff1-10d2-49cb-8ca3-e8d48979001b	33	primary_reverse	3	36	dummyII	1	T	true	T	20	23	182
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	28	-1	0
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	29	-1	0
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	30	-1	0
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	32	-1	0
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	43	-1	77
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	44	-1	0
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					T	true	T	3	-1	221
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					T	true	T	8	-1	242
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					T	true	T	27	-1	0
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					T	true	T	39	-1	47
a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					T	true	T	47	-1	239
```

## Simulate mod bam

If you are a developer who needs BAM files with defined single-molecule modification patterns
to help develop/test your tool, nanalogue can also help you create BAM files from scratch
using artificial data created using parameters defined by you.

```python
import pynanalogue as pn
print(pn.simulate_mod_bam.__doc__)
```

### Example

This example generates two contigs with random DNA sequences with the given properties and stores
them in a fasta file. Then, it generates a BAM file with the modification pattern and other properties
as shown. In this example, the reads are methylated, with 5Cs with a probability drawn randomly in the
range 30-70% and three Cs with a probability in the range 10%-50%. This 5C, 3C pattern repeats throughout
the read. You can then analyze this pattern with your tool and test its functionality.
You can set up multiple modifications etc. Please have a look at the documentation
[here](https://docs.rs/nanalogue/latest/nanalogue_core/simulate_mod_bam/index.html) for the options
available in the json configuration.

```python
import pynanalogue

json_config = '''
{
"contigs": {
    "number": 2,
    "len_range": [100, 200],
    "repeated_seq": "ACGTACGT"
},
"reads": [
    {
        "number": 10,
        "mapq_range": [10, 30],
        "base_qual_range": [20, 40],
        "len_range": [0.1, 0.9],
        "barcode": "ACGTAA",
        "mods": [{
            "base": "C",
            "is_strand_plus": true,
            "mod_code": "m",
            "win": [5, 3],
            "mod_range": [[0.3, 0.7], [0.1, 0.5]]
        }]
    }
]
}
'''

pynanalogue.simulate_mod_bam(
json_config=json_config,
bam_path="output.bam",
fasta_path="output.fasta"
)
```

# Further documentation

In addition to this repository, we are developing a
companion cookbook [here](https://www.nanalogue.com).

# Changelog

For a detailed list of changes in each version, please see CHANGELOG.md in the repository.

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
