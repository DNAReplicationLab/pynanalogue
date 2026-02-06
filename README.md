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
- [Functions](#functions)
  - [Peek](#peek)
    - [Documentation](#documentation)
    - [Sample input and output](#sample-input-and-output)
  - [Read info](#read-info)
    - [Documentation](#documentation-1)
    - [Sample input and output](#sample-input-and-output-1)
  - [Window reads](#window-reads)
    - [Documentation](#documentation-2)
    - [Sample input and output](#sample-input-and-output-2)
    - [Gradient mode](#gradient-mode)
  - [Seq table](#seq-table)
    - [Documentation](#documentation-3)
    - [Sample input and output](#sample-input-and-output-3)
  - [Polars bam mods](#polars-bam-mods)
    - [Documentation](#documentation-4)
    - [Sample input and output](#sample-input-and-output-4)
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

## Peek

Quickly extract BAM file metadata without processing all records.
This function returns information about contigs (reference sequences) and
modifications present in the BAM file, making it useful for understanding
the structure of your data before running more intensive analyses.

### Documentation

```python
import pynanalogue as pn
print(pn.peek.__doc__)
```

### Sample input and output

A sample execution and output follows.

<!-- TEST CODE: START peek -->
```python
import pynanalogue as pn
metadata = pn.peek("tests/data/examples/example_1.bam")
print(metadata)
```
<!-- TEST CODE: END peek -->

The output is a dictionary with two keys: `contigs` and `modifications`.

<!-- TEST OUTPUT: START peek -->
```python
{'contigs': {'dummyI': 22, 'dummyII': 48, 'dummyIII': 76}, 'modifications': [['G', '-', '7200'], ['T', '+', 'T']]}
```
<!-- TEST OUTPUT: END peek -->

The `contigs` dictionary maps contig names to their lengths.
The `modifications` list contains modification information as `[base, strand, code]` where
`+` indicates the basecalled strand and `-` indicates its complement.

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

A sample execution and output follows.
You will get one record per alignment.

<!-- TEST CODE: START read_info -->
```python
import pynanalogue as pn
import json
result_bytes = pn.read_info("tests/data/examples/example_1.bam")
decoded_output = json.loads(result_bytes)
print(json.dumps(decoded_output, indent=2))
```
<!-- TEST CODE: END read_info -->

<!-- TEST OUTPUT: START read_info -->
```json
[
  {
    "read_id": "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
    "sequence_length": 8,
    "contig": "dummyI",
    "reference_start": 9,
    "reference_end": 17,
    "alignment_length": 8,
    "alignment_type": "primary_forward",
    "mod_count": "T+T:0;(probabilities >= 0.5020, PHRED base qual >= 0)"
  },
  {
    "read_id": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
    "sequence_length": 48,
    "contig": "dummyIII",
    "reference_start": 23,
    "reference_end": 71,
    "alignment_length": 48,
    "alignment_type": "primary_forward",
    "mod_count": "T+T:3;(probabilities >= 0.5020, PHRED base qual >= 0)"
  },
  {
    "read_id": "fffffff1-10d2-49cb-8ca3-e8d48979001b",
    "sequence_length": 33,
    "contig": "dummyII",
    "reference_start": 3,
    "reference_end": 36,
    "alignment_length": 33,
    "alignment_type": "primary_reverse",
    "mod_count": "T+T:1;(probabilities >= 0.5020, PHRED base qual >= 0)"
  },
  {
    "read_id": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
    "sequence_length": 48,
    "alignment_type": "unmapped",
    "mod_count": "G-7200:0;T+T:3;(probabilities >= 0.5020, PHRED base qual >= 0)"
  }
]
```
<!-- TEST OUTPUT: END read_info -->

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

A sample execution and output follows.

<!-- TEST CODE: START window_reads -->
```python
import pynanalogue as pn
df = pn.window_reads("tests/data/examples/example_1.bam", win=2, step=1)
print(df.write_csv(separator='\t'), end='')
```
<!-- TEST CODE: END window_reads -->

The output is a polars dataframe printed as TSV.
(This was generated from a file without basecalling quality information, which is why we show 255s here).

<!-- TEST OUTPUT: START window_reads -->
```text
contig	ref_win_start	ref_win_end	read_id	win_val	strand	base	mod_strand	mod_type	win_start	win_end	basecall_qual
dummyI	9	13	5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	0.0	+	T	+	T	0	4	255
dummyI	12	14	5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	0.0	+	T	+	T	3	5	255
dummyI	13	17	5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	0.0	+	T	+	T	4	8	255
dummyIII	26	32	a4f36092-b4d5-47a9-813e-c22c3b477a0c	1.0	+	T	+	T	3	9	255
dummyIII	31	51	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	+	T	+	T	8	28	255
dummyIII	50	63	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	T	+	T	27	40	255
dummyIII	62	71	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	+	T	+	T	39	48	255
dummyII	15	17	fffffff1-10d2-49cb-8ca3-e8d48979001b	0.0	-	T	+	T	12	14	255
dummyII	16	20	fffffff1-10d2-49cb-8ca3-e8d48979001b	0.0	-	T	+	T	13	17	255
dummyII	19	23	fffffff1-10d2-49cb-8ca3-e8d48979001b	0.0	-	T	+	T	16	20	255
dummyII	22	24	fffffff1-10d2-49cb-8ca3-e8d48979001b	0.5	-	T	+	T	19	21	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	.	G	-	7200	28	30	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	.	G	-	7200	29	31	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	.	G	-	7200	30	33	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	.	G	-	7200	32	44	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	.	G	-	7200	43	45	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	1.0	.	T	+	T	3	9	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	.	T	+	T	8	28	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	.	T	+	T	27	40	255
.	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	.	T	+	T	39	48	255
```
<!-- TEST OUTPUT: END window_reads -->

### Gradient mode

The `window_reads` function supports a `win_op` parameter that controls the windowing operation.
By default, `win_op="density"` reports the modification density within each window.
Setting `win_op="grad_density"` instead reports the gradient (slope) of modification density
within each window, which can be useful for detecting transitions in modification patterns.

<!-- TEST CODE: START gradient_mode -->
```python
import pynanalogue as pn
df = pn.window_reads(
    "tests/data/examples/example_10.bam",
    win=10,
    step=1,
    win_op="grad_density"
)
print(df.write_csv(separator='\t'), end='')
```
<!-- TEST CODE: END gradient_mode -->

The output format is identical to the density mode, but the `win_val` column now contains
the gradient value instead of the density value.

<!-- TEST OUTPUT: START gradient_mode -->
```text
contig	ref_win_start	ref_win_end	read_id	win_val	strand	base	mod_strand	mod_type	win_start	win_end	basecall_qual
dummyIII	23	33	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	0	10	255
dummyIII	24	34	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	1	11	255
dummyIII	25	35	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	2	12	255
dummyIII	26	36	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	3	13	255
dummyIII	27	37	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	4	14	255
dummyIII	28	38	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	5	15	255
dummyIII	29	39	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	6	16	255
dummyIII	30	40	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	7	17	255
dummyIII	31	41	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	8	18	255
dummyIII	32	42	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	9	19	255
dummyIII	33	43	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	10	20	255
dummyIII	34	44	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	11	21	255
dummyIII	35	45	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	12	22	255
dummyIII	36	46	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	13	23	255
dummyIII	37	47	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	14	24	255
dummyIII	38	48	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	15	25	255
dummyIII	39	49	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	16	26	255
dummyIII	40	50	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.054545455	+	N	+	N	17	27	255
dummyIII	41	51	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.096969694	+	N	+	N	18	28	255
dummyIII	42	52	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.12727273	+	N	+	N	19	29	255
dummyIII	43	53	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.14545454	+	N	+	N	20	30	255
dummyIII	44	54	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.15151516	+	N	+	N	21	31	255
dummyIII	45	55	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.14545454	+	N	+	N	22	32	255
dummyIII	46	56	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.12727273	+	N	+	N	23	33	255
dummyIII	47	57	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.096969694	+	N	+	N	24	34	255
dummyIII	48	58	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.054545455	+	N	+	N	25	35	255
dummyIII	49	59	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	26	36	255
dummyIII	50	60	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	27	37	255
dummyIII	51	61	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	28	38	255
dummyIII	52	62	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	29	39	255
dummyIII	53	63	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	30	40	255
dummyIII	54	64	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	31	41	255
dummyIII	55	65	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	32	42	255
dummyIII	56	66	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	33	43	255
dummyIII	57	67	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	34	44	255
dummyIII	58	68	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	35	45	255
dummyIII	59	69	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	36	46	255
dummyIII	60	70	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	37	47	255
dummyIII	61	71	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	38	48	255
dummyIII	62	71	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	39	49	255
dummyIII	63	71	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	40	50	255
dummyIII	64	71	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	41	51	255
dummyIII	65	71	c4f36092-b4d5-47a9-813e-c22c3b477a0c	0.0	+	N	+	N	42	52	255
```
<!-- TEST OUTPUT: END gradient_mode -->

## Seq table

Extract read sequences and base qualities for a genomic region as a Polars DataFrame.
This function is useful for retrieving sequences from a particular region
with modification information overlaid. If mods are not present or the modification
probabilities are not high enough, then they are not shown.
Insertions are shown as lowercase letters, deletions as periods (`.`),
and modified bases as `Z` (or `z` for modifications in insertions).

### Documentation

```python
import pynanalogue as pn
print(pn.seq_table.__doc__)
```

### Sample input and output

A sample execution follows. Note that the `region` parameter is required.

<!-- TEST CODE: START seq_table -->
```python
import pynanalogue as pn
df = pn.seq_table(
    "tests/data/examples/example_pynanalogue_1.bam",
    region="contig_00000:0-10"
).sort("read_id")
print(df.write_csv(separator='\t'), end='')
```
<!-- TEST CODE: END seq_table -->

The output is a Polars DataFrame with three columns: `read_id`, `sequence`, and `qualities`.

<!-- TEST OUTPUT: START seq_table -->
```text
read_id	sequence	qualities
0.dc09ae0d-6b6e-4cb2-b092-078f251a778e	AZGTAZGTAZ	20.20.20.20.20.20.20.20.20.20
1.cb098e1d-26d6-4e14-b979-b089e492c068	ACGTACGTAC	30.30.30.30.30.30.30.30.30.30
```
<!-- TEST OUTPUT: END seq_table -->

Sequence column conventions:
- Uppercase letters: bases aligned to reference
- Lowercase letters: inserted bases
- `.` (period): deleted bases
- `Z`: modified base on the reference
- `z`: modified base in an insertion

The qualities column contains period-separated base quality scores (0-255),
with 255 indicating a deleted position or unknown quality.

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

A sample execution and output follows.

<!-- TEST CODE: START polars_bam_mods -->
```python
import pynanalogue as pn
df = pn.polars_bam_mods("tests/data/examples/example_1.bam")
print(df.write_csv(separator='\t'), end='')
```
<!-- TEST CODE: END polars_bam_mods -->

The output is a polars dataframe printed as TSV.
Mod quality is a probability represented as a number between 0 and 255,
where 0 means not modified and 255 means modified with certainty.
This is how modification data is stored in the mod BAM format.

<!-- TEST OUTPUT: START polars_bam_mods -->
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
<!-- TEST OUTPUT: END polars_bam_mods -->

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

<!-- TEST CODE: NOOUTPUT simulate_mod_bam -->
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
<!-- TEST CODE: END simulate_mod_bam -->

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

# README example testing

Code examples in this README are automatically tested by `tests/test_readme_examples.py`.
HTML comment markers identify which code blocks to test and what output to expect.

The following marker types are used (shown without angle brackets to avoid parser interference;
in practice, wrap each marker in standard HTML comment delimiters i.e. `<` + `!--` ... `--` + `>`):

- `!-- TEST CODE: START my_example --` / `!-- TEST CODE: END my_example --` wraps a testable code block.
- `!-- TEST OUTPUT: START my_example --` / `!-- TEST OUTPUT: END my_example --` wraps the expected stdout.
- `!-- TEST CODE: NOOUTPUT my_example --` / `!-- TEST CODE: END my_example --` wraps code that is executed
  but has no expected output (e.g. it just verifies the code runs without error).

The marker name (e.g. `my_example` above) is a plain identifier that links a code block
to its output block. Each tested code block must include `print()` calls that produce
exactly the text shown in the corresponding output block.

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
