# Simulation Configuration Files

This directory contains JSON configuration files for generating test BAM files using `simulate_mod_bam()`.

These configurations omit `seed`, so generated sequences and modification calls vary
between runs. Fixed workload dimensions make benchmarks comparable, not identical.
For repeatable simulated data with the same simulator version and configuration,
add a top-level integer `seed` to the JSON. Mapping-quality ranges below apply to
mapped reads; simulated unmapped reads have MAPQ 0.

## simple_bam.json

**Purpose**: Basic test fixture for most test cases

**Specifications**:
- **Contigs**: 2 contigs, each 10,000 bp long
- **Reads**: 1,000 reads total
  - Mapping quality: 10-20
  - Base quality: 10-20
  - Initial length: 50% of contig length (5,000 bp)
  - Insert sequence: "ATCG" in the middle of each read, giving a final length of 5,004 bp
- **Modifications**: Thymine (T) modifications
  - Base: T on plus strand
  - Modification code: T
  - Window size: 40T
  - Modification probability ranges: 10-20% and 70-80%

**Use cases**:
- Testing basic BAM reading functionality
- Validating modification extraction
- Testing filtering by sequence length, mapping quality, etc.
- Quick integration tests

## benchmark_bam.json

**Purpose**: Performance benchmarking with fixed workload dimensions

**Specifications**:
- **Contigs**: 2 contigs, each 1,000,000 bp (1 Mbp) long
- **Reads**: 5,000 reads total
  - Mapping quality: 10-20
  - Base quality: 10-20
  - Initial length: 1% of contig length (10,000 bp)
  - Insert sequence: "ATCG" in the middle of each read, giving a final length of 10,004 bp
- **Modifications**: Thymine (T) modifications
  - Base: T on plus strand
  - Modification code: T
  - Window size: 40T
  - Modification probability ranges: 10-20% and 70-80%

**Use cases**:
- Performance benchmarking with pytest-benchmark
- Timing comparisons across different hardware
- Regression testing for performance
- Demonstrating package speed to users
- Testing performance with realistic data volumes

## two_mods_bam.json

**Purpose**: Test fixture for BAM files containing two different modification types

**Specifications**:
- **Contigs**: 2 contigs, each 10,000 bp long
- **Reads**: 100 reads total
  - Mapping quality: 20-30
  - Base quality: 20-30
  - Initial length: 50% of contig length (5,000 bp)
  - Insert sequence: "ATCG" in the middle of each read, giving a final length of 5,004 bp
- **Modifications**: Two different modification types
  1. Thymine (T) modifications on minus strand
     - Base: T, Strand: minus
     - Modification code: T
     - Window size: 40T
     - Modification probability ranges: 10-20% and 70-80%
  2. Cytosine (C) modifications on plus strand
     - Base: C, Strand: plus
     - Modification code: 76792
     - Window size: 30C
     - Modification probability ranges: 10-20% and 40-80%

**Use cases**:
- Testing peek functionality with multiple modification types
- Validating modification detection across different strands
- Testing filtering and grouping by modification type
