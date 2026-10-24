# HIV-1 Gag-Nucleocapsid Mutation Analysis

## Overview

This repository contains the R script used to process HIV-1
Gag-nucleocapsid (NC) consensus sequences for the study:

> **Non-canonical resistance-associated mutations in HIV-1
> Gag-nucleocapsid are associated with dolutegravir pressure**

The script reads nucleotide consensus sequences in FASTA format, translates
codons containing ambiguous IUPAC nucleotide codes, and identifies amino-acid
substitutions across positions 1–55 of the HIV-1 NC region.

Detected amino-acid substitutions are classified into three mutually exclusive
categories:

1. **Non-polymorphic NC resistance-associated mutations (`np.NC-RAMs`)**:
   the predefined non-polymorphic NC mutations included as NC RAMs in the
   associated manuscript, based on evidence of phenotypic drug resistance.

2. **Polymorphic NC mutations (`pNC`)**:
   amino-acid substitutions classified as polymorphic in HIV-1 subtype C.

3. **Other mutations (`other_mutations`)**:
   non-polymorphic amino-acid substitutions for which there is no conclusive
   evidence of phenotypic drug resistance.


## Scientific context

Integrase strand transfer inhibitor (INSTI) resistance is usually assessed
using mutations in the HIV-1 polymerase gene. However, mutations in the HIV-1
Gag-NC zinc-finger and adjacent regions have also been associated with reduced
susceptibility to INSTIs, including dolutegravir.

This script was developed to identify and classify NC amino-acid substitutions
in consensus sequences from individuals with virological non-suppression while
receiving dolutegravir-based antiretroviral therapy.


## Mutation categories

### Non-polymorphic NC resistance-associated mutations (`np.NC-RAMs`)

The script identifies the following ten predefined non-polymorphic NC
resistance-associated mutations included in the associated manuscript:

- G22R
- T24V
- N27S
- G40R
- K41G
- M46I
- M46L
- D48E
- T50S
- E51G

In the associated study, NC RAMs were defined using two criteria:

1. a prevalence below 5% among 1,240 HIV-1 subtype C NC sequences from
   South Africa obtained from the Los Alamos National Laboratory sequence
   database between 1997 and 2017, before the introduction of dolutegravir; and
2. evidence of a statistically significant reduction in drug susceptibility
   in phenotypic testing.

These mutations are reported in the `np.NC-RAMs` output column.

### Polymorphic NC mutations (`pNC`)

Amino-acid substitutions classified as polymorphic in HIV-1 subtype C are
reported in the `pNC` output column.

These variants are distinguished from non-polymorphic mutations and are not
classified as NC resistance-associated mutations by this analysis.

### Other mutations (`other_mutations`)

Non-polymorphic amino-acid substitutions that are not among the predefined NC
RAMs are reported as `other_mutations`.

This category therefore represents non-polymorphic mutations without
conclusive evidence of phenotypic drug resistance.


## Repository contents

```text
.
├── README.md
├── NC_fasta_process.R
├── CITATION.cff
└── LICENSE
```

The repository does not contain participant-level sequence data.


## Requirements

The script requires:

- R
- Biostrings
- tidyverse

`Biostrings` is distributed through Bioconductor.

### Installing the required packages

Run the following commands in R if the packages are not already installed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("Biostrings")
install.packages("tidyverse")
```

The analysis script checks whether the required packages are available but
does not install packages automatically.


## Input data

The input must be a nucleotide FASTA file containing HIV-1 consensus
sequences that are aligned consistently with the NC reference positions used
by the script.

Example:

```text
>sequence_001
ATGCAGAGAGGCAATTTTAGGAATCAGAGAC...
>sequence_002
ATGCAGAGAGGCAATTTTAGGAATCAGAGAC...
```

FASTA header names are retained as sequence identifiers in the output.


## Sequence offset

The `offset_to_pos1` argument specifies the number of nucleotides preceding
NC amino-acid position 1 in each input sequence.

If the input sequences begin exactly at NC position 1, use:

```r
offset_to_pos1 = 0
```

If the input contains upstream sequence, the offset must be adjusted to the
number of nucleotides preceding NC position 1.

For example, if NC begins after 377 amino acids of upstream Gag sequence:

```r
offset_to_pos1 = 377 * 3
```

which is equivalent to:

```r
offset_to_pos1 = 1131
```

Incorrect specification of the offset will cause codons to be assigned to
incorrect NC positions.


## Running the analysis

Open `NC_fasta_process.R` in R or RStudio and run the script.

The following line opens a file-selection window:

```r
my_fasta_file <- file.choose()
```

After selecting the FASTA file, the analysis is run using:

```r
final_results <- process_nc_mutations(
  fasta_file = my_fasta_file,
  ref_data = reference_data,
  offset_to_pos1 = 0
)
```

Change `offset_to_pos1` if the sequences do not begin at NC position 1.

The main function can also be called directly:

```r
final_results <- process_nc_mutations(
  fasta_file = "path/to/NC_consensus_sequences.fasta",
  ref_data = reference_data,
  offset_to_pos1 = 0
)
```


## Handling ambiguous nucleotide codes

The script supports standard IUPAC ambiguous nucleotide codes.

For each ambiguous codon, it:

1. expands the ambiguous nucleotide code into all possible codons;
2. translates each possible codon;
3. retains all unique possible amino-acid translations; and
4. evaluates each possible amino acid against the reference definitions.

Codons containing an alignment gap are represented by `-`.

Codons containing `N` or `?` are represented by `X` and are not reported as
mutations. Stop codons are also excluded from the mutation output.


## Mutation classification

For each NC position, the script compares the translated amino acid with the
HXB2 reference amino acid.

A detected non-reference amino acid is classified using the following
hierarchy:

1. **`np.NC-RAMs`** – predefined non-polymorphic NC mutations
