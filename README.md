# Comparative_Genomic_Analysis_Identifying_and_Characterizing_-SNPs-_in_DNA_Sequences

This project focuses on identifying and analyzing Single Nucleotide Polymorphisms (SNPs) between two DNA sequences using Biopython. SNPs are important for understanding genetic diversity, disease susceptibility, and evolutionary relationships.

## Features

- **Sequence Reading**: 
  - Load DNA sequences from FASTA files using `SeqIO` from Biopython.

- **Sequence Alignment**: 
  - Align two DNA sequences with `PairwiseAligner` to ensure corresponding nucleotides are matched.

- **SNP Identification**: 
  - Detect and record positions where nucleotides differ between aligned sequences, indicating SNPs.

- **SNP Analysis**: 
  - **SNP Distribution Plot**: Visualize SNP distribution along the sequence.
  - **SNP Frequency Calculation**: Calculate the frequency of SNP changes (e.g., A->G, C->T).
  - **Planned Analyses**: 
    - **SNP Density Calculation**: Measure SNP density along sequences.
    - **Transition/Transversion Ratio**: Compute the ratio of transition to transversion mutations.
    - **Motif Identification**: Analyze sequence motifs around SNP regions.

## Installation

To install the necessary dependencies, run:
```bash
pip install biopython
```

## Usage

## 1. Load DNA Sequences:
```python
from Bio import SeqIO
seq1 = SeqIO.read("sequence1.fasta", "fasta")
seq2 = SeqIO.read("sequence2.fasta", "fasta")
```

## 2. Align Sequences:

```python
from Bio.Align import PairwiseAligner
aligner = PairwiseAligner()
alignment = aligner.align(seq1.seq, seq2.seq)
```

## 3. Identify SNPs:
```python
snps = [(i, alignment[0][i], alignment[1][i]) for i in range(len(alignment[0])) if alignment[0][i] != alignment[1][i]]
```

## 4. Analyze SNPs:
Plot SNP Distribution: Use Matplotlib or another plotting library to visualize SNPs along the sequence.
Calculate SNP Frequency: Count occurrences of each SNP type.

## Applications:
Genetic Diversity Studies: Identify variations between species or populations.
Disease Susceptibility: Understand genetic predispositions to certain diseases.
Evolutionary Research: Explore evolutionary relationships through genetic differences.

## License:
This project is licensed under the MIT License.

## Contact
For questions, contact [Yasmine Elnasharty](elnashartyasmine@gmail.com).
