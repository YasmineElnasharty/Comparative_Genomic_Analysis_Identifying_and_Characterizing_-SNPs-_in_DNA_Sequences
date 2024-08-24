# Comparative_Genomic_Analysis_Identifying_and_Characterizing_-SNPs-_in_DNA_Sequences

This project identifies and analyzes Single Nucleotide Polymorphisms (SNPs) between two DNA sequences using Biopython.

**Features**
Sequence Alignment: Align DNA sequences to match nucleotides.
SNP Detection: Identify differences (SNPs) between aligned sequences.
SNP Analysis:
   Plot SNP distribution
   Calculate SNP frequency
   (Planned) Analyze SNP density and transition/transversion ratios

**Requirements**
Python 3.x
Biopython (pip install biopython)
Quick Start
Load Sequences:

**python**
from Bio import SeqIO
seq1 = SeqIO.read("seq1.fasta", "fasta")
seq2 = SeqIO.read("seq2.fasta", "fasta")
Align Sequences:

**python**
from Bio.Align import PairwiseAligner
aligner = PairwiseAligner()
alignment = aligner.align(seq1.seq, seq2.seq)
Identify SNPs:

**python**
snps = [(i, alignment[0][i], alignment[1][i]) for i in range(len(alignment[0])) if alignment[0][i] != alignment[1][i]]
Applications
Genetic diversity studies
Disease susceptibility analysis
Evolutionary research

**License**
MIT License
