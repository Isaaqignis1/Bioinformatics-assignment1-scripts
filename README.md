# Bioinformatics assignment1 scripts

This repository contains Python scripts developed for the MSc Bioinformatics assignment on variant identification and annotation in candidate genes related to Tetralogy of Fallot. The scripts automate extraction of coding sequences (CDS), pairwise comparison of cDNA variants against reference transcripts, and annotation of coding changes in HGVS format.

## Requirements
- **Python:** Version 3.8 or higher  
- **Packages:**  
  - [Biopython](https://biopython.org/)  

## Scripts

### Trim sequences to CDS & generate name.py
Extracts the coding sequence (CDS) region from a provided mRNA FASTA file using RefSeq coordinates.

**Input:** mRNA FASTA file  
**Output:** Trimmed CDS FASTA (`*_CDS.fa`)  
**Main libraries:** Bio.SeqIO, os  

### Nucleotide and Protein Difference Identifier
Compares two nucleotide sequences, identifies aminoo acid change amd outputs information in HGVS format.

**Input:** two mRNA FASTA files 
**Output:** Terminal information on the amino acid change in HGVS notation
