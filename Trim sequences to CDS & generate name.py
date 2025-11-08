#install seqio if not already installed as well as biopython

from Bio import SeqIO
import os

# === USER INPUTS ===
# Navigate to Assignment 1 folder from Tools folder
assignment1_folder = os.path.join("C:\\Users\\Crach\\OneDrive - The University of Manchester\\Bioinformatics\\Bioinformatics unit\\Assignment 1")
input_fasta = os.path.join(assignment1_folder, "Homo sapiens BEN domain containing 7 (BEND7), transcript variant 3, mRNA.fasta")   # your FASTA filename
cds_start = 482                 # start coordinate (1-based)
cds_end = 1723                   # end coordinate (1-based)
output_fasta = ""  # leave empty for automatic naming, or specify a filename
if not output_fasta:
    # Automatically generate output filename from input filename
    input_name = os.path.splitext(os.path.basename(input_fasta))[0]
    output_fasta = os.path.join(assignment1_folder, f"{input_name}_CDS.fa")

# === PROCESS ===
record = SeqIO.read(input_fasta, "fasta")

# Convert to 0-based indexing for Python slicing
cds_seq = record.seq[cds_start - 1: cds_end]

# === WRITE NEW FASTA ===
record.seq = cds_seq
record.id = record.id + "_CDS"
record.description = f"CDS from {cds_start}..{cds_end}"
SeqIO.write(record, output_fasta, "fasta")

print(f"CDS extracted: {len(cds_seq)} bp written to {output_fasta}")
print(f"First 30 bases: {cds_seq[:30]}")
print(f"Last 30 bases:  {cds_seq[-30:]}")
