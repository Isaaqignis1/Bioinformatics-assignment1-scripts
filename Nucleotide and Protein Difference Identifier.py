#input: two fasta files containing the CDS sequences of the reference and alternative alleles
#output: list of nucleotide differences and their corresponding protein changes (if any)
#as well as predicted HGVS nomenclature for each difference

from Bio import SeqIO  # Import SeqIO module from Biopython for reading sequence files
import os  # Import os module for operating system interface functions

# Change current working directory to the assignment folder
os.chdir("C:\\Users\\Crach\\OneDrive - The University of Manchester\\Bioinformatics\\Bioinformatics unit\\Assignment 1")

# ====================================
# user inputs
ref = "BEND7"
alt = "cDNA_3"
#===================================

REF_FA = ref + "_CDS.fa"  
ALT_FA = alt + "_CDS.fa"

# ==========================================
genetic_code = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def load_seq(path):
    rec = SeqIO.read(path, "fasta")
    return str(rec.seq).upper().replace("\n","").replace("\r","")

ref = load_seq(REF_FA)
alt = load_seq(ALT_FA)

print(f"[INFO] Loaded REF {REF_FA} length={len(ref)} bp")
print(f"[INFO] Loaded ALT {ALT_FA} length={len(alt)} bp")

if len(ref) != len(alt):
    print("[WARN] Lengths differ; you likely have an indel. Use your global alignment to localize the gap(s) and report c.start_enddel/ins.")
    # You can stop here if you want to handle indels manually.
    # Or proceed to a naive position-by-position diff on the common length:
    L = min(len(ref), len(alt))
else:
    L = len(ref)

diffs = [(i, ref[i], alt[i]) for i in range(L) if ref[i] != alt[i]]

print(f"[RESULT] Nucleotide mismatches found: {len(diffs)}")
if not diffs:
    print("[CALL] No coding difference (CDS identical).")
else:
    for i, r, a in diffs:
        c_pos = i + 1  # 1-based within CDS
        codon_start = (i // 3) * 3
        ref_codon = ref[codon_start:codon_start+3]
        alt_codon = alt[codon_start:codon_start+3]
        # Only translate if we still have complete codons
        if len(ref_codon) == 3 and len(alt_codon) == 3:
            ref_aa = genetic_code.get(ref_codon, '?')
            alt_aa = genetic_code.get(alt_codon, '?')
            aa_pos = (codon_start // 3) + 1
            if ref_aa == alt_aa:
                p_change = f"p.{ref_aa}{aa_pos}="
            else:
                p_change = f"p.{ref_aa}{aa_pos}{alt_aa}"
        else:
            ref_aa = alt_aa = '?'
            aa_pos = (codon_start // 3) + 1
            p_change = "p.(codon context incomplete)"

        print(f"[VAR] c.{c_pos}{r}>{a}  | codon {aa_pos}: {ref_codon}->{alt_codon}  | {p_change}")

#create HGVS nomenclature for each difference
        hgvs_c = f"c.{c_pos}{r}>{a}"
        hgvs_p = p_change
        print(f"     HGVS: {hgvs_c}, {hgvs_p}")


#print("\n[NOTE] Validate the final HGVS with Mutalyzer Name Checker using your transcript version (e.g., NM_017617.5:c.<pos><ref>><alt>).")
