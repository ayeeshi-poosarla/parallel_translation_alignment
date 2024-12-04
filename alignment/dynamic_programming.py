import sys
from typing import Dict

if len(sys.argv) < 3:
    print("Usage: python3 script.py <reference> <seq1> <seq2> ... ")
    sys.exit(1)

def read_Fastq(filename):
    sequences = []
    with open(filename, 'r') as file:
        while True:
            file.readline()
            seq = file.readline().strip()
            file.readline()
            file.readline()
            if len(seq) == 0:
                break
            sequences.append(seq)
    return sequences

def translate_dna_to_protein(dna_sequence: str) -> str:
    # Genetic code dictionary
    codon_table: Dict[str, str] = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

    # Ensure the sequence length is a multiple of 3
    n = len(dna_sequence)
    if n % 3 != 0:
        dna_sequence = dna_sequence[:n - (n % 3)]
    
    # Translate the sequence
    protein_sequence = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        if codon in codon_table:
            amino_acid = codon_table[codon]
            if amino_acid == '*':  # Stop codon
                break
            protein_sequence.append(amino_acid)
    
    return ''.join(protein_sequence)

def global_alignment(seq1, seq2, match=1, mismatch=-1, gap=-2):
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize DP table
    for i in range(1, m + 1):
        dp[i][0] = i * gap
    for j in range(1, n + 1):
        dp[0][j] = j * gap

    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                score = match
            else:
                score = mismatch
            dp[i][j] = max(dp[i - 1][j - 1] + score, dp[i - 1][j] + gap, dp[i][j - 1] + gap)

    # Final similarity score
    return dp[m][n]

reference = read_Fastq(sys.argv[1])
sequences = [read_Fastq(filename) for filename in sys.argv[2:]]