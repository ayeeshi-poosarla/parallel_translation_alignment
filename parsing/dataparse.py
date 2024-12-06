import sys
from typing import Dict

if len(sys.argv) < 3:
    print("Usage: python3 script.py <sequence> <protein1> <protein2> ... <best_match.txt>")
    sys.exit(1)

def translate_dna_to_protein(dna_sequence: str) -> str:
    """
    Translates a DNA sequence into a protein sequence based on the standard genetic code.
    
    Args:
        dna_sequence (str): A string representing a DNA sequence composed of the bases
                            'A', 'T', 'C', and 'G'. The sequence length should ideally
                            be a multiple of 3 (codons).
    
    Returns:
        str: A string representing the translated protein sequence, where each codon
             is converted to its corresponding amino acid. Translation halts at the first
             stop codon ('TAA', 'TAG', 'TGA').
    
    Notes:
        - Non-multiple-of-three sequences are truncated from the end to ensure codon alignment.
        - Stop codons ('*') are not included in the returned protein sequence.
    """
    
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