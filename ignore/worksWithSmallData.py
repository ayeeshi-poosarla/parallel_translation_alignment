from Bio.Seq import Seq
from Bio import pairwise2
from typing import List, Tuple, Dict
import numpy as np

# Function to translate a genomic DNA sequence into a protein sequence
def translate_genomic_sequence(dna_sequence: str) -> str:
    dna_seq = Seq(dna_sequence)
    return str(dna_seq.translate())

# Function to compare two protein sequences using global pairwise alignment
def compare_protein_sequences(protein1: str, protein2: str) -> int:
    alignments = pairwise2.align.globalxx(protein1, protein2)
    best_alignment = max(alignments, key=lambda x: x[2])
    return best_alignment[2]

# Function to find the best matching protein name from a protein database
def find_best_matching_protein(genomic_sequence: str, protein_database: List[Tuple[str, str]]) -> str:
    translated_protein = translate_genomic_sequence(genomic_sequence)
    best_score, best_match_name = -1, None
    for header, protein in protein_database:
        score = compare_protein_sequences(translated_protein, protein)
        if score > best_score:
            best_score, best_match_name = score, header
    return best_match_name

# Function to parse genomic/protein sequences from a FASTA file
def get_genomic_sequence(filename: str) -> List[Tuple[str, str]]:
    sequences = []
    with open(filename, 'r') as file:
        header = None
        sequence = []
        
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # FASTA header line
                if header:
                    sequences.append((header, ''.join(sequence)))
                header = line[1:]  # Remove the '>' character
                sequence = []
            else:
                sequence.append(line)
        
        # Append the last sequence
        if header:
            sequences.append((header, ''.join(sequence)))
    
    return sequences

# Function to parse protein sequences (alias for get_genomic_sequence)
def get_proteins(filename: str) -> List[Tuple[str, str]]:
    return get_genomic_sequence(filename)

# Example usage
if __name__ == "__main__":
    genomic_file = "gene.fna"  # Replace with your genomic FASTA file
    protein_file = "shot_uniprot.fasta"  # Replace with your protein FASTA file
    
    # Get genomic sequences
    genomic_sequences = get_genomic_sequence(genomic_file)
    print("Genomic Sequences:")
    for header, seq in genomic_sequences:
        print(f"Header: {header}, Sequence: {seq[:50]}...")  # Truncate for readability
    
    # Get protein sequences
    protein_database = get_proteins(protein_file)
    print("\nProtein Sequences:")
    for header, seq in protein_database:
        print(f"Header: {header}, Sequence: {seq[:50]}...")  # Truncate for readability

    # Find the best matching protein for each genomic sequence
    print("\nBest Matching Proteins:")
    for header, sequence in genomic_sequences:
        best_match = find_best_matching_protein(sequence, protein_database)
        print(f"Genomic Header: {header}, Best Protein Match: {best_match}")
