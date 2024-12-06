from Bio.Seq import Seq
from Bio import pairwise2
from pyspark.sql import SparkSession
from pyspark.broadcast import Broadcast
from typing import List, Tuple
import numpy as np

# Function to translate a genomic DNA sequence into a protein sequence
def translate_genomic_sequence(dna_sequence: str) -> str:
    dna_seq = Seq(dna_sequence)
    return str(dna_seq.translate())

# Function to generate k-mers from a sequence
def generate_kmers(sequence: str, k: int) -> List[str]:
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

# Build a hashmap of k-mers in the reference sequence to the positions they match in the protein sequence
def build_kmer_index(reference: str, protein: str, k: int) -> dict:
    kmer_map = {}
    for i in range(len(reference) - k + 1):
        kmer = reference[i:i + k]
        match_positions = [
            j for j in range(len(protein) - k + 1) if protein[j:j + k] == kmer
        ]
        if match_positions:
            kmer_map[kmer] = match_positions
    print(f"Kmer Map for reference {reference}: {kmer_map}")  # Debugging line
    return kmer_map

# Get the top k-mer with the highest number of matches
def get_top_kmer_match(kmer_map: dict) -> Tuple[str, int]:
    if not kmer_map:
        return None, 0
    top_kmer = max(kmer_map, key=lambda k: len(kmer_map[k]))
    return top_kmer, len(kmer_map[top_kmer])

# Rank proteins by the number of k-mer matches with the reference
def rank_proteins_by_matches(reference: str, proteins: List[str], k: int) -> List[str]:
    protein_match_data = []
    for protein in proteins:
        kmer_map = build_kmer_index(reference, protein, k)
        _, max_matches = get_top_kmer_match(kmer_map)
        protein_match_data.append((protein, max_matches))

    # Sort proteins by the number of matches in descending order
    protein_match_data.sort(key=lambda x: x[1], reverse=True)
    return [protein for protein, _ in protein_match_data]

# BLOSUM62 matrix creation
def create_blosum62():
    blosum62_str = """
        A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
    A   4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
    R  -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
    N  -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
    D  -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
    C   0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
    Q  -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
    E  -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
    G   0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
    H  -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
    I  -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
    L  -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
    K  -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
    M  -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
    F  -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
    P  -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
    S   1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
    T   0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
    W  -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
    Y  -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
    V   0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
    """
    lines = blosum62_str.strip().split("\n")
    keys = lines[0].split()
    matrix = {}
    for i, line in enumerate(lines[1:]):
        values = list(map(int, line.split()[1:]))
        matrix[keys[i]] = dict(zip(keys, values))
    return matrix

# Compare two protein sequences using local alignment (Smith-Waterman)
def compare_protein_sequences(seq1, seq2) -> int:
    match_score = 1
    mismatch_score = 0
    gap_penalty = -1  # Optionally, add a gap penalty
    
    rows, cols = len(seq1) + 1, len(seq2) + 1
    score_matrix = np.zeros((rows, cols), dtype=int)
    max_score = 0
    max_i, max_j = 0, 0

    # Fill scoring matrix (dynamic programming approach)
    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i-1] == seq2[j-1]:  # Match case
                score = match_score
            else:  # Mismatch case
                score = mismatch_score

            # Calculate diagonal, up, and left scores (gap penalties can be applied here)
            diagonal = score_matrix[i-1][j-1] + score
            up = score_matrix[i-1][j] + gap_penalty
            left = score_matrix[i][j-1] + gap_penalty
            
            score_matrix[i][j] = max(0, diagonal, up, left)  # The max score should never go below 0

            if score_matrix[i][j] >= max_score:
                max_score = score_matrix[i][j]
                max_i, max_j = i, j

    # Reconstruct the optimal alignment from the scoring matrix
    aligned_seq1, aligned_seq2 = [], []
    i, j = max_i, max_j

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        if seq1[i-1] == seq2[j-1]:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif score_matrix[i-1][j] >= score_matrix[i][j-1]:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append("-")
            i -= 1
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j-1])
            j -= 1

    # Return the aligned sequences
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), max_score

# Setup PySpark session
spark = SparkSession.builder.appName("GenomicAnalysis").getOrCreate()

# Main function to process and rank proteins
def rank_proteins_for_reference(reference: str, protein_sequences: List[str], k: int) -> List[str]:
    """
    Given a reference protein sequence and a list of other protein sequences, 
    ranks proteins by the number of k-mer matches with the reference sequence.
    """
    ranked_proteins = rank_proteins_by_matches(reference, protein_sequences, k)
    return ranked_proteins

# Example usage:
reference_sequence = "MKTAYIAKQRQISFVKSHFSRQHFLASVHTKLNQ"
protein_sequences = [
    "MKTAYIAKQRQISFVKSHFSRQHFLASVHTKLNQ",
    "MKTAYIAKQRQISFVKHFSRQHFLASVHTKLNQ",  # Slight modification
    "MKTAHIQKRQISFVKSHFSRQHFLASVHTKLNQ"
]
k = 3  # K-mer length

# Rank proteins
ranked_proteins = rank_proteins_for_reference(reference_sequence, protein_sequences, k)
print("Ranked Proteins: ", ranked_proteins)
