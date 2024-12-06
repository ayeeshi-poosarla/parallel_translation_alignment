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

# Smith-Waterman alignment implementation
def smith_waterman(reference, query, match=2, mismatch=-1, gap_penalty=-2):
    len_ref, len_query = len(reference), len(query)
    scoring_matrix = np.zeros((len_ref + 1, len_query + 1))

    max_score = 0
    max_pos = (0, 0)
    
    for i in range(1, len_ref + 1):
        for j in range(1, len_query + 1):
            if reference[i - 1] == query[j - 1]:
                score = scoring_matrix[i - 1, j - 1] + match
            else:
                score = scoring_matrix[i - 1, j - 1] + mismatch
            
            score = max(score, scoring_matrix[i - 1, j] + gap_penalty, scoring_matrix[i, j - 1] + gap_penalty, 0)
            scoring_matrix[i, j] = score
            
            if score > max_score:
                max_score = score
                max_pos = (i, j)

    aligned_ref = []
    aligned_query = []
    
    i, j = max_pos
    while scoring_matrix[i, j] != 0:
        if scoring_matrix[i, j] == scoring_matrix[i - 1, j] + gap_penalty:
            aligned_ref.append(reference[i - 1])
            aligned_query.append('-')
            i -= 1
        elif scoring_matrix[i, j] == scoring_matrix[i, j - 1] + gap_penalty:
            aligned_ref.append('-')
            aligned_query.append(query[j - 1])
            j -= 1
        else:
            aligned_ref.append(reference[i - 1])
            aligned_query.append(query[j - 1])
            i -= 1
            j -= 1
    
    aligned_ref = ''.join(reversed(aligned_ref))
    aligned_query = ''.join(reversed(aligned_query))
    
    return aligned_ref, aligned_query, max_score

# Setup PySpark session
spark = SparkSession.builder.appName("GenomicAnalysis").getOrCreate()

# Main function to process and rank proteins
def rank_proteins_for_reference(reference: str, protein_sequences: List[str], k: int) -> List[str]:
    best_score = float('-inf')
    best_alignment = None
    best_query_sequence = None

    # Translate the reference genome into protein sequence
    translated_reference = translate_genomic_sequence(reference)
    print(f"Translated Reference Protein Sequence: {translated_reference}")

    ranked_proteins = rank_proteins_by_matches(translated_reference, protein_sequences, k)
    top_2 = ranked_proteins[0:3]

    for two in top_2:
        aligned_ref, aligned_query, alignment_score = smith_waterman(translated_reference, two)
        
        print(f"\nAlignment for {two} (Query):")
        print("Reference:", aligned_ref)
        print("Query:    ", aligned_query)
        print(f"Alignment Score: {alignment_score}")

        # Update the best protein if necessary
        if alignment_score > best_score:
            best_score = alignment_score
            best_alignment = (aligned_ref, aligned_query)
            best_query_sequence = two

    # Returning the top-ranked protein based on the alignment score
    return best_query_sequence, best_score, best_alignment

# Example of how to use the function
reference_genome = "ATGGCTGCTGCTGCCGCCGCTGCTGCCGCTGCTGCTGCTGCTGCTGCCGCTGCCGCTGCTTAA"
protein_list = ["MAAAAAAAAAAAAAAAAAAA*", "MAAAAAAATMPLAAAAAAAA*"]
k = 4

best_protein, best_score, alignment = rank_proteins_for_reference(reference_genome, protein_list, k)
print(f"\nBest protein: {best_protein}, Score: {best_score}")
