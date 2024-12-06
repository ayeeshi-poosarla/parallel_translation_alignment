from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from pyspark.sql import SparkSession
from pyspark.broadcast import Broadcast
from typing import List, Tuple, Dict
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

# Get the top 3 k-mers with the highest number of matches
def get_top_kmers(kmer_map: dict, top_n=3) -> List[Tuple[str, int]]:
    if not kmer_map:
        return []
    # Sort k-mers by the number of matches (descending order)
    sorted_kmers = sorted(kmer_map.items(), key=lambda x: len(x[1]), reverse=True)
    return sorted_kmers[:top_n]

# Rank proteins by the number of k-mer matches with the reference
def rank_proteins_by_matches(reference: str, proteins: List[str], k: int) -> List[str]:
    protein_match_data = []
    for protein in proteins:
        kmer_map = build_kmer_index(reference, protein, k)
        top_kmers = get_top_kmers(kmer_map)  # Get top 3 k-mers
        max_matches = sum(len(positions) for _, positions in top_kmers)  # Sum of top k-mer match counts
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

# Function to read FASTA/FASTQ file and return sequences as dictionary
def read_sequences_from_file(file_path: str, file_format: str = "fasta") -> Dict[str, str]:
    sequences = {}
    for record in SeqIO.parse(file_path, file_format):
        sequences[record.id] = str(record.seq)
    return sequences

# Main function to process and rank proteins against multiple genomes
def rank_proteins_for_multiple_references(reference_genomes: List[str], reference_names: List[str], protein_sequences: List[str], protein_names: List[str], k: int) -> Dict[str, Tuple[str, int, Tuple[str, str]]]:
    best_scores = {}
    best_alignments = {}
    best_query_sequences = {}

    # Broadcasting the protein sequences to workers
    protein_sequences_broadcast: Broadcast[List[str]] = spark.sparkContext.broadcast(protein_sequences)

    # Define a function to process each reference genome
    def process_reference(reference_genome, reference_name):
        best_score = float('-inf')
        best_alignment = None
        best_query_sequence = None

        # Translate the reference genome into protein sequence
        translated_reference = translate_genomic_sequence(reference_genome)
        print(f"Translated Reference Protein Sequence for {reference_name}: {translated_reference}")

        # Rank proteins by matches for this reference genome
        ranked_proteins = rank_proteins_by_matches(translated_reference, protein_sequences_broadcast.value, k)
        top_2 = ranked_proteins[0:3]  # Take the top 3 proteins

        # For each of the top proteins, align them with the reference
        for protein in top_2:
            aligned_ref, aligned_query, alignment_score = smith_waterman(translated_reference, protein)
            
            print(f"\nAlignment for {protein} (Query):")
            print("Reference:", aligned_ref)
            print("Query:    ", aligned_query)
            print(f"Alignment Score: {alignment_score}")

            # Update the best protein if necessary
            if alignment_score > best_score:
                best_score = alignment_score
                best_alignment = (aligned_ref, aligned_query)
                best_query_sequence = protein

        return reference_name, best_query_sequence, best_score, best_alignment

    # Use PySpark to process each reference genome in parallel
    reference_rdd = spark.sparkContext.parallelize(zip(reference_genomes, reference_names))
    results = reference_rdd.map(lambda x: process_reference(x[0], x[1])).collect()

    # Store the results for each reference genome
    for reference_name, best_protein, best_score, alignment in results:
        best_scores[reference_name] = best_score
        best_alignments[reference_name] = alignment
        best_query_sequences[reference_name] = best_protein

    return best_query_sequences, best_scores, best_alignments

# Example of how to use the function with multiple reference genomes
reference_genomes_file = "C:/Users/ishit/CompGenom/project/parallel_translation_alignment/distributed_computing/official_tests_keep/multiple_match/reference.fastq"  # Path to reference genomes file
protein_file = "C:/Users/ishit/CompGenom/project/parallel_translation_alignment/distributed_computing/official_tests_keep/multiple_match/proteins.fasta"  # Path to protein sequences file

# Read the reference genomes and protein sequences from files
reference_genomes = read_sequences_from_file(reference_genomes_file, "fasta")
reference_names = list(reference_genomes.keys())

protein_sequences = read_sequences_from_file(protein_file, "fasta")
protein_names = list(protein_sequences.keys())

k = 4

best_proteins, best_scores, best_alignments = rank_proteins_for_multiple_references(
    list(reference_genomes.values()), reference_names, list(protein_sequences.values()), protein_names, k
)

# Print out results for each reference genome
for genome_name in best_proteins.keys():
    print(f"Best protein for reference genome {genome_name}: {best_proteins[genome_name]}, Score: {best_scores[genome_name]}")
    print(f"Alignment: {best_alignments[genome_name]}")

#throwing memory issues due to spark/pyspark 