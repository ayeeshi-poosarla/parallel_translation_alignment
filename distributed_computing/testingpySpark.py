from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from pyspark.sql import SparkSession
from pyspark.broadcast import Broadcast
from typing import List, Tuple
import numpy as np
import csv

# Function to translate a genomic DNA sequence into a protein sequence
def translate_genomic_sequence(dna_sequence: str) -> str:
    dna_seq = Seq(dna_sequence)
    return str(dna_seq.translate())

# Function to generate k-mers from a sequence
def generate_kmers(sequence: str, k: int) -> List[str]:
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

# Build a hashmap of k-mers in the reference sequence to the positions they match in the protein sequence
def build_kmer_index(reference: str, protein: str, k: int) -> dict:
    reference_kmers = generate_kmers(reference, k)
    protein_kmers = generate_kmers(protein, k)
    
    kmer_map = {}
    for ref_kmer in reference_kmers:
        match_positions = [i for i, protein_kmer in enumerate(protein_kmers) if protein_kmer == ref_kmer]
        if match_positions:
            kmer_map[ref_kmer] = match_positions
    
    return kmer_map

# Get the top k-mer with the highest number of matches
def get_top_kmer_match(kmer_map: dict) -> Tuple[str, int]:
    if not kmer_map:
        return None, 0
    top_kmer = max(kmer_map, key=lambda k: len(kmer_map[k]))
    return top_kmer, len(kmer_map[top_kmer])

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

# Read reference genomes from a FASTQ file and return gene names and sequences
def read_reference_genomes(fastq_file: str) -> List[Tuple[str, str]]:
    references = []
    for record in SeqIO.parse(fastq_file, "fastq"):
        references.append((record.id, str(record.seq)))  # Store gene name (id) and sequence
    return references

# Read protein sequences from a FASTA file and return gene names and sequences
def read_protein_sequences(fasta_file: str) -> List[Tuple[str, str]]:
    proteins = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        proteins.append((record.id, str(record.seq)))  # Store gene name (id) and sequence
    return proteins

# Main function to process and rank proteins using parallelism
def rank_proteins_for_reference(reference: str, reference_name: str, protein_sequences: List[Tuple[str, str]], k: int) -> List[str]:
    best_score = float('-inf')
    best_alignment = None
    best_query_sequence = None
    best_query_name = None

    # Translate the reference genome into protein sequence
    translated_reference = translate_genomic_sequence(reference)
    print(f"Translated Reference Protein Sequence: {translated_reference}")

    # Broadcast the translated reference to all workers in the cluster
    broadcast_ref = spark.sparkContext.broadcast(translated_reference)

    # Parallelize protein sequences to distribute computation across workers
    protein_rdd = spark.sparkContext.parallelize(protein_sequences)

    # Define a function for ranking proteins on each worker
    def process_protein(protein: Tuple[str, str]) -> Tuple[str, int, str, Tuple[str, str]]:
        protein_name, protein_seq = protein
        # Rank proteins by the number of k-mer matches with the reference
        kmer_map = build_kmer_index(broadcast_ref.value, protein_seq, k)
        _, max_matches = get_top_kmer_match(kmer_map)
        
        # Perform Smith-Waterman alignment on the top k-mer protein sequence
        aligned_ref, aligned_query, alignment_score = smith_waterman(broadcast_ref.value, protein_seq)
        
        return protein_name, alignment_score, reference_name, (aligned_ref, aligned_query)

    # Apply the function to each protein in parallel
    results = protein_rdd.map(process_protein).collect()

    # Find the best protein with the highest alignment score
    for protein_name, alignment_score, ref_name, (aligned_ref, aligned_query) in results:
        print(f"\nAlignment for {protein_name} (Query):")
        print("Reference:", aligned_ref)
        print("Query:    ", aligned_query)
        print(f"Alignment Score: {alignment_score}")

        # Update the best protein if necessary
        if alignment_score > best_score:
            best_score = alignment_score
            best_alignment = (aligned_ref, aligned_query)
            best_query_sequence = protein_name
            best_query_name = ref_name

    # Returning the top-ranked protein based on the alignment score
    return best_query_name, best_query_sequence, best_score, best_alignment

# Output results to a CSV file
def output_to_csv(results: List[Tuple[str, int, str, Tuple[str, str]]], output_file: str):
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Reference_Gene_Name", "Protein_Gene_Name", "Alignment_Score"])
        for ref_gene, protein_gene, score, (aligned_ref, aligned_query) in results:
            writer.writerow([ref_gene, protein_gene, score])

# Example usage
reference_fastq = "C:/Users/ishit/CompGenom/project/parallel_translation_alignment/distributed_computing/official_tests_keep/multiple_match/reference.fastq"
protein_fasta = "C:/Users/ishit/CompGenom/project/parallel_translation_alignment/distributed_computing/official_tests_keep/multiple_match/proteins.fasta"
output_csv = "C:/Users/ishit/CompGenom/project/parallel_translation_alignment/distributed_computing/official_tests_keep/multiple_match/alignment2_results.csv"
k = 4

# Read data
reference_genomes = read_reference_genomes(reference_fastq)
protein_sequences = read_protein_sequences(protein_fasta)

# Process each reference genome
all_results = []
for reference_name, reference in reference_genomes:
    best_ref_name, best_protein_name, best_score, alignment = rank_proteins_for_reference(reference, reference_name, protein_sequences, k)
    all_results.append((best_ref_name, best_protein_name, best_score, alignment))

# Output results to CSV
output_to_csv(all_results, output_csv)
print(f"Results have been saved to {output_csv}")
