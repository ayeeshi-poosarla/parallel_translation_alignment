from Bio.Seq import Seq
from Bio import pairwise2
from pyspark.sql import SparkSession
from pyspark.broadcast import Broadcast
from typing import List, Tuple
import numpy as np
import time
from sklearn.metrics import f1_score

# Function to translate a genomic DNA sequence into a protein sequence
def translate_genomic_sequence(dna_sequence: str) -> str:
    """
    Translates a genomic DNA sequence into a protein sequence.
    
    Parameters:
        dna_sequence (str): The genomic DNA sequence to be translated.
        
    Returns:
        str: The resulting protein sequence.
    """
    dna_seq = Seq(dna_sequence)
    return str(dna_seq.translate())

# Function to generate k-mers from a sequence
def generate_kmers(sequence: str, k: int) -> List[str]:
    """
    Translates a genomic DNA sequence into a protein sequence.
    
    Parameters:
        dna_sequence (str): The genomic DNA sequence to be translated.
        
    Returns:
        str: The resulting protein sequence.
    """
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

# Build a hashmap of k-mers in the reference sequence to the positions they match in the protein sequence
def build_kmer_index(reference: str, protein: str, k: int) -> dict:
    """
    Builds a hashmap of k-mers in the reference sequence and the positions they match in the protein sequence.
    
    Parameters:
        reference (str): The reference sequence.
        protein (str): The protein sequence to be indexed.
        k (int): The length of the k-mer.
        
    Returns:
        dict: A dictionary with k-mers as keys and their match positions as values.
    """
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
    """
    Gets the top k-mers with the highest number of matches from a k-mer map.
    
    Parameters:
        kmer_map (dict): A dictionary of k-mers and their match positions.
        top_n (int): The number of top k-mers to return.
        
    Returns:
        List[Tuple[str, int]]: A list of the top k-mers and their match counts.
    """
    if not kmer_map:
        return []
    # Sort k-mers by the number of matches (descending order)
    sorted_kmers = sorted(kmer_map.items(), key=lambda x: len(x[1]), reverse=True)
    return sorted_kmers[:top_n]

# Rank proteins by the number of k-mer matches with the reference
def rank_proteins_by_matches(reference: str, proteins: List[str], k: int) -> List[str]:
    """
    Ranks proteins based on the number of k-mer matches with the reference sequence.
    
    Parameters:
        reference (str): The reference sequence.
        proteins (List[str]): A list of protein sequences to be ranked.
        k (int): The length of the k-mer.
        
    Returns:
        List[str]: A list of proteins ranked by the number of k-mer matches.
    """
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
    """
    Performs Smith-Waterman alignment between the reference and query sequences.
    
    Parameters:
        reference (str): The reference sequence.
        query (str): The query sequence.
        match (int): The score for a match (default 2).
        mismatch (int): The score for a mismatch (default -1).
        gap_penalty (int): The penalty for gaps (default -2).
        
    Returns:
        Tuple[str, str, int]: The aligned reference and query sequences, and the alignment score.
    """
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
            
            score = max(
                score,
                scoring_matrix[i - 1, j] + gap_penalty,
                scoring_matrix[i, j - 1] + gap_penalty,
                0
            )
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
spark = SparkSession.builder \
    .appName("GenomicAnalysis") \
    .config("spark.driver.memory", "8g") \
    .config("spark.executor.memory", "8g") \
    .getOrCreate()


# Function to calculate F1 score (Assuming ground truth is available)
def calculate_f1_score(predicted: List[str], true: List[str]) -> float:
    """
    Calculates the F1 score by comparing the predicted protein sequences with the true protein sequences.
    
    Parameters:
        predicted (List[str]): The list of predicted protein sequences.
        true (List[str]): The list of true protein sequences.
        
    Returns:
        float: The F1 score.
    """
    # Here, we are using a simplistic assumption of comparing predicted proteins to true proteins
    return f1_score(true, predicted, average='micro')

# Main function to process and rank proteins against multiple genomes
def rank_proteins_for_multiple_references(reference_genomes: List[str], protein_sequences: List[str], k: int) -> List[str]:
    """
    Processes and ranks proteins against multiple reference genomes using k-mer matching and Smith-Waterman alignment.
    
    Parameters:
        reference_genomes (List[str]): A list of reference genome sequences.
        protein_sequences (List[str]): A list of protein sequences to rank.
        k (int): The length of the k-mer.
        
    Returns:
        Tuple[dict, dict, dict, float, float, float]: Best query sequences, best scores, best alignments, sequential processing time, parallel processing time, and F1 score.
    """
    best_scores = {}
    best_alignments = {}
    best_query_sequences = {}

    # Broadcasting the protein sequences to workers
    protein_sequences_broadcast: Broadcast[List[str]] = spark.sparkContext.broadcast(protein_sequences)

    # Define a function to process each reference genome
    def process_reference(reference_genome):
        
        best_score = float('-inf')
        best_alignment = None
        best_query_sequence = None

        # Translate the reference genome into protein sequence
        translated_reference = translate_genomic_sequence(reference_genome)
        print(f"Translated Reference Protein Sequence: {translated_reference}")

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

        return reference_genome, best_query_sequence, best_score, best_alignment

    # Sequential processing (for timing comparison)
    start_time = time.time()
    results_sequential = [process_reference(genome) for genome in reference_genomes]
    sequential_time = time.time() - start_time

    # Use PySpark to process each reference genome in parallel
    reference_rdd = spark.sparkContext.parallelize(reference_genomes)
    start_time = time.time()
    results_parallel = reference_rdd.map(process_reference).collect()
    parallel_time = time.time() - start_time

    # Store the results for each reference genome
    for reference_genome, best_protein, best_score, alignment in results_parallel:
        best_scores[reference_genome] = best_score
        best_alignments[reference_genome] = alignment
        best_query_sequences[reference_genome] = best_protein

    # Assume you have true labels for F1 score calculation (these would be known in a real use case)
    true_labels = [
        "DTPGTDYECETLFSWNVTRR" ] 
    
    predicted_proteins = list(best_query_sequences.values())
    predicted_proteinsf1 = predicted_proteins[:len(true_labels)]
    f1 = calculate_f1_score(predicted_proteinsf1, true_labels)

    # Return the results along with timing and F1 score
    return best_query_sequences, best_scores, best_alignments, sequential_time, parallel_time, f1

# Copy from test_case 
reference_genomes = [
    "ATGGCCGCTGCTGCTGCTGCCGCCGCCGCTGCCGCCGCTGCCGCCGCCGCTGCCGCTGCTTAA"
]
    
protein_list = [
        "MPTGACDAKLVTRVQHGTSPY", "DTPGTDYECETLFSWNVTRR","MPTGTDYECSYTLFSWNVTRR" ] 


k = 4

best_proteins, best_scores, best_alignments, sequential_time, parallel_time, f1 = rank_proteins_for_multiple_references(reference_genomes, protein_list, k)

# Print out results for each reference genome
for genome, best_protein in best_proteins.items():
    print(f"Best protein for reference genome {genome}: {best_protein}, Score: {best_scores[genome]}")

# Print timing and F1 score
print(f"Parallel Processing Time: {sequential_time} seconds")
print(f"Sequential Processing Time: {parallel_time} seconds")
print(f"F1 Score: {f1}")
