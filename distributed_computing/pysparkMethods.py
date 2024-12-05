from Bio.Seq import Seq
from Bio import pairwise2
from pyspark.sql import SparkSession
from pyspark.broadcast import Broadcast
from typing import List, Tuple

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
def parse_fasta(filename: str) -> List[Tuple[str, str]]:
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

# Initialize PySpark
spark = SparkSession.builder \
    .appName("GenomicSequenceMatching") \
    .getOrCreate()

sc = spark.sparkContext

# Paths to files
genomic_file = "gene.fna"  # Replace with your genomic FASTA file
protein_file = "short_uniprot.fasta"  # Replace with your protein FASTA file

# Load and broadcast protein database
protein_database = parse_fasta(protein_file)
broadcast_protein_db = sc.broadcast(protein_database)

# Load genomic sequences and parallelize
genomic_sequences = parse_fasta(genomic_file)
genomic_rdd = sc.parallelize(genomic_sequences)

# Map genomic sequences to their best matching protein
def find_best_match(genomic_entry: Tuple[str, str], broadcast_db: Broadcast) -> Tuple[str, str]:
    header, sequence = genomic_entry
    best_match = find_best_matching_protein(sequence, broadcast_db.value)
    return header, best_match

results_rdd = genomic_rdd.map(lambda x: find_best_match(x, broadcast_protein_db))

# Collect and print results
results = results_rdd.collect()
print("\nBest Matching Proteins:")
for genomic_header, protein_match in results:
    print(f"Genomic Header: {genomic_header}, Best Protein Match: {protein_match}")
