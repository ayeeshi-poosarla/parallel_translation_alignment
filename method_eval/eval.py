import time
import psutil
from sklearn.metrics import precision_score, recall_score, f1_score
from Bio import SeqIO
from collections import defaultdict

# Function to calculate mismatch rate and indel rate
def calculate_alignment_errors(predicted_alignment, reference_alignment):
    mismatches = 0
    indels = 0
    for p, r in zip(predicted_alignment, reference_alignment):
        if p != r:
            if p == '-' or r == '-':  
                indels += 1
            else: 
                mismatches += 1
    mismatch_rate = mismatches / len(reference_alignment)
    indel_rate = indels / len(reference_alignment)
    return mismatch_rate, indel_rate

# Function to evaluate accuracy using precision, recall, and F1 score
def evaluate_accuracy(true_labels, predicted_labels):
    precision = precision_score(true_labels, predicted_labels)
    recall = recall_score(true_labels, predicted_labels)
    f1 = f1_score(true_labels, predicted_labels)
    return precision, recall, f1

# Function to evaluate runtime and memory usage
def measure_performance(function, *args, **kwargs):
    start_time = time.time()
    
    process = psutil.Process()
    memory_before = process.memory_info().rss  # in bytes
    
    result = function(*args, **kwargs)
    
    end_time = time.time()
    
    memory_after = process.memory_info().rss  # in bytes
    
    execution_time = end_time - start_time
    memory_usage = (memory_after - memory_before) / (1024 * 1024)  # Convert to MB
    
    return result, execution_time, memory_usage

# Example function to process FASTQ/FASTA files and evaluate them
def process_genomic_data(fasta_file, reference_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    predictions = []  
    for seq in sequences:
        # analysis output here
        predictions.append(seq.id) 
    
    # Simulate reading reference data
    reference = defaultdict(str)
    for seq in SeqIO.parse(reference_file, "fasta"):
        reference[seq.id] = str(seq.seq)
    
    # Compare predictions with reference
    true_labels = [reference.get(pred, '') for pred in predictions]
    precision, recall, f1 = evaluate_accuracy(true_labels, predictions)
    
    return predictions, precision, recall, f1

# Function to benchmark on synthetic data (you can use real datasets here)
def benchmark_method_on_synthetic_data():
    synthetic_fasta = "synthetic_data.fasta"  
    reference_fasta = "reference_data.fasta" 
    
    # Measure performance (runtime and memory usage)
    predictions, precision, recall, f1 = measure_performance(process_genomic_data, synthetic_fasta, reference_fasta)
    
    # Print evaluation results
    print(f"Predictions: {predictions[:5]}...")  
    print(f"Precision: {precision}")
    print(f"Recall: {recall}")
    print(f"F1 Score: {f1}")
    print(f"Execution Time: {execution_time:.2f} seconds")
    print(f"Memory Usage: {memory_usage:.2f} MB")

benchmark_method_on_synthetic_data()
