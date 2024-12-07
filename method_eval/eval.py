import time
import psutil
from sklearn.metrics import precision_score, recall_score, f1_score
from Bio import SeqIO
from collections import defaultdict

def calculate_alignment_errors(predicted_alignment, reference_alignment):
    """
    Calculate mismatch rate and indel rate between predicted and reference alignments.

    Parameters:
    predicted_alignment (str): The predicted alignment sequence.
    reference_alignment (str): The reference alignment sequence.

    Returns:
    tuple: A tuple containing the mismatch rate and the indel rate.
    """
    mismatches = 0
    indels = 0
    for p, r in zip(predicted_alignment, reference_alignment):
        if p != r:
            if p == '-' or r == '-':  
                indels += 1  # Indel detected if there's a gap ('-')
            else: 
                mismatches += 1  # Mismatch detected
    mismatch_rate = mismatches / len(reference_alignment)
    indel_rate = indels / len(reference_alignment)
    return mismatch_rate, indel_rate

def evaluate_accuracy(true_labels, predicted_labels):
    """
    Evaluate the accuracy of the predictions using precision, recall, and F1 score.

    Parameters:
    true_labels (list): The true labels for the data.
    predicted_labels (list): The predicted labels for the data.

    Returns:
    tuple: A tuple containing the precision, recall, and F1 score.
    """
    precision = precision_score(true_labels, predicted_labels)
    recall = recall_score(true_labels, predicted_labels)
    f1 = f1_score(true_labels, predicted_labels)
    return precision, recall, f1

def measure_performance(function, *args, **kwargs):
    """
    Measure the performance of a function, including execution time and memory usage.

    Parameters:
    function (function): The function to benchmark.
    *args: Arguments for the function.
    **kwargs: Keyword arguments for the function.

    Returns:
    tuple: A tuple containing the result of the function, execution time, and memory usage.
    """
    start_time = time.time()
    
    process = psutil.Process()
    memory_before = process.memory_info().rss  # Memory usage before execution
    
    result = function(*args, **kwargs)  # Call the function
    
    end_time = time.time()
    
    memory_after = process.memory_info().rss  # Memory usage after execution
    
    execution_time = end_time - start_time
    memory_usage = (memory_after - memory_before) / (1024 * 1024)  # Convert to MB
    
    return result, execution_time, memory_usage

def process_genomic_data(fasta_file, reference_file):
    """
    Process genomic data, compare predictions with a reference, and evaluate performance.

    Parameters:
    fasta_file (str): Path to the FASTA file containing the sequences to predict.
    reference_file (str): Path to the FASTA file containing the reference data.

    Returns:
    tuple: A tuple containing the predictions, precision, recall, and F1 score.
    """
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    predictions = [seq.id for seq in sequences]  # Extract predicted sequence IDs
    
    reference = defaultdict(str)
    for seq in SeqIO.parse(reference_file, "fasta"):
        reference[seq.id] = str(seq.seq)  # Read the reference data
    
    true_labels = [reference.get(pred, '') for pred in predictions]
    precision, recall, f1 = evaluate_accuracy(true_labels, predictions)
    
    return predictions, precision, recall, f1

def benchmark_method_on_synthetic_data():
    """
    Benchmark the genomic data processing method on synthetic data.

    Returns:
    None
    """
    synthetic_fasta = "synthetic_data.fasta"  
    reference_fasta = "reference_data.fasta" 
    
    predictions, precision, recall, f1, execution_time, memory_usage = measure_performance(process_genomic_data, synthetic_fasta, reference_fasta)
    
    # Print evaluation results
    print(f"Predictions: {predictions[:5]}...")  
    print(f"Precision: {precision}")
    print(f"Recall: {recall}")
    print(f"F1 Score: {f1}")
    print(f"Execution Time: {execution_time:.2f} seconds")
    print(f"Memory Usage: {memory_usage:.2f} MB")

benchmark_method_on_synthetic_data()
