import numpy as np

# Smith-Waterman algorithm for local alignment
def smith_waterman(reference, query, match=2, mismatch=-1, gap_penalty=-2):
    # Create the scoring matrix
    len_ref, len_query = len(reference), len(query)
    scoring_matrix = np.zeros((len_ref + 1, len_query + 1))

    # Fill the matrix
    max_score = 0
    max_pos = (0, 0)
    
    for i in range(1, len_ref + 1):
        for j in range(1, len_query + 1):
            if reference[i - 1] == query[j - 1]:
                score = scoring_matrix[i - 1, j - 1] + match
            else:
                score = scoring_matrix[i - 1, j - 1] + mismatch
            
            # Determine the score with the gap penalty
            score = max(score, scoring_matrix[i - 1, j] + gap_penalty, scoring_matrix[i, j - 1] + gap_penalty, 0)
            
            scoring_matrix[i, j] = score
            
            # Track the maximum score and its position
            if score > max_score:
                max_score = score
                max_pos = (i, j)

    # Traceback to get the aligned sequences
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
    
    # Reverse the sequences
    aligned_ref = ''.join(reversed(aligned_ref))
    aligned_query = ''.join(reversed(aligned_query))
    
    return aligned_ref, aligned_query, max_score

# Sample k-mer index (you'd have generated this from your previous code)
k_mer_index = {
    'seq1': 'MKTAYIAKQRQISFVKSHFSRQHFLASVHTKLNQ',
    'seq2': 'MKTAYIAKQRQISFVKHFSRQHFLASVHTKLNQ',
    'seq3': 'MKTAFIAKQRQISFVKSHFSRQHFLASVHTKLNQ',
}

# Top 2 sequences based on k-mer matching
top_sequences = ['seq1', 'seq2']

# Reference sequence from the k-mer index (for alignment)
reference_sequence = k_mer_index['seq3']

# Function to compare the top 2 sequences using Smith-Waterman
def compare_top_kmers(k_mer_index, top_sequences, reference_sequence):
    best_score = float('-inf')
    best_alignment = None
    best_query_sequence = None
    
    for seq_id in top_sequences:
        query_sequence = k_mer_index[seq_id]
        
        # Perform Smith-Waterman alignment
        aligned_ref, aligned_query, alignment_score = smith_waterman(reference_sequence, query_sequence)
        
        print(f"\nAlignment for {seq_id} (Query):")
        print("Reference:", aligned_ref)
        print("Query:    ", aligned_query)
        print("Alignment Score:", alignment_score)
        
        # Track the best alignment score
        if alignment_score > best_score:
            best_score = alignment_score
            best_alignment = (aligned_ref, aligned_query)
            best_query_sequence = query_sequence
    
    return best_query_sequence, best_alignment, best_score

# Run comparison
best_query_sequence, best_alignment, best_score = compare_top_kmers(k_mer_index, top_sequences, reference_sequence)

print("\nBest Match Based on Alignment Score:")
print(f"Best Query Sequence: {best_query_sequence}")
print("Best Alignment:")
print("Reference:", best_alignment[0])
print("Query:    ", best_alignment[1])
print("Best Alignment Score:", best_score)
