import numpy as np

def create_blosum62():
    """
    Create the BLOSUM62 substitution matrix as a dictionary of dictionaries.
    """
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

def local_alignment_no_traceback(seq1, seq2, gap_penalty=-4):
    """
    Perform local alignment of two amino acid sequences using the Smith-Waterman algorithm 
    with a custom BLOSUM62 matrix, but without traceback.
    
    Args:
        seq1 (str): First amino acid sequence.
        seq2 (str): Second amino acid sequence.
        gap_penalty (int): Penalty for introducing a gap.

    Returns:
        int: The maximum alignment score.
    """
    # Create the BLOSUM62 matrix
    blosum62 = create_blosum62()

    # Initialize the scoring matrix
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score_matrix = np.zeros((rows, cols), dtype=int)

    # Track the maximum score
    max_score = 0

    # Fill the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            # Get substitution score from BLOSUM62
            match_score = blosum62[seq1[i-1]].get(seq2[j-1], gap_penalty)

            # Compute scores for match/mismatch, gap in seq1, gap in seq2
            diagonal = score_matrix[i-1][j-1] + match_score
            up = score_matrix[i-1][j] + gap_penalty
            left = score_matrix[i][j-1] + gap_penalty

            # Select the best score or reset to 0 (local alignment)
            best_score = max(0, diagonal, up, left)
            score_matrix[i][j] = best_score

            # Track the maximum score
            if best_score > max_score:
                max_score = best_score

    return max_score