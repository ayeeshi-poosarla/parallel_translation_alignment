from collections import Counter
from typing import List

def global_alignment(seq1, seq2, match=1, mismatch=-1, gap=-2):
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize DP table
    for i in range(1, m + 1):
        dp[i][0] = i * gap
    for j in range(1, n + 1):
        dp[0][j] = j * gap

    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                score = match
            else:
                score = mismatch
            dp[i][j] = max(dp[i - 1][j - 1] + score, dp[i - 1][j] + gap, dp[i][j - 1] + gap)

    # Final similarity score
    return dp[m][n]