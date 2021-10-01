import numpy as np
import os
import parameters as params


# Needleman-Wunsch algorithm
def global_alignment(seq1, seq2, scoring_matrix):

    gap_penalty = -5

    n = len(seq1)
    m = len(seq2)

    dp = np.zeros((n + 1, m + 1), dtype=int)

    for i in range(n + 1):
        dp[i][0] = gap_penalty * i
    for i in range(m + 1):
        dp[0][i] = gap_penalty * i
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            dp[i][j] = max(
                            dp[i-1][j-1] + scoring_matrix[seq1[i-1]][seq2[j-1]], 
                            dp[i][j-1] + gap_penalty, 
                            dp[i-1][j] + gap_penalty
                        )

    gen_seq1 = ""
    gen_seq2 = ""

    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i-1][j-1] == dp[i][j] - scoring_matrix[seq1[i-1]][seq2[j-1]]:
            gen_seq1 += seq1[i-1]
            gen_seq2 += seq2[j-1]
            i -= 1
            j -= 1
        elif j > 0 and dp[i][j-1] == dp[i][j] - gap_penalty:
            gen_seq1 += '_'
            gen_seq2 += seq2[j-1]
            j -= 1
        elif i > 0 and dp[i-1][j] == dp[i][j] - gap_penalty:
            gen_seq1 += seq1[i-1]
            gen_seq2 += '_'
            i -= 1

    gen_seq1 = gen_seq1[::-1]
    gen_seq2 = gen_seq2[::-1]

    print('\nGlobal Alignment:\n\n', 'Score:', dp[n][m], '\n\n', 'Sequence 1:', gen_seq1, '\n', 'Sequence 2:', gen_seq2, '\n', sep='')


# Smith-Waterman algorithm
def local_alignment(seq1, seq2, scoring_matrix):

    gap_penalty = -5
    
    n = len(seq1)
    m = len(seq2)
    
    dp = np.zeros((n + 1, m + 1), dtype=int)
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            dp[i][j] = max(
                            dp[i-1][j-1] + scoring_matrix[seq1[i-1]][seq2[j-1]], 
                            dp[i][j-1] + gap_penalty, 
                            dp[i-1][j] + gap_penalty, 
                            0
                        )
    
    gen_seq1 = ""
    gen_seq2 = ""
    
    i, j = 0, 0
    
    score = 0
    for k in range(n+1):
        for l in range(m+1):
            if score <= dp[k][l]:
                i, j = k, l
                score = dp[k][l]

    while (i > 0 or j > 0) and dp[i][j] > 0:
        if i > 0 and j > 0 and dp[i-1][j-1] == dp[i][j] - scoring_matrix[seq1[i-1]][seq2[j-1]]:
            gen_seq1 += seq1[i-1]
            gen_seq2 += seq2[j-1]
            i -= 1
            j -= 1
        elif j > 0 and dp[i][j-1] == dp[i][j] - gap_penalty:
            gen_seq1 += '_'
            gen_seq2 += seq2[j-1]
            j -= 1
        elif i > 0 and dp[i-1][j] == dp[i][j] - gap_penalty:
            gen_seq1 += seq1[i-1]
            gen_seq2 += '_'
            i -= 1
    
    gen_seq1 = gen_seq1[::-1]
    gen_seq2 = gen_seq2[::-1]
    
    print('\nLocal Alignment:\n\n', 'Score:', score, '\n\n', 'Sequence 1:', gen_seq1, '\n', 'Sequence 2:', gen_seq2, '\n', sep='')


# Main loop
if __name__ == "__main__":


    seq_type = input("For nucleotide, press 0 and for protein, press 1:\n")
    align_type = input("\nFor global alignment, press 0 and for local alignment, press 1:\n")
    seq1_path = input("\nEnter path to first sequence:\n")
    seq2_path = input("\nEnter path to second sequence:\n")

    if os.path.isfile(seq1_path) and os.path.isfile(seq2_path):
        with open(seq1_path, 'r') as file:
            seq1 = file.read().replace('\n', '')
        with open(seq2_path, 'r') as file:
            seq2 = file.read().replace('\n', '')

        scoring_matrix = params.DNAFull
        if seq_type == "1":
            scoring_matrix = params.BLOSUM62

        if align_type == "1":
            local_alignment(seq1, seq2, scoring_matrix)
        else:
            global_alignment(seq1, seq2, scoring_matrix)
    
    else:
        print("Invalid input path")