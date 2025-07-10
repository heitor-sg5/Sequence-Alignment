import numpy as np

def build_score_matrix(pmm, p_gap):
    bases = ['A', 'C', 'G', 'T', '-']
    score_matrix = {}
    for b1 in bases:
        for b2 in bases:
            if b1 == '-' or b2 == '-':
                score = -p_gap
            elif b1 == b2:
                score = 1
            else:
                score = -pmm
            score_matrix[(b1, b2)] = score
    return score_matrix

def local_alignment(v, w, pmm, p_gap):
    score_matrix = build_score_matrix(pmm, p_gap)
    n, m = len(v), len(w)

    s = np.zeros((n+1, m+1), dtype=int)
    backtrack = np.full((n+1, m+1), '', dtype=object)

    max_score = 0
    max_pos = (0, 0) 

    for i in range(1, n+1):
        for j in range(1, m+1):
            match = s[i-1][j-1] + score_matrix[(v[i-1], w[j-1])]
            delete = s[i-1][j] + score_matrix[(v[i-1], '-')]
            insert = s[i][j-1] + score_matrix[('-', w[j-1])]

            s[i][j] = max(match, delete, insert, 0)

            if s[i][j] == 0:
                backtrack[i][j] = '0'
            elif s[i][j] == match:
                backtrack[i][j] = 'D'
            elif s[i][j] == delete:
                backtrack[i][j] = 'V'
            else:
                backtrack[i][j] = 'H'

            if s[i][j] > max_score:
                max_score = s[i][j]
                max_pos = (i, j)

    aligned_v = []
    aligned_w = []
    i, j = max_pos

    while s[i][j] != 0: 
        if backtrack[i][j] == 'D':
            aligned_v.append(v[i-1])
            aligned_w.append(w[j-1])
            i -= 1
            j -= 1
        elif backtrack[i][j] == 'V':
            aligned_v.append(v[i-1])
            aligned_w.append('-')
            i -= 1
        elif backtrack[i][j] == 'H':
            aligned_v.append('-')
            aligned_w.append(w[j-1])
            j -= 1

    return ''.join(reversed(aligned_v)), ''.join(reversed(aligned_w)), max_score

v = "GCCCAGTCTATGTCAGGGGGCACGAGCATGCACA"
w = "GCCGCCGTCGTTTTCAGCAGTTATGTTCAGAT"

aligned_v, aligned_w, score = local_alignment(v, w, pmm=1, p_gap=1)
print("Aligned v:", aligned_v)
print("Aligned w:", aligned_w)
print("Score:", score)
