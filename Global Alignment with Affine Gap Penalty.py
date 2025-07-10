import numpy as np

def build_score_matrix(pmm):
    bases = ['A', 'C', 'G', 'T']
    score_matrix = {}
    for b1 in bases:
        for b2 in bases:
            if b1 == b2:
                score = 1
            else:
                score = -pmm
            score_matrix[(b1, b2)] = score
    return score_matrix

def global_alignment_affine(v, w, pmm, gap_open, gap_extend):
    score_matrix = build_score_matrix(pmm)
    n, m = len(v), len(w)

    M = np.full((n+1, m+1), -np.inf)
    Ix = np.full((n+1, m+1), -np.inf)
    Iy = np.full((n+1, m+1), -np.inf)
    backtrack = np.full((n+1, m+1), '', dtype=object)

    M[0][0] = 0
    for i in range(1, n+1):
        Iy[i][0] = -gap_open - (i - 1) * gap_extend
        M[i][0] = Iy[i][0]
        backtrack[i][0] = 'V'
    for j in range(1, m+1):
        Ix[0][j] = -gap_open - (j - 1) * gap_extend
        M[0][j] = Ix[0][j]
        backtrack[0][j] = 'H'

    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score_matrix[(v[i-1], w[j-1])]
            M[i][j] = max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]) + match
            Ix[i][j] = max(M[i][j-1] - gap_open, Ix[i][j-1] - gap_extend)
            Iy[i][j] = max(M[i-1][j] - gap_open, Iy[i-1][j] - gap_extend)

            if M[i][j] >= Ix[i][j] and M[i][j] >= Iy[i][j]:
                backtrack[i][j] = 'D'
            elif Iy[i][j] >= Ix[i][j]:
                backtrack[i][j] = 'V'
            else:
                backtrack[i][j] = 'H'

    aligned_v = []
    aligned_w = []
    i, j = n, m
    while i > 0 or j > 0:
        if backtrack[i][j] == 'D':
            aligned_v.append(v[i-1])
            aligned_w.append(w[j-1])
            i -= 1
            j -= 1
        elif backtrack[i][j] == 'V':
            aligned_v.append(v[i-1])
            aligned_w.append('-')
            i -= 1
        else:
            aligned_v.append('-')
            aligned_w.append(w[j-1])
            j -= 1

    return ''.join(reversed(aligned_v)), ''.join(reversed(aligned_w)), int(M[n][m])

v = "GCCCAGTCTATGTCAGGGGGCACGAGCATGCACA"
w = "GCCGCCGTCGTTTTCAGCAGTTATGTTCAGAT"

aligned_v, aligned_w, score = global_alignment_affine(v, w, pmm=1, gap_open=1, gap_extend=1)
print("Aligned v:", aligned_v)
print("Aligned w:", aligned_w)
print("Score:", score)
