from Bio.Align import substitution_matrices
import numpy as np

pam250 = substitution_matrices.load("PAM250")

def global_alignment(v, w, score_matrix, gap_penalty):
    n, m = len(v), len(w)
    s = np.zeros((n+1, m+1), dtype=int)
    backtrack = np.full((n+1, m+1), '', dtype=object)

    for i in range(1, n+1):
        s[i][0] = s[i-1][0] + gap_penalty
        backtrack[i][0] = 'V'
    for j in range(1, m+1):
        s[0][j] = s[0][j-1] + gap_penalty
        backtrack[0][j] = 'H'

    for i in range(1, n+1):
        for j in range(1, m+1):
            match_score = score_matrix.get((v[i-1], w[j-1]), -float('inf'))
            match = s[i-1][j-1] + match_score
            delete = s[i-1][j] + gap_penalty
            insert = s[i][j-1] + gap_penalty
            s[i][j] = max(match, delete, insert)

            if s[i][j] == match:
                backtrack[i][j] = 'D'
            elif s[i][j] == delete:
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

    return ''.join(reversed(aligned_v)), ''.join(reversed(aligned_w)), s[n][m]

v = "HEAGAWGHEE"
w = "PAWHEAE"

aligned_v, aligned_w, score = global_alignment(v, w, pam250, gap_penalty=-1)

print("Aligned v:", aligned_v)
print("Aligned w:", aligned_w)
print("Score:", score)
