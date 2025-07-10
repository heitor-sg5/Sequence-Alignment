import numpy as np
from itertools import product
from collections import Counter
import math

def entropy_score(chars, p_gap):
    counts = Counter(chars)
    total = len(chars)

    entropy = 0
    for char, count in counts.items():
        p = count / total
        if char != '-':
            entropy -= p * math.log2(p)
        else:
            entropy += p_gap * count 

    return -entropy  

def t_dimensional_global_alignment(sequences, p_gap):
    t = len(sequences)
    lengths = [len(seq) for seq in sequences]

    s = np.full(tuple(l+1 for l in lengths), float('-inf'))
    backtrack = np.empty(tuple(l+1 for l in lengths), dtype=object)
    s[(0,) * t] = 0

    moves = [move for move in product([0, 1], repeat=t) if any(move)]
    it = np.ndindex(*tuple(l+1 for l in lengths))

    for idx in it:
        if idx == (0,) * t:
            continue
        max_score = float('-inf')
        best_move = None
        for move in moves:
            prev_idx = tuple(idx[i] - move[i] for i in range(t))
            if any(x < 0 for x in prev_idx):
                continue

            chars = []
            for i, m in enumerate(move):
                if m == 1:
                    chars.append(sequences[i][prev_idx[i]])
                else:
                    chars.append('-')

            prev_score = s[prev_idx]
            if prev_score == float('-inf'):
                continue

            curr_score = prev_score + entropy_score(chars, p_gap)
            if curr_score > max_score:
                max_score = curr_score
                best_move = move
        s[idx] = max_score
        backtrack[idx] = best_move

    aligned = [''] * t
    idx = tuple(lengths)
    while idx != (0,) * t:
        move = backtrack[idx]
        for i in range(t):
            if move[i] == 1:
                aligned[i] = sequences[i][idx[i]-1] + aligned[i]
            else:
                aligned[i] = '-' + aligned[i]
        idx = tuple(idx[i] - move[i] for i in range(t))

    return aligned, s[tuple(lengths)]

seqs = [
    "GCCCAGTCTATGTCAGGGGGCACGAGCATGCACA",
    "GCCGCCGTCGTTTTCAGCAGTTATGTTCAGAT",
    "GCCAGTCTATGTCAGGGGGCACGAGCAT"
]

aligned_seqs, score = t_dimensional_global_alignment(seqs, p_gap=1)

for i, aln in enumerate(aligned_seqs):
    print(f"Aligned {i+1}: {aln}")
print(f"Entropy Score: {score:.0f}")
