# Sequence Alignment Algorithms

This repository contains implementations of fundamental sequence alignment algorithms widely used in bioinformatics. It covers a range of techniques including global, local, and fitting alignments, as well as advanced methods like affine gap penalties, PAM scoring, and linear-space optimizations. Multiple sequence alignment and entropy-based scoring are also included, supporting diverse applications in comparative genomics and evolutionary analysis.

---

## 🧬 What is Sequence Alignment?

Sequence alignment is a method used to arrange DNA, RNA, or protein sequences to identify regions of similarity that may indicate functional, structural, or evolutionary relationships. It helps detect conserved sequences and understand biological processes across different organisms.

Mutations such as substitutions, insertions, and deletions cause variations in sequences. Alignment algorithms account for these mutations by introducing gaps or mismatches to best match related sequences, enabling researchers to study genetic variation, evolutionary divergence, and identify important functional sites.

---

## 📁 Files in This Repository

- `needleman-wunsch_global_alignment.py`: Aligns two sequences end-to-end using a 2D dynamic programming matrix.
- `smith-waterman_local_alignment.py`: Finds the best matching subsequences between two sequences via local alignment.
- `fitting_alignment.py`: Aligns one (shorter) sequence completely to a substring of another (longer) sequence.
- `global_alignment_with_affine_gap_penalty.py`: Performs global alignment with separate penalties for gap opening (high) and gap extension (lower).
- `global_alignment_with_pam_scoring.py`: Uses PAM substitution matrices to score global alignment of protein sequences.
- `global_alignment_with_hirschberg_linear_space.py`: Implements a memory-efficient global alignment using Hirschberg’s divide-and-conquer method.
- `global_msa.py`: Produces a global multiple sequence alignment for a collection of sequences using t-way score matrices.
- `global_msa_with_entropy_scoring.py`: Performs multiple sequence alignment using Shannon's entropy-based scoring for columns.

---

## ⚙️ How to Use

### 1. Prepare Input

The inputs vary depending on the alignment algorithm used:

- For most algorithms, input two sequences (v, w) of DNA, RNA, or protein (on-letter) codes, with similar lengths. Except for `fitting_alignment.py` which requires w to be shorter.
- For `global_alignment_with_pam_scoring.py`, the inputs are two protein sequences.
- For the `global_msa.py` and `global_msa_with_entropy_scoring.py` algorithms, the input are multiple sequences (more than two) as strings of DNA, RNA, or proteins.

### 2. Run the Algorithms

Each algorithm produces:

- Two aligned sequences showing insertions/deletions (gaps as '-') and matches/mismatches.
- A score indicating the quality of the alignment; the highest positive value for the best alignment or the most negative value, if using `global_msa_with_entropy_scoring.py`.

---

#### Needleman-Wunsch Global Alignment

  bash
```needleman-wunsch_global_alignment.py```

#### Smith-Waterman Local Alignment

  bash
```smith-waterman_local_alignment.py```

#### Fitting Alignment

  bash
```fitting_alignment.py```

#### Global Alignment with Affine Gap

  bash
```global_alignment_with_affine_gap_penalty.py```

#### Global Alignment with PAM

  bash
```global_alignment_with_pam_scoring.py```

#### Hirschberg's Global Alignment

  bash
```global_alignment_with_hirschberg_linear_space.py```

#### Global Multiple Sequence Alignment (MSA)

  bash
```global_msa.py```

#### Global MSA with Entropy

  bash
```global_msa_with_entropy_scoring.py```

The gap penalty (p_gap) and mismatch penalty (pmm) are configurable parameters for the pairwise alignments. Additionally, `global_alignment_with_affine_gap_penalty.py` has gap start (gap_open) and extension (gap_extend) penalty values.

---

## 🧠 Algorithm Overviews

### Needleman-Wunsch Global Alignment

- Aligns two sequences end-to-end using dynamic programming, optimizing a score based on matches, mismatches, and gap penalties.
- Builds a scoring matrix where matches score positively, mismatches and gaps are penalized with specified parameters (pmm and p_gap).
- Uses backtracking through the matrix to reconstruct the highest-scoring alignment and returns aligned sequences with the total alignment score.

### Smith-Waterman Local Alignment

- Identifies the highest-scoring matching subsequences between two sequences using dynamic programming with match, mismatch, and gap penalties.
- Scores are reset to zero when negative, allowing identification of optimal local regions rather than full-length alignments.
- Backtracks from the highest scoring cell to reconstruct the best local alignment and returns the aligned subsequences with their score.

### Fitting Alignment

- Aligns a shorter sequence completely to a substring of a longer sequence, useful when one sequence is expected to fit within another.
- Uses dynamic programming with match, mismatch, and gap penalties, but allows free gaps at the start of the longer sequence.
- Finds the highest score at the end of the shorter sequence and backtracks to reconstruct the optimal fitting alignment.

### Global Alignment with Affine Gap

- Implements global alignment considering different penalties for gap opening (high) and gap extension (lower), reflecting biological reality more accurately.
- Uses three dynamic programming matrices to track match/mismatch and gaps separately, improving gap scoring precision.
- Outputs the optimal end-to-end alignment with the highest score accounting for both substitution and affine gap penalties.

### Global Alignment with PAM

- Performs global alignment of protein sequences using (Biopython's) PAM substitution matrix to score amino acid matches and mismatches.
- Incorporates a uniform gap penalty for insertions and deletions while leveraging biologically informed substitution scores.
- Returns the highest-scoring end-to-end alignment that reflects evolutionary likelihoods based on PAM scores.

### Hirschberg's Global Alignment

- Implements a space-efficient global alignment method using a divide-and-conquer strategy, reducing memory from quadratic to linear in sequence length.
- Recursively splits one sequence and computes alignment scores from both directions to find the optimal split point.
- Uses Needleman-Wunsch for base cases, reconstructing the full optimal alignment while only storing linear space score arrays.

### Global MSA

- Extends pairwise global alignment to multiple sequences using dynamic programming on a t-dimensional scoring matrix, where t is the number of sequences.
- Considers all possible moves (combinations of gaps and matches across sequences) at each step and scores them by summing pairwise substitution and gap penalties.
- Uses backtracking from the matrix’s endpoint to reconstruct the optimal multiple alignment, producing aligned sequences with gaps to maximize the overall score.

### Global MSA with Entropy

- Extends multiple sequence alignment by scoring alignments using an entropy-based function that rewards conservation and penalizes gaps, reflecting sequence variability.
- Uses dynamic programming over a t-dimensional matrix where each cell considers all possible insertion/deletion/match combinations across sequences.
- Recursively computes scores by combining previous scores with entropy of the aligned characters, then backtracks to produce the optimal multiple alignment and its overall entropy score.

---

## 🧪 Example Output

Aligned v: GCC-CAGTC-TATGTCAGGGGGCACGAGCATG--CACA-

Aligned w: GCCGCCGTCGT-TTTCA----GCA-G-TTATGTTCAGAT

Score: 5

---

## 👤 Author

Heitor Gelain do Nascimento
Email: heitorgelain@outlook.com
GitHub: @heitor-sg5

---

## 📚 References

Bioinformatics Algorithms: An Active Learning Approach (Chapter 5) by
Phillip Compeau & Pavel Pevzner
https://bioinformaticsalgorithms.com
