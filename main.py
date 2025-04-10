from Bio import SeqIO

# Scoring scheme
match_score = 1
mismatch_penalty = -1
gap_penalty = -2

def read_fasta(file_path):
    sequences = list(SeqIO.parse(file_path, "fasta"))
    return str(sequences[0].seq), str(sequences[1].seq)

def needleman_wunsch(seq1, seq2):
    m, n = len(seq1), len(seq2)
    score = [[0 for j in range(n+1)] for i in range(m+1)]

    # Initialization
    for i in range(m+1):
        score[i][0] = i * gap_penalty
    for j in range(n+1):
        score[0][j] = j * gap_penalty

    # Scoring
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = score[i-1][j] + gap_penalty
            insert = score[i][j-1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    # Traceback
    align1, align2 = "", ""
    i, j = m, n
    while i > 0 and j > 0:
        current = score[i][j]
        diagonal = score[i-1][j-1]
        up = score[i-1][j]
        left = score[i][j-1]

        if current == diagonal + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif current == up + gap_penalty:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1

    # Final gaps
    while i > 0:
        align1 = seq1[i-1] + align1
        align2 = "-" + align2
        i -= 1
    while j > 0:
        align1 = "-" + align1
        align2 = seq2[j-1] + align2
        j -= 1

    return align1, align2, score[m][n]

if __name__ == "__main__":
    seq1, seq2 = read_fasta("sequences.fasta")
    aligned_seq1, aligned_seq2, score = needleman_wunsch(seq1, seq2)
    
    print("Sequence 1: ", aligned_seq1)
    print("Sequence 2: ", aligned_seq2)
    print("Alignment Score:", score)
