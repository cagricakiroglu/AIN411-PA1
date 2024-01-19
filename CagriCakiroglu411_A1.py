import sys
import blosum as bl

def initialize_matrix(rows, cols, gap_open, gap_extend):
    """ Initialize the scoring matrix with gap penalties. """
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    for i in range(1, rows):
        matrix[i][0] = gap_open + gap_extend * (i - 1)
    for j in range(1, cols):
        matrix[0][j] = gap_open + gap_extend * (j - 1)
    return matrix

def global_alignment_with_identity(seq1, seq2, scoring_matrix, gap_open, gap_extend):
    """ Perform global alignment including the identity calculation. """
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = initialize_matrix(rows, cols, gap_open, gap_extend)
    # Track the source of each score for the traceback
    traceback = [[None for _ in range(cols)] for _ in range(rows)]

    # Filling the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            match = matrix[i-1][j-1] + scoring_matrix[seq1[i-1]][seq2[j-1]]
            delete = matrix[i-1][j] + (gap_extend if traceback[i-1][j] == 'D' else gap_open)
            insert = matrix[i][j-1] + (gap_extend if traceback[i][j-1] == 'I' else gap_open)
            if match >= delete and match >= insert:
                matrix[i][j] = match
                traceback[i][j] = 'M'
            elif delete > insert:
                matrix[i][j] = delete
                traceback[i][j] = 'D'
            else:
                matrix[i][j] = insert
                traceback[i][j] = 'I'

    # Traceback
    aligned_seq1, aligned_seq2 = '', ''
    i, j = len(seq1), len(seq2)
    while i > 0 and j > 0:
        if traceback[i][j] == 'M':
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback[i][j] == 'D':
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    while i > 0:
        aligned_seq1 = seq1[i-1] + aligned_seq1
        aligned_seq2 = '-' + aligned_seq2
        i -= 1
    while j > 0:
        aligned_seq1 = '-' + aligned_seq1
        aligned_seq2 = seq2[j-1] + aligned_seq2
        j -= 1

    # Calculate identity and match representation
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    identity_percentage = (matches / len(aligned_seq1)) * 100
    match_representation = ''.join('|' if a == b else ' ' for a, b in zip(aligned_seq1, aligned_seq2))

    final_score = matrix[-1][-1]
    return aligned_seq1, match_representation, aligned_seq2, final_score, matches, len(aligned_seq1), identity_percentage

# Main function to read file and perform alignment
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python ass1.py input.txt gap_open gap_extend blosum_param1 blosum_param2")
        sys.exit(1)

    filename = sys.argv[1]
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            if len(lines) < 2:
                print("Input file must contain two lines with sequences.")
                sys.exit(1)
            seq1, seq2 = lines[0].strip(), lines[1].strip()
    except FileNotFoundError:
        print("File " + filename + " not found.")
        sys.exit(1)

    try:
        gap_open = int(sys.argv[2])
        gap_extend = int(sys.argv[3])
        blosum_param1 = int(sys.argv[4])
        blosum_param2 = int(sys.argv[5])
    except ValueError:
        print("Gap penalties and BLOSUM parameters must be integers.")
        sys.exit(1)

    # BLOSUM parameter 1
    scoring_matrix = bl.BLOSUM(blosum_param1, default=0)
    print("\n---------------BLOSUM " + str(blosum_param1) + " ---------------------- \n ")
    aligned_seq1, match_representation, aligned_seq2, final_score, matches, alignment_length, identity_percentage = global_alignment_with_identity(seq1, seq2, scoring_matrix, gap_open, gap_extend)
    print("Aligned Sequence 1:  " + aligned_seq1)
    print("Match Representation:" + match_representation)
    print("Aligned Sequence 2:  " + aligned_seq2)
    print("Alignment score: " + str(final_score))
    print("Identity value: " + str(matches) + "/" + str(alignment_length) + " (" + str(round(identity_percentage, 1)) + "%)")

    # BLOSUM parameter 2
    scoring_matrix = bl.BLOSUM(blosum_param2, default=0)
    print("\n---------------BLOSUM " + str(blosum_param2) + " ---------------------- ")
    aligned_seq1, match_representation, aligned_seq2, final_score, matches, alignment_length, identity_percentage = global_alignment_with_identity(seq1, seq2, scoring_matrix, gap_open, gap_extend)
    print("Aligned Sequence 1:  " + aligned_seq1)
    print("Match Representation:" + match_representation)
    print("Aligned Sequence 2:  " + aligned_seq2)
    print("Alignment score: " + str(final_score))
    print("Identity value: " + str(matches) + "/" + str(alignment_length) + " (" + str(round(identity_percentage, 1)) + "%)")