from Bio import pairwise2
from Bio.Seq import Seq

# Define two sequences
from Bio import SeqIO
from Bio import AlignIO
# Open and read the sequence


# Access the sequence and other information
##print(f"Sequence ID: {seq1.id}")
#print(f"Sequence: {seq1.seq}")
#print(f"Description: {seq1.description}")


file_path = "output.sto"

# Read the alignment
alignment = AlignIO.read(file_path, "stockholm")

# Loop through the alignment and extract sequences
for seq1 in alignment:
  for seq2 in alignment:
    print(f"Sequence ID: {seq1.id}")
    #print(f"Sequence: {record.seq}")
    #print(f"Description: {record.description}")
    #print()  # Blank line between sequences

    # Align the sequences using global alignment
    #seq2 = Seq(seq1.seq.replace('-', ''))

    #print(seq1.seq)
    #print(seq2)

    alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)
    best_alignment = alignments[0]
    from icecream import ic;import sys;ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr), includeContext=True,contextAbsPath=True);ic.configureOutput(prefix='')
    ic(best_alignment)
    # Extract the aligned sequences
    aligned_seq1, aligned_seq2 = best_alignment[0], best_alignment[1]

    # Calculate similarity (number of matches / length of alignment)
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    alignment_length = len(aligned_seq1)
    similarity = (matches / alignment_length) * 100


    # Print the alignment and similarity score
    #print(f"Alignment:\n{aligned_seq1}\n{aligned_seq2}")
    print(f"{seq1.id.strip()} Similarity: {similarity:.2f}%")

    #exit()
# Take the best alignment
