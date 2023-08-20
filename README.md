# Sequence_Local_Alignment_GapPenalty

This is a set of python scripts for affine gap local alignment of gene sequences and comparison of different gap penalty strategies, including affine gap, linear gap, constant gap and convex gap.

To use affine gap alignment, set this as working directory, type in: 

python BA_ICA_Code_AffineGap.py dna1 dna2 1

It will print out the result of sequence dna1 and dna2.
 
To compare gap penalties, set Comparison_gap_penalties_code as working directory, type in:

python comparison.py

It print out the alignment result of dna1&2, dna3&4, dna5&6, dna7&8, dna9&10, as different type of sequences, including short or long sequences with high or low similarities, and sequences with large gap in alignment.

For plotting the result, type in:

python stats_plot.py

It generate Align_Performance.png and Time_Performance.png files.
