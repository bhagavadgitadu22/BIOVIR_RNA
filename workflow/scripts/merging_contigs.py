import os
import sys

def merge_contigs_with_gaps(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    merged_sequence = ''
    first = True
    for line in lines:
        if line.startswith('>'):
            if not(first):
                merged_sequence += 'N' * 1000  # Add 1000 Ns between contigs
        else:
            merged_sequence += line.strip()
        first = False

    name_contig = os.path.basename(input_file).replace(".fa", "")
    with open(output_file, 'w') as f:
        f.write('>'+name_contig+'\n')  # Create a header for the merged contig
        # Write the merged sequence with contigs separated by 1000 Ns
        f.write(merged_sequence)

input_file = 'input.fasta'  # Replace 'input.fasta' with the path to your input FASTA file
output_file = 'merged_contig.fasta'  # Replace 'merged_contig.fasta' with the desired output file path

merge_contigs_with_gaps(sys.argv[1], sys.argv[2])
