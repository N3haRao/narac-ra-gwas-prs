# Python code to obtain RA_pcs.txt from ra_pruned.evec

# Read the input file
input_file = "ra_pruned.evec"
output_file = "RA_pcs.txt"

with open(input_file, 'r') as f:
    lines = f.readlines()

# Remove the first line
lines.pop(0)

# Split the first column into two columns and remove the last column
formatted_lines = []
for line in lines:
    split_line = line.split()
    split_column = split_line[0].split(":")
    formatted_line = split_column + split_line[1:-1]
    formatted_lines.append(formatted_line)

# Add header with consistent spacing
header = "FID      IID      PC1       PC2       PC3       PC4       PC5       PC6       PC7       PC8       PC9       PC10\n"

# Write to the output file
with open(output_file, 'w') as f:
    f.write(header)
    for line in formatted_lines[1:]:
        f.write(' '.join(line) + '\n')
