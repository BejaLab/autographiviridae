from pandas import read_excel, notna

xlsx_file = str(snakemake.input)
fasta_file = str(snakemake.output)

df = read_excel(xlsx_file)

# Filter rows based on the 'exclude' column
df = df[df['whole_genome'].isna()]

# Create a function to format each row as a FASTA record
def format_fasta(row, gene):
    sequence = row[f"{gene}_cds"]
    if notna(sequence) and sequence.strip() != '':
        return f">{row['name']}_{gene} {row['host']}\n{sequence}\n"
    else:
        return ''

# Write the FASTA data to a file
with open(fasta_file, 'w') as fasta:
    # Apply the function to each row and join the results
    for record in df.apply(format_fasta, axis = 1, gene = 'pol'):
        fasta.write(record)
    for record in df.apply(format_fasta, axis = 1, gene = 'mcp'):
        fasta.write(record)
