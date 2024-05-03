from Bio import SeqIO
from collections import defaultdict

input_file = str(snakemake.input)
output_file = str(snakemake.output)

with open(output_file, 'w') as file:
    for record in SeqIO.parse(input_file, 'uniprot-xml'):
        for xref in record.dbxrefs:
            db, acc = xref.split(':', maxsplit = 1)
            if db == "GO":
                file.write(f'{acc}\t{record.id}\n')
