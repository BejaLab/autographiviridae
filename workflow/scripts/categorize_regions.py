from Bio import SeqIO
from pandas import read_excel
import warnings

phages_file = snakemake.input['phages']
gff_file = snakemake.input['gff']
fasta_file = snakemake.input['fasta']
saf_file = str(snakemake.output)

fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

with warnings.catch_warnings(record=True):
    warnings.simplefilter("always")
    phages = read_excel(phages_file)
clades = dict(zip(phages['Accession'], phages['Clade']))

with open(gff_file) as gff:
    with open(saf_file, 'w') as saf:
        for line in gff:
            contig, source, gene, start, end, score, strand, frame, category = line.split()
            if contig in fasta:
                if gene == "exonuclease":
                    clade = clades[category]
                    exo_end = int(end)
                else:
                    start = max(exo_end + 1, int(start))
                saf.write(f"{gene}_{clade}\t{contig}\t{start}\t{end}\t+\n")
