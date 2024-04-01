from Bio import SeqIO
import warnings
import json

jtree_file = snakemake.input['jtree']
blast_file = snakemake.input['blast']
gff_file = snakemake.input['gff']
fasta_file = snakemake.input['fasta']
saf_file = str(snakemake.output)

identity = snakemake.params['identity']

fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

genome_clades = {}
with open(jtree_file) as file:
    jtree = json.load(file)
    for edge in jtree['data']:
        if 'Accession' in edge and 'Clade_assigned' in edge:
            genome_clades[edge['Accession']] = edge['Clade_assigned']

clades = {}
contigs = {}
with open(blast_file) as blast:
    for line in blast:
        contig, genome, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
        if contig not in contigs and float(pident) >= identity and genome in genome_clades:
            clades[contig] = genome_clades[genome]
        contigs[contig] = 1

with open(gff_file) as gff:
    with open(saf_file, 'w') as saf:
        for line in gff:
            contig, source, gene, start, end, score, strand, frame, attribute = line.split()
            if contig in fasta and contig in clades:
                if gene == "exonuclease":
                    clade = clades[contig]
                    exo_end = int(end)
                else:
                    start = max(exo_end + 1, int(start))
                saf.write(f"{gene}_{clade}\t{contig}\t{start}\t{end}\t+\n")
