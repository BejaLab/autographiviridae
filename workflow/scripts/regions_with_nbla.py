from Bio import SeqIO
from pandas import read_excel
import jsonlines

blast_file = snakemake.input['blast']
fasta_file = snakemake.input['fasta']
jsonl_file = snakemake.input['jsonl']
phages_file = snakemake.input['phages']

phages = read_excel(phages_file)
clades = dict(zip(phages['Phage'], phages['Subclade']))

out_file = str(snakemake.output)

fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

records = {}
with jsonlines.open(jsonl_file) as jsonl:
    for record in jsonl:
        contig = record["seq_name"].rsplit('_', maxsplit = 1)[0]
        if contig not in records or records[contig]["full"]["score"] < record["full"]["score"]:
            records[contig] = record

with open(out_file, 'w') as out:
    with open(blast_file) as blast:
        for line in blast:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
            if sseqid in fasta:
                score = evalue = -1
                description = ''
                if sseqid in records:
                    score = records[sseqid]['full']['score']
                    evalue = records[sseqid]['full']['evalue']
                    description = records[sseqid]['description']
                out.write(f'{sseqid}\t{qseqid}\t{clades[qseqid]}\t{pident}\t{length}\t{description}\t{score}\t{evalue}\n')
