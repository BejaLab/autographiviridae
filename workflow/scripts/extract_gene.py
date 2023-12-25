from Bio import SeqIO

tsv_file = snakemake.input['tsv']
fna_file = snakemake.input['fna']
out_file = str(snakemake.output)

max_pad = snakemake.params['max_pad']
min_pad = snakemake.params['min_pad']

fasta = SeqIO.to_dict(SeqIO.parse(fna_file, 'fasta'))

def coords(start, end):
    c1 = int(start)
    c2 = int(end)
    if c1 < c2:
        strand = 1
        s = c1 - 1
        e = c2
    else:
        strand = -1
        s = c2 - 1
        e = c1
    return s, e, strand

with open(tsv_file) as tsv:
    with open(out_file, 'w') as out:
        for line in tsv:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
            record = fasta[sseqid]
            start, stop, strand = coords(sstart, send)
            if strand > 0:
                record.seq = record.seq[start:stop + max_pad]
            else:
                record.seq = record.seq[max(start-max_pad,0):stop]
                record.seq = record.seq.reverse_complement()
            pad = len(record.seq) - (stop - start)
            if pad >= min_pad:
                record.description = f'{start}:{stop}({"+" if strand > 0 else "-"})'
                SeqIO.write(record, out, 'fasta')
