from Bio import SeqIO

bed_file = snakemake.input['bed']
fna_file = snakemake.input['fna']
out_file = str(snakemake.output)

max_upstream = snakemake.params['max_upstream']
max_downstream = snakemake.params['max_downstream']
min_downstream = snakemake.params['min_downstream']

fasta = SeqIO.to_dict(SeqIO.parse(fna_file, 'fasta'))

with open(bed_file) as bed:
    with open(out_file, 'w') as out:
        for line in bed:
            seqid, start, stop, _, _, strand, *rest = line.split()
            start, stop = int(start), int(stop)
            record = fasta[seqid]
            seq_len = len(record.seq)
            if strand == '+':
                reg_start = max(start - max_upstream, 0)
                reg_stop  = min(stop + max_downstream, seq_len)
                record.seq = record.seq[reg_start:reg_stop]
                upstream = start - reg_start
                downstream = reg_stop - stop
            else:
                reg_start = max(start - max_downstream, 0)
                reg_stop  = min(stop + max_upstream, seq_len)
                record.seq = record.seq[reg_start:reg_stop]
                record.seq = record.seq.reverse_complement()
                downstream = start - reg_start
                upstream = reg_stop - stop
            feature_len = stop - start
            if downstream >= min_downstream and upstream >= 0:
                record.description = f"{reg_start}:{reg_stop}({strand})|{upstream}:{upstream + feature_len}"
                SeqIO.write(record, out, 'fasta')
