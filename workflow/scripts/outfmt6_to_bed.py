
outfmt6_file = str(snakemake.input)
bed_file = str(snakemake.output)

with open(outfmt6_file) as outfmt6:
    with open(bed_file, 'w') as bed:
        for line in outfmt6:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
            sstart = int(sstart)
            send = int(send)
            strand = '+' if sstart < send else '-'
            sstart, send = min(sstart, send), max(sstart, send)
            line = f"{sseqid}\t{sstart-1}\t{send}\t{qseqid}\t{bitscore}\t{strand}\t{pident}\n"
            bed.write(line)
