from Bio import SeqIO
import jsonlines
import re

bed_file = snakemake.input['bed']
fasta_file = snakemake.input['fasta']
jsonl_file = snakemake.input['jsonl']

out_file = str(snakemake.output)

fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

hmm_records = {}
with jsonlines.open(jsonl_file) as jsonl:
    for record in jsonl:
        contig = record["seq_name"].rsplit('_', maxsplit = 1)[0]
        if contig not in hmm_records or hmm_records[contig]["full"]["score"] < record["full"]["score"]:
            hmm_records[contig] = record

orf_re = re.compile('\[(\d+) - (\d+)\]')
exo_re = re.compile('\|(\d+):(\d+)')

with open(out_file, 'w') as out:
    with open(bed_file) as bed:
        for bed_line in bed:
            sseqid, soffset, send, qseqids, best_score, strand, scores, pidents = bed_line.split()
            if sseqid in fasta:
                fasta_record = fasta[sseqid]
                exo_stop = int(send) - int(soffset)
                best_qseqid = ''
                best_pident = ''
                for qseqid, score, pident in zip(qseqids.split(','), scores.split(','), pidents.split(',')):
                    if score == best_score:
                        best_qseqid = qseqid
                        best_pident = pident
                exo_start, exo_stop = exo_re.findall(fasta_record.description)[0]
                exo_start, exo_stop = int(exo_start), int(exo_stop)
                out.write(f"{sseqid}\ttblastn\texonuclease\t{exo_start + 1}\t{exo_stop}\t{best_score}\t+\t.\t{best_qseqid}\n")
                if sseqid in hmm_records:
                    record = hmm_records[sseqid]
                    hmm_score = record['full']['score']
                    hmm_evalue = record['full']['evalue']
                    description = record['description']
                    env_from = record['domains'][0]['env']['from']
                    env_to   = record['domains'][0]['env']['to']
                    orf_start, orf_end = orf_re.findall(record['description'])[0]
                    orf_start, orf_end = int(orf_start), int(orf_end)
                    env_start = orf_start + (env_from - 1) * 3
                    env_stop = orf_start + env_to * 3
                    out.write(f"{sseqid}\thmmsearch\tnblA\t{env_start}\t{env_stop}\t{hmm_score}\t+\t.\tnblA\n")
