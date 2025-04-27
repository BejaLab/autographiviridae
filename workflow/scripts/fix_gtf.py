import re

tr_re = re.compile('transcript_id "(\w+)"')

with open(str(snakemake.input)) as fh1:
    with open(str(snakemake.output), 'w') as fh2:
        prev_seqname = ''
        for line in fh1:
            seqname, source, feature, *rest = line.split()
            if seqname != prev_seqname:
                tr = 1
            elif feature == 'transcript':
                tr += 1
            line = tr_re.sub('transcript_id "{seqname}_{tr:04d}"'.format(seqname = seqname, tr = tr), line)
            fh2.write(line)
            prev_seqname = seqname
