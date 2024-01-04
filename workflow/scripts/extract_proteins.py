from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation
import re
import jsonlines

jsonl_file = snakemake.input['jsonl']
fna_file = snakemake.input['fna']
output_file = str(snakemake.output)

coords_re = re.compile(r'\[([0-9]+) *- *([0-9]+)\]')
def extract_coords(desc):
    coord1, coord2 = coords_re.match(desc).groups()
    coord1 = int(coord1)
    coord2 = int(coord2)
    if coord1 < coord2:
        start = coord1 - 1
        stop = coord2
        strand = 1
    else:
        stop = coord1
        start = coord2 - 1
        strand = -1
    return start, stop, strand

records = SeqIO.to_dict(SeqIO.parse(fna_file, 'fasta'))
with jsonlines.open(jsonl_file) as hits:
    with open(output_file, 'w') as out:
        for hit in hits:
            name, description = hit["seq_name"], hit["description"]
            scaffold, orf_num = name.rsplit('_', 1)
            record = records[scaffold]
            start, stop, strand = extract_coords(description)
            loc = SimpleLocation(start, stop, strand = strand)
            feature = SeqFeature(location = loc, type = "CDS")
            translation = feature.translate(record, cds = False).seq.rstrip('*')
            strand_str = '+' if strand > 0 else '-'
            out.write(f">{name} cds:{scaffold}:{start+1}-{stop}({strand_str})\n{translation}\n")
