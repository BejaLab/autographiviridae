from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation
import csv
import re

tsv_file = snakemake.input['tsv']
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
names = {}
with open(tsv_file) as tsv:
    reader = csv.DictReader(tsv, delimiter = '\t')
    with open(output_file, 'w') as out:
        for row in reader:
            if row["name"] not in names:
                scaffold, orf_num = row["name"].rsplit('_', 1)
                record = records[scaffold]
                start, stop, strand = extract_coords(row["description"])
                loc = SimpleLocation(start, stop, strand = strand)
                feature = SeqFeature(location = loc, type = "CDS")
                translation = feature.translate(record, cds = False).seq.rstrip('*')
                desc = coords_re.sub('', row["description"])
                strand_str = '+' if strand > 0 else '-'
                out.write(f">{row['name']} cds:{scaffold}:{start+1}-{stop}({strand_str}) {desc}\n{translation}\n")
                names[row["name"]] = True
