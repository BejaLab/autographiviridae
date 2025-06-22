from pandas import read_csv, DataFrame
import re
import textwrap
from collections import defaultdict

nexus_file = snakemake.input['nexus']
data_file = snakemake.input['data']
clstr_file = snakemake.input['clstr']

output_file = str(snakemake.output)

def read_taxa(nexus_file):
    taxa = []
    with open(nexus_file) as file:
        in_tax_list = False
        ntax = -1
        for line in file:
            if re.match('taxlabels', line, re.IGNORECASE):
                in_tax_list = True
                break
            if match := re.match('dimensions ntax=(\\d+)', line, re.IGNORECASE):
                ntax = int(match.group(1))
        assert in_tax_list and ntax > 0, "Taxon labels not found"
        for i in range(ntax):
            line = next(file)
            match = re.match('\\[(\\d+)\\] \'(.+)\'', line)
            taxa.append(match.group(2))
    return taxa

def read_clstr(clstr_file):
    clusters = defaultdict(list)
    with open(clstr_file) as file:
        for line in file:
            if line.startswith('>Cluster'):
                cluster = line
            else:
                match = re.match('(\\d+)\t(\\d+)aa, >(.+)[.][.][.] (.+)', line)
                num, seq_len, seq_name, identity = match.groups()
                clusters[cluster].append(seq_name)
    data = []
    for cluster, labels in clusters.items():
        for label in labels:
            data.append({ 'label': label, 'cluster_members': ','.join(labels) })
    return DataFrame(data).set_index('label')

taxa = read_taxa(nexus_file)
clstr = read_clstr(clstr_file)
data = read_csv(data_file, index_col = 'label').join(clstr).loc[taxa]

with open(output_file, 'w') as out_nex:
    with open(nexus_file) as file:
        for line in file:
            out_nex.write(line)

    n_traits = len(data.columns)
    trait_labels = ' '.join(data.columns)
    matrix = data.to_csv(sep = "\t", header = False, index = False, na_rep = '?').replace(' ', '_')
    out_nex.write(textwrap.dedent(f"""
        BEGIN Traits;
        DIMENSIONS nTRAITS={n_traits};
        FORMAT LABELS=no MISSING=? SEPARATOR=TAB;
        TRAITLABELS {trait_labels};
    """))
    out_nex.write(f"MATRIX\n{matrix};\n")
    out_nex.write("END;\n")
