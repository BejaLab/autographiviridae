
def parse_header(line, sequence_id, comment):
    comment = (comment + 1) % 3
    if comment == 0:
        sequence_id = line.split()[1]
    elif comment == 1: 
        assert line.startswith('# gc'), f"Unexpected 'gc' comment: {line}"
    elif comment == 2: 
        assert line.startswith('# self'), f"Unexpected 'self' comment: {line}"
    return sequence_id, comment

def parse_annotation(line, sequence_id):
    fields = line.split()
    gene_id, start, end, strand, frame, completeness, score, model, rbs_start, rbs_end, rbs_score = fields
    attributes = f"ID={gene_id};Name={gene_id};completeness={completeness};model={model};rbs_start={rbs_start};rbs_end={rbs_end};rbs_score={rbs_score}"
    source = 'MetaGeneAnnotator'
    feature = 'CDS'
    if frame != '0':
        if strand == '+':
            start = str(int(start) + int(frame))
        else:
            end = str(int(end) - int(frame))
    phase = '.'
    gff = [ sequence_id, source, feature, start, end, score, strand, phase, attributes]
    return gff

def meta_gene_annotator_to_gff(input_file, output_file):
    with open(input_file, 'r') as file:
        comment = -1
        sequence_id = None
        for line in file:
            if line.startswith('#'):
                sequence_id, comment = parse_header(line, sequence_id, comment)
                yield line.strip()
            else:
                gff = parse_annotation(line, sequence_id)
                yield '\t'.join(gff)

input_file = str(snakemake.input)
output_file = str(snakemake.output)

with open(output_file, 'w') as f:
    f.write("##gff-version 3\n")
    for gff_line in meta_gene_annotator_to_gff(input_file, output_file):
        f.write(gff_line + '\n')
