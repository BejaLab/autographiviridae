import jsonlines
import csv

jsonl_file = str(snakemake.input)
out_file = str(snakemake.output)

min_dom_score = snakemake.params['dom_score']
overlap = snakemake.params['overlap']

all_data = {}
with jsonlines.open(jsonl_file) as reader:
    for data in reader:
        seq_name = data["seq_name"]
        all_data[seq_name] = data

header = [ "name", "description", "profile", "ali_from", "ali_to", "hmm_from", "hmm_to", "hmm_to_end" ]
with open(out_file, 'w') as file:
    tsv = csv.writer(file, delimiter = "\t")
    tsv.writerow(header)
    for seq_name, data in all_data.items():
        profile = data['profile']['name']
        blocks = []
        for domain in data['domains']:
            if domain['score'] >= min_dom_score:
                blocks.append(( domain['hmm']['from'], domain['hmm']['to'], domain['ali']['from'], domain['ali']['to'] ))
        blocks.sort(key = lambda x: x[0])
        groups = []
        for i in range(len(blocks)):
            hmm1_from, hmm1_to, ali1_from, ali1_to = blocks[i]
            connected = False
            for group in groups:
                if i in group:
                    connected = True
            if not connected:
                for j in range(i + 1, len(blocks)):
                    hmm2_from, hmm2_to, ali2_from, ali2_to = blocks[j]
                    if hmm1_to - hmm2_from < overlap and ali1_to - ali2_from < overlap:
                        chained = False
                        for group in groups:
                            if i in group or j in group:
                                group.add(i)
                                group.add(j)
                                chained = True
                        if not chained:
                            groups.append({i, j})
                        connected = True
                if not connected:
                    groups.append({i})

        ranges = []
        for group in groups:
            hmm_from = min([ blocks[i][0] for i in group ])
            hmm_to   = max([ blocks[i][1] for i in group ])
            ali_from = min([ blocks[i][2] for i in group ])
            ali_to   = max([ blocks[i][3] for i in group ])
            ranges.append((hmm_from, hmm_to, ali_from, ali_to))

        ranges.sort(key = lambda x: x[2] - x[3])
        ranges_done = []
        for hmm_from, hmm_to, ali_from, ali_to in ranges:
            nested = False
            for hmm0_from, hmm0_to, ali0_from, ali0_to in ranges_done:
                if hmm0_from <= hmm_from <= hmm0_to and hmm0_from <= hmm_to <= hmm0_to:
                    nested = True
            if not nested:
                hmm_to_end = data['profile']['len'] - hmm_to
                line = [ seq_name, data['description'], profile, ali_from, ali_to, hmm_from, hmm_to, hmm_to_end ]
                tsv.writerow(line)
                ranges_done.append((hmm_from, hmm_to, ali_from, ali_to))
