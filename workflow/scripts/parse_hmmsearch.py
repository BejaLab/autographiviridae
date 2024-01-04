import re
import jsonlines

input_file = str(snakemake.input)
output_file = str(snakemake.output)
threshold = snakemake.params['threshold']

domains = {}

hmm_len = None

def parse_hmm_name_len(line):
    """Parse the profile length line."""
    hmm_name, hmm_len = re.findall(r'([^ ]+) +\[M=(\d+)\]', line)[0]
    return hmm_name, int(hmm_len)

def parse_seq_name_desc(line):
    """Parse sequence name and description."""
    seq_name, seq_desc = re.findall(r'>>\s*(\S+)\s*(.*)', line)[0]
    return seq_name, seq_desc.strip()

def parse_complete_data(fh):
    """Parse the full scores block.

    The block ends with an empty line.
    """
    values = {}

    inclusion = True
    while (line := next(fh).strip()) != "":
        if 'inclusion threshold' in line:
            inclusion = False
        if not line.startswith('--') and not line.startswith('E-value'):
            full_E_value, full_score, full_bias, best_E_value, best_score, best_bias, exp, N, Sequence, *rest = line.split()
            values[Sequence] = {
            "full": {
                    "evalue": float(full_E_value),
                    "score": float(full_score),
                    "bias": float(full_bias)
                },
                "best": {
                    "evalue": float(best_E_value),
                    "score": float(best_score),
                    "bias": float(best_bias)
                },
                "exp": float(exp),
                "num_domains": int(N),
                "included": inclusion
            }
    return values

def completeness(compl):
    return [ compl[0] == '[', compl[1] == ']' ]

def parse_domain_data(fh):
    """Parse a block for one input sequence.
    
    The block starts with ">>" and
    ends with last domain's last alignment block.
    """
    domains = []
    no_domains = False

    # The scores section starts immediately and ends with an empty line
    while (line := next(fh).strip()) != '' and not no_domains:
        no_domains = line.startswith("[No individual domains")
        if not line.startswith('#') and not line.startswith('--') and not no_domains:
            num, passed, score, bias, c_evalue, i_evalue, hmm_from, hmm_to, hmm_compl, ali_from, ali_to, ali_compl, env_from, env_to, env_compl, acc = line.split()
            vals = {
                "num": int(num),
                "passed": passed == "!",
                "score": float(score),
                "bias": float(bias),
                "c_evalue": float(c_evalue),
                "i_evalue": float(i_evalue),
                "hmm": {
                    "from": int(hmm_from),
                    "to":   int(hmm_to),
                    "complete": completeness(hmm_compl),
                    "seq": ""
                },
                "ali": {
                    "from": int(ali_from),
                    "to":   int(ali_to),
                    "complete": completeness(ali_compl),
                    "seq": ""
                },
                "env": {
                    "from": int(env_from),
                    "to":   int(env_to),
                    "complete": completeness(env_compl)
                },
                "acc": float(acc),
                "RF": "",
                "PP": "",
                "matches": ""
            }
            domains.append(vals)

    alignments_for_each_domain = next(fh) # skip this line

    i = -1
    while i < len(domains) - 1:
        line = next(fh).strip()
        if line.startswith('=='):
            i += 1
            hmm_right = -1
            while hmm_right < domains[i]['hmm']['to']:
                block = parse_aln_block(fh)
                domains[i]["RF"] += block['RF']
                domains[i]["PP"] += block['PP']
                domains[i]["hmm"]["seq"] += block['hmm_seq']
                domains[i]["ali"]["seq"] += block['ali_seq']
                domains[i]["matches"] += block['matches']
                if block['hmm_right'] != "-":
                    hmm_right = int(block['hmm_right'])
    return(domains)

def parse_aln_block(fh):
    """Parse profile-query alignment (sub)block.

    The block ends with an empty line.
    """
    vals = {}

    line = next(fh)
    fields = re.split(' +', line.strip(), 2)
    vals['RF'] = fields[0]

    line = next(fh)
    fields = re.split(" +", line.strip())
    vals['hmm_left']  = fields[-3]
    vals['hmm_seq']   = fields[-2]
    vals['hmm_right'] = fields[-1]

    line = next(fh)
    vals['matches'] = line.rstrip('\n\r')[-len(vals['hmm_seq']):]

    line = next(fh)
    fields = re.split(" +", line.strip())
    vals['ali_left']  = fields[-3]
    vals['ali_seq']   = fields[-2]
    vals['ali_right'] = fields[-1]
    
    line = next(fh)
    fields = re.split(' +', line.strip(), 2)
    vals['PP'] = fields[0]

    line = next(fh).strip()
    assert line == '', "An empty line expected, got %s" % line

    return vals

# Parse the hmmsearch output
with jsonlines.open(output_file, 'w') as writer:
    with open(input_file) as file:
        full_scores = {}
        for line in file:
            if line.startswith('Query:'):
                hmm_name, hmm_len = parse_hmm_name_len(line)
            elif line.startswith('Scores'):
                complete_data = parse_complete_data(file)
            elif line.startswith('>>'):
                seq_name, seq_desc = parse_seq_name_desc(line)
                assert seq_name in complete_data, f"Sequence does not appear in the complete sequences table: {seq_name}"
                data = complete_data[seq_name]
                if data["full"]["score"] >= threshold:
                    data["description"] = seq_desc
                    data["profile"] = {
                        "name": hmm_name,
                        "len": hmm_len
                    }
                    data["domains"] = parse_domain_data(file)
                    data["seq_name"] = seq_name
                    writer.write(data)
