
checkpoint cyanorak_organisms:
    output:
        "analysis/cyanorak/organisms.csv"
    params:
        url = "https://cyanorak.sb-roscoff.fr/cyanorak/organisms.html",
        attrs = { 'data-toggle': 'table' }
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/html_table.py"

def cyanorak_files(suffix):
    csv_file = str(checkpoints.cyanorak_organisms.get().output)
    df = read_csv(csv_file)
    files = [ f"analysis/cyanorak/{org['Name']}{suffix}" for i, org in df.iterrows() ]
    return files

rule download_organism:
    output:
        "analysis/cyanorak/{organism}.fna"
    params:
        url = "https://cyanorak.sb-roscoff.fr/cyanorak/svc/export/organism/contigs/fna/{organism}"
    shell:
        "wget -qO {output} {params.url}"

rule cyanorak_hmmsearch:
    input:
        faa = "analysis/cyanorak/{organism}.gffread.faa",
        hmm = "analysis/profile/{profile}.hmm"
    output:
        "analysis/cyanorak/{organism}_{profile}.txt"
    conda:
        "envs/search.yaml"
    shell:
        "hmmsearch --max --cpu 1 -o {output} {input.hmm} {input.faa}"

rule cyanorak_hmmsearch_parse:
    input:
        "analysis/cyanorak/{organism}_{profile}.txt"
    output:
        "analysis/cyanorak/{organism}_{profile}.jsonl"
    params:
        threshold = config["nblA"]["threshold"]
    conda:
        "envs/python.yaml"
    script:
        "scripts/parse_hmmsearch.py"

rule cyanorak_proteins:
    input:
        jsonl = "analysis/cyanorak/{organism}_{profile}.jsonl",
        faa = "analysis/cyanorak/{organism}.gffread.faa"
    output:
        "analysis/cyanorak/{organism}_{profile}.faa"
    conda:
        "envs/gene_tools.yaml"
    shell:
        "jq -r .seq_name {input.jsonl} | seqkit grep -f- -o {output} {input.faa}"

rule cyanorak_nbla_align:
    input:
        lambda w: cyanorak_files('_nblA.faa')
    output:
        "analysis/profile/nblA-cyanorak.fas"
    conda:
        "envs/mafft.yaml"
    shell:
        "cat {input} | mafft --auto --reorder - > {output}"

rule cyanorak_fas_a2m:
    input:
        "analysis/profile/nblA-cyanorak.fas"
    output:
        "analysis/profile/nblA-cyanorak.a2m"
    conda:
        "envs/reformat.yaml"
    shell:
        """
        export PATH=$CONDA_PREFIX/scripts:$PATH
        reformat.pl fas a2m {input} {output} -M 50
        """

rule cyanorak_build_hmm:
    input:
        "analysis/profile/nblA-cyanorak.a2m"
    output:
        "analysis/profile/nblA-cyanorak.hmm"
    params:
        profile = "nblA-cyanorak"
    conda:
        "envs/search.yaml"
    shell:
        "hmmbuild -n {params.profile} {output} {input}"

rule all_cyanorak_nblAs:
    input:
        lambda w: cyanorak_files('_nblA-cyanorak.faa')
    output:
        "analysis/nbla/nblA_cyanorak.faa"
    shell:
        "cat {input} > {output}"


