
img_df = read_csv("input/IMG_organisms.csv")

checkpoint cyanorak_organisms:
    output:
        "analysis/picocyano/cyanorak/organisms.csv"
    params:
        url = "https://cyanorak.sb-roscoff.fr/cyanorak/organisms.html",
        attrs = { 'data-toggle': 'table' }
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/html_table.py"

def cyanorak_files(suffix):
    csv_file = checkpoints.cyanorak_organisms.get().output[0]
    df = read_csv(csv_file)
    files = [ f"analysis/picocyano/cyanorak/{org['Name']}{suffix}" for i, org in df.iterrows() ]
    return files

def img_files(suffix):
    df = img_df
    files = [ f"analysis/picocyano/IMG/{org['Name']}{suffix}" for i, org in df.iterrows() ]
    return files

rule download_cyanorak_fasta:
    output:
        "analysis/picocyano/cyanorak/{organism}.fna"
    params:
        url = "https://cyanorak.sb-roscoff.fr/cyanorak/svc/export/organism/contigs/fna/{organism}"
    shell:
        "wget -qO {output} {params.url}"

rule picocyano_hmmsearch:
    input:
        faa = "analysis/picocyano/{source}/{organism}.gffread.faa",
        hmm = "analysis/profile/{profile}.hmm"
    output:
        "analysis/picocyano/{source}/{organism}_{profile}.txt"
    conda:
        "envs/search.yaml"
    shell:
        "hmmsearch --max --cpu 1 -o {output} {input.hmm} {input.faa}"

rule picocyano_hmmsearch_parse:
    input:
        "analysis/picocyano/{source}/{organism}_{profile}.txt"
    output:
        "analysis/picocyano/{source}/{organism}_{profile}.jsonl"
    params:
        threshold = config["nblA"]["threshold"]
    conda:
        "envs/python.yaml"
    script:
        "scripts/parse_hmmsearch.py"

rule picocyano_proteins:
    input:
        jsonl = "analysis/picocyano/{source}/{organism}_{profile}.jsonl",
        faa = "analysis/picocyano/{source}/{organism}.gffread.faa"
    output:
        "analysis/picocyano/{source}/{organism}_{profile}.faa"
    conda:
        "envs/gene_tools.yaml"
    shell:
        "jq -r .seq_name {input.jsonl} | seqkit grep -f- {input.faa} | seqkit rmdup -s -o {output}"

rule picocyano_nbla_align:
    input:
        lambda w: cyanorak_files('_nblA.faa') + img_files('_nblA.faa')
    output:
        "analysis/profile/nblA-picocyano.fas"
    conda:
        "envs/mafft.yaml"
    shell:
        "cat {input} | mafft --auto --reorder - > {output}"

rule picocyano_fas_a2m:
    input:
        "analysis/profile/nblA-picocyano.fas"
    output:
        "analysis/profile/nblA-picocyano.a2m"
    conda:
        "envs/reformat.yaml"
    shell:
        """
        export PATH=$CONDA_PREFIX/scripts:$PATH
        reformat.pl fas a2m {input} {output} -M 50
        """

rule picocyano_build_hmm:
    input:
        "analysis/profile/nblA-picocyano.a2m"
    output:
        "analysis/profile/nblA-picocyano.hmm"
    params:
        profile = "nblA-cyanorak"
    conda:
        "envs/search.yaml"
    shell:
        "hmmbuild -n {params.profile} {output} {input}"

rule img_rename:
    input:
        lambda w: "input/IMG/{accession}.fna".format(accession = img_df[img_df['Name'] == w.organism]['Accession'].item())
    output:
        "analysis/picocyano/IMG/{organism}.fna"
    conda:
        "envs/gene_tools.yaml"
    shell:
        "seqkit replace -p ^ -r {wildcards.organism}_ -o {output} {input}"

rule all_picocyano_nblAs:
    input:
        lambda w: cyanorak_files("_nblA-picocyano.faa") + img_files("_nblA-picocyano.faa")
    output:
        "analysis/nbla/nblA_picocyano.faa"
    shell:
        "cat {input} > {output}"
