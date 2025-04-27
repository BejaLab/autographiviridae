rule metagene:
    input:
        "{prefix}.fna"
    output:
        "{prefix}.mga.txt"
    params:
        mode = 's'
    conda:
        "envs/metagene.yaml"
    shell:
        "mga -{params.mode} {input} > {output}"

rule prodigal_single:
    input:
        "{prefix}.fna"
    output:
        "{prefix}.prodigal.gff"
    params:
        mode = 'single'
    shadow:
        "minimal"
    conda:
        "envs/gene_tools.yaml"
    shell:
        "prodigal -i {input} -m -p {params.mode} -f gff -o {output}"

rule metagene_gff:
    input:
        "{prefix}.mga.txt"
    output:
        "{prefix}.mga.gff"
    conda:
        "envs/python.yaml"
    script:
        "scripts/metagene_gff.py"

# NB: genemarks not under conda
rule genemarks_prok:
    input:
        "{prefix}.fna"
    output:
        "{prefix}.genemarks.gff"
    conda:
        "envs/gene_tools.yaml"
    shadow:
        "minimal"
    shell:
        "seqkit seq -i {input} -o seq.fna && gms2.pl --format gff3 --seq seq.fna --genome-type auto --out {output}"

rule gffcompare:
    input:
        "{prefix}.genemarks.gff",
        "{prefix}.prodigal.gff"
    output:
        "{prefix}.gffcmp-raw.gtf"
    shadow:
        "minimal"
    conda:
        "envs/gene_tools.yaml"
    shell:
        "gffcompare -CTo gffcmp {input} && mv gffcmp.combined.gtf {output}"

rule fix_gtf:
    input:
        "{prefix}.gffcmp-raw.gtf"
    output:
        "{prefix}.gffcmp.gtf"
    conda:
        "envs/python.yaml"
    script:
        "scripts/fix_gtf.py"

rule gffread:
    input:
        gtf = "{prefix}.gffcmp.gtf",
        fna = "{prefix}.fna",
        fai = "{prefix}.fna.fai"
    output:
        gff = "{prefix}.gffread.gff",
        faa = "{prefix}.gffread.faa"
    conda:
        "envs/gene_tools.yaml"
    shell:
        "gffread -g {input.fna} -w - -o {output.gff} {input.gtf} | seqkit translate | seqkit replace -sp '[*]$' -o {output.faa}"
