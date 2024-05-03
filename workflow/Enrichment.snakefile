
proteome = 'UP000002780'

rule download_uniprot:
    output:
        "analysis/enrichment/uniprot.xml"
    params:
        url = f"https://rest.uniprot.org/uniprotkb/stream?download=true&format=xml&query=proteome:{proteome}"
    shell:
        "wget -O {output} '{params.url}'"

rule extract_go_terms:
    input:
        "analysis/enrichment/uniprot.xml"
    output:
        "analysis/enrichment/go-terms.tsv"
    conda:
        "envs/python.yaml"
    script:
        "scripts/extract_go_terms.py"

rule dload_go:
    output:
        "analysis/enrichment/go-basic.obo"
    params:
        url = "https://current.geneontology.org/ontology/go-basic.obo"
    shell:
        "wget -O {output} {params.url}"

rule enrichment_gomap:
    input:
        tsv = "analysis/enrichment/go-terms.tsv",
        obo = "analysis/enrichment/go-basic.obo"
    output:
        BP = "analysis/enrichment/go-terms-BP.tsv",
        MF = "analysis/enrichment/go-terms-MF.tsv",
        CC = "analysis/enrichment/go-terms-CC.tsv"
    conda:
        "envs/enrichment.yaml"
    script:
        "scripts/enrichment-gomap.R"

rule enrichment:
    input:
        results_file = "input/Day2_NblAOE_For_Go Enrichment.xlsx",
        annot_file = "analysis/enrichment/go-terms-{db}.tsv"
    output:
        "output/enrichment-{db}.pdf"
    params:
        sheet = "SemiN Ratio <-1 only NblA",
        alpha = 0.05
    conda:
        "envs/enrichment.yaml"
    script:
        "scripts/enrichment.R"
