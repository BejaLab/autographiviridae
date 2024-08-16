from pandas import ExcelFile

proteome = 'UP000002780'

results_day1 = ExcelFile('input/Day1_NblAOE_For_Go Enrichment.xlsx')
results_day2 = ExcelFile('input/Day2_NblAOE_For_Go Enrichment.xlsx')

rule all_enrichment:
    input:
        expand("output/enrichment-{sheet}.svg", sheet = "SemiN Ratio <-1 only NblA")

rule download_uniprot:
    output:
        "analysis/enrichment/uniprot.xml"
    params:
        url = f"https://rest.uniprot.org/uniprotkb/stream?download=true&format=xml&query=proteome:{proteome}"
    shell:
        "wget -O {output} '{params.url}'"

rule uniprot_orgdb:
    input:
        "analysis/enrichment/uniprot.xml"
    output:
        directory("analysis/orgDb/org.SWH8109.eg.db")
    conda:
        "envs/orgdb.yaml"
    script:
        "scripts/orgdb.R"

rule enrichment:
    input:
        results = "input/{day}_NblAOE_For_Go Enrichment.xlsx",
        orgdb = "analysis/orgDb/org.SWH8109.eg.db"
    output:
        "analysis/enrichment/{day}-{sheet}.csv"
    wildcard_constraints:
        day = 'Day\d+'
    params:
        cutoff = 0.5,
        alpha = 0.05
    conda:
        "envs/enrichment.yaml"
    script:
        "scripts/enrichment.R"

rule enrichment_plot:
    input:
        expand("analysis/enrichment/{day}-{{sheet}}.csv", day = [ "Day1", "Day2" ])
    output:
        "output/enrichment-{sheet}.svg",
    params:
        alpha = 0.05
    conda:
        "envs/enrichment.yaml"
    script:
        "scripts/enrichment_plot.R"
