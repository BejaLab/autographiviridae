
metag_genes = {
    "terL": 1.8,
    "nbla_ref1": 0.75,
    "nbla_ref2": 1.0
}

sras ,= glob_wildcards("viral_contig_databases/GOV2/fastq/{sra}_1.fastq.gz")

rule metagenomics:
    input:
        expand("analysis/gov2/bwa/{gene}/{sra}.txt", gene = metag_genes.keys(), sra = sras)
        # expand("analysis/OM-RGC/markers/{marker}.domtblout_{type}.tsv", marker = marker_scores.keys(), type = [ "metaG", "metaT" ]),
        # expand("analysis/OM-RGC/nbla/{ref}.domtblout_{type}.tsv", ref = nbla_scores.keys(), type = [ "metaG", "metaT" ])

rule select_ref:
    input:
        expand("input/target_genes/{gene}.faa", gene = metag_genes.keys())
    output:
        "analysis/gov2/genes/ref.faa"
    params:
        ref = ref_virus
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit grep -rp {params.ref} -o {output} {input}"

rule copy_gene_gov2:
    input:
        "input/target_genes/{gene}.faa"
    output:
        "analysis/gov2/genes/{gene}.faa"
    shell:
        "cp {input} {output}"

rule makedb_gov2:
    input:
        "analysis/gov2/genes/{gene}.faa"
    output:
        "analysis/gov2/genes/{gene}.faa.dmnd"
    conda:
        "envs/diamond.yaml"
    shell:
        "diamond makedb --in {input} -d {input}"

rule makeblastdb_gov2:
    input:
        "analysis/gov2/genes/{gene}.faa"
    output:
        "analysis/gov2/genes/{gene}.faa.pdb"
    conda:
        "envs/tools.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"

rule diamond_blastx_gov2:
    input:
        scafs = "viral_contig_databases/GOV2/scaffolds.fasta",
        dmnd = "analysis/gov2/genes/{gene}.faa.dmnd",
        db = "analysis/gov2/genes/{gene}.faa"
    output:
        "analysis/gov2/diamond_blastx/{gene}.outfmt6"
    conda:
        "envs/tools.yaml"
    threads:
        10
    shell:
        "diamond blastx --sensitive -d {input.db} -q {input.scafs} -o {output} -f 6 -p {threads}"

rule blastx_gov2:
    input:
        scafs = "viral_contig_databases/GOV2/scaffolds.fasta",
        pdb = "analysis/gov2/genes/{gene}.faa.pdb",
        db = "analysis/gov2/genes/{gene}.faa"
    output:
        "analysis/gov2/blastx/{gene}.outfmt6"
    conda:
        "envs/tools.yaml"
    threads:
        10
    shell:
        "blastx -db {input.db} -query {input.scafs} -out {output} -outfmt 6 -num_threads {threads} -max_target_seqs 1 -evalue 0.01"

rule gov2_tblastn:
    input:
        db = "viral_contig_databases/GOV2/scaffolds.fasta",
        ref = "analysis/gov2/genes/ref.faa"
    output:
        "analysis/gov2/tblastn/ref.outfmt6"
    params:
        max_seq = 1000 # 1000000000
    conda:
        "envs/tools.yaml"
    threads:
        workflow.cores
    shell:
        "tblastn -db {input.db} -query {input.ref} -out {output} -outfmt 6 -num_threads {threads} -max_target_seqs {params.max_seq} -evalue 0.01"

rule gov2_stat:
    input:
        "viral_contig_databases/GOV2/fastq/{sra}_1.fastq.gz"
    output:
        "analysis/gov2/fastq/{sra}_1.stat"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit stat -o {output} {input}"

rule gov2_collect:
    input:
        db = "viral_contig_databases/GOV2/scaffolds.fasta",
        tblastn = "analysis/gov2/tblastn/ref.outfmt6"
    output:
        "analysis/gov2/contigs.fna"
    conda:
        "envs/tools.yaml"
    shell:
        "cut -f2 {input.tblastn} | sort -u | blastdbcmd -db {input.db} -entry_batch - > {output}"

rule gov2_hmmsearch:
    input:
        faa = "analysis/gov2/contigs.gffread.faa",
        hmm = "input/target_genes/{gene}.hmm"
    output:
        "analysis/gov2/hmmsearch/{gene}.domtblout"
    conda:
        "envs/tools.yaml"
    threads:
        2
    shell:
        "hmmsearch --cpu {threads} -o /dev/null --domtblout {output} {input.hmm} {input.faa}"

rule gov2_filter_marker:
    input:
        "analysis/gov2/hmmsearch/{gene}.domtblout"
    output:
        "analysis/gov2/hmmsearch/{gene}.domtblout.txt"
    params:
        threshold = lambda w: metag_genes[w.gene]
    conda:
        "envs/r.yaml"
    script:
        "scripts/filter_domtblout.R"

rule gov2_subset_gff:
    input:
        txt = "analysis/gov2/hmmsearch/{gene}.domtblout.txt",
        gff = "analysis/gov2/contigs.gffread.gff"
    output:
        "analysis/gov2/hmmsearch/{gene}.gff"
    shell:
        "grep -f {input.txt} {input.gff} > {output}"

rule gov2_gff_slop:
    input:
        gff = "analysis/gov2/hmmsearch/{gene}.gff",
        fasta = "analysis/gov2/contigs.fna"
    output:
        gff = "analysis/gov2/slop/{gene}.gff",
        fasta = "analysis/gov2/slop/{gene}.fasta"
    params:
        min_contig_len = 1000,
        extension = 100,
        type = "transcript"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/slop_features.py"

rule gov2_gff_to_gtf:
    input:
        "analysis/gov2/slop/{gene}.gff"
    output:
        "analysis/gov2/slop/{gene}.gtf"
    conda:
        "envs/tools.yaml"
    shell:
        "gffread --gtf {input} > {output}"

rule gov2_dload_metadata:
    output:
        "analysis/gov2/metadata.xlsx"
    params:
        url = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867419303411-mmc3.xlsx"
    shell:
        "wget -O {output} {params.url}"

rule tara_dload_tax_profile:
    output:
        "analysis/gov2/taxonomic_profiles.tsv"
    params:
        url = "http://ocean-microbiome.embl.de/data/miTAG.taxonomic.profiles.release.tsv.gz"
    shell:
        "wget -O- {params.url} | gzip -cd > {output}"

rule filter_marker:
    input:
        "analysis/OM-RGC/markers/{marker}.domtblout"
    output:
        "analysis/OM-RGC/markers/{marker}.domtblout.txt"
    params:
        threshold = lambda w: marker_scores[w.marker]
    conda:
        "envs/r.yaml"
    script:
        "scripts/filter_domtblout.R"

rule hmmsearch_nblA:
    input:
        faa = "databases/OM-RGC/OM-RGC_v2_orfs.faa",
        hmm = "analysis/nbla_refs/nbla_{ref}.hmm"
    output:
        "analysis/OM-RGC/nbla/{ref}.domtblout"
    conda:
        "envs/tools.yaml"
    threads:
        2
    shell:
        "hmmsearch --cpu {threads} -o /dev/null --domtblout {output} {input.hmm} {input.faa}"

rule filter_nbla:
    input:
        "analysis/OM-RGC/nbla/{ref}.domtblout"
    output:
        "analysis/OM-RGC/nbla/{ref}.domtblout.txt"
    params:
        threshold = lambda w: nbla_scores[w.ref]
    conda:
        "envs/r.yaml"
    script:
        "scripts/filter_domtblout.R"

rule bwa_index:
    input:
        "{prefix}"
    output:
        "{prefix}.bwt"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa index {input}"

rule bwa_mem:
    input:
        fna = "analysis/gov2/slop/{gene}.fasta",
        bwa = "analysis/gov2/slop/{gene}.fasta.bwt",
        fq1 = "viral_contig_databases/GOV2/fastq/{sra}_1.fastq.gz",
        fq2 = "viral_contig_databases/GOV2/fastq/{sra}_2.fastq.gz"
    output:
        "analysis/gov2/bwa/{gene}/{sra}.bam"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa mem {input.fna} {input.fq1} {input.fq2} | samtools view -F4 -bS - | samtools sort -o {output}"

rule get_profile:
    input:
        tsv = "databases/OM-RGC/OM-RGC_v2_gene_profile_{type}.tsv",
        txt = "analysis/OM-RGC/{prefix}.domtblout.txt"
    output:
        "analysis/OM-RGC/{prefix}.domtblout_{type}.tsv"
    shell:
        "(head -n1 {input.tsv}; grep -f {input.txt} {input.tsv}) > {output}"

rule feature_counts:
    input:
        bam = "analysis/gov2/bwa/{gene}/{sra}.bam",
        gtf = "analysis/gov2/slop/{gene}.gtf",
        fasta = "analysis/gov2/slop/{gene}.fasta"
    output:
        "analysis/gov2/bwa/{gene}/{sra}.txt"
    params:
        min_overlap = 10
    conda:
        "envs/subread.yaml"
    shell:
        """
        reads=$(samtools view {input.bam} | wc -l)
        if [ "$reads" -gt 0 ]; then
            featureCounts -p --minOverlap {params.min_overlap} --countReadPairs --fraction -O -M -a {input.gtf} -o {output} -R CORE {input.bam}
        else
            touch {output}
        fi
        """

rule dload_gtdb_ssu:
    output:
        "analysis/gov2/ssu/ssu_all_r214.fasta"
    params:
        url = "https://data.gtdb.ecogenomic.org/releases/release214/214.1/genomic_files_all/ssu_all_r214.tar.gz"
    shell:
        "wget -O- {params.url} | tar -xOzf - > {output}"

rule gtdb_ssu_dada2:
    input:
        fasta = "analysis/gov2/ssu/ssu_all_r214.fasta",
        clades = "metadata/clades.xlsx"
    output:
        "analysis/gov2/ssu/ssu_all_r214_{taxon}.fa.gz"
    params:
        taxon = "p__Cyanobacteriota",
        min_len = 800
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/gtdb_dada2.py"

rule dload_tara_ssu_cyanobacteria:
    output:
        "analysis/gov2/ssu/16S.OTU.SILVA.reference.sequences_{taxon}.fna"
    params:
        url = "http://ocean-microbiome.embl.de/data/16S.OTU.SILVA.reference.sequences.fna.gz",
        taxon = "Cyanobacteria"
    conda:
        "envs/tools.yaml"
    shell:
        "wget -O- {params.url} | seqkit grep -rnp '{wildcards.taxon};' -o {output}"

rule annotate_cyanobacteria:
    input:
        query = "analysis/gov2/ssu/16S.OTU.SILVA.reference.sequences_Cyanobacteria.fna",
        db = "analysis/gov2/ssu/ssu_all_r214_p__Cyanobacteriota.fa.gz"
    output:
        "analysis/gov2/ssu/16S.OTU.SILVA.reference.sequences_Cyanobacteria_dada2.csv"
    params:
        min_boot = 80
    conda:
        "envs/dada2.yaml"
    threads:
        20
    script:
        "scripts/dada2.R"

rule estus:
    input:
        "metadata/ESTUs.xlsx"
    output:
        "analysis/gov2/estus.csv"
    conda:
        "envs/tidyverse.yaml"
    script:
        "scripts/estus.R"

rule plot_gov2:
    input:
        txt = expand("analysis/gov2/bwa/{gene}/{sra}.txt", gene = metag_genes.keys(), sra = sras),
        stat = expand("analysis/gov2/fastq/{sra}_1.stat", sra = sras),
        metadata = "analysis/gov2/metadata.xlsx",
        tax_profiles = "analysis/gov2/taxonomic_profiles.tsv",
        estus = "analysis/gov2/estus.csv"
    output:
        "output/gov2.pdf"
    conda:
        "envs/r.yaml"
    script:
        "scripts/plot_gov2.R"
