
rule nbla_align:
    input:
        "input/profile/nblA.fasta",
        "input/profile/nblA_more.fasta"
    output:
        "analysis/profile/nblA.fas"
    conda:
        "envs/mafft.yaml"
    shell:
        "cat {input} | mafft --auto --reorder - > {output}"

rule fas_a2m:
    input:
        "analysis/profile/nblA.fas"
    output:
        "analysis/profile/nblA.a2m"
    conda:
        "envs/reformat.yaml"
    shell:
        """
        export PATH=$CONDA_PREFIX/scripts:$PATH
        reformat.pl fas a2m {input} {output} -M 50
        """

rule nbla_mafft:
    input:
        "analysis/nbla/nblA_{set}.faa"
    output:
        "analysis/nbla/nblA_{set}.mafft"
    conda:
        "envs/phylophlan.yaml"
    shell:
        "mafft --auto --reorder {input} > {output}"

rule nbla_trimal:
    input:
        "analysis/nbla/nblA_{set}.mafft"
    output:
        "analysis/nbla/nblA_{set}.trimal"
    conda:
        "envs/phylophlan.yaml"
    shell:
        "trimal -in {input} -out {output} -automated1"

rule nbla_tree:
    input:
        "analysis/nbla/nblA_{set}.trimal"
    output:
        "analysis/nbla/nblA_{set}.treefile"
    params:
        seed = 123,
        B = 200,
        prefix = "analysis/nbla/nblA_{set}"
    conda:
        "envs/phylophlan.yaml"
    shell:
        "iqtree -s {input} --prefix {params.prefix} -b {params.B} --seed {params.seed} -redo"

rule plot_nbla:
    input:
        tree = "analysis/nbla/nblA_{set}.treefile",
        clusters = "input/host_clusters.xlsx",
        jtree = "output/phylophlan_{set}.jtree",
        fasta = "analysis/genomes/all_genomes_nblA.faa"
    output:
        plot = "output/nbla_{set}.svg",
        jtree = "output/nbla_{set}.jtree"
    conda:
        "envs/r.yaml"
    script:
        "scripts/plot_nbla.R"

rule dload_profile:
    output:
        "analysis/profile/{profile}.hmm"
    wildcard_constraints:
        profile = "TIGR[0-9.]+"
    params:
        url = lambda w: f"https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/{w.profile}.HMM"
    shell:
        "wget -O {output} {params.url}"

rule dload_pfam:
    output:
        "analysis/profile/{profile}.hmm"
    wildcard_constraints:
        profile = "PF[0-9.]+"
    params:
        url = lambda w: f"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{w.profile}?annotation=hmm"
    shell:
        "wget -O- {params.url} | gzip -cd > {output}"

rule dload_cyanobacteria:
    output:
        "analysis/nbla/nblA_cyanobacteria.faa"
    params:
        # NblAs from Cyanobacteriota[1117] without "other cyanobacteria"[1983111], "Prochlorococcaceae"[1890426] and "Synechococcaceae"[2881426] except Synechococcus elongatus[32046]
        query = "xref:pfam-PF04485 AND taxonomy_id:1117 ((NOT taxonomy_id:1890426 NOT taxonomy_id:2881426) OR taxonomy_id:32046) NOT taxonomy_id:1983111"
    shell:
        "wget -qO {output} 'https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query={params.query}'"

rule nbla_metadata:
    input:
        jtree = "output/phylophlan_markers9.jtree",
        ingroup = "analysis/nbla/nblA_markers9.faa",
        pro_syn = "analysis/nbla/nblA_picocyano.faa",
        cyanorak = "analysis/picocyano/cyanorak/organisms.csv",
        img = "input/IMG_organisms.csv",
        refs = "input/profile/nblA_refs.faa",
        cyanobacteria = "analysis/nbla/nblA_cyanobacteria.faa"
    output:
        "analysis/nbla/nblA_metadata.csv"
    params:
        show_genomes = [ 'Pro_SS120' ]
    conda:
        "envs/r.yaml"
    script:
        "scripts/nbla_metadata.R"

rule logo:
    input:
        aln = "analysis/nbla/nblA_all_cdhit.a3m",
        metadata = "analysis/nbla/nblA_metadata.csv"
    output:
        plot = "output/logo_{set}.pdf",
        data = "output/logo_{set}.csv"
    params:
        min_seqs = 10
    conda:
        "envs/logo.yaml"
    script:
        "scripts/logo.R"

rule build_hmm:
    input:
        "analysis/profile/nblA.a2m"
    output:
        "analysis/profile/nblA.hmm"
    params:
        profile = "nblA"
    conda:
        "envs/search.yaml"
    shell:
        "hmmbuild -n {params.profile} {output} {input}"

rule cat_fasta:
    input:
        "input/profile/nblA_refs.faa",
        "analysis/nbla/nblA_markers9.faa",
        "analysis/nbla/nblA_picocyano.faa",
        "analysis/nbla/nblA_cyanobacteria.faa"
    output:
        "analysis/nbla/nblA_all.faa"
    params:
        min_len = 40
    conda:
        "envs/search.yaml"
    shell:
        "seqkit rmdup {input} | seqkit seq -gm {params.min_len} -o {output}"

rule all_nblas_cdhit:
    input:
        "analysis/nbla/nblA_all.faa"
    output:
        fasta = "analysis/nbla/nblA_all_cdhit.faa",
        clstr = "analysis/nbla/nblA_all_cdhit.faa.clstr"
    params:
        c = 0.90
    conda:
        "envs/cdhit.yaml"
    shell:
        "cd-hit -c {params.c} -d 0 -i {input} -o {output.fasta}"

rule all_nblas_mafft:
    input:
        "analysis/nbla/nblA_all_cdhit.faa"
    output:
        "analysis/nbla/nblA_all_cdhit.mafft"
    shell:
        "mafft --reorder --localpair --maxiterate 1000 {input} > {output}"

rule all_nblas_trim:
    input:
        "analysis/nbla/nblA_all_cdhit.mafft"
    output:
        fasta = "analysis/nbla/nblA_all_cdhit.trim.faa",
        cols = "analysis/nbla/nblA_all_cdhit.trim.txt"
    params:
        mode = "strict",
        block = 15
    conda:
        "envs/phylophlan.yaml"
    shell:
        "trimal -in {input} -{params.mode} -block {params.block} -out {output.fasta} -colnumbering | grep -o '[[:digit:]]*' > {output.cols}"

rule all_nblas_strip:
    input:
        fasta = "analysis/nbla/nblA_all_cdhit.mafft",
        cols = "analysis/nbla/nblA_all_cdhit.trim.txt"
    output:
        "analysis/nbla/nblA_all_cdhit.strip.faa"
    conda:
        "envs/phylophlan.yaml"
    shell:
        """
        readarray COLS < {input.cols}
        trimal -selectcols {{ "${{COLS[0]}}-${{COLS[-1]}}" }} -complementary -in {input.fasta} -out {output}
        """

rule add_ss:
    input:
        "analysis/nbla/nblA_all_cdhit.strip.faa"
    output:
        "analysis/nbla/nblA_all_cdhit.a3m"
    conda:
        "envs/reformat.yaml"
    params:
        M = 50
    shell:
        """
        export PATH=$CONDA_PREFIX/scripts:$PATH
        reformat.pl fas a3m {input} - -v 0 -M {params.M} | hhconsensus -i /dev/stdin -o /dev/stdout | addss.pl -v 0 stdin stdout | grep -v '^#' > {output}
        """

rule splitstree_exec_file:
    input:
        "workflow/resources/splitstree.nex"
    output:
        "analysis/nbla/nblA_all_cdhit_splitstree_exec.nex"
    params:
        char_transform = 'Uncorrected_P'
    shell:
        "char_transform='{params.char_transform}' envsubst < {input} > {output}"

# not in conda
rule splitstree_run:
    input:
        fasta = "analysis/nbla/nblA_all_cdhit.trim.faa",
        nexus = "analysis/nbla/nblA_all_cdhit_splitstree_exec.nex"
    output:
        "analysis/nbla/nbla_all_cdhit.nex" # case is not preserved
    shell:
        "xvfb-run -a SplitsTreeCMD -x 'IMPORT FILE={input.fasta} DATATYPE=PROTEIN; EXECUTE FILE={input.nexus}; UPDATE; SAVE FILE={output}; QUIT;'"

rule network_add_traits:
    input:
        nexus = "analysis/nbla/nbla_all_cdhit.nex",
        data = "analysis/nbla/nblA_metadata.csv",
        clstr = "analysis/nbla/nblA_all_cdhit.faa.clstr"
    output:
        "output/NblA_network.nex"
    conda:
        "envs/python.yaml"
    script:
        "scripts/nexus_add_traits.py"

rule plot_network:
    input:
        network = "analysis/nbla/nbla_all_cdhit.nex",
        metadata = "analysis/nbla/nblA_metadata.csv",
        clstr = "analysis/nbla/nblA_all_cdhit.faa.clstr"
    output:
        plot = "output/NblA_network.svg",
        data = "output/NblA_network.csv"
    conda:
        "envs/r_plot_network.yaml"
    script:
        "scripts/plot_network.R"
