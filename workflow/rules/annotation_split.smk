scattergather:
    split = config["annotation"]["splits"]

localrules:
    split_fasta,
    pfam_scan_gather

rule split_fasta:
    input:
        "results/annotation/{assembly}/final_contigs.faa"
    output:
        expand("results/annotation/{{assembly}}/splits/split_{n}-of-{N}.faa",
        n = list(range(1, config["annotation"]["splits"]+1)),
        N =  config["annotation"]["splits"])
    params:
        n_files = lambda wildcards: config["annotation"]["splits"],
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    log: "results/annotation/{assembly}/splits/splits.log"
    shell:
        "python workflow/scripts/split_fasta_file.py {input} -n {params.n_files} -o {params.outdir} > {log} 2>&1"

rule split_assembly:
    input:
        "results/assembly/{assembly}/final_contigs.fa"
    output:
        expand("results/assembly/{{assembly}}/splits/split_{n}-of-{N}.fa",
        n = list(range(1, config["annotation"]["assembly_splits"]+1)),
        N = config["annotation"]["assembly_splits"])
    params:
        n_files = config["annotation"]["assembly_splits"],
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    log: "results/assembly/{assembly}/splits/splits.log"
    shell:
        "python workflow/scripts/split_fasta_file.py {input} -n {params.n_files} -o {params.outdir} > {log} 2>&1"

rule contigtax_search_split:
    input:
        fasta="results/assembly/{assembly}/splits/split_{scatteritem}.faa",
        db=expand("resources/{db}/diamond.dmnd",db=config["taxonomy"]["database"])
    output:
        "results/annotation/{assembly}/taxonomy/splits/split_{scatteritem}.final_contigs.{db}.tsv.gz"
    log:
        "results/annotation/{assembly}/taxonomy/splits/split_{scatteritem}.{db}.contigtax.log"
    params:
        tmpdir=temppath,
        min_len=config["taxonomy"]["min_len"],
        settings=config["taxonomy"]["search_params"]
    threads: 20
    resources:
        runtime= lambda wildcards,attempt: attempt ** 2 * 60 * 10
    conda:
        "../envs/taxonomy.yml"
    shell:
        """
        contigtax search {params.settings} -p {threads} \
            --tmpdir {params.tmpdir} -l {params.min_len} \
            {input.fasta} {input.db} {output} >{log} 2>&1
        """

rule contigtax_assign_split:
    input:
        tsv="results/annotation/{assembly}/taxonomy/splits/split_{scatteritem}.final_contigs.{db}.tsv.gz",
        sql=ancient("resources/taxonomy/taxonomy.sqlite")
    output:
        "results/annotation/{assembly}/taxonomy/splits/split_{scatteritem}.contigtax.taxonomy.{db}.tsv"
    log:
        "results/annotation/{assembly}/taxonomy/splits/split_{scatteritem}.contigtax.taxonomy.{db}.log"
    params:
        taxonomy_ranks=" ".join(config["taxonomy"]["ranks"]),
        taxdir="resources/taxonomy",
        settings=config["taxonomy"]["assign_params"]
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*6
    conda:
        "../envs/taxonomy.yml"
    shell:
         """
         contigtax assign {params.settings} -p {threads} -m rank_lca \
            --reportranks {params.taxonomy_ranks} -t {params.taxdir} \
            {input.tsv} {output} > {log} 2>&1
         """

rule contigtax_gather:
    input:
        gather.split("results/annotation/{{assembly}}/taxonomy/splits/split_{scatteritem}.contigtax.taxonomy.{{db}}.tsv")
    output:
        touch("results/annotation/{assembly}/taxonomy/{assembly}.contigtax.{db}.gathered")
    shell:
        """
        egrep "^#" {input[0]} > {output[0]}
        cat {input} | egrep -v "^#" >> {output[0]}
        """

rule pfam_scan_split:
    input:
        "results/annotation/{assembly}/splits/split_{scatteritem}.faa",
        expand("resources/pfam/Pfam-A.hmm.h3{suffix}",
               suffix=["f", "i", "m", "p"])
    output:
        "results/annotation/{assembly}/splits/split_{scatteritem}.pfam.out"
    log:
        "results/annotation/{assembly}/splits/split_{scatteritem}.pfam.log"
    conda:
        "../envs/annotation.yml"
    params:
        dir="resources/pfam",
        tmp_out=temppath+"/{assembly}.pfam.out"
    threads: 2
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        pfam_scan.pl -fasta {input[0]} -dir {params.dir} -cpu {threads} \
            -outfile {params.tmp_out} >{log} 2>&1
        mv {params.tmp_out} {output[0]}
        """

rule pfam_scan_gather:
    input:
        gather.split("results/annotation/{{assembly}}/splits/split_{scatteritem}.pfam.out")
    output:
        "results/annotation/{assembly}/{assembly}.pfam.out",
        touch("results/annotation/{assembly}/{assembly}.pfam.gathered")
    shell:
        """
        egrep "^#" {input[0]} > {output[0]}
        cat {input} | egrep -v "^#" >> {output[0]}
        """
