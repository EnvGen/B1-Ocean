scattergather:
    assemsplit = config["annotation"]["assembly_splits"]

localrules:
    split_assembly,
    contigtax_gather

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
        "python workflow/scripts/split_fasta_file.py {input} -n {params.n_files} -o {params.outdir} --suffix fa > {log} 2>&1"

rule contigtax_search_split:
    input:
        fasta="results/assembly/{assembly}/splits/split_{scatteritem}.fa",
        db=expand("resources/{db}/diamond.dmnd",db=config["taxonomy"]["database"])
    output:
        "results/annotation/{assembly}/taxonomy/splits/split_{scatteritem}.final_contigs.tsv.gz"
    log:
        "results/annotation/{assembly}/taxonomy/splits/split_{scatteritem}.contigtax.log"
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
        tsv="results/annotation/{assembly}/taxonomy/splits/split_{scatteritem}.final_contigs.tsv.gz",
        sql=ancient("resources/taxonomy/taxonomy.sqlite")
    output:
        "results/annotation/{assembly}/taxonomy/splits/split_{scatteritem}.contigtax.taxonomy.tsv"
    log:
        "results/annotation/{assembly}/taxonomy/splits/split_{scatteritem}.contigtax.taxonomy.log"
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
        gather.assemsplit("results/annotation/{{assembly}}/taxonomy/splits/split_{scatteritem}.contigtax.taxonomy.tsv")
    output:
        results+"/annotation/{assembly}/taxonomy/contigtax.taxonomy.tsv",
        touch("results/annotation/{assembly}/taxonomy/{assembly}.contigtax.gathered")
    log:
        results+"/annotation/{assembly}/taxonomy/contigtax_gather.log"
    shell:
        """
        exec &>{log}
        head -1 {input[0]} > {output[0]}
        cat {input} | egrep -v "^query" >> {output[0]}
        """
