rule merge_tcc:
    input:
        mtxList = expand("transcripts_quant/{sample}/eq_counts/output.mtx", sample = samples),
        ecToGeneList = expand("transcripts_quant/{sample}/eq_counts/ec-to-gene.txt", sample = samples),
        ecList = expand("transcripts_quant/{sample}/eq_counts/output.ec.txt", sample = samples),
        bcList = expand("transcripts_quant/{sample}/eq_counts/output.barcodes.txt", sample = samples)
    output:
        mtx = "transcripts_quant/tcc_merged.mtx",
        ECmap = "transcripts_quant/tcc_merged.ECs.txt",
        bc = "transcripts_quant/tcc_merged.barcodes.txt"
    params:
        rscript = os.path.join(workflow.basedir, "tools", "merge_tcc_mtx.R"),
        ecToGeneList = ",".join(expand("transcripts_quant/{sample}/eq_counts/ec-to-gene.txt", sample = samples)),
        mtxList = ",".join(expand("transcripts_quant/{sample}/eq_counts/output.mtx", sample = samples)),
        bcList = ",".join(expand("transcripts_quant/{sample}/eq_counts/output.barcodes.txt", sample = samples)),
        ecList = ",".join(expand("transcripts_quant/{sample}/eq_counts/output.ec.txt", sample = samples)),
        samples = ",".join(expand("{sample}", sample = samples)),
        mergeBy = "TxSet"
    log: "logs/merge_TCCs.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "Rscript {params.rscript} {params.mtxList} {params.ecToGeneList} {params.ecList} {params.bcList} {params.samples} {params.mergeBy} \
         {output.mtx} {output.ECmap} {output.bc} 2> {log} 2>&1"

rule merge_genes:
    input:
        mtx = expand("transcripts_quant/{sample}/gene_counts/output.mtx", sample = samples),
        barcodes = expand("transcripts_quant/{sample}/gene_counts/output.barcodes.txt", sample = samples),
        genes = expand("transcripts_quant/{sample}/gene_counts/output.genes.txt", sample = samples)
    output:
        merged_mtx = "transcripts_quant/gene_merged.mtx",
        merged_bc = "transcripts_quant/gene_merged.barcodes.txt",
        merged_genes = "transcripts_quant/gene_merged.genes.txt"
    params:
        pyscript = os.path.join(workflow.basedir, "tools", "merge_genes_wrapper.py"),
        out_dir =  "transcripts_quant/",
        samples = ",".join(expand("{sample}", sample = samples))
    log: "logs/merge_genes.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "{params.pyscript} -s {params.samples} -o {params.out_dir} > {log} 2>&1"

rule cluster_tcc:
    input:
        mtx = "transcripts_quant/tcc_merged.mtx",
        ECmap = "transcripts_quant/tcc_merged.ECs.txt",
        bc = "transcripts_quant/tcc_merged.barcodes.txt"
    output:
        adata = "clustering_tcc/anndata.loom",
        cluster = "clustering_tcc/cluster.mtx",
        cl_bc = "clustering_tcc/barcode_cluster.tsv",
        cl_var = "clustering_tcc/var_cluster.tsv",
        preprocessed_fig = "clustering_tcc/preprocess_stats.pdf",
        cluster_fig = "clustering_tcc/clustering.pdf"
    params:
        clustering = os.path.join(workflow.basedir, "tools", "clustering_wrapper.py"),
        out_dir = "clustering_tcc",
        annotation = annotation
    log: "logs/cluster_tcc.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "{params.clustering} -o {params.out_dir} -s {input.mtx} -b {input.bc} -v {input.ECmap} \
        -t ECs -hv 5000 -a {params.annotation} > {log} 2>&1"

rule cluster_genes:
    input:
        mtx = "transcripts_quant/gene_merged.mtx",
        bc = "transcripts_quant/gene_merged.barcodes.txt",
        genes = "transcripts_quant/gene_merged.genes.txt"
    output:
        adata = "clustering_genes/anndata.loom",
        cluster = "clustering_genes/cluster.mtx",
        cl_bc = "clustering_genes/barcode_cluster.tsv",
        cl_var = "clustering_genes/var_cluster.tsv",
        preprocessed_fig = "clustering_genes/preprocess_stats.pdf",
        cluster_fig = "clustering_genes/clustering.pdf"
    params:
        clustering = os.path.join(workflow.basedir, "tools", "clustering_wrapper.py"),
        out_dir = "clustering_genes",
        annotation = annotation
    log: "logs/cluster_genes.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "{params.clustering} -o {params.out_dir} -s {input.mtx} -b {input.bc} -v {input.genes} \
         -t genes -hv 5000 -a {params.annotation} > {log} 2>&1"
