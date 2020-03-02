
#rule get_ec_geneMap:
#    input:
#        tr2g = "annotations/tr2g.tsv",
#        txlist = "transcripts_quant/{sample}/transcripts.txt",
#        ecToTx = "transcripts_quant/{sample}/eq_counts/output.ec.txt"
#    output:
#        temp("transcripts_quant/{sample}/eq_counts/ec-to-gene.txt")
#    params:
#        rscript = os.path.join(workflow.basedir, "tools", "get_ec_geneMap.R")
#    log: "logs/get_ec_geneMap.{sample}.out"
#    threads: 1
#    conda: CONDA_SHARED_ENV
#    shell:
#        "Rscript {params.rscript} {input.tr2g} {input.txlist} {input.ecToTx} {output} 2> {log} 2>&1"

## Use awk to subset the EC matrix and get the transcripts
# awk -v FS="\t" 'NR==FNR{rows[$1]++;next} ($1 in rows)' \
# output_ESC.txt ../transcripts_quant/ESC_merged/eq_counts/tcc.ec.txt > tcc_subset.txt

## paste matching EC records
# paste output_ESC.txt tcc_subset.txt | \
# awk 'OFS="\t" {if($3 == $1) {print $1,$2,$4 > "ESC_EC_subset_merged.txt"} else {print $0 > "ESC_EC_subset_unmerged.txt"}}'

#rule merge_ECmap:
#    input:
#        ecToTx = "transcripts_quant/{sample}/matrix.ec",
#        ecToGene = "transcripts_quant/{sample}/eq_counts/ec-to-gene.txt"
#    output:
#        "transcripts_quant/{sample}/eq_counts/ECtoGene_map.txt"
#    params:
#        tmpfile = 'transcripts_quant/{sample}/ec_subset.tmp',
#        mismatchFile = 'transcripts_quant/{sample}/ec_mismatch.txt'
#    log: "logs/merge_ECmap.{sample}.out"
#    threads: 1
#    conda: CONDA_SHARED_ENV
#    shell:
#        """
#        awk -v FS="\\t" 'NR==FNR {{ rows[$1]++;next }} ($1 in rows)' \
#        {input.ecToGene} {input.ecToTx} > {params.tmpfile} &&
        ##
#        paste {input.ecToGene} {params.tmpfile} | \
#        awk 'OFS="\\t" {{ if($3 == $1) {{ print $1,$2,$4 > "{output}" }} else {{ print $0 > "{params.mismatchFile}" }} }}' \
#        2> {log} 2>&1
#        """

rule merge_TCCs:
    input:
        mtxList = expand("transcripts_quant/{sample}/eq_counts/output.mtx", sample = samples),
        ecToGeneList = expand("transcripts_quant/{sample}/eq_counts/ec-to-gene.txt", sample = samples),
        ecList = expand("transcripts_quant/{sample}/eq_counts/output.ec.txt", sample = samples),
        bcList = expand("transcripts_quant/{sample}/eq_counts/output.barcodes.txt", sample = samples)
    output:
        mtx = "transcripts_quant/TCCs_filtered_merged.mtx",
        ECmap = "transcripts_quant/ECs_filtered_merged.txt",
        bc = "transcripts_quant/barcodes_merged.txt"
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

rule cluster_tcc:
    input:
        mtx = "transcripts_quant/TCCs_filtered_merged.mtx",
        ECmap = "transcripts_quant/ECs_filtered_merged.txt",
        bc = "transcripts_quant/barcodes_merged.txt"
    output:
        #preprocessed = "clustering_tcc/preprocessed.tsv",
        adata = "clustering_tcc/anndata_filtered.loom",
        cluster = "clustering_tcc/cluster.mtx",
        cl_bc = "clustering_tcc/barcode_cluster.tsv",
        cl_var = "clustering_tcc/var_cluster.tsv",
        preprocessed_fig = "clustering_tcc/preprocessed.pdf",
        cluster_fig = "clustering_tcc/clustering.pdf"
    params:
        clustering = os.path.join(workflow.basedir, "tools", "clustering_wrapper.py"),
        out_dir = "clustering_tcc",
        extendedVar = extendedVar
    log: "logs/cluster_tcc.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "{params.clustering} -o {params.out_dir} -s {input.mtx} -b {input.bc} -v {input.ECmap} \
        --cells 5 --count 1000 --genes 100 --normalize 1e4 -g {col_groups} \
        -t ECs -hv -ev {params.extendedVar} > {log} 2>&1"

rule merge_genes:
    input:
        mtx = expand("transcripts_quant/{sample}/gene_counts/output.mtx", sample = samples),
        barcodes = expand("transcripts_quant/{sample}/gene_counts/output.barcodes.txt", sample = samples),
        genes = expand("transcripts_quant/{sample}/gene_counts/output.genes.txt", sample = samples)
    output:
        merged_mtx = "transcripts_quant/gene_merged.mtx",
        merged_bc = "transcripts_quant/barcodes_gene_merged.txt",
        merged_genes = "transcripts_quant/genes_gene_merged.txt"
    params:
        pyscript = os.path.join(workflow.basedir, "tools", "merge_genes_wrapper.py"),
        out_dir =  "transcripts_quant/",
        samples = ",".join(expand("{sample}", sample = samples))
    log: "logs/merge_genes.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "{params.pyscript} -s {params.samples} -o {params.out_dir} > {log} 2>&1"

rule cluster_genes:
    input:
        mtx = "transcripts_quant/gene_merged.mtx",
        bc = "transcripts_quant/barcodes_gene_merged.txt",
        genes = "transcripts_quant/genes_gene_merged.txt"
    output:
        #preprocessed = "clustering_genes/preprocessed.tsv",
        adata = "clustering_genes/anndata_filtered.loom",
        cluster = "clustering_genes/cluster.mtx",
        cl_bc = "clustering_genes/barcode_cluster.tsv",
        cl_var = "clustering_genes/var_cluster.tsv",
        preprocessed_fig = "clustering_genes/preprocessed.pdf",
        cluster_fig = "clustering_genes/clustering.pdf"
    params:
        clustering = os.path.join(workflow.basedir, "tools", "clustering_wrapper.py"),
        out_dir = "clustering_genes",
        extendedVar = extendedVar
    log: "logs/cluster_genes.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "{params.clustering} -o {params.out_dir} -s {input.mtx} -b {input.bc} -v {input.genes} \
        --cells 5 --count 1000 --genes 100 --dispersity 0.5 --normalize 1e4 -g {col_groups} \
        -t genes -hv -ev {params.extendedVar} > {log} 2>&1"
