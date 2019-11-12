### velocit indexing (prepare files)

rule velocity_index:
    input:
        "annotations/cDNA_introns.fa"
    output:
        "annotations/cDNA_introns.idx"
    params:
        kallisto = os.path.join(workflow.basedir, "tools", "kallisto")
    log:
        out = "logs/velocity_index.out",
        err = "logs/velocity_index.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "{params.kallisto} index -i {output} -k 31 {input} > {log.out} 2> {log.err}"

### Velocity mapping (to combined cDNA-intron index)
rule velocity_map:
    input:
        R1 = "FASTQ_trimmed/{sample}"+reads[0]+".fastq.gz" if trim else "FASTQ/{sample}"+reads[0]+".fastq.gz",
        R2 = "FASTQ_trimmed/{sample}"+reads[1]+".fastq.gz" if trim else "FASTQ/{sample}"+reads[1]+".fastq.gz",
        idx = "annotations/cDNA_introns.idx"
    output:
        bus = "velocity_quant/{sample}/output.bus",
        matrix = "velocity_quant/{sample}/matrix.ec",
        transcripts = "velocity_quant/{sample}/transcripts.txt"
    params:
        outdir = "velocity_quant/{sample}",
        kallisto = os.path.join(workflow.basedir, "tools", "kallisto"),
        protocol = 'VASASeq'
    log:
        out = "logs/velocity_map.{sample}.out",
        err = "logs/velocity_map.{sample}.err"
    threads: 20
    conda: CONDA_SHARED_ENV
    shell:
        "{params.kallisto} bus -i {input.idx} -x {params.protocol} -t {threads} -o {params.outdir} {input.R1} {input.R2} > {log.out} 2> {log.err}"

rule velocyto_correct_sort:
    input:
        whitelist = "whitelist_barcodes.txt",
        busfile = "velocity_quant/{sample}/output.bus"
    output:
        "velocity_quant/{sample}/output.correct.sort.bus"
    params:
        outdir = 'velocity_quant/{sample}'
    log: "logs/correct_sort_velocyto_{sample}.out"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        """
        mkdir -p {params.outdir};
        bustools correct -w "{input.whitelist}" \
        -o {params.outdir}/output.correct.bus {input.busfile} > {log} 2>&1;
        bustools sort -t {threads} -o {output} {params.outdir}/output.correct.bus >> {log};
        rm {params.outdir}/output.correct.bus
        """
## bustools v 0.39.3
## for spliced counts, we subset the bus file to get the complement of the "introns" transcript set
rule spliced_capture:
    input:
        mtx = "velocity_quant/{sample}/matrix.ec",
        introns = "annotations/introns_tx_to_capture.txt",
        busfile = "velocity_quant/{sample}/output.correct.sort.bus",
        transcripts = "velocity_quant/{sample}/transcripts.txt"
    output: "velocity_quant/{sample}/spliced.bus"
    log: "logs/spliced_capture_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bustools capture -s -x -o {output} -c {input.introns} \
        -e {input.mtx} -t {input.transcripts} {input.busfile} > {log} 2>&1"

## for unspliced counts, we subset the bus file to get the complement of the captured set from "cDNA" list
rule unspliced_capture:
    input:
        mtx = "velocity_quant/{sample}/matrix.ec",
        cdna = "annotations/cDNA_tx_to_capture.txt",
        busfile = "velocity_quant/{sample}/output.correct.sort.bus",
        transcripts = "velocity_quant/{sample}/transcripts.txt"
    output: "velocity_quant/{sample}/unspliced.bus"
    log: "logs/unspliced_capture_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bustools capture -s -x -o {output} -c {input.cdna} \
        -e {input.mtx} -t {input.transcripts} {input.busfile} > {log} 2>&1"

rule spliced_counts:
    input:
        t2g = "annotations/tr2g.tsv",
        ec = "velocity_quant/{sample}/matrix.ec",
        bus = "velocity_quant/{sample}/spliced.bus",
        transcripts = "velocity_quant/{sample}/transcripts.txt"
    output: "velocity_quant/{sample}/spliced_counts/spliced.mtx"
    params:
        prefix = "velocity_quant/{sample}/spliced_counts/spliced"
    log: "logs/spliced_counts_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bustools count -o {params.prefix} -g {input.t2g} -e {input.ec} \
        -t {input.transcripts} --genecounts {input.bus} > {log} 2>&1"


rule unspliced_counts:
    input:
        t2g = "annotations/tr2g.tsv",
        ec = "velocity_quant/{sample}/matrix.ec",
        bus = "velocity_quant/{sample}/unspliced.bus",
        transcripts = "velocity_quant/{sample}/transcripts.txt"
    output: "velocity_quant/{sample}/unspliced_counts/unspliced.mtx"
    params:
        prefix = "velocity_quant/{sample}/unspliced_counts/unspliced"
    log: "logs/unspliced_counts_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bustools count -o {params.prefix} -g {input.t2g} -e {input.ec} \
        -t {input.transcripts} --genecounts {input.bus} > {log} 2>&1"

rule velocyto:
    input:
        unspliced = expand("velocity_quant/{sample}/unspliced_counts/unspliced.mtx", sample = samples),
        spliced = expand("velocity_quant/{sample}/spliced_counts/spliced.mtx", sample = samples),
        t2g = "annotations/tr2g.tsv"
    output:
        adata_all = "velocity_output/anndata.loom",
        adata_filt = "velocity_output/anndata_filtered.loom",
        qc_metrics = "velocity_output/qc-metrics.csv",
        velo_fig1 = "velocity_output/velocity-grid_louvain.png",
        velo_fig2 = "velocity_output/velocity-grid_samples.png"
    params:
        scvelo = os.path.join(workflow.basedir, "tools", "scVelo_wrapper.py"),
        samples = " ".join(samples),
        outdir = "velocity_output"
    log: "logs/velocity_report.out"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell:
        "{params.scvelo} -s {params.samples} -o {params.outdir} -t {input.t2g} > {log} 2>&1"



#rule velocyto:
#    input:
#        unspliced = expand("velocity_quant/{sample}/unspliced_counts/unspliced.mtx", sample = samples),
#        spliced = expand("velocity_quant/{sample}/spliced_counts/spliced.mtx", sample = samples)
#    output: "velocity_report.html"
#    params:
#        rscript = os.path.join(workflow.basedir, "tools", "velocity_report.R"),
#        rmd = os.path.join(workflow.basedir, "tools", "velocity_report.Rmd")
#    log:
#        out = "logs/velocity_report.out",
#        err = "logs/velocity_report.err"
#    threads: 2
#    conda: CONDA_SHARED_ENV
#    shell:
#        "cat {params.rscript} | R --vanilla --quiet --args {params.rmd} > {log.out} 2> {log.err}"
