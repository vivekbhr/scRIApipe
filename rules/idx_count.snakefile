
rule create_whitelist:
    input: barcode_tsv
    output: temp("whitelist_barcodes.txt")
    shell: "cut -f2 {input} > {output}"

# the tr2g file created here would be needed everywhere, other files are needed only
# for the velocity
rule prep_annotation:
    input: input_gtf
    output:
        fa = "annotations/cDNA_introns.fa",
        cdna = "annotations/cDNA_tx_to_capture.txt",
        introns = "annotations/introns_tx_to_capture.txt",
        tr2g = "annotations/tr2g.tsv",
        txdb = "annotations/gtf.txdb"
    params:
        genome = genome_id,
        readLength = read_length,
        script=os.path.join(workflow.basedir, "tools", "velocyto_indexing.R"),
        outdir = "annotations"
    log:
        out = "logs/prep_annotation.out",
        err = "logs/prep_annotation.err"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell: "{params.script} {params.genome} {input} {params.readLength} {params.outdir} > {log.out} 2> {log.err}"

rule transcript_index:
    input: cdna_fasta
    output: "annotations/cDNA.all.idx"
#    params:
#        kallisto = os.path.join(workflow.basedir, "tools", "kallisto")
    log:
        out = "logs/transcript_index.out",
        err = "logs/transcript_index.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "kallisto index -i {output} -k 31 {input} > {log.out} 2> {log.err}"

rule transcript_map:
    input:
        R1 = "FASTQ_trimmed/{sample}"+reads[0]+".fastq.gz" if trim else "FASTQ/{sample}"+reads[0]+".fastq.gz",
        R2 = "FASTQ_trimmed/{sample}"+reads[1]+".fastq.gz" if trim else "FASTQ/{sample}"+reads[1]+".fastq.gz",
        idx = "annotations/cDNA.all.idx"
    output:
        bus = "transcripts_quant/{sample}/output.bus",
        matrix = "transcripts_quant/{sample}/matrix.ec",
        transcripts = "transcripts_quant/{sample}/transcripts.txt"
    params:
        outdir = "transcripts_quant/{sample}",
        protocol = lambda wildcards: "0,6,14:0,0,6:1,0,0" if protocol == 'VASASeq' else protocol
    log:
        out = "logs/transcript_map.{sample}.out",
        err = "logs/transcript_map.{sample}.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell: "kallisto bus -i {input.idx} -x {params.protocol} -t {threads} -o {params.outdir} {input.R1} {input.R2} > {log.out} 2> {log.err}"

rule correct_sort:
    input:
        whitelist = "whitelist_barcodes.txt",
        busfile = "transcripts_quant/{sample}/output.bus"
    output: "transcripts_quant/{sample}/output.correct.sort.bus"
    params:
        outdir = 'transcripts_quant/{sample}'
    log:
        out = "logs/correct_sort.{sample}.out",
        err = "logs/correct_sort.{sample}.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        """
        mkdir -p {params.outdir};
        bustools correct -w {input.whitelist} \
        -o {params.outdir}/output.correct.bus {input.busfile} > {log.out} 2> {log.err};
        bustools sort -t {threads} -o {output} {params.outdir}/output.correct.bus >> {log.out} 2>> {log.err};
        rm {params.outdir}/output.correct.bus
        """

rule get_tcc:
    input:
        t2g = "annotations/tr2g.tsv",
        mtx = "transcripts_quant/{sample}/matrix.ec",
        transcripts = "transcripts_quant/{sample}/transcripts.txt",
        busfile = "transcripts_quant/{sample}/output.correct.sort.bus"
    output:
         mtx = "transcripts_quant/{sample}/eq_counts/output.mtx",
         txt = "transcripts_quant/{sample}/eq_counts/output.ec.txt",
         bc = "transcripts_quant/{sample}/eq_counts/output.barcodes.txt"
    params:
        out = "transcripts_quant/{sample}/eq_counts/",
        bustools = os.path.join(workflow.basedir, "tools", "bustools")
    log: "logs/get_tcc_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:  "{params.bustools} count -o {params.out} -g {input.t2g} -e {input.mtx} -t {input.transcripts} {input.busfile} > {log} 2>&1"

rule get_geneCounts:
    input:
        t2g = "annotations/tr2g.tsv",
        mtx = "transcripts_quant/{sample}/matrix.ec",
        transcripts = "transcripts_quant/{sample}/transcripts.txt",
        busfile = "transcripts_quant/{sample}/output.correct.sort.bus"
    output:
         mtx = "transcripts_quant/{sample}/gene_counts/output.mtx",
         txt = "transcripts_quant/{sample}/gene_counts/output.genes.txt",
         bc = "transcripts_quant/{sample}/gene_counts/output.barcodes.txt"
    params:
        out = "transcripts_quant/{sample}/gene_counts/",
        bustools = os.path.join(workflow.basedir, "tools", "bustools")
    log: "logs/get_geneCounts.{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:  "{params.bustools} count --genecounts -o {params.out} -g {input.t2g} -e {input.mtx} -t {input.transcripts} {input.busfile} > {log} 2>&1"

rule get_counts_txt:
    input:
        "transcripts_quant/{sample}/output.correct.sort.bus"
    output:
        "transcripts_quant/{sample}/output.txt"
    log: "logs/get_counts.{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "bustools text -o {output} {input} > {log}"
