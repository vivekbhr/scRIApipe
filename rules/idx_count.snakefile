
# wget ftp://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz


rule create_whitelist:
    input: barcode_tsv #/hpc/hub_oudenaarden/vbhardwaj/annotations/celseq2_barcodes.csv
    output: temp("whitelist_barcodes.txt")
    shell: "cut -f2 {input} > {output}"

rule prep_velocity_files:
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
        out = "logs/prep_velocyto.out",
        err = "logs/prep_velocyto.err"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell: "{params.script} {params.genome} {input} {params.readLength} {params.outdir} > {log.out} 2> {log.err}"

rule transcript_index:
    input: cdna_fasta #"Mus_musculus.GRCm38.cdna.all.fa"
    output:
        "annotations/cdna.all.idx"
    params:
        kallisto = os.path.join(workflow.basedir, "tools", "kallisto")
    log:
        out = "logs/transcript_index.out",
        err = "logs/transcript_index.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "{params.kallisto} index -i {output} -k 31 {input} > {log.out} 2> {log.err}"

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

rule transcript_map:
    input:
        R1 = "FASTQ_trimmed/{sample}"+reads[0]+".fastq.gz" if trim else "FASTQ/{sample}"+reads[0]+".fastq.gz",
        R2 = "FASTQ_trimmed/{sample}"+reads[1]+".fastq.gz" if trim else "FASTQ/{sample}"+reads[1]+".fastq.gz",
        idx = "annotations/cdna.all.idx"
    output:
        bus = "transcripts_quant/{sample}/output.bus",
        matrix = "transcripts_quant/{sample}/matrix.ec",
        transcripts = "transcripts_quant/{sample}/transcripts.txt"
    params:
        outdir = "transcripts_quant/{sample}",
        kallisto = os.path.join(workflow.basedir, "tools", "kallisto"),
        protocol = 'VASASeq'
    log:
        out = "logs/transcript_map.{sample}.out",
        err = "logs/transcript_map.{sample}.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell: "{params.kallisto} bus -i {input.idx} -x {params.protocol} -t {threads} -o {params.outdir} {input.R1} {input.R2} > {log.out} 2> {log.err}"

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
        kallisto = os.path.join(workflow.basedir, "tools", "kallisto")
    log:
        out = "logs/velocity_map.{sample}.out",
        err = "logs/velocity_map.{sample}.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "{params.kallisto} bus -i {input.idx} -x VASASeq -t {threads} -o {params.outdir} {input.R1} {input.R2} > {log.out} 2> {log.err}"

rule correct_sort:
    input:
        whitelist = "whitelist_barcodes.txt",
        busfile = "transcripts_quant/{sample}/output.bus"
    output: "transcripts_quant/{sample}/output.correct.sort.bus"
    params:
        outdir = 'transcripts_quant/{sample}'
    log:
        out = "logs/correct_sort_{sample}.out",
        err = "logs/correct_sort_{sample}.err"
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

rule velocyto_correct_sort:
    input:
        whitelist = "whitelist_barcodes.txt",
        busfile = "velocity_quant/{sample}/output.bus"
    output:
        "velocity_quant/{sample}/output.correct.sort.bus"
    params:
        outdir = 'velocity_quant/{sample}'
    log:
        out = "logs/correct_sort_velocyto_{sample}.out",
        err = "logs/correct_sort_velocyto_{sample}.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        """
        mkdir -p {params.outdir};
        bustools correct -w "{input.whitelist}" \
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
         mtx = "transcripts_quant/{sample}/eq_counts/tcc.mtx",
         txt = "transcripts_quant/{sample}/eq_counts/tcc.ec.txt"
    params:
        out = "transcripts_quant/{sample}/eq_counts/tcc"
    log:
        out = "logs/get_tcc_{sample}.out",
        err = "logs/get_tcc_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:  "bustools count -o {params.out} -g {input.t2g} -e {input.mtx} -t {input.transcripts} {input.busfile} > {log.out} 2> {log.err}"

rule get_geneCounts:
    input:
        t2g = "annotations/tr2g.tsv",
        mtx = "transcripts_quant/{sample}/matrix.ec",
        transcripts = "transcripts_quant/{sample}/transcripts.txt",
        busfile = "transcripts_quant/{sample}/output.correct.sort.bus"
    output:
         mtx = "transcripts_quant/{sample}/gene_counts/gene.mtx"
    params:
        out = "transcripts_quant/{sample}/gene_counts/gene"
    log:
        out = "logs/get_geneCounts_{sample}.out",
        err = "logs/get_geneCounts_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:  "bustools count --genecounts -o {params.out} -g {input.t2g} -e {input.mtx} -t {input.transcripts} {input.busfile} > {log.out} 2> {log.err}"

rule get_counts_txt:
    input:
        "transcripts_quant/{sample}/output.correct.sort.bus"
    output:
        "transcripts_quant/{sample}/output.txt"
    log: "logs/get_counts_{sample}.out",
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "bustools text -o {output} {input} > {log}"
