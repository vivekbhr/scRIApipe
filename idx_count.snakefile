
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
        outdir = "annotations"
    log:
        out = "logs/transcript_index.out",
        err = "logs/transcript_index.err"
    threads: 2
    conda: CONDA_scRIA_ENV
    shell: "./velocyto_indexing.R {params.genome} {input} {params.readLength} {params.outdir} > {log.out} 2> {log.err}"

rule transcript_index:
    input: cdna_fasta #"Mus_musculus.GRCm38.cdna.all.fa"
    output:
        "annotations/cdna.all.idx"
    log:
        out = "logs/transcript_index.out",
        err = "logs/transcript_index.err"
    threads: 2
    conda: CONDA_scRIA_ENV
    shell: "kallisto index -i {output} -k 31 {input} > {log.out} 2> {log.err}"

rule velocity_index:
    input:
        "annotations/cDNA_introns.fa"
    output:
        "annotations/cDNA_introns.idx"
    log:
        out = "logs/velocity_index.out",
        err = "logs/velocity_index.err"
    threads: 2
    conda: CONDA_scRIA_ENV
    shell: "kallisto index -i {output} -k 31 {input} > {log.out} 2> {log.err}"

rule transcript_map:
    input:
        R1 = "trimmed_fastq/{sample}_R1.fastq.gz",
        R2 = "trimmed_fastq/{sample}_R2.fastq.gz",
        idx = "annotations/cdna.all.idx"
    output:
        bus = "transcripts_quant/{sample}/output.bus",
        matrix = "transcripts_quant/{sample}/matrix.ec",
        transcripts = "transcripts_quant/{sample}/transcripts.txt"
    params:
        outdir = "transcripts_quant/{sample}"
    log:
        out = "logs/transcript_map.{sample}.out",
        err = "logs/transcript_map.{sample}.err"
    threads: 2
    conda: CONDA_scRIA_ENV
    shell:
    "kallisto bus -i {input.idx} -x CELSeq -t {threads} -o {params.outdir} {input.R1} {input.R2} > {log.out} 2> {log.err}"

rule velocity_map:
    input:
        R1 = "trimmed_fastq/{sample}_R1.fastq.gz",
        R2 = "trimmed_fastq/{sample}_R2.fastq.gz",
        idx = "annotations/cDNA_introns.idx"
    output:
        bus = "velocity_quant/{sample}/output.bus",
        matrix = "velocity_quant/{sample}/matrix.ec",
        transcripts = "velocity_quant/{sample}/transcripts.txt"
    params:
        outdir = "velocity_quant/{sample}"
    log:
        out = "logs/velocity_map.{sample}.out",
        err = "logs/velocity_map.{sample}.err"
    threads: 2
    conda: CONDA_scRIA_ENV
    shell:
    "kallisto bus -i {input.idx} -x CELSeq -t {threads} -o {params.outdir} {input.R1} {input.R2} > {log.out} 2> {log.err}"

rule correct_sort:
    input:
        whitelist = "whitelist_barcodes.txt",
        busfile = "transcripts_quant/{sample}/output.bus"
    output: "transcripts_quant/{sample}/output.correct.sort.bus"
    log:
        out = "logs/correct_sort_{sample}.out",
        err = "logs/correct_sort_{sample}.err"
    threads: 2
    conda: CONDA_scRIA_ENV
    shell:
        """
        bustools correct -w {input.whitelist} \
        -o transcripts_quant/{sample}/output.correct.bus {input.busfile} > {log.out} 2> {log.err};
        bustools sort -t 10 -o {output} transcripts_quant/{sample}/output.correct.bus >> {log.out} 2>> {log.err};
        rm transcripts_quant/{sample}/output.correct.bus
        """

rule velocyto_correct_sort:
    input:
        whitelist = "whitelist_barcodes.txt",
        busfile = "velocity_quant/{sample}/output.bus"
    output:
        "velocity_quant/{sample}/output.correct.sort.bus"
    log:
        out = "logs/correct_sort_velocyto_{sample}.out",
        err = "logs/correct_sort_velocyto_{sample}.err"
    threads: 2
    conda: CONDA_scRIA_ENV
    shell:
        """
        bustools correct -w {input.whitelist} \
        -o velocity_quant/{sample}/output.correct.bus {input.busfile} > {log.out} 2> {log.err};
        bustools sort -t 10 -o {output} velocity_quant/{sample}/output.correct.bus >> {log.out} 2>> {log.err};
        rm velocity_quant/{sample}/output.correct.bus
        """

rule get_tcc:
    input:
        t2g = "annotations/tr2g.tsv",
        mtx = "transcripts_quant/{sample}/matrix.ec",
        transcripts = "annotations/cDNA_tx_to_capture.txt",
        busfile = "transcripts_quant/{sample}/output.correct.sort.bus"
    output:
         mtx = "transcripts_quant/{sample}/eq_counts/tcc.mtx",
         txt = "transcripts_quant/{sample}/eq_counts/tcc.ec.txt"
    params:
        out = "transcripts_quant/{sample}/eq_counts/tcc"
    log:
        out = "logs/get_tcc.out",
        err = "logs/get_tcc.err"
    threads: 2
    conda: CONDA_scRIA_ENV
    shell:  "bustools count -o {params.out} -g {input.t2g} -e {input.mtx} -t {input.transcripts} {input.busfile} > {log.out} 2> {log.err}"

rule get_geneCounts:
    input:
        t2g = "annotations/tr2g.tsv",
        mtx = "transcripts_quant/{sample}/matrix.ec",
        transcripts = "annotations/cDNA_tx_to_capture.txt",
        busfile = "transcripts_quant/{sample}/output.correct.sort.bus"
    output:
         mtx = "transcripts_quant/{sample}/gene_counts/tcc.mtx",
         txt = "transcripts_quant/{sample}/gene_counts/tcc.ec.txt"
    params:
        out = "transcripts_quant/{sample}/gene_counts/gene"
    log:
        out = "logs/get_geneCounts.out",
        err = "logs/get_geneCounts.err"
    threads: 2
    conda: CONDA_scRIA_ENV
    shell:  "bustools count --genecounts -o {params.out} -g {input.t2g} -e {input.mtx} -t {input.transcripts} {input.busfile} > {log.out} 2> {log.err}"

rule get_counts_txt:
    input:
        "annotations/cDNA_introns.fa"
    output:
        "annotations/cDNA_introns.idx"
    log:
        out = "logs/velocity_index.out",
        err = "logs/velocity_index.err"
    threads: 2
    conda: CONDA_scRIA_ENV
    shell: "kallisto index -i {output} -k 31 {input} > {log.out} 2> {log.err}"
