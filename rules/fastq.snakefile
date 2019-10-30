rule FASTQ1:
      input:
          indir+"/{sample}"+reads[0]+ext
      output:
          "FASTQ/{sample}"+reads[0]+".fastq.gz"
      shell:
          "( [ -f {output} ] || ln -s -r {input} {output} )"

rule FASTQ2:
      input:
          indir+"/{sample}"+reads[1]+ext
      output:
          "FASTQ/{sample}"+reads[1]+".fastq.gz"
      shell:
          "( [ -f {output} ] || ln -s -r {input} {output} )"

rule FastQC:
    input:
        "FASTQ/{sample}{read}.fastq.gz"
    output:
        "FASTQ/FastQC/{sample}{read}_fastqc.html"
    log:
        out = "logs/FastQC.{sample}{read}.out",
        err = "logs/FastQC.{sample}{read}.err"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell: "fastqc -o FastQC {input} > {log.out} 2> {log.err}"

if downsample:
    rule FASTQdownsample:
        input:
            r1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = "FASTQ/downsample_{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ/downsample_{sample}"+reads[1]+".fastq.gz"
        params:
            num_reads = downsample
        threads: 10
        conda: CONDA_SHARED_ENV
        shell:
            """
            seqtk sample -s 100 {input.r1} {params.num_reads} | pigz -p {threads} -9 > {output.r1}
            seqtk sample -s 100 {input.r2} {params.num_reads} | pigz -p {threads} -9 > {output.r2}
            """

if trim:
    rule cutadapt:
        input:
            R1 = "FASTQ/downsample_{sample}"+reads[0]+".fastq.gz" if downsample else "FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2 = "FASTQ/downsample_{sample}"+reads[1]+".fastq.gz" if downsample else "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            R1 = "FASTQ_trimmed/{sample}"+reads[0]+".fastq.gz",
            R2 = "FASTQ_trimmed/{sample}"+reads[1]+".fastq.gz"
        params:
            opts = str(trimmerOptions or '')
        log:
            out = "logs/Cutadapt.{sample}.out",
            err = "logs/Cutadapt.{sample}.err"
        threads: 8
        conda: CONDA_SHARED_ENV
        shell:
            """
            cutadapt -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 10 \
            -a AGATCGGAAGAGC -A AGATCGGAAGAGC --nextseq-trim=16 \
            -b TGGAATTCTCGGGTGCCAAGG -B TGGAATTCTCGGGTGCCAAGG \
            -b ATCTCGTATGCCGTCTTCTGCTTG -B ATCTCGTATGCCGTCTTCTGCTTG \
            -b GTTCAGAGTTCTACAGTCCGACGATC -B GTTCAGAGTTCTACAGTCCGACGATC \
            --match-read-wildcards {params.opts} \
            -o "{output.R1}" -p "{output.R2}" "{input.R1}" "{input.R2}" > {log.out} 2> {log.err}
            """

    rule FastQC_trimmed:
        input:
            "FASTQ_trimmed/{sample}{read}.fastq.gz"
        output:
            "FASTQ_trimmed/FastQC/{sample}{read}_fastqc.html"
        params:
            outdir = "FASTQ_trimmed/FastQC"
        log:
            out = "logs/FastQC_trimmed.{sample}{read}.out",
            err = "logs/FastQC_trimmed.{sample}{read}.err"
        threads: 2
        conda: CONDA_SHARED_ENV
        shell: "fastqc -o {params.outdir} {input} > {log.out} 2> {log.err}"
