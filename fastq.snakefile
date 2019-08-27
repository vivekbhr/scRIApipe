downsample = 10000
trimmerOptions = None
reads = ['', '']


rule origFASTQ1:
      input:
          indir+"/{sample}"+reads[0]+ext
      output:
          "FASTQ/{sample}"+reads[0]+".fastq.gz"
      shell:
          "( [ -f {output} ] || ln -s -r {input} {output} )"

rule origFASTQ2:
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
        "FastQC/{sample}{read}_fastqc.html"
    log:
        out = "FastQC/logs/FastQC.{sample}{read}.out",
        err = "FastQC/logs/FastQC.{sample}{read}.err"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell: "fastqc -o FastQC {input} > {log.out} 2> {log.err}"

if downsample:
    rule FASTQdownsample:
        input:
            r1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = temp("FASTQ/downsample_{sample}"+reads[0]+".fastq.gz"),
            r2 = temp("FASTQ/downsample_{sample}"+reads[1]+".fastq.gz")
        params:
            num_reads = downsample
        threads: 10
        conda: CONDA_SHARED_ENV
        shell:
            """
            seqtk sample -s 100 {input.r1} {params.num_reads} | pigz -p {threads} -9 > {output.r1}
            seqtk sample -s 100 {input.r2} {params.num_reads} | pigz -p {threads} -9 > {output.r2}
            """

## TODO: trimming (cutadapt or fastp?, single or paired?)

rule cutadapt:
    input:
        r1 = "FASTQ/downsample_{sample}"+reads[0]+".fastq.gz" if downsample else "FASTQ/{sample}"+reads[0]+".fastq.gz"
    output:
        "FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz"
    params:
        opts = str(trimmerOptions or '')
    log:
        out = "FASTQ_Cutadapt/logs/Cutadapt.{sample}.out",
        err = "FASTQ_Cutadapt/logs/Cutadapt.{sample}.err"
    benchmark:
        "FASTQ_Cutadapt/.benchmark/Cutadapt.{sample}.benchmark"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        """
        cutadapt -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 25 -a AGATCGGAAGAGC {params.opts} \
            -o "{output}" "{input.r1}" > {log.out} 2> {log.err}
        """
