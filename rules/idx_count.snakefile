
rule create_whitelist:
    input: barcode_tsv
    output: temp("whitelist_barcodes.txt")
    shell: "cut -f2 {input} > {output}"

# the tr2g file created here would be needed everywhere, other files are needed only
# for the velocity
rule prep_annotation:
    input: input_gtf
    output:
        velocity_fa = "annotations/cDNA_introns.fa",
        cDNA_fa = "annotations/cDNA.fa",
        cdna = "annotations/cDNA_tx_to_capture.txt",
        introns = "annotations/introns_tx_to_capture.txt",
        tr2g = "annotations/tr2g.tsv",
        txdb = "annotations/gtf.txdb"
    params:
        genome = genome_id,
        readLength = read_length,
        script=os.path.join(workflow.basedir, "tools", "prep_annotation.R"),
        outdir = "annotations"
    log:
        out = "logs/prep_annotation.out",
        err = "logs/prep_annotation.err"
    threads: 2
    conda: CONDA_BUS_ENV
    shell: "{params.script} {params.genome} {input} {params.readLength} {params.outdir} > {log.out} 2> {log.err}"

rule transcript_index:
    input: "annotations/cDNA.fa"
    output: "annotations/cDNA.idx"
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
        idx = "annotations/cDNA.idx"
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

# instead of counting using bustools. I count using the text files myself
rule get_counts_txt:
    input:
        "transcripts_quant/{sample}/output.correct.sort.bus"
    output:
        "transcripts_quant/{sample}/output.correct.sort.txt"
    log: "logs/get_counts.{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "bustools text -o {output} {input} > {log}"

rule get_counts:
    input:
        "transcripts_quant/{sample}/output.correct.sort.txt"
    output:
         mtx = "transcripts_quant/{sample}/eq_counts/output.mtx",
         ec = "transcripts_quant/{sample}/eq_counts/output.ec.txt",
         bc = "transcripts_quant/{sample}/eq_counts/output.barcodes.txt"
    params:
        bashScript = os.path.join(workflow.basedir, "tools", "count.sh"),
        bc_index = temp("transcripts_quant/{sample}/eq_counts/bc_index.txt"),
        ec_index = temp("transcripts_quant/{sample}/eq_counts/ec_index.txt"),
        header = temp("transcripts_quant/{sample}/eq_counts/mtx_header.txt"),
        tmp1 = temp("transcripts_quant/{sample}/eq_counts/tmp1.txt"),
        tmp2 = temp("transcripts_quant/{sample}/eq_counts/tmp2.txt")
    log: "logs/get_counts_{sample}.out"
    conda: CONDA_SHARED_ENV
    shell:
        "{params.bashScript} {input} {output.mtx} {output.ec} {output.bc} \
        {params.bc_index} {params.ec_index} {params.header} {params.tmp1} {params.tmp2} > {log} 2>&1"

        # """
        # bustools text -p {input} 2> {log}| \
        # awk -F "\t" 'BEGIN{{OFS=FS}} {{print($1, $3)}}' - | sort | uniq -c | \
        # awk -F " " 'BEGIN{{OFS="\t"}} {{print($2, $3, $1)}}' > {params.tmp1} &&\
        # echo "still going 1" >> {log} &&\
        # cut -f 1 {params.tmp1} | sort | uniq | awk 'BEGIN{{OFS="\t"}} {{print($1, NR)}}' > {params.bc_index} &&\
        # echo "still going 2" >> {log} &&\
        # awk 'FNR==NR{{a[$1]=$2; next}} {{print(a[$1], $2, $3)}}' {params.bc_index} {params.tmp1} > {params.tmp2} &&\
        # echo "still going 3" >> {log} &&\
        # cut -f 1 {params.bc_index} > {output.bc} &&\
        # echo "still going 4" >> {log} &&\
        # cut -d " " -f 2 {params.tmp2} | sort -n | uniq | awk 'BEGIN{{OFS="\t"}} {{print($1, NR)}}' > {params.ec_index} &&\
        # echo "still going 5" >> {log} &&\
        # awk 'FNR==NR{{a[$1]=$2; next}} {{print($1, a[$2], $3)}}' {params.ec_index} {params.tmp2} > {params.tmp1} &&\
        # tail -n 1 {params.bc_index} | cut -f 2 > {params.header} &&\
        # echo "still going 6" >> {log} &&\
        # cut -d " " -f 2 {params.tmp1} | sort -n -r | head -n 1 >> {log} &&\
        # cut -d " " -f 2 {params.tmp1} | sort -n -r | head -n 1 >> {params.header} &&\
        # awk 'END{{print(NR)}}' {params.tmp1} >> {log} &&\
        # echo "still going 7" >> {log} &&\
        # awk 'END{{print(NR)}}' {params.tmp1} >> {params.header} &&\
        # echo 'still going 8' >> {log} &&\
        # echo '%%MatrixMarket matrix coordinate real general' > {output.mtx} &&\
        # awk 'BEGIN{{ORS=" "}}1' {params.header} >> {output.mtx} &&\
        # echo "still going 9" >> {log} &&\
        # echo '' >> {output.mtx} &&\
        # echo "still going 10" >> {log} &&\
        # cat {params.tmp1} >> {output.mtx} &&\
        # echo "still going 11" >> {log} &&\
        # cut -f 1 {params.ec_index} > {output.ec} &&\
        # echo 'done' >> {log}
        # """
        #"{params.bashScript} {input} {params.out} 2> {log}" # {output.mtx} {output.bc} {output.ec}

## this rule gets TCCs and also creates EC to gene map
## NOTE: the TCC matrix output is not filtered for multi-genic ECs, but the
## ECtoGene map is, the rule "merge TCCs" would filter all samples for the multigenic
## ECs before merging
rule get_tcc:
    input:
        tr2g = "annotations/tr2g.tsv",
        ecToTr = "transcripts_quant/{sample}/matrix.ec",
        transcripts = "transcripts_quant/{sample}/transcripts.txt",
        ecList = "transcripts_quant/{sample}/eq_counts/output.ec.txt"
    output:
         ecToGene = "transcripts_quant/{sample}/eq_counts/ec-to-gene.txt"
    params:
        out = "transcripts_quant/{sample}/eq_counts/",
        rscript = os.path.join(workflow.basedir, "tools", "get_ec_geneMap.R")
    log: "logs/get_tcc_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
    #"{params.bustools} count -o {params.out} -g {input.t2g} -e {input.mtx} -t {input.transcripts} {input.busfile} > {log} 2>&1"
        "Rscript {params.rscript} {input.tr2g} {input.transcripts} {input.ecList} {input.ecToTr} {params.out} 2> {log} 2>&1"

rule get_geneCounts:
    input:
        ecToGene = "transcripts_quant/{sample}/eq_counts/ec-to-gene.txt",
        mtx = "transcripts_quant/{sample}/eq_counts/output.mtx",
        ec = "transcripts_quant/{sample}/eq_counts/output.ec.txt",
        bc = "transcripts_quant/{sample}/eq_counts/output.barcodes.txt"
    output:
         mtx = "transcripts_quant/{sample}/gene_counts/output.mtx",
         bc = "transcripts_quant/{sample}/gene_counts/output.barcodes.txt",
         genes = "transcripts_quant/{sample}/gene_counts/output.genes.txt"
    params:
        out = "transcripts_quant/{sample}/gene_counts/",
        rscript = os.path.join(workflow.basedir, "tools", "get_geneCounts.R")
    log: "logs/get_geneCounts.{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:  #"{params.bustools} count --genecounts -o {params.out} -g {input.t2g} -e {input.mtx} -t {input.transcripts} {input.busfile} > {log} 2>&1"
        "Rscript {params.rscript} {input.ecToGene} {input.mtx} {input.ec} {input.bc} {params.out} > {log} 2>&1"
