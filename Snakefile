import os
import glob
import yaml
### snakemake_workflows initialization ########################################
maindir = os.path.dirname(os.path.dirname(workflow.basedir))
#workflow_rscripts=os.path.join(maindir, "shared", "rscripts")

## some internal functions  ###################################################
def load_configfile(configfile):
   with open(configfile, "r") as f:
       config = yaml.load(f, Loader=yaml.FullLoader)
   return(config)

def set_condaEnv():
    return{'CONDA_SHARED_ENV': 'env.yaml',
            'CONDA_BUS_ENV' : 'bus_env.yaml'}

def get_sample_names(infiles, ext, reads):
    """
    Get sample names without file extensions
    """
    s = set()
    lext = len(ext)
    l0 = len(reads[0])
    l1 = len(reads[1])
    for x in infiles:
        x = os.path.basename(x)[:-lext]
        if x.endswith(reads[0]):
            x = x[:-l0]
        elif x.endswith(reads[1]):
            x = x[:-l1]
        else:
            continue
        s.add(x)
    return sorted(list(s))

# update envs
globals().update(set_condaEnv())
# load config file
globals().update(load_configfile(workflow.overwrite_configfiles[0]))

## load samples
infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = get_sample_names(infiles,ext,reads)

### include modules of other snakefiles ########################################
################################################################################
include: os.path.join(workflow.basedir, "rules", "fastq.snakefile")
include: os.path.join(workflow.basedir, "rules", "idx_count.snakefile")
include: os.path.join(workflow.basedir, "rules", "tcc_clustering.snakefile")
include: os.path.join(workflow.basedir, "rules", "velocity.snakefile")

### conditional/optional rules #################################################
################################################################################
def runTrimming(trim):
    if trim:
        file_list = [
        expand("FASTQ_trimmed/{sample}{read}.fastq.gz", sample = samples, read = reads),
        expand("FASTQ_trimmed/FastQC/{sample}{read}_fastqc.html", sample = samples, read = reads)
        ]
        return(file_list)
    else:
        return([])

def runVelocity():
    if velocity:
        file_list = [
        expand("velocity_quant/{sample}/output.correct.sort.bus", sample = samples),
        expand("velocity_quant/{sample}/spliced.bus", sample = samples),
        expand("velocity_quant/{sample}/unspliced.bus", sample = samples),
        expand("velocity_quant/{sample}/spliced.txt", sample = samples),
        expand("velocity_quant/{sample}/unspliced.txt", sample = samples),
        expand("velocity_quant/{sample}/TCC_spliced/output.mtx", sample = samples),
        expand("velocity_quant/{sample}/TCC_spliced/ec-to-gene.txt", sample = samples),
        expand("velocity_quant/{sample}/TCC_unspliced/output.mtx", sample = samples),
        expand("velocity_quant/{sample}/TCC_unspliced/ec-to-gene.txt", sample = samples),
        "velocity_output/anndata.loom"
        ]
        return(file_list)
    else:
        return([])

def getIdx(idxOnly):
    file_list = [
    expand("FASTQ/FastQC/{sample}{read}_fastqc.html", sample = samples, read=reads),
    "annotations/cDNA_tx_to_capture.txt",
    "annotations/tr2g.tsv",
    "annotations/gtf.txdb",
    "annotations/cDNA.idx",
    "annotations/cDNA_introns.fa",
    "annotations/cDNA_introns.idx"
    ]
    if not idxOnly:
        file_list += [
        expand("transcripts_quant/{sample}/output.correct.sort.bus", sample = samples),
        expand("transcripts_quant/{sample}/eq_counts/output.mtx", sample = samples),
        expand("transcripts_quant/{sample}/gene_counts/output.mtx", sample = samples),
        expand("transcripts_quant/{sample}/output.correct.sort.txt", sample = samples),
        expand("transcripts_quant/{sample}/eq_counts/ec-to-gene.txt", sample = samples),
        "transcripts_quant/gene_merged.mtx",
        "transcripts_quant/tcc_merged.mtx",
        "clustering_genes/anndata.loom",
        "clustering_tcc/anndata.loom"
        ]
    return(file_list)
### main rule ##################################################################
################################################################################
localrules: FASTQ1, FASTQ2
rule all:
    input:
        runTrimming(trim),
        getIdx(idxOnly),
        runVelocity()

### execute after workflow finished ############################################
################################################################################
onsuccess:
    if "verbose" in config and config["verbose"]:
        print("\n--- scRIA workflow finished successfully! --------------------------------\n")
onerror:
    print("\n !!! ERROR in scRIA workflow! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
