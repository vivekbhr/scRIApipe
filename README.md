# scRIApipe

single-cell **R**NA **I**soform **A**nalysis **Pipe**line

**Author: @vivekbhr**

## DAG

![](./workflow_dag.png)

## Quick Installation

From Conda:

```
conda create -n scria -c vivekbhr -c bioconda -c conda-forge scriapipe
```

From GitHub:

```
conda create -n scria python pip
conda activate scria && pip install snakemake
pip install git+https://github.com/vivekbhr/scRIApipe.git@master
```

## Configure paths

The semi-permanent parameters to the workflow could be configured with a yaml file (config.yaml). Download [here](./scRIApipe/config.yaml)

After conda install, move to the folder where you wan to run the workflow, and prepare the config.yaml

```
cd <your output dir>
vim config.yaml ## now modify the paths as per your requirements
```

The workflow needs
1) path to a GTF file (preferably pre-filter the GTF to remove pseudogene and low confidence annotations)
2) UCSC ID of the genome


cDNA fasta and GTF can be downloaded [here](https://www.ensembl.org/info/data/ftp/index.html)
UCSC ID is, for example "mm10" (mouse) or "hg38" (human)

Additional parameters for the workflow can be accessed via --help

```
scRIA --help
```

## Test-drive the workflow

```
## inside the output dir
scRIA -i <fastq_folder> -o . -c <your>config.yaml -j <jobs> -s ' -np'
```

## 4. Submission parameters

#### Running on HPC Cluster
  - In the workflow command above, **j** is the number of parallel jobs you want to run, **-cl** means submit to cluster (default is to run locally). Therefore if you wish to run the workflow on a cluster, simply use the workflow with the -cl command on the submission node.

  - cluster configuration, such as memory and cluster submission command are placed in [cluster_config.yaml](./cluster_config.yaml), and can be modified to suite the users internal infrastructure.

#### Dry-run
In order to just test what the workflow would do, use the command `-s ' -np' `

#### memory errors
Index builing needs >40G of memory, if the workflow fails and the *logs/velocity_index.err* says something like `std::badalloc`, increase memory in the file `cluster_config.yaml` in the scRIA folder.


### Other technical Notes

  - After running the pipeline, **LOG** file are stored in the **<output>/log/** directory and the workflow top-level log is in scRIA.log file.

  - Currently the -o option is not very flexible and and pipeline works only when it's executed in the output directory.

  - Use the -t argument to specify a local directory for temp files. Default is to use the /tmp/ folder, which might have low space on cluster (unless tmpspace is specified in cluster_config.yaml)

  - **Manual interruption of the workflow**: Simple Ctrl+C is enough to cancel/inturrupt the workflow. However, in some cases re-running the workflow after inturruption might fail with message "Locked working directory". In that case, please run the workflow with `-s ' --unlock'` once.

## Output

Major outputs of the workflow are:

  - Transcript compatibility counts (TCC) in folder `<outdir>/transcripts_quant/<sample>/eq_counts/tcc.mtx`
  - Gene counts in folder `<outdir>/transcripts_quant/<sample>/gene_counts/gene.mtx`
  - RNA velocity output in folder `<outdir>/velocity_output` (normal/filtered loom files, velocity plots)
