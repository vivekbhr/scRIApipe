#!/usr/bin/env python3

__description__ = """
scVelo wrapper for velocity analysis on kallisto output
usage example:
    scVelo_wrapper -i input-dir -o output-dir mm10
"""

# for wrapper
import argparse
import os
import sys
import textwrap

# for scVelo
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as scp
import sklearn
import sys
import loompy
import scipy.optimize
import scvelo as scv
import glob
import pickle
import anndata

# for plots
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

from collections import Counter
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
scv.settings.set_figure_params('scvelo')

def get_matrix(sampleName, spliced_dir, unspliced_dir, prefix):
    spliced_dir = "velocity_quant/" + sampleName + "/" + spliced_dir + "/" #spliced_counts
    unspliced_dir = "velocity_quant/" + sampleName + "/" + unspliced_dir + "/" #unspliced_counts

    s = sc.read_mtx(spliced_dir + prefix + "_isect.mtx")
    u = sc.read_mtx(unspliced_dir + prefix + "_isect.mtx")

    s_bcs = pd.read_csv(spliced_dir + prefix + ".barcodes.txt", header=None)
    u_bcs = pd.read_csv(unspliced_dir + prefix + ".barcodes.txt", header=None)

    s_genes = pd.read_csv(spliced_dir + prefix + ".genes_isect.txt", header=None)
    u_genes = pd.read_csv(unspliced_dir + prefix + ".genes_isect.txt", header=None)

    s.obs.index = s_bcs[0].values
    u.obs.index = u_bcs[0].values

    s.var.index = s_genes[0].values
    u.var.index = u_genes[0].values


    s_bcs["sample"] = sampleName
    u_bcs["sample"] = sampleName

    s_bcs.columns = ["bcs", "sample"]
    u_bcs.columns = ["bcs", "sample"]

    s_bcs.index = s_bcs["bcs"] + "." + s_bcs["sample"]
    u_bcs.index = u_bcs["bcs"] + "." + u_bcs["sample"]

    out = {'s': s,
           'u': u,
           's_bcs': s_bcs,
           'u_bcs': u_bcs,
           'genes': s_genes}
    return(out)

## to plot a given gene list (not exposed to the wrapper yet)
def plotGenesOnVelocity(geneList, outPdf):
    pdf = matplotlib.backends.backend_pdf.PdfPages(outPdf)
    for gene in geneList:
        try:
            fig = scv.pl.velocity_embedding_grid(adata, basis='umap',
                                                 color= gene, color_map = 'Reds',
                                                dpi = 300, use_raw = True, alpha=0.3, show=False)
            pdf.savefig(fig.get_figure())
        except ValueError:
            print(gene+" not in velocity genes.")
            continue
    pdf.close()

    return(None)


#samples = ['BM_all_S3merged', 'BM_all_S4merged', 'BM_progenitors_S1merged', 'BM_progenitors_S2merged']

def parse_args(defaults=None):
    """
    Parse arguments from the command line.
    """
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__description__),
        add_help=False
    )

    # Workflow options
    parser.add_argument("-s", "--sampleNames",
                          dest="samples",
                          metavar="STR",
                          help="sample names to combine for the analysis. The spliced matrices \
                                should exist as velocity_quant/<sampleName>/spliced_counts/spliced.mtx, \
                                (and vice-versa for unspliced counts).",
                          type=str,
                          nargs='+',
                          required=True,
                          default=None)

    parser.add_argument("-t", "--t2geneFile",
                          dest="tr2gene_file",
                          metavar="STR",
                          help="transcript to gene conversion file for kallisto",
                          type=str,
                          required=True,
                          default=None)

    parser.add_argument("-o", "--outdir",
                          dest="outdir",
                          required=True,
                          help="output directory for plots and anndata object",
                          default=None)

    parser.add_argument("-g", "--minGenes",
                         dest="min_genes",
                         metavar="INT",
                         help="minimum number of genes required to keep the cell (default: '%(default)s')",
                         type=int,
                         default=100)

    parser.add_argument("-c", "--minCells",
                         dest="min_cells",
                         metavar="INT",
                         help="minimum number of cells required to keep the gene (default: '%(default)s')",
                         type=int,
                         default=10)

    parser.add_argument("-sb", "--samplesAsBatches",
                         dest="samples_as_batches",
                         action="store_true",
                         default=False,
                         help="treat samples as batches and correct for them using combat (default: '%(default)s')")

    parser.add_argument("-cs", "--cellCycleGenes",
                          dest="cell_cycle_genes_tsv",
                          metavar="STR",
                          help="TSV file containing cell cycle genes. If provided, the wrapper will make a PCA plot based on those genes",
                          type=str,
                          required=False,
                          default=None)

    parser.add_argument("-h", "--help",
                         action="help",
                         help="show this help message and exit")


    return parser


def main():

    ## parse args
    parser = parse_args()
    args = parser.parse_args()
    outdir = args.outdir
    scv.settings.figdir = outdir
    sc.settings.figdir = outdir
    scv.settings.autoshow = False
    sc.settings.autoshow = False

    ## prepare matrices
    print("Preparing matrices")
    quant = [get_matrix(x, 'geneCounts_spliced', 'geneCounts_unspliced', 'output') for x in args.samples]
    # merge samples based on geneIDs
    s = quant[0]['s'].concatenate([sample['s'] for sample in quant[1:]], join = 'inner',
                              batch_key= 'sample', batch_categories=args.samples)
    u = quant[0]['u'].concatenate([sample['u'] for sample in quant[1:]], join = 'inner',
                              batch_key= 'sample', batch_categories=args.samples)
    #s_bcs = pd.concat([sample['s_bcs'] for sample in quant],ignore_index=False)
    #u_bcs = pd.concat([sample['u_bcs'] for sample in quant],ignore_index=False)

    ## get common indicies bw spliced and unspliced
    #intersect = u_bcs.index.intersection(s_bcs.index)
    #u_bcs = u_bcs.loc[intersect]
    #s_bcs = s_bcs.loc[intersect]
    #u = u[u_bcs.index.get_indexer_for(intersect)]
    #s = s[s_bcs.index.get_indexer_for(intersect)]

    ## make anndata
    print("Making anndata object")
    #genes = quant[0]['genes']
    #genes.columns=["gid"]
    #sadata = anndata.AnnData(X=s, obs=s_bcs, var=genes)
    #uadata = anndata.AnnData(X=u, obs=u_bcs, var=genes)

    adata = s.copy()
    adata.layers["spliced"] = s.X
    adata.layers["unspliced"] = u.X
    adata.layers["ambiguous"] = scp.sparse.csr_matrix(np.zeros(adata.X.shape))
    adata.obs = s.obs

    ## add gene names (from t2g)
    #t2g = pd.read_csv(args.tr2gene_file, header=None, sep="\t", names=["tid", "gid", "gene"], index_col=False)
    #t2g = t2g.drop_duplicates(["gid", "gene"])
    #t2g = t2g.set_index("gid")
    #adata.var["Gene"] = adata.var.index.map(t2g["gene"])
    #adata.var["Transcript"] = adata.var.index.map(t2g["tid"])

    ## show spliced/unspliced props
    scv.pp.show_proportions(adata)

    print("Writing QC metrics\n")
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.obs.to_csv(outdir+"/qc-metrics.csv")
    ## write out the unfiltered anndata
    adata.write_loom(outdir+"/anndata.loom")

    ## first, scanpy filter
    print("Filtering cells with min expressed genes >="+str(args.min_genes))
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    print("Cells left: {}".format(len(adata.obs.index)))

    print("Filtering genes with min expressed cells >="+str(args.min_cells))
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print("Genes left: {}".format(len(adata.var_names)))

    # batch correct?
    if args.samples_as_batches:
        sc.pp.combat(adata, key='sample')

    if args.cell_cycle_genes_tsv:
        # need to make it more generic
        gene_column = "mouse_gene_name"
        # check for cell cycle effects if cell cycle file given
        cell_cycle_genes = pd.read_csv(args.cell_cycle_genes_tsv, sep="\t")
        cell_cycle_genes = cell_cycle_genes[cell_cycle_genes[gene_column].isin (adata.var_names)]
        s_genes = cell_cycle_genes[cell_cycle_genes["stage"] == "S_phase"][gene_column]
        g2m_genes = cell_cycle_genes[cell_cycle_genes["stage"] == "G2M_phase"][gene_column]
        adata_copy = adata.copy()
        sc.pp.normalize_per_cell(adata_copy, counts_per_cell_after=1e4)
        sc.pp.log1p(adata_copy)
        sc.pp.scale(adata_copy)
        sc.tl.score_genes_cell_cycle(adata_copy, s_genes=s_genes, g2m_genes=g2m_genes)
        adata_cc_genes = adata_copy[:, cell_cycle_genes['mouse_gene_name']]
        sc.tl.pca(adata_cc_genes)
        sc.pl.pca_scatter(adata_cc_genes, color='phase', save = "/pca_cell_cycle_genes.png")

    ## then, velocyto filter
    print("Velocity filtering")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=5000)
    print("Cells left: {}".format(len(adata.obs.index)))
    print("Genes left: {}".format(len(adata.var_names)))
    ## plot
    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.louvain(adata)
    sc.tl.umap(adata, min_dist=0.1, spread=4)
    fig = sc.pl.umap(adata, color=['louvain', 'sample'], return_fig=True)
    fig.savefig(outdir+"/allsamples_UMAP.png")


    ## Velocity
    print("Computing velocity")
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata)
    # kept velocity genes
    nveloGenes = adata.var['velocity_genes'].sum()
    print("Number of velocity genes kept after default filtering: {}".format(nveloGenes))
    if nveloGenes < 10:
        print("Reducing the default threshold to get more velocity genes. New threshold: -0.05")
        scv.tl.velocity_genes(adata, min_r2=-0.05)
        if nveloGenes < 10:
            print("Still less than 10 velocity genes left. Output might not be meaningful")

    #adata.layers['velocity']
    # velocity Rsq histogram
    hs = scv.pl.hist(adata.var['velocity_r2'], show = False)
    histfig = hs.get_figure()
    histfig.savefig(outdir+"/velocity_Rsq-histogram.png")

    scv.tl.velocity_graph(adata)
    ## plots
    print("plotting and saving")
    fig1 = scv.pl.velocity_embedding_stream(adata, basis='umap', dpi = 300, use_raw = True, show = False)
    fig2 = scv.pl.velocity_embedding_grid(adata, basis='umap', color= 'louvain',dpi = 300, use_raw = True, show = False)
    fig3 = scv.pl.velocity_embedding_grid(adata, basis='umap', color= 'sample',dpi = 300, use_raw = True, show = False)

    fig1.get_figure().savefig(outdir+'/velocity-stream_louvain.png')
    fig2.get_figure().savefig(outdir+'/velocity-grid_louvain.png')
    fig3.get_figure().savefig(outdir+'/velocity-grid_samples.png')

    ## write the filtered anndata
    adata.write_loom(outdir+"/anndata_filtered.loom", write_obsm_varm=True)

if __name__ == "__main__":
    main()
