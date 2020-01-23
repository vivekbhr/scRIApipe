#!/usr/bin/env python3

__description__ = """
Wrapper for clustering on equivalance classes, transcripts or genes from
kallisto output.
usage example:
    ecClustering_wrapper -s tcc.mtx -b barcodes.csv -o output-dir
"""

# for wrapper
import argparse
import os
import sys
import textwrap

# for clustering
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as scp

# for plots
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

# load matrix as AnnData object
def get_matrix(countFile, outdir, barcodes, varlist, groups):
    # load adata
    adata = sc.read_mtx(countFile)

    # set obs.index and add other obs
    obs = pd.read_csv(barcodes, header=None)
    adata.obs.index = obs[0].values
    obs_new = obs[0].str.split("_", n=2, expand=True)
    n_col = obs_new.shape[1]
    adata.obs['sample'] = obs_new.iloc[:, 0:n_col-1].apply(lambda x: '_'.join(x), axis=1).astype('category').values
    adata.obs['barcode'] = obs_new[2].astype('category').values

    # set index
    adata.var.index = pd.read_csv(varlist, sep='\t', header=None)[2].astype('category').values

    # add predivined group
    group_keys = None
    if groups is not None:
        obs_groups = pd.read_csv(groups, sep="\t")
        obs_groups.set_index(obs_groups.keys()[0], inplace=True)
        group_keys = obs_groups.keys()
        obs_groups = obs_groups[group_keys].astype('category') # assumed category based but maybe make this optional for continues scales
        adata.obs = adata.obs.merge(obs_groups, how='left', left_index=True, right_index=True)

    return(adata, group_keys)

# preprocess tcc/gene matrix
# preprocess tcc/gene matrix
def preprocess(adata, outdir, mcells, mgenes, mcounts, md, ntotal, highlyvar):
    # open pdf file for plots
    pdf = matplotlib.backends.backend_pdf.PdfPages(outdir+"/preprocessed.pdf")
    plt.figure(figsize=(10,5))

    # basic filtering
    sc.pp.filter_cells(adata, min_genes= mgenes)
    sc.pp.filter_cells(adata, min_counts= mcounts)
    sc.pp.filter_genes(adata, min_cells = mcells)

    # plot to see count distribution and gene distrubution
    fig1 = sc.pl.violin(adata, 'n_counts', jitter=0.4, return_fig=True)
    fig2 = sc.pl.violin(adata, 'n_genes', jitter=0.4, return_fig=True)

    # normalize and log transform (should we use normalize per cell?)
    sc.pp.normalize_total(adata, target_sum = ntotal,
                          exclude_highly_expressed = True)
    sc.pp.log1p(adata)

    if highlyvar == True:
        # set raw data before highly variable genes selection
        adata.raw = adata

        # look at highly variable genes/ec
        sc.pp.highly_variable_genes(adata, min_mean=0.25, max_mean=6,
                                    min_disp=0.5, inplace=True)
        fig3 = sc.pl.highly_variable_genes(adata, return_fig=True)
        #actually filter for highly variable genes/ec
        sc.pp.highly_variable_genes(adata, min_mean=0.25, max_mean=6,
                                    min_disp=0.5, inplace=True, subset=True)

    # store adata object
    adata.to_df().to_csv(outdir+"/preprocessed.tsv", sep='\t', header=True)

    # save figures close pdf file
    pdf.savefig(fig1)
    pdf.savefig(fig2)
    if highlyvar == True:
        pdf.savefig(fig3)
    pdf.close()

    return(adata)

# cluster cells use louvain clustering
# cluster cells use louvain clustering
def clustercells(adata, outdir, group_keys):
    # open pdf file for plots
    pdf = matplotlib.backends.backend_pdf.PdfPages(outdir+"/clustering.pdf")
    plt.figure(figsize=(10,5))

    # do the clustering
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata)

    # plot louvain on pca
    fig1 = sc.pl.pca(adata, color="louvain", return_fig=True)
    fig2 = sc.pl.pca_variance_ratio(adata)

    # plot louvain on umap
    sc.tl.umap(adata, min_dist=0.1, spread=0.5)
    fig3 = sc.pl.umap(adata, color="louvain", return_fig=True)

    # plot sample on umap
    fig4 = sc.pl.umap(adata, color='sample', return_fig=True)

    # plot also on groups if indicated
    if group_keys is not None:
        for key in group_keys:
            fig4 = sc.pl.pca(adata, color=key, return_fig=True)
            sc.tl.umap(adata, min_dist=0.1, spread=0.5)
            fig5 = sc.pl.umap(adata, color=key, return_fig=True)

    # save barcode cluster tsv
    if group_keys is None:
        barcode_cluster = adata.obs.loc[:,'louvain']
    else:
        barcode_cluster = adata.obs.loc[:,['louvain'] + group_keys.to_list()]
    barcode_cluster.to_csv(outdir+"/barcode_cluster.tsv", sep="\t")

    # store adata object
    adata.to_df().to_csv(outdir+"/cluster.tsv", sep='\t', header=True)

    # save figures and close pdf
    pdf.savefig(fig1)
    pdf.savefig(fig2)
    pdf.savefig(fig3)
    if group_keys is not None:
        pdf.savefig(fig4)
        pdf.savefig(fig5)
    pdf.close()

    return(adata)


# parse arguments from commandline
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
    parser.add_argument("-o", "--outdir",
                          dest="outdir",
                          required=True,
                          help="output directory for plots and anndata object",
                          default=None,
                          type=str)

    parser.add_argument("-s", "--sample",
                          dest="sample",
                          required=True,
                          help="ec/gene count matrix path/<file>",
                          default=None,
                          type=str)

    parser.add_argument("-b", "--barcodes",
                          dest="barcodes",
                          required=True,
                          help="file with barcode list matching the tcc.mtx",
                          default=None,
                          type=str)

    parser.add_argument("-v", "--varlist",
                          dest="varlist",
                          required=True,
                          help="provide ec/gene name file",
                          default=None,
                          type=str)

    parser.add_argument("-g", "--groups",
                          dest="col_groups",
                          required=False,
                          help="file indicating cell types matching tcc.mtx",
                          default=None)

    parser.add_argument("--cells",
                         dest="cells",
                         required=False,
                         help="minimum number cells containing gene/TCC - used in filtering",
                         default=5,
                         type=int)

    parser.add_argument("--count",
                         dest="count",
                         required=False,
                         help="minimum number of counts per cell - used in filtering",
                         default=20,
                         type=int)

    parser.add_argument("--genes",
                         dest="genes",
                         required=False,
                         help="minimum number of genes/TCCs per cell - used in filtering",
                         default=100,
                         type=int)

    parser.add_argument("--dispersity",
                         dest="dispersity",
                         required=False,
                         help="minimum amount of dispersity - used in filtering",
                         default=0.5,
                         type=float)

    parser.add_argument("--normalize",
                        dest="normalize",
                        required=False,
                        help="normalize cell reads to indicated amount",
                        default=1e6,
                        type=float)

    parser.add_argument("-hv", "--highly-variable",
                         dest="highlyvar",
                         required=False,
                         help="to turn on filtering for highly variable genes/TCCs",
                         default=False)

    return parser

def main():

    # parse arguments
    parser = parse_args()
    args = parser.parse_args()
    outdir = args.outdir
    sc.settings.figdir = outdir
    sc.settings.autoshow = False

    # get matrix
    print("get matrix")
    adata, group_keys = get_matrix(args.sample, outdir, args.barcodes, args.varlist, args.col_groups)

    print("start preprocessing")
    adata = preprocess(adata, outdir, args.cells, args.genes,
                                      args.count, args.dispersity,
                                      args.normalize, args.highlyvar)
    print("start clustering")
    barcode_cluster = clustercells(adata, outdir, group_keys)

if __name__ == "__main__":
    main()


# def diff_expression_ec():
#
#     return()

# def diff_expression_gene(adata, outdir, groups, n_top_genes):
#     # open pdf file for plots
#     pdf = matplotlib.backends.backend_pdf.PdfPages(outdir+"/diff_expression_gene")
#
#     # use wilcoxon instead of t-test
#     sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
#     fig1 = sc.pl.rank_genes_groups(adata, n_genes=n_top_genes, sharey=False, return_fig=True)
#
#     # do logistic regression analyses mostly usefull for datasets with large cell populations
#     sc.tl.rank_genes_groups(adata, 'louvain', method='logreg')
#     fig2 = sc.pl.rank_genes_groups(adata, n_genes=n_top_genes, sharey=False, return_fig=True)
#
#     # save figures
#     pdf.savefig(fig1.get_figure())
#     pdf.savefig(fig2.get_figure())
#
#     if groups is not None:
#         # use wilcoxon instead of t-test
#         sc.tl.rank_genes_groups(adata, groups, method='wilcoxon')
#         fig3 = sc.pl.rank_genes_groups(adata, n_genes=n_top_genes, sharey=False, return_fig=True)
#
#         # do logistic regression analyses mostly usefull for datasets with large cell populations
#         sc.tl.rank_genes_groups(adata, groups, method='logreg')
#         fig4 = sc.pl.rank_genes_groups(adata, n_genes=n_top_genes, sharey=False, return_fig=True)
#
#         # save figures
#         pdf.savefig(fig3.get_figure())
#         pdf.savefig(fig4.get_figure())
#
#     # close pdf
#     pdf.close()
#
#     return()
