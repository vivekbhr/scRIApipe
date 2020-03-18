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

# for clustering and pp
import pandas as pd
import scanpy as sc
import numpy as np
import scipy.io
import loompy
import math

# for plots
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

# load matrix as AnnData object
def get_matrix(countFile, outdir, barcodes, varlist, groups, extendedVar, type):
    # load adata
    adata = sc.read_mtx(countFile)

    # set obs.index and add other obs
    obs = pd.read_csv(barcodes, header=None)
    adata.obs.index = obs[0].values
    obs_new = obs[0].str.split("_", n=2, expand=True)
    n_col = obs_new.shape[1]
    if n_col > 1:
        adata.obs['sample'] = obs_new.iloc[:, 0:n_col-1].apply(lambda x: '_'.join(x), axis=1).astype('category').values
        adata.obs['barcode'] = obs_new[2].astype('category').values

    # set index
    if type == 'ECs':
        var = pd.read_csv(varlist, sep='\t', header=None, index_col=2)
        var.index = var.index.astype('str').values
        n_col = var.shape[1]
        if n_col > 2:
            var.drop([1]+list(range(3,n_col+1)), axis=1, inplace=True)
        else:
            var.drop([1], axis=1, inplace=True)
        var.columns = ['geneID']
        adata.var = var
    elif type == 'genes':
        var = pd.read_csv(varlist, sep='\t', header=None, index_col=0)
        adata.var.index = var.index.values

    # extend variable df assumes same index
    if extendedVar is not None and type == 'genes':
        exVar = pd.read_csv(extendedVar, sep='\t', header=None, index_col=0)
        adata.var = pd.merge(adata.var, exVar, left_index=True, right_index=True, how='left')
    elif extendedVar is not None and type == 'ECs':
        exVar = pd.read_csv(extendedVar, sep='\t', header=None, index_col=0)
        adata.var = pd.merge(adata.var, exVar, left_on='geneID', right_index=True, how='left')

    # add predivined group
    group_keys = None
    if groups is not None:
        obs_groups = pd.read_csv(groups, sep="\t")
        obs_groups.set_index(obs_groups.keys()[0], inplace=True)
        group_keys = obs_groups.keys()
        obs_groups = obs_groups[group_keys].astype('category') # assumed category based but maybe make this optional for continues scales
        adata.obs = adata.obs.merge(obs_groups, how='left', left_index=True, right_index=True)

    # color by gene
    print("\nData before filtering:")
    print(adata)
    return(adata, group_keys)

# preprocess tcc/gene matrix
def preprocess(adata, outdir,  highlyvar):
    # open pdf file for plots
    pdf = matplotlib.backends.backend_pdf.PdfPages(outdir+"/preprocessed.pdf")
    plt.figure(figsize=(10,5))

    # make figure for basic filtering
    qc = sc.pp.calculate_qc_metrics(adata)

    # plot total counts per cell histogram
    fig_counts_cell = plt.figure(figsize=(10,5))
    fig_counts_cell = qc[0]['total_counts'].plot(kind='hist', logy=True, logx=True, bins=200, title='total counts per cell').get_figure()
    print('\nMedian counts per cell: {}'.format(qc[0]['total_counts'].sort_values(ascending=False).median()))
    print('\nAdata before filtering cells:\n{}'.format(adata))
    sc.pp.filter_cells(adata, min_genes= 100)
    sc.pp.filter_cells(adata, min_counts= 1000)
    print('\nAdata after filtering cells:\n{}'.format(adata))
    pdf.savefig(fig_counts_cell)

    # plot counts per gene
    fig_counts_gene = plt.figure(figsize=(10,5))
    fig_counts_gene = qc[1]['total_counts'].sort_values(ascending=False).plot(kind='line', logy=True, title='total counts per gene').get_figure()
    print('\nMedian counts per gene: {}'.format(qc[1]['total_counts'].sort_values(ascending=False).median()))
    print('\nAdata before filtering genes:\n{}'.format(adata))
    sc.pp.filter_genes(adata, min_cells = 3)
    # sc.pp.filter_genes(adata, min_counts = 6)
    print('\nAdata after filtering genes:\n{}'.format(adata))
    pdf.savefig(fig_counts_gene)

    # plot to see count distribution and gene distrubution
    fig1 = sc.pl.violin(adata, 'n_counts', jitter=0.4, return_fig=True)
    fig2 = sc.pl.violin(adata, 'n_genes', jitter=0.4, return_fig=True)

    # normalize and log transform (should we use normalize per cell?)
    qc = sc.pp.calculate_qc_metrics(adata)
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    ntotal = math.ceil(np.quantile(qc[0]['total_counts'].values, .10)/100)*100 #round 90th quantile up to the nearest 100
    sc.pp.normalize_total(adata, target_sum = ntotal,
                          exclude_highly_expressed = True,
                          inplace = True)
    sc.pp.log1p(adata)

    if highlyvar == True:
        # set raw data before highly variable genes selection
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, n_top_genes=5000, inplace=True, subset=True)
        sc.pp.regress_out(adata, ['n_counts'])
        sc.pp.scale(adata, max_value=10)
        print("\nAfter selection highly variable genes/TCCs:")
        print(adata)

    # save figures close pdf file
    pdf.savefig(fig1)
    pdf.savefig(fig2)

    pdf.close()

    return(adata)

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
    sc.tl.umap(adata, min_dist=0.1, spread=5)
    fig3 = sc.pl.umap(adata, color="louvain", return_fig=True, palette='tab20')

    # plot sample on umap
    fig4 = sc.pl.umap(adata, color='sample', return_fig=True, palette='tab20')

    # plot also on groups if indicated
    if group_keys is not None:
        for key in group_keys:
            fig4 = sc.pl.pca(adata, color=key, return_fig=True)
            sc.tl.umap(adata, min_dist=0.1, spread=5)
            fig5 = sc.pl.umap(adata, color=key, return_fig=True, palette='tab20')

    # save barcode cluster tsv
    # if group_keys is None:
    #     barcode_cluster = adata.obs.loc[:,'louvain']
    # else:
    #     barcode_cluster = adata.obs.loc[:,['louvain'] + group_keys.to_list()]
    # barcode_cluster.to_csv(outdir+"/barcode_cluster.tsv", sep="\t")
    adata.obs.to_csv(outdir+"/barcode_cluster.tsv", sep="\t")

    # store genes
    adata.var.to_csv(outdir+"/var_cluster.tsv", sep="\t")

    # store adata object
    scipy.io.mmwrite(outdir+"/cluster", adata.X)
    #adata.to_df().to_csv(outdir+"/cluster.tsv", sep='\t', header=True)

    ## write the filtered anndata
    adata.write_loom(outdir+"/anndata_filtered.loom", write_obsm_varm=True)

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
                          help="provide file with ordered variables",
                          default=None,
                          type=str)

    parser.add_argument("-t", "--type",
                          dest="type",
                          required=True,
                          help="provide type of data ECs or genes",
                          type=str)

    parser.add_argument("-g", "--groups",
                          dest="col_groups",
                          required=False,
                          help="file indicating cell types matching tcc.mtx",
                          default=None)

    parser.add_argument("-ev", "--extendedVar",
                          dest="extendedVar",
                          required=False,
                          help="tsv file containing additional information to load into AnnData var e.g. gene name, chromosome or biotype",
                          default=None)

    parser.add_argument("-hv", "--highly-variable",
                         dest="highlyvar",
                         action = 'store_true',
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

    # correct for passing 'None' by
    if args.col_groups == 'None':
        args.col_groups = None
    if args.extendedVar == 'None':
        args.extendedVar = None

    print("parameters that were used:")
    print("clustering_wrapper: \n--outdir {} \n--sample {} \n--barcodes {} \n--varlist {} \
            \n--type {} \n--highly-variable {}".format(args.outdir,
            args.sample, args.barcodes, args.varlist, args.type, args.highlyvar))

    # get matrix
    print("\nget matrix")
    adata, group_keys = get_matrix(args.sample, outdir, args.barcodes, args.varlist,
                                   args.col_groups, args.extendedVar, args.type)

    print("\nstart preprocessing")
    adata = preprocess(adata, outdir, args.highlyvar)
    #print(adata.X)
    print("\nstart clustering")
    barcode_cluster = clustercells(adata, outdir, group_keys)

if __name__ == "__main__":
    main()
