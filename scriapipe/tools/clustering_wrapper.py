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
def get_matrix(countFile, outdir, barcodes, varlist, groups, annotation, type):
    # load adata
    adata = sc.read_mtx(countFile)

    # set obs.index and add other obs
    obs = pd.read_csv(barcodes, header=None)
    adata.obs.index = obs[0].values
    obs_new = obs[0].str.split("_", n=2, expand=True)
    n_col = obs_new.shape[1]
    if n_col > 1:
        adata.obs['sample'] = obs_new.iloc[:, 0:n_col-1].apply(lambda x: '_'.join(x),
                                                                axis=1).astype('category').values
        adata.obs['barcode'] = obs_new[2].astype('category').values
    else:
        adata.obs['sample'] = 'one_sample'

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
    if annotation is not None and type == 'genes':
        an = pd.read_csv(annotation, sep='\t', header=None, index_col=0)
        adata.var = pd.merge(adata.var, an, left_index=True, right_index=True,
                            how='left')
    elif annotation is not None and type == 'ECs':
        an = pd.read_csv(annotation, sep='\t', header=None, index_col=0)
        adata.var = pd.merge(adata.var, an, left_on='geneID', right_index=True,
                            how='left')

    # add predivined groups to obs, assumes categorical data and same index and header
    group_keys = None
    if groups is not None:
        obs_groups = pd.read_csv(groups, sep="\t")
        obs_groups.set_index(obs_groups.keys()[0], inplace=True)
        group_keys = obs_groups.keys()
        obs_groups = obs_groups[group_keys].astype('category') # assumed category based but maybe make this optional for continues scales
        adata.obs = adata.obs.merge(obs_groups, how='left', left_index=True,
                                    right_index=True)

    # color by gene
    print("\nData before filtering:")
    print(adata)
    return(adata, group_keys)

# preprocess tcc/gene matrix
def preprocess(adata, outdir,  highlyvar):
    # open pdf file for plots
    pdf = matplotlib.backends.backend_pdf.PdfPages(outdir+"/preprocess_stats.pdf")
    plt.figure(figsize=(10,5))

    # variables might be nice to expose these variables
    min_counts_cell = 1000
    min_genes_cell = 200
    min_cells_gene = 3

    print('\nMinumum counts per cell: {}'.format(min_counts_cell))
    print('Minimum genes per cell: {}'.format(min_genes_cell))
    print('Minimum cells expressing a gene: {}'.format(min_cells_gene))

    samples = set(adata.obs['sample'].tolist())

    # make figure for basic filtering
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.obs['n_genes'] = (adata.X>0).sum(axis=1).A1

    # plot total counts per cell histogram
    fig, axs = plt.subplots(1,3, figsize=(20,5))
    colors = ['b', 'g', 'm', 'c', 'y']
    max_counts = max(adata.obs['total_counts'].tolist())
    xmin = 8e-1
    xmax = max_counts
    y0min = 8e-1
    y0max = max(adata.obs['n_genes'].tolist())
    y1min = 8e-1
    y1max = max_counts

    for i, sample in enumerate(samples):
        # setup in loop variables
        bc = adata.obs.index[adata.obs['sample'] == sample].tolist()
        color = colors[i % 5]

        # plot total counts against number of genes per barcode
        x = adata.obs.loc[bc, 'total_counts'].tolist()
        y = adata.obs.loc[bc, 'n_genes'].tolist()
        axs[0].scatter(x, y, c=color, s=2)

        # plot total counts per barcode
        x = range(len(bc))
        y = adata.obs.loc[bc, 'total_counts'].sort_values(ascending=False).tolist()
        axs[1].plot(x, y, c=color)

        # plot total counts as histogram
        data = adata.obs.loc[bc, 'total_counts'].tolist()
        axs[2].hist(data, bins=200, color=color)


    axs[0].plot([min_counts_cell, min_counts_cell], [y0min, y0max], linewidth=1,
                linestyle='--', color='r')
    axs[0].plot([xmin, xmax], [min_genes_cell, min_genes_cell], linewidth=1,
                linestyle='--', color='r')
    axs[0].set_xlabel('total counts'); axs[0].set_ylabel('number of genes')
    axs[0].set_ylim([y0min, y0max]); axs[0].set_xlim([xmin, xmax])
    axs[0].set_yscale('log'); axs[0].set_xscale('log')
    axs[0].set_title('Total counts and number of genes correlation')

    axs[1].plot([xmin, xmax], [min_counts_cell, min_counts_cell], linewidth=1,
                linestyle='--', color='r')
    axs[1].set_xlabel('barcodes'); axs[1].set_ylabel('total counts')
    axs[1].set_ylim([y1min, y1max]); axs[1].set_xlim([xmin, xmax])
    axs[1].set_yscale('log'); axs[1].set_xscale('log')
    axs[1].set_title('Total counts sorted')

    # axs[2].set_ylim([1e0, 5e3]); axs[2].set_xlim([1e0, 5e4])
    y2min, y2max = axs[2].get_ylim()
    axs[2].plot([min_counts_cell, min_counts_cell], [y2min, y2max], linewidth=1,
                linestyle='--', color='r')
    axs[2].set_xlabel('Total counts'); axs[2].set_ylabel('Number of barcodes')
    axs[2].set_yscale('log'); axs[2].set_xscale('linear')
    axs[2].set_title('Total counts sorted')

    # save the figure to preprocess_stats.pdf
    pdf.savefig(fig)

    #filter minimum number of genes per cell
    print('\nAdata dimensions before filtering cells:\n{}'.format(adata.shape))
    if adata.shape[1] <= min_genes_cell:
        print('WARNING: few genes only {}'.format(adata.shape[1]))
        sc.pp.filter_cells(adata, min_genes = math.floor(adata.shape[1]/2))
    else:
        sc.pp.filter_cells(adata, min_genes= min_genes_cell)
    # filter minimum number of counts per cell, check if amount of cells with at least
    if math.ceil(np.quantile(adata.obs['total_counts'].values, .90)) <= min_counts_cell:
        # if 10% or less of cells are kept at threshold set threshold to keep 50%
        print('WARNING: few counts per cell keep top 50% with cut off at: {} counts'.format(math.ceil(np.quantile(qc[0]['total_counts'].values, .50))))
        sc.pp.filter_cells(adata, min_counts = math.ceil(np.quantile(adata.obs['total_counts'].values, .50)))
    else:
        sc.pp.filter_cells(adata, min_counts= min_counts_cell)
    print('\nAdata dimensions after filtering cells:\n{}'.format(adata.shape))

    for sample in samples:
        bc = adata.obs.index[adata.obs['sample'] == sample].tolist()
        print('\nMedian counts sample {} per cell: {}'.format(sample, adata.obs.loc[bc, 'total_counts'].sort_values(ascending=False).median()))

    # plot counts per gene
    fig, axs = plt.subplots(1,2, figsize=(20,8))

    # plot total counts against number of barcodes per gene
    x = adata.var['total_counts'].tolist()
    y = adata.var['n_cells_by_counts'].tolist()
    axs[0].scatter(x, y, c='g', s=2)
    xmin, xmax = axs[0].get_xlim()
    ymin, ymax = axs[0].get_ylim()
    # plot from below 1, 0 not possible due to log axis
    xmin = 8e-1; ymin = 8e-1
    axs[0].set_ylim([ymin, ymax]); axs[0].set_xlim([xmin, xmax])
    axs[0].plot([xmin, xmax], [min_cells_gene, min_cells_gene], linewidth=1,
                linestyle='--', color='r')
    axs[0].set_xlabel('total counts'); axs[0].set_ylabel('number of cells')
    axs[0].set_yscale('log'); axs[0].set_xscale('log')
    axs[0].set_title('Total counts and number of cells correlation')

    # plot counts per gene normalized to precence in number of cells
    adata.var['mean_count_per_expressing_cell'] = adata.var['total_counts']/adata.var['n_cells_by_counts']
    axs[1].hist(adata.var['mean_count_per_expressing_cell'].sort_values(ascending=False),
                bins=50, color='g')
    axs[1].set_xlabel('Mean counts per expressing cell'); axs[1].set_ylabel('Number of genes')
    axs[1].set_yscale('log'); axs[1].set_xscale('linear')
    axs[1].set_title('Mean counts per gene per cell')

    # save the figure to preprocess_stats.pdf
    pdf.savefig(fig)

    # filter min of cells per gene
    sc.pp.filter_genes(adata, min_cells = min_cells_gene)
    print('\nAdata dimensions after filtering genes:\n{}'.format(adata.shape))

    # normalize and log transform (should we use normalize per cell?)
    ntotal = math.ceil(np.quantile(adata.obs['total_counts'].values, .90)/100)*100 #round 90th quantile up to the nearest 100
    sc.pp.normalize_total(adata, target_sum = ntotal,
                          exclude_highly_expressed = True,
                          inplace = True)
    print('\nNormalize to value: {}'.format(ntotal))
    sc.pp.log1p(adata)

    if highlyvar > 0 and highlyvar < adata.shape[1]:
        print('\nSelection of highly variable genes')
        # set raw data before highly variable genes selection
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, n_top_genes=highlyvar, inplace=True,
                                    subset=True)
        sc.pp.regress_out(adata, ['n_counts'])
        sc.pp.scale(adata, max_value=10)
        print("Adata dimensions after selection highly variable genes/TCCs: {}".format(adata.shape))
    else:
        print('No selection for highly variable genes')

    # close preprocess_stats.pdf
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

    adata.obs.to_csv(outdir+"/barcode_cluster.tsv", sep="\t")

    # store genes
    adata.var.to_csv(outdir+"/var_cluster.tsv", sep="\t")

    # store adata object
    scipy.io.mmwrite(outdir+"/cluster", adata.X)

    ## write the filtered anndata
    adata.write_loom(outdir+"/anndata.loom", write_obsm_varm=True)

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

    parser.add_argument("-a", "--annotation",
                          dest="annotation",
                          required=False,
                          help="tsv file containing additional information to load into AnnData var e.g. gene name, chromosome or biotype",
                          default=None)

    parser.add_argument("-hv", "--highly-variable",
                         dest="highlyvar",
                         required=False,
                         help="to turn on filtering for highly variable genes/TCCs",
                         type=int,
                         default=0)

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
    if args.annotation == 'None':
        args.annotation = None

    print("parameters that were used:")
    print("clustering_wrapper: \n--outdir {} \n--sample {} \n--barcodes {} \
            \n--varlist {} \n--type {} \n--highly-variable {}".format(args.outdir,
            args.sample, args.barcodes, args.varlist, args.type, args.highlyvar))

    # get matrix
    print("\nget matrix")
    adata, group_keys = get_matrix(args.sample, outdir, args.barcodes, args.varlist,
                                   args.col_groups, args.annotation, args.type)

    print("\nstart preprocessing")
    adata = preprocess(adata, outdir, args.highlyvar)
    #print(adata.X)
    print("\nstart clustering")
    barcode_cluster = clustercells(adata, outdir, group_keys)

if __name__ == "__main__":
    main()
