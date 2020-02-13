#!/usr/bin/env python3

__description__ = """
Wrapper for merging gene.mtx files.
usage example:
    gene_merged_wrapper -s samples -o output_path
"""

# import for wrapper
import argparse
import os
import sys
import textwrap

# import for merging
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io

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
    parser.add_argument("-s", "--samples",
                          dest="samples",
                          required=True,
                          help="samples that need to be merged",
                          default=None,
                          type=str)

    parser.add_argument("-o", "--out_dir",
                          dest="out_dir",
                          required=True,
                          help="output directory for gene_merged.mtx, barcodes_gene_merged.txt and genes_gene_merged.txt",
                          default=None,
                          type=str)

    return parser


def merge_files(samples, out_dir):
    # open each file as anndata
    adata_list = []

    # read mtx file
    for i in range(len(samples)):
        adata_list += [sc.read_mtx('transcripts_quant/'+samples[i]+'/gene_counts/output.mtx')]

        # read barcode file and set as obs
        barcodes = pd.read_csv('transcripts_quant/'+samples[i]+'/gene_counts/output.barcodes.txt', header=None)[0].values
        barcodes = samples[i] + '_' + barcodes # could also do this at the concatination step
        genes = pd.read_csv('transcripts_quant/'+samples[i]+'/gene_counts/output.genes.txt', header=None)[0].values

        # store obs and var
        adata_list[i].obs.index = barcodes
        adata_list[i].var.index = genes

    # concatenate all anndata objects
    if len(samples) > 0:
        con_adata = adata_list[0].concatenate(adata_list[1:len(samples)], join='outer', index_unique=None)
    else:
        con_adata = adata_list[0]

    # store result
    scipy.io.mmwrite(out_dir +'gene_merged', con_adata.X.astype(int))
    con_adata.obs.drop('batch', axis=1).to_csv(out_dir + 'barcodes_gene_merged.txt', header=None)
    con_adata.var.to_csv(out_dir + 'genes_gene_merged.txt', header=None)

def main():
    print("parsing arguments")
    # parse arguments
    parser = parse_args()
    args = parser.parse_args()
    samples = args.samples.split(',')

    print("start merging")
    # actually do the merging
    merge_files(samples, args.out_dir)

if __name__ == "__main__":
    main()
