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
    # read mtx file
    print("reading first mtx file")
    mtx = scipy.io.mmread(out_dir+samples[0]+'/gene_counts/output.mtx')
    X = mtx.todense()

    # read barcode file and set as obs
    print("reading barcode file")
    barcodes = pd.read_csv(out_dir+samples[0]+'/gene_counts/output.barcodes.txt', header=None, sep=' ')[0].values
    barcodes = samples[0] + '_' + barcodes
    # read gene file and set as var
    print("reading gene file")
    genes = pd.read_csv(out_dir+samples[0]+'/gene_counts/output.genes.txt', header=None, sep=' ')[0].values

    # make df
    print("making first df")
    df_old = pd.DataFrame(X, columns=genes, index=barcodes).astype(int)

    for sample in samples[1::]:
        print("sample: {}".format(sample))
        # read mtx file
        mtx = scipy.io.mmread(out_dir+sample+'/gene_counts/output.mtx')
        X = mtx.todense()

        # read barcode file and set as obs
        barcodes = pd.read_csv(out_dir+sample+'/gene_counts/output.barcodes.txt', header=None, sep=' ')[0].values
        barcodes = sample + '_' + barcodes
        # read gene file and set as var
        genes = pd.read_csv(out_dir+sample+'/gene_counts/output.genes.txt', header=None, sep=' ')[0].values

        # make df
        df_new = pd.DataFrame(X, columns=genes, index=barcodes).astype(int)

        # concatenate new to total
        df_old = df_old.append(df_new)

    # store final result
    scipy.io.mmwrite(out_dir+"gene_merged", scipy.sparse.csr_matrix(df_old))

    # store barcodes in format: sample_barcode
    barcode_file = open(out_dir+'barcodes_gene_merged.txt', 'w')
    barcode_file.write('\n'.join(df_old.index.values))
    barcode_file.close()

    # store genes
    genes_file = open(out_dir+'genes_gene_merged.txt', 'w')
    genes_file.write('\n'.join(df_old.columns.values))
    genes_file.close()

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
