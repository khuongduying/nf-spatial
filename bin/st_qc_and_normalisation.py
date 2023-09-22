#!/usr/bin/env python

import argparse
import scanpy as sc
import pandas as pd
import re
import warnings

def main(args):
    warnings.filterwarnings("ignore")
    sc.settings.verbosity = 0

    # Data Loading and Initial QC Metrics
    st_adata = sc.read(args.raw_adata)

    # Mark mitochondrial genes
    st_adata.var['mt'] = st_adata.var_names.str.startswith('MT-')

    # Mark hemoglobin genes based on their naming patterns
    hemoglobin_pattern = re.compile(r'^HBA|HBB')
    st_adata.var['hb'] = st_adata.var_names.str.match(hemoglobin_pattern)

    # Calculate QC metrics with both mt and hb genes as QC variables
    sc.pp.calculate_qc_metrics(st_adata, qc_vars=["mt", "hb"], inplace=True)

    # Filter spots based on tissue presence and QC parameters
    st_adata = st_adata[st_adata.obs["in_tissue"] == 1]
    sc.pp.filter_cells(st_adata, min_counts=args.min_counts)
    sc.pp.filter_cells(st_adata, min_genes=args.min_genes)

    # Normalisation and variable genes detection
    st_adata.write(args.data_plain_name)
    sc.pp.normalize_total(st_adata, inplace=True)
    sc.pp.log1p(st_adata)
    sc.pp.highly_variable_genes(st_adata, flavor="seurat", n_top_genes=2000)
    st_adata.write(args.data_norm_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform QC and normalization on spatial transcriptomics data.')
    parser.add_argument('--raw_adata', required=True, help='Path to the raw data')
    parser.add_argument('--min_counts', type=int, default=500, help='Minimum counts')
    parser.add_argument('--min_genes', type=int, default=250, help='Minimum genes')
    parser.add_argument('--min_cells', type=int, default=1, help='Minimum cells')
    parser.add_argument('--data_plain_name', default="st_adata_plain.h5ad", help='Name of the raw data save file')
    parser.add_argument('--data_norm_name', default="st_adata_norm.h5ad", help='Name of the normalized data save file')

    args = parser.parse_args()
    main(args)