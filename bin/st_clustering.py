import argparse
import scanpy as sc
import numpy as np
import pandas as pd
from umap import UMAP

def main(args):
    sc.settings.verbosity = 0
    sc.set_figure_params(dpi_save=300, facecolor="white")

    # Read in the Spatial Transcriptomics data
    st_adata = sc.read(args.file_name_st)

    # PCA, Neighbors, UMAP, and Leiden clustering
    sc.pp.pca(st_adata)
    sc.pp.neighbors(st_adata)
    sc.tl.umap(st_adata)
    sc.tl.leiden(st_adata, key_added="clusters", resolution=args.resolution)

    # Fixing Scanpy log1p issue
    st_adata.uns['log1p']['base'] = None

    # Save processed data
    if args.save_file_st:
        st_adata.write(args.save_file_st)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Cluster and visualize spatial transcriptomics data.')
    parser.add_argument('--file_name_st', required=True, help='Path to the ST data file')
    parser.add_argument('--resolution', type=float, default=1.0, help='Resolution for Leiden clustering')
    parser.add_argument('--save_file_st', default=None, help='Save path for processed ST data')

    args = parser.parse_args()
    main(args)