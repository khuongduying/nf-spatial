#!/usr/bin/env python

import argparse
import scanpy as sc
import pandas as pd
import SpatialDE

def main(args):
    # Read spatial transcriptomics data
    st_adata = sc.read(args.file_name_st)

    # Spatially variable genes
    results = SpatialDE.run(st_adata.obsm["spatial"], st_adata.to_df())
    results.set_index("g", inplace=True)
    results = results.loc[~results.index.duplicated(keep="first")]
    st_adata.var = pd.concat([st_adata.var, results.loc[st_adata.var.index.values, :]], axis=1)

    sorted_results = st_adata.var.sort_values("qval", ascending=True)
    sorted_results.to_csv(args.save_spatial_de_file_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform spatial differential gene expression.')
    parser.add_argument('--file_name_st', required=True, help='Path to the ST data file')
    parser.add_argument('--save_spatial_de_file_name', required=True, help='Save path for spatial DE results')
    
    args = parser.parse_args()
    main(args)
