#!/usr/bin/env python

# Import necessary libraries
import argparse
import json
import pandas as pd
from anndata import AnnData
from matplotlib.image import imread
from pathlib import Path
from scanpy import read_10x_h5
from typing import Union, Optional


def read_visium_mtx(
    path: Union[str, Path],
    *,
    load_images: bool = True,
    library_id: Optional[str] = None
) -> AnnData:
    """
    Read data from Visium spaceranger output and return an AnnData object.
    Parameters
    ----------
    path : Union[str, Path]
        Path to the spaceranger output directory.
    load_images : bool, optional
        If True, images will be loaded into the AnnData object. Defaults to True.
    library_id : str, optional
        Identifier for the visium library. Can be used when concatenating multiple adata objects.
        Defaults to None.
    Returns
    -------
    AnnData
        Annotated Data object containing matrix, barcode, features, and optional image data.
    """

    path = Path(path)
    adata = read_10x_h5(path / "raw_feature_bc_matrix.h5")
    adata.var["gene_symbol"] = adata.var_names
    adata.var.set_index("gene_ids", inplace=True)

    adata.uns["spatial"] = dict()

    if library_id is None:
        library_id = "library_id"

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        # Define paths for the required files
        files = {
            'tissue_positions_file': path / "tissue_positions.csv",
            'scalefactors_json_file': path / "scalefactors_json.json",
            'hires_image': path / "tissue_hires_image.png",
            'lowres_image': path / "tissue_lowres_image.png"
        }

        # Check if files exist; continue if images are missing
        for f in files.values():
            if not f.exists():
                if "image" in str(f):
                    print(f"You seem to be missing an image file: {f}")
                else:
                    raise OSError(f"Could not find '{f}'")

        # Load scale factors for the spatial data
        with open(files["scalefactors_json_file"], 'r') as f:
            scalefactors = json.load(f)
        adata.uns["spatial"][library_id]["scalefactors"] = scalefactors

        # Load spatial data
        positions = pd.read_csv(files["tissue_positions_file"], index_col="barcode", dtype={"in_tissue": bool})
        adata.obs = adata.obs.join(positions, how="left")
        adata.obsm["spatial"] = adata.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
        adata.obs.drop(columns=["pxl_row_in_fullres", "pxl_col_in_fullres"], inplace=True)

        # Load images
        adata.uns["spatial"][library_id]["images"] = dict()
        for res in ["hires", "lowres"]:
            if files[f"{res}_image"].exists():
                adata.uns["spatial"][library_id]["images"][res] = imread(files[f"{res}_image"])

    return adata


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Load spatial transcriptomics data from MTX matrices and aligned images.")
    parser.add_argument("--SRCountDir", metavar="SRCountDir", type=str, required=True, help="Input directory with Spaceranger data.")
    parser.add_argument("--outAnnData", metavar="outAnnData", type=str, required=True, help="Output h5ad file path.")
    args = parser.parse_args()

    # Read Visium data
    st_adata = read_visium_mtx(args.SRCountDir, library_id=None, load_images=True)

    # Write raw AnnData to file
    st_adata.write(args.outAnnData)
