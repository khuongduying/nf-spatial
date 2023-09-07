# Import necessary libraries
import argparse
import json
import pandas as pd
from anndata import AnnData
from matplotlib.image import imread
from pathlib import Path
import scanpy as sc
from typing import Union, Optional


def read_visium_mtx(
        path: Union[str, Path], 
        *, 
        load_images: bool = True, 
        library_id: Optional[str] = None
        ) -> AnnData:
    """
    Function to read data from Visium spaceranger output and return an AnnData object.

    Parameters:
    ----------
    path : Union[str, Path]
        Path to the spaceranger output directory.
    load_images : bool, optional
        If True, images will be loaded into the AnnData object. Defaults to True.
    library_id : str, optional
        Identifier for the visium library. Can be used when concatenating multiple adata objects.
        Defaults to None.

    Returns:
    -------
    AnnData
        Annotated Data object containing matrix, barcode, features, and optional image data.
    """

    # Convert path to Path object for easy manipulations
    path = Path(path)

    # Define paths to matrix, barcode, and features
    matrix_file = path / "filtered_feature_bc_matrix" / "matrix.mtx"
    barcodes_file = path / "filter_feature_bc_matrix" / "barcodes.tsv"
    genes_file = path / "filter_feature_bc_matrix" / "features.tsv"

    # Load matrix and transpose to get cells in rows
    adata = sc.read_mtx(matrix_file).T
    adata.obs_names = pd.read_csv(barcodes_file, header=None)[0].values

    # Load gene features
    features = pd.read_csv(genes_file, sep="\t", header=None)
    adata.var_names = features[1].values
    adata.var["gene_ids"] = features[0].values

    # Load scale factors for the spatial data
    scalefactors_file_path = path / 'spatial' / 'scalefactors_json.json'
    with open(scalefactors_file_path, 'r') as f:
        scalefactors = json.load(f)
    adata.uns['spatial'][spatial_key]['scalefactors'] = scalefactors

    # Read spatial data and match with adata barcodes
    positions = pd.read_csv(path / 'spatial' / 'tissue_positions_list.csv', index_col="barcode", dtype={"in_tissue": bool})
    matched_spatial_data = positions.loc[adata.obs_names]
    matched_spatial_data.dropna(inplace=True)  # Ensure no NaNs
    adata = adata[matched_spatial_data.index]
    adata.obs = adata.obs.join(matched_spatial_data, how="left")

    # Load images and decide which to use: high-resolution or low-resolution
    hires_image_path = path / 'spatial' / 'tissue_hires_image.png'
    lowres_image_path = path / 'spatial' / 'tissue_lowres_image.png'
    hires_image = None if not hires_image_path.exists() else plt.imread(hires_image_path)
    lowres_image = None if not lowres_image_path.exists() else plt.imread(lowres_image_path)

    # Preferentially select high-res, else low-res
    image = None
    image_key = None
    if hires_image is not None:
        image_key = "hires"
        image = hires_image
    elif lowres_image is not None:
        image_key = "lowres"
        image = lowres_image

    # Assign spatial data and images to adata object
    adata.obsm['spatial'] = matched_spatial_data[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values
    adata.obs.drop(columns=["pxl_row_in_fullres", "pxl_col_in_fullres"], inplace=True)  # Drop spatial columns
    if image is not None:
        adata.uns['spatial'][spatial_key]['images'] = {image_key: image}

    return adata


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Load spatial transcriptomics data from MTX matrices and aligned images."
    )
    parser.add_argument(
        "--SRCountDir", required=True, help="Input directory with Spaceranger data."
    )
    parser.add_argument(
        "--outAnnData", required=True, help="Output h5ad file path."
    )
    args = parser.parse_args()

    # Read Visium data
    st_adata = read_visium_mtx(args.SRCountDir, library_id=None, load_images=True)

    # Write raw anndata to file
    with st_adata.file:
        st_adata.write(args.outAnnData)


if __name__ == "__main__":
    main()