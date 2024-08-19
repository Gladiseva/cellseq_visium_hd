import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

def load_dataframe(file_path):
    """Load a DataFrame from a Parquet file."""
    return pd.read_parquet(file_path)


def read_microns_per_pixel(json_path):
    """Read microns per pixel from a JSON file."""
    import json
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data['tissue_hires_scalef']


def calculate_one_micron(microns_per_pixel):
    """Calculate pixel equivalent of one micron."""
    return 1 / microns_per_pixel

def get_barcode_coordinates(df, one_micron):
    """Extract barcode coordinates and calculate the bounding box."""
    first_barcode_pxl_row = df['pxl_row_in_fullres'].min()
    first_barcode_pxl_col = df['pxl_col_in_fullres'].min()
    last_barcode_pxl_row = df['pxl_row_in_fullres'].max() + one_micron
    last_barcode_pxl_col = df['pxl_col_in_fullres'].max() + one_micron
    return first_barcode_pxl_row, first_barcode_pxl_col, last_barcode_pxl_row, last_barcode_pxl_col


def create_geodataframe(df_tissue_positions):
    """Convert DataFrame to GeoDataFrame using spatial coordinates."""
    geometry = [Point(xy) for xy in zip(df_tissue_positions['pxl_col_in_fullres'], df_tissue_positions['pxl_row_in_fullres'])]
    return gpd.GeoDataFrame(df_tissue_positions, geometry=geometry)
    
def process_ann_data(adata, df_tissue_positions):
    """
    Process AnnData object with tissue positions.
    
    Parameters:
    - adata: AnnData object to process.
    - df_tissue_positions: DataFrame with tissue positions.
    
    Returns:
    - Processed AnnData object.
    """
    # Example processing logic
    df_tissue_positions = df_tissue_positions.set_index('barcode')
    adata.obs = pd.merge(adata.obs, df_tissue_positions, left_index=True, right_index=True)
    return adata
   
def filter_adata_by_area_and_counts(adata, gdf, area_cutoff=1000, count_cutoff=100):
    """
    Filter AnnData object based on area and counts.
    
    Parameters:
    - adata: AnnData object to filter.
    - gdf: GeoDataFrame with nuclei areas.
    - area_cutoff: Minimum area of nuclei to include.
    - count_cutoff: Minimum total counts to include.
    
    Returns:
    - Filtered AnnData object.
    """
    gdf['area'] = gdf['geometry'].area
    mask_area = adata.obs['id'].isin(gdf[gdf['area'] >= area_cutoff]['id'])
    mask_count = adata.obs['total_counts'] > count_cutoff
    filtered_adata = adata[mask_area & mask_count, :]
    return filtered_adata

def filter_barcodes(result_spatial_join, adata):
    """Filter barcodes that are within non-overlapping nuclei."""
    barcodes_in_overlaping_polygons = pd.unique(result_spatial_join[result_spatial_join.duplicated(subset=['index'])]['index'])
    result_spatial_join['is_not_in_an_polygon_overlap'] = ~result_spatial_join['index'].isin(barcodes_in_overlaping_polygons)
    barcodes_in_one_polygon = result_spatial_join[result_spatial_join['is_within_polygon'] & result_spatial_join['is_not_in_an_polygon_overlap']]
        
    filtered_obs_mask = adata.obs_names.isin(barcodes_in_one_polygon['index'])
    return adata[filtered_obs_mask, :]
