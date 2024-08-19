import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon

class DataFiltering:
    @staticmethod
    def filter_barcodes(result_spatial_join, adata):
        """
        Filter barcodes to keep only those that are within non-overlapping nuclei.
        
        Parameters:
        - result_spatial_join: A GeoDataFrame resulting from a spatial join operation.
        - adata: An AnnData object containing the spatial transcriptomics data.
        
        Returns:
        - A filtered AnnData object.
        """
        # Identify barcodes in overlapping polygons
        barcodes_in_overlapping_polygons = pd.unique(result_spatial_join[result_spatial_join.duplicated(subset=['index'])]['index'])
        
        # Mark barcodes that are in overlapping regions
        result_spatial_join['is_not_in_an_polygon_overlap'] = ~result_spatial_join['index'].isin(barcodes_in_overlapping_polygons)
        
        # Filter to get barcodes in non-overlapping polygons
        barcodes_in_one_polygon = result_spatial_join[
            result_spatial_join['is_within_polygon'] & 
            result_spatial_join['is_not_in_an_polygon_overlap']
        ]
        
        # Create a mask to filter the AnnData object
        filtered_obs_mask = adata.obs_names.isin(barcodes_in_one_polygon['index'])
        return adata[filtered_obs_mask, :]

    @staticmethod
    def filter_by_area_and_counts(adata, gdf, area_cutoff=1000, count_cutoff=100):
        """
        Further filter the AnnData object based on nuclei area and total counts.
        
        Parameters:
        - adata: An AnnData object containing the spatial transcriptomics data.
        - gdf: A GeoDataFrame containing nuclei geometries and areas.
        - area_cutoff: Minimum area of nuclei to include.
        - count_cutoff: Minimum total counts to include.
        
        Returns:
        - A filtered AnnData object.
        """
        # Store the area of each nucleus in the GeoDataFrame
        gdf['area'] = gdf['geometry'].area
        
        # Create masks based on the area and total counts
        mask_area = adata.obs['id'].isin(gdf[gdf['area'] >= area_cutoff]['id'])
        mask_count = adata.obs['total_counts'] > count_cutoff
        
        # Apply both masks to filter the AnnData object
        filtered_adata = adata[mask_area & mask_count, :]
        return filtered_adata

    @staticmethod
    def filter_low_quality_data(adata, quality_metrics):
        """
        Filter the AnnData object based on quality control metrics.
        
        Parameters:
        - adata: An AnnData object containing the spatial transcriptomics data.
        - quality_metrics: Dictionary containing quality control metrics to filter.
        
        Returns:
        - A filtered AnnData object.
        """
        for metric, cutoff in quality_metrics.items():
            if metric in adata.obs.columns:
                adata = adata[adata.obs[metric] > cutoff, :]
        return adata
