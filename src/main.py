import sys
import os
import importlib

from utils.image_processing import load_image, normalize_image, save_image
from utils.data_processing import filter_adata_by_area_and_counts, process_ann_data
from utils.clustering import perform_leiden_clustering
from utils.visualization import plot_mask_and_save_image, plot_clusters_and_save_image
from models.nuclei_segmentation import NucleiSegmentation
from config import Config

def main():
    # Step 1: Load and process the image
    print("Start: loading image")
    img = load_image()
    img = normalize_image(img, Config.MIN_PERCENTILE, Config.MAX_PERCENTILE)
    save_image(img)
    print("Finish: image loaded")

    # Step 2: Perform nuclei segmentation
    print("Start: performing nuclei segmentation")
    segmentation_model = NucleiSegmentation()
    gdf = segmentation_model.segment_nuclei(img)
    gdf['area'] = gdf['geometry'].area
    print("Finish: nuclei segmented")

    # Step 3: Load and filter AnnData
    print("Start: processing and filtering andata")
    adata = process_ann_data()
    filtered_adata = filter_adata_by_area_and_counts(adata, gdf)
    print("Finish: pprocessed andata")

    # Step 4: Perform clustering
    print("Start: performing clustering")
    perform_leiden_clustering(filtered_adata, resolution=Config.CLUSTER_RESOLUTION)
    print("Finish: clustering performed")

    # Step 5: Visualize and save results
    print("Start: visualising results")
    plot_mask_and_save_image("Region of Interest 1", gdf, img, adata)
    plot_clusters_and_save_image("Region of Interest 1", gdf, img, filtered_adata)
    print("Finish: pipeline finished")

if __name__ == '__main__':
    main()