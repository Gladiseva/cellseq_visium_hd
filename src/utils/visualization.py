import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import geopandas as gpd

def plot_mask_and_save_image(title, gdf, bbox, cmap, img, output_name):
    """Plot and save the masked image with the region of interest."""
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(img[bbox[0]:bbox[2], bbox[1]:bbox[3]], cmap='gray')
    gdf.plot(ax=ax, facecolor="none", edgecolor='red', linewidth=0.5)
    plt.title(title)
    plt.savefig(output_name)
    plt.close()

def plot_clusters_and_save_image(title, gdf, img, adata, bbox, color_by_obs, output_name):
    """Plot clusters and save the image."""
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(img[bbox[0]:bbox[2], bbox[1]:bbox[3]], cmap='gray')
        
    # Color map for clusters
    cmap = ListedColormap(['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00'])
        
    clusters = adata.obs[color_by_obs]
    gdf['cluster'] = clusters
    gdf.plot(column='cluster', ax=ax, cmap=cmap, edgecolor='black', linewidth=0.5, legend=True)
        
    plt.title(title)
    plt.savefig(output_name)
    plt.close()


def plot_nuclei_area_distribution(gdf, area_cutoff):
    """Plot the distribution of nuclei areas before and after filtering."""
    plt.figure(figsize=(8, 6))
    plt.hist(gdf['area'], bins=50, color='skyblue', edgecolor='black')
    plt.axvline(area_cutoff, color='red', linestyle='dashed', linewidth=2)
    plt.title('Nuclei Area Distribution')
    plt.xlabel('Area')
    plt.ylabel('Frequency')
    plt.show()

def plot_total_umi_distribution(adata, umi_cutoff):
    """Plot the distribution of total UMI counts."""
    plt.figure(figsize=(8, 6))
    plt.hist(adata.obs['total_counts'], bins=50, color='skyblue', edgecolor='black')
    plt.axvline(umi_cutoff, color='red', linestyle='dashed', linewidth=2)
    plt.title('Total UMI Counts Distribution')
    plt.xlabel('Total UMI Counts')
    plt.ylabel('Frequency')
    plt.show()
