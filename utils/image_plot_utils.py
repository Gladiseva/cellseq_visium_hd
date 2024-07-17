# utils/image_plot_utils.py

import matplotlib.pyplot as plt
from shapely.geometry import Polygon

def plot_mask_and_save_image2(title, gdf, img, output_name=None, bbox=None):
    if bbox is not None:
        # Crop the image to the bounding box
        cropped_img = img[bbox[1]:bbox[3], bbox[0]:bbox[2]]
    else:
        cropped_img = img

    # Plot options
    fig, ax = plt.subplots(figsize=(12, 12))

    # Plot the cropped image
    ax.imshow(cropped_img, cmap='gray', origin='lower')
    ax.set_title(title)
    ax.axis('off')

    # Create filtering polygon
    if bbox is not None:
        bbox_polygon = Polygon([(bbox[0], bbox[1]), (bbox[2], bbox[1]), (bbox[2], bbox[3]), (bbox[0], bbox[3])])
        # Filter for polygons in the box
        intersects_bbox = gdf['geometry'].intersects(bbox_polygon)
        filtered_gdf = gdf[intersects_bbox]
    else:
        filtered_gdf = gdf

    # Overlay the filtered polygons on the image
    for poly in filtered_gdf['geometry']:
        x, y = poly.exterior.xy
        ax.plot(x, y, color='red', linewidth=1)  # Change color and linewidth as needed

    # Save the plot if output_name is provided
    if output_name is not None:
        plt.savefig(output_name, bbox_inches='tight')
    else:
        plt.show()
