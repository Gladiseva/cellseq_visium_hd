# utils/geodataframe_creation.py

import geopandas as gpd
from shapely.geometry import Polygon

def create_geodataframe(polys):
    # Creating a list to store Polygon geometries
    geometries = []
    
    # Iterating through each nuclei in the 'polys' DataFrame
    for nuclei in range(len(polys['coord'])):
        # Extracting coordinates for the current nuclei and converting them to (y, x) format
        coords = [(y, x) for x, y in zip(polys['coord'][nuclei][0], polys['coord'][nuclei][1])]
        
        # Creating a Polygon geometry from the coordinates
        geometries.append(Polygon(coords))
    
    # Creating a GeoDataFrame using the Polygon geometries
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf['id'] = [f"ID_{i+1}" for i, _ in enumerate(gdf.index)]
    
    return gdf
