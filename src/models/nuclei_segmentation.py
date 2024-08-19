from stardist.models import StarDist2D
import geopandas as gpd
from shapely.geometry import Polygon

class NucleiSegmentation:
    def __init__(self, model_name='2D_versatile_he'):
        self.model = StarDist2D.from_pretrained(model_name)
    
    def segment_nuclei(self, img, **kwargs):
        labels, polys = self.model.predict_instances_big(img, **kwargs)
        geometries = [Polygon([(y, x) for x, y in zip(p[0], p[1])]) for p in polys['coord']]
        return gpd.GeoDataFrame(geometry=geometries)

