from tifffile import imread, imwrite
from csbdeep.utils import normalize
import os
from config import Config

def load_image(image_path=None):
    if image_path is None:
        image_path = Config.get_image_path()
    return imread(image_path)

def normalize_image(image, min_percentile, max_percentile):
    return normalize(image, min_percentile, max_percentile)

def save_image(image, output_path=None):
    if output_path is None:
        output_path = Config.get_output_image_path()
    imwrite(output_path, image)
