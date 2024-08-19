import os

# Configuration parameters
class Config:
    BASE_DIR = '../data/raw/p5/'
    PROCESSED_DIR = '../data/processed/p5_stardist_whole/'
    IMAGE_FILENAME = 'Visium_HD_Human_Colon_Cancer_P5_tissue_image.btf'
    SCALING_FACTORS_JSON = 'spatial/scalefactors_json.json'
    TISSUE_POSITIONS_FILE = 'spatial/tissue_positions.parquet'
    OUTPUT_IMAGE_NAME = 'normalized_image.tiff'
    CLUSTER_RESOLUTION = 0.2
    MIN_PERCENTILE = 5
    MAX_PERCENTILE = 95

    @staticmethod
    def get_image_path():
        return os.path.join(Config.BASE_DIR, Config.IMAGE_FILENAME)

    @staticmethod
    def get_output_image_path():
        return os.path.join(Config.PROCESSED_DIR, Config.OUTPUT_IMAGE_NAME)

    @staticmethod
    def get_tissue_positions_path():
        return os.path.join(Config.BASE_DIR, Config.TISSUE_POSITIONS_FILE)

    @staticmethod
    def get_scaling_factors_path():
        return os.path.join(Config.BASE_DIR, Config.SCALING_FACTORS_JSON)
