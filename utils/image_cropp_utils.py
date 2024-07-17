# utils/image_cropp_utils.py

import pandas as pd
import json
from skimage.io import imread, imsave
import matplotlib.pyplot as plt
import tifffile as tiff

def load_dataframe(parquet_path):
    df = pd.read_parquet(parquet_path)
    return df

def read_microns_per_pixel(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data['microns_per_pixel']

def calculate_one_micron(microns_per_pixel):
    return 1 / microns_per_pixel

def get_barcode_coordinates(df, one_micron):
    first_barcode_pxl_row = df.loc[0, 'pxl_row_in_fullres'] - one_micron
    first_barcode_pxl_col = df.loc[0, 'pxl_col_in_fullres'] - one_micron

    num_rows = len(df)
    last_barcode_pxl_row = df.loc[num_rows - 1, 'pxl_row_in_fullres'] + one_micron
    last_barcode_pxl_col = df.loc[num_rows - 1, 'pxl_col_in_fullres'] + one_micron

    return (int(first_barcode_pxl_row), int(first_barcode_pxl_col), 
            int(last_barcode_pxl_row), int(last_barcode_pxl_col))

def read_image(image_path):
    img = imread(image_path)
    return img

def crop_image(img, first_barcode_pxl_row, first_barcode_pxl_col, last_barcode_pxl_row, last_barcode_pxl_col):
    cropped_img = img[last_barcode_pxl_row:first_barcode_pxl_row, first_barcode_pxl_col:last_barcode_pxl_col]
    return cropped_img

def save_image_as_btf(cropped_img, output_path):
    tiff.imwrite(output_path, cropped_img, bigtiff=True)
    print(f"Cropped image saved as {output_path}")
