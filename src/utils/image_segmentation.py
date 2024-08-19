# utils/image_segmentation.py

from csbdeep.utils import normalize
from stardist.models import StarDist2D

def segment_image(cropped_img):
    # Load the pretrained StarDist2D model
    model = StarDist2D.from_pretrained('2D_versatile_he')
    
    # Normalize the image
    min_percentile = 5
    max_percentile = 95
    img = normalize(cropped_img, min_percentile, max_percentile)
    
    # Predict instances using the model
    labels, polys = model.predict_instances_big(img,
                                                axes='YXC',
                                                block_size=1000,
                                                prob_thresh=0.01,
                                                nms_thresh=0.001,
                                                min_overlap=128,
                                                context=128,
                                                normalizer=None,
                                                n_tiles=(4,4,1))
    
    return polys
