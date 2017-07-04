#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

#
# QPATH2.ANNOT.OPERATORS - various operators for annotation: apply,...
#

__all__ = ['poly_annot_inside']

import vigra
import numpy as np
from ..compgeom import simple_polygon_intersection, rect_inside_polygon
from ..masks import add_region, apply_mask


##-
def poly_annot_inside(img, roi, ann, outside_value=0):
    """Keep only the pixels inside annotated regions, all the rest being set to
    a given value (default 0). Changes are operated in situ, in the specified ROI,
    to all channels.

    Args:
        img (vigra.VigraArray): an input image (ROI)
        roi (tuple): (x0, y0, x1, y1) coordinates of the ROI within the original image
        ann (dict): a dictionary of polygons - one polygon per annotated region
        outside_value (img.dtype): a value to set outside the regions.

    Returns:
        vigra.VigraArray
    """
    x0, y0, x1, y1 = roi
    r = np.array([[x0, y0], [x1, y1]])
    roi_mask = vigra.VigraArray((img.width, img.height), dtype=img.dtype, axistags=vigra.AxisTags('xy'))
    roi_mask.fill(0)

    for a in ann:
        P = np.array(ann[a])

        if rect_inside_polygon([x0, y0], [x1, y1], P):
            # the ROI is inside the annotated region, nothing to change
            return img

        intr = simple_polygon_intersection(r, P)
        if len(intr) > 0:
            # add the intersection to the ROI mask
            #-translate the intersections such that (x0,y0) -> (0,0) (from ROI)
            #-update the mask
            for q in intr:  # for each part of the intersection
                q -= [x0, y0]
                add_region(roi_mask, q)
        # else nothing to do
    apply_mask(img, roi_mask)  # sets to 0 all pixels outside the annotated region

    roi_mask = 1 - roi_mask  # invert the mask
    roi_mask *= outside_value
    ch_iter = img.channelIter()
    for c in ch_iter:
        c += roi_mask  # pixels previously set to 0 now get +outside_value

    return img
##-