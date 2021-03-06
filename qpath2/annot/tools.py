#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

#
# QPATH2.ANNOT.TOOLS - tools for handling various proprietary annotations.
#

__all__ = ['ndpa2xy', 'ndpa_read_single', 'ndpa_read']

import numpy as np
from qpath2.core import Error
import xml.etree.ElementTree as ET

##-
def ndpa2xy(ndpa_pts, level, wsi_params):
    """Convert a set of points from Hamamatsu's annotation file (.ndpa)
    to (x,y) image coordinates.

    Args:
        ndpa_pts (list): a list of (x,y)  wsi coordinates
        level (int): magnification level (0: maximum magnification, 1: half of
            the maximum magnification, etc.)
        wsi_params (dict): a structure holding the parameters describing the
            whole slide image (WSI)

    Returns:
        a numpy.array (n x 2), with (x, y) coordinates by rows
    """
    if wsi_params['vendor'] != 'hamamatsu':
        raise Error('vendor mismatch')
    if level > wsi_params['level_count']:
        raise Error('level out of bounds')

    d = 2**level
    xy_coords = list()
    for p in ndpa_pts:
        x, y = p
        x -= wsi_params['x_offset']
        y -= wsi_params['y_offset']
        x /= (1000 * wsi_params['x_mpp'])
        y /= (1000 * wsi_params['y_mpp'])
        x = long((x + wsi_params['levels'][0]['x_size'] / 2) / d)  # in pixels, relative to UL corner
        y = long((y + wsi_params['levels'][0]['y_size'] / 2) / d)

        xy_coords.append([x, y])

    xy = np.array(xy_coords, dtype=np.int64)

    if np.any(xy < 0):
        raise Error('negative coordinates')

    return xy
##-


##-
def ndpa_read_single(ndpa_file, ann_title):
    """Read a single annotation object from the NDPA file. Note that an
    annotation object may actually have several components, each with the
    same title. All these components are collected in a list of lists.

    Args:
        ndpa_file (str): filename
        ann_title (str): name of the annotation object

    Returns:
        -a list of lists [(x,y),...] with the coordinates of annotation
        points in slide coordinate system
        -None if the annotation object was not found

    See also:
        ndap_read
    """

    xml_file = ET.parse(ndpa_file)
    xml_root = xml_file.getroot()

    xy_coords = []

    for ann in list(xml_root):
        name = ann.find('title').text
        if not name == ann_title:
            continue
        p = ann.find('annotation')
        if p is None:
            continue

        p = p.find('pointlist')
        if p is None:
            continue

        xy_coords.append([(long(pts.find('x').text), long(pts.find('y').text)) for pts in list(p)])

    if len(xy_coords) == 0:
        xy_coords = None

    return xy_coords
##-


##-
def ndpa_read(ndpa_file):
    """Read all annotations.

    Args:
        ndpa_file (str): annotation file name

    Returns:
        -a dictionary with keys corresponding to annotation object
        names and with values the corresponding lists of points

    See also:
        ndpa_read_single
    """

    xml_file = ET.parse(ndpa_file)
    xml_root = xml_file.getroot()

    annot = dict()

    for ann in list(xml_root):
        name = ann.find('title').text
        annot[name] = []

        p = ann.find('annotation')
        if p is None:
            continue

        p = p.find('pointlist')
        if p is None:
            continue

        annot[name].append([(long(pts.find('x').text), long(pts.find('y').text)) for pts in list(p)])

    return annot
##-