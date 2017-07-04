#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

import os
import re
import numpy
import vigra
from lazyflow.graph import Operator, InputSlot, OutputSlot
from lazyflow.utility.helpers import get_default_axisordering
import openslide as osl
from qpath2.core import WSIInfo, MRI

##-
class OpNdpiReader(Operator):
    """Operator for reading Hamamatsu's NDPI files, using OpenSlide."""

    name = "OpNdpiReader"

    FilePath = InputSlot(stype='filestring')
    Level = InputSlot(stype='int')
    Output = OutputSlot()

    class NdpiReadError(Exception):
        pass

    def __init__(self, *args, **kwargs):
        super(OpNdpiReader, self).__init__(*args, **kwargs)
        self._filepath = None
        self._level = None
        self._wsi = None
        self._mri = None

    def cleanUp(self):
        self._filepath = None
        self._level = None
        self._wsi = None
        self._mri = None
        super(OpNdpiReader, self).cleanUp()

    def setupOutputs(self):
        self._level = self.Level.value
        self._filepath = self.FilePath.value
        filename = os.path.split(self._filepath)[1]
        self._wsi = WSIInfo(self._filepath)
        self._mri = MRI(self._wsi)

        self.Output.meta.dtype = numpy.dtype('uint8')
        shape = (self._wsi.info['levels'][self._level]['x_size'],
                 self._wsi.info['levels'][self._level]['y_size'],
                 4)
        self.Output.meta.shape = shape
        self.Output.meta.axistags = vigra.AxisTags('xyc')
        self.Output.meta.x_resolution = self._wsi.info['x_mpp']
        self.Output.meta.y_resolution = self._wsi.info['y_mpp']
        self.Output.meta.downsample_factor = self._wsi.info['levels'][self._level]['downsample_factor']

        return

    def execute(self, slot, subindex, roi, result):
        res = self._mri.get_region_px(roi.start[0], roi.start[1],
                                      roi.stop[0]-roi.start[0], roi.stop[1]-roi.start[1],
                                      self._level, as_type=numpy.dtype('uint8'))
        result[...] = res[:, :, roi.start[2]:roi.stop[2]]
        return result

    def propagateDirty(self, slot, subindex, roi):
        if slot == self.FilePath:
            self.Output.setDirty( slice(None) )
        elif slot == self.Level:
            self.Output.setDirty(slice(None))
##-

if __name__ == "__main__":
    from lazyflow.graph import Graph

    graph = Graph()
    opReader = OpNdpiReader(graph=graph)
    opReader.FilePath.setValue('/data/hpath/crc/radboudmc/batch_1/BP7680-1.ndpi')
    opReader.Level.setValue(3)

    print('AxisTags: ' + str(opReader.Output.meta.axistags))
    print('Shape: ' + str(opReader.Output.meta.shape))
    print('Dtype: ' + str(opReader.Output.meta.dtype))
    print('Data [200:210, 1000:1010, 1]' + str(opReader.Output[200:210,1000:1010,1].wait()))
# end
