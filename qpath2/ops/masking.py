#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

#
# OPS.MASKING: various operators for masking operations, including annotation-based masking
#

import lazyflow.operator as lzop
from lazyflow.utility import OrderedSignal


##-
class OpAnnotMasking(lzop.Operator):
    """Defines an masking operator based on annotation."""

    OutsideValue = lzop.InputSlot()   # value to set for pixels "outside" the annotated region
    Annotation = lzop.InputSlot()     # annotation as a dict of lists of contour point coordinates
    InImage = lzop.InputSlot()        # input image to be filtered
    OutImage = lzop.OutputSlot()      # result

    def __init__(self, *args, **kwargs):
        super(OpAnnotMasking, self).__init__(*args, **kwargs)
        self.progressSignal = OrderedSignal()

    def setupOutputs(self):
        assert self.InImage.meta.getAxisKeys() == list('xyc')
        self.OutImage.meta.assignFrom(self.InImage.meta)

    def execute(self, slot, subindex, roi, result):
        key = roi.toSlice()
        if slot.name == 'OutImage':
            self.InImage[key].writeInto(result).wait()
#############
## TODO: complete the operator
        return result

    def propagateDirty(self, dirtySlot, subindex, roi):
        self.OutImage.setDirty(roi)
##-


##-
