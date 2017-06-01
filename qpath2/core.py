#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

"""QPATH2.CORE: Core classes and functions.

Defines exception classes and other basic classes.
"""

class Error(Exception):
    """Basic error exception for QPATH2.

    Args:
        msg (str): Human-readable string describing the exception.
        code (:obj:`int`, optional): Error code.

    Attributes:
        msg (str): Human-readable string describing the exception.
        code (int): Error code.
    """

    def __init__(self, msg, code=1, *args):
        self.message = "QPATH2: " + msg
        self.code = code
        super(Error, self).__init__(msg, code, *args)

