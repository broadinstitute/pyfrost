"""
Module with pyfrost specific exceptions
"""

import networkx


class PyfrostError(Exception):
    pass


class PyfrostInvalidNodeError(PyfrostError, networkx.NodeNotFound):
    pass


class PyfrostInvalidPathError(PyfrostError):
    pass


class PyfrostLinkMismatchError(PyfrostError):
    pass
