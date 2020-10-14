from networkx.classes.coreviews import AtlasView


class PyfrostAtlasView(AtlasView):
    """
    The Pyfrost _succ and _pred objects have a specialized __contains__. This class adds support for that function in
    NetworkX AtlasView.
    """

    def __contains__(self, item):
        return item in self._atlas


class PyfrostAdjacencyView(PyfrostAtlasView):
    __slots__ = ()

    def __getitem__(self, name):
        return PyfrostAtlasView(self._atlas[name])

    def copy(self):
        return {n: self[n].copy() for n in self._atlas}
