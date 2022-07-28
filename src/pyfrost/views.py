from networkx.classes.coreviews import AtlasView
from networkx.classes.reportviews import NodeView, NodeDataView


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


class PyfrostNodeDataView(NodeDataView):
    """
    A data view class for Pyfrost's colored compacted de Bruijn graph.

    This NodeDataView subclass adds the option to only iterate over nodes in its forward strand version,
    and thus excludes the reverse complements of nodes.
    """

    def __init__(self, nodedict, data=False, default=None, with_rev_compl=True):
        super().__init__(nodedict, data, default)

        self.with_rev_compl = with_rev_compl

    def __len__(self):
        return super().__len__() if self.with_rev_compl else super().__len__() // 2

    def __iter__(self):
        data = self._data

        if self.with_rev_compl:
            return super().__iter__()
        else:
            if data is False:
                return self._nodes.iter_no_rev_compl()

            if data is True:
                return ((n, self._nodes[n]) for n in self._nodes.iter_no_rev_compl())
            else:
                return ((n, self._nodes[n].get(data, self._default)) for n in self._nodes.iter_no_rev_compl())


class PyfrostNodeView(NodeView):
    def __call__(self, data=False, default=None, with_rev_compl=True):
        if data is False and with_rev_compl is True:
            return self

        return PyfrostNodeDataView(self._nodes, data, default, with_rev_compl)

    def data(self, data=False, default=None, with_rev_compl=True):
        if data is False and with_rev_compl is True:
            return self

        return PyfrostNodeDataView(self._nodes, data, default, with_rev_compl)
