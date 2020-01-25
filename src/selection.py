import numpy as np

from matplotlib.widgets import LassoSelector
from matplotlib.path import Path


class Lasso:
    """Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to `alpha_other`.
    """

    def __init__(self, ax, button=3):
        self.canvas = ax.figure.canvas

        self.Npts = len(self.xys)

        self.lasso = LassoSelector(ax, onselect=self.onselect, button=[button])
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.canvas.draw_idle()
