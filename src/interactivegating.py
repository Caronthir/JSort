import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.widgets import RectangleSelector, LassoSelector
from matplotlib.path import Path
import numpy as np
from pathlib import Path
import turbo


class GatePicker:
    def __init__(self, matrix, x, y):
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.values = np.asarray(matrix)

        self.patches = []
        self.extent = []

        fig, self.ax = plt.subplots(figsize=(10, 10))
        self.plot_matrix()
        # cid = fig.canvas.mpl_connect("button_press_event", self.click_factory())
        # cid = fig.canvas.mpl_connect("key_press_event", self.press_factory())
        #cid = fig.canvas.mpl_connect("scroll_event", self.handle_scroll)

        def format_coord(x, y):
            xarr = self.x
            yarr = self.y
            if ((x > xarr.min()) & (x <= xarr.max())
               & (y > yarr.min()) & (y <= yarr.max())):
                col = np.searchsorted(xarr, x)-1
                row = np.searchsorted(yarr, y)-1
                z = self.values.T[row, col]
                return f'x={x:1.2f}, y={y:1.2f}, z={z:1.2E}'
                # return f'x={x:1.0f}, y={y:1.0f}, z={z:1.3f}   [{row},{col}]'
            else:
                return f'x={x:1.0f}, y={y:1.0f}'
        self.ax.format_coord = format_coord
        self.rs = RectangleSelector(self.ax, self.gate_factory(),
                                    drawtype='box', useblit=True, button=[1],
                                    spancoords='pixels', interactive=True,
                                    rectprops={"facecolor":"none"})
        self.ls = LassoSelector(self.ax, self.onselect_factory(), button=[3])

    def plot_matrix(self):
        q = self.ax.pcolormesh(self.x, self.y, self.values.T,
                          cmap="turbo", norm=LogNorm())
        self.quad_mesh = q

    def gate_factory(self):
        def handle_gate(event_click, event_release):
            xmin, xmax, ymin, ymax = self.rs.extents
            self.ax.figure.canvas.draw_idle()
            plt.close()
            #print(self.rs.extents)
            self.extent = [xmin, xmax, ymin, ymax]

        return handle_gate

    def onselect_factory(self):
        def handle_onselect(vertices):
            path = Path(vertices)
            values, (x0, x1), (y0, y1) = self.visible()
            x, y = np.meshgrid(self.x[x0:x1], self.y[y0:y1])
            x, y = x.ravel(), y.ravel()
            xy = np.stack([x, y]).T
            indices = np.nonzero(path.contains_points(xy))[0]
            x = indices // values.shape[0]
            y = indices %  values.shape[0]
            x, y = y, x
            # fig, ax = plt.subplots()
            # m = np.zeros_like(values)
            # m[x, y] = values[x, y]
            #xm, ym = centerofmass(values, x, y)
            #xm = int(xm + x0)
            #ym = int(ym + y0)
            self.add_rectangles(x, y)
            #scat = self.ax.scatter(self.x[xm], self.y[ym], color="r", marker="x")
            #self.patches[-1].append(scat)
            self.ax.figure.canvas.draw_idle()
            # ax.pcolormesh(self.x[x0:x1], self.y[y0:y1], m.T, norm=LogNorm(), cmap='turbo')
            # plt.show()

        return handle_onselect

    def add_rectangles(self, x, y):
        values, (x0, x1), (y0, y1) = self.visible()
        m = np.zeros_like(values, dtype=bool)
        m[x, y] = True
        m = m.T
        x_, y_ = np.meshgrid(self.x[x0:x1], self.y[y0:y1])
        poss = np.c_[x_[m], y_[m]]
        h = self.x[1] - self.x[0]
        w = self.y[1] - self.y[0]
        patches = []
        for pos in poss:
            r = plt.Rectangle(pos-0.5, h, w, facecolor='none', edgecolor='k', linewidth=1, alpha=0.1)
            self.ax.add_patch(r)
            patches.append(r)
        self.patches.append(patches)


def read(fname: Path):
    fname = Path(fname)
    if fname.suffix == ".bin":
        data = np.fromfile(fname, dtype="float32")
        return data.reshape((-1, 2))
    elif fname.suffix == ".csv":
        return np.loadtxt(fname, delimiter=",", skiprows=1)


if __name__ == "__main__":
    data = read("../si/edef1.bin")
    hist, xedges, yedges = np.histogram2d(data[:, 0], data[:, 1], bins=1000)
    # x = np.linspace(0, 2, 10)
    # y = np.linspace(1, 4.5, 10)
    # vals = np.random.random((10, 10))
    click = GatePicker(hist, xedges, yedges)
    plt.show()
