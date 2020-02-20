import matplotlib.pyplot as plt
import numpy as np
import turbo
import time
from matplotlib.colors import LogNorm, Normalize
from matplotlib.widgets import RectangleSelector, LassoSelector
from matplotlib.path import Path
from queue import Queue


class Clicky2D():
    def __init__(self, matrix, x, y):
        self.x = x
        self.y = y
        self.values = matrix

        self.patches = []
        self.peaks = []
        self.peak_numbers = []

        fig, self.ax = plt.subplots(figsize=(10, 10))
        self.plot_matrix()
        # cid = fig.canvas.mpl_connect("button_press_event", self.click_factory())
        cid = fig.canvas.mpl_connect("key_press_event", self.press_factory())
        cid = fig.canvas.mpl_connect("scroll_event", self.handle_scroll)

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
        self.rs = RectangleSelector(self.ax, self.zoom_factory(),
                                    drawtype='box', useblit=True, button=[1],
                                    spancoords='pixels', interactive=True,
                                    rectprops={"facecolor":"none"})
        self.ls = LassoSelector(self.ax, self.onselect_factory(), button=[3])

        # For communication
        self.queue = Queue()

    def plot_matrix(self):
        q = self.ax.pcolormesh(self.x, self.y, self.values.T,
                          cmap="turbo", norm=LogNorm())
        self.quad_mesh = q

    def calibrate(self, calibrate_e, calibrate_Δe):
        raise NotImplementedError()
        self.x = calibrate_e(self.x)
        self.y = calibrate_Δe(self.y)

        # Move peaks
        peaks = []
        for i, peak in enumerate(self.peaks):
            self.peaks[i] = [calibrate_e(peak[0]),
                             calibrate_Δe(peak[1])]
        self.quad_mesh.remove()
        self.plot_matrix()
        self.ax.figure.canvas.draw_idle()

    def add_peak(self, peak):
        self.peaks.append(peak)

    def add_peak_number(self, num: int) -> None:
        numpeaks = len(self.peaks)
        numnums = len(self.peak_numbers)
        if numpeaks == 0:
            print("Select a peak before assigning number")
            return
        if numpeaks <= numnums:
            print("All numbers assigned, overwriting last")
            self.peak_numbers.pop()

        if num in self.peak_numbers:
            print(f"Warning: {num} is already assigned a peak")
        self.peak_numbers.append(num)
        self.queue.put((self.peak_numbers,
                        self.peaks[:numnums+1]))

    def remove_peak(self):
        if self.patches:
            for p in self.patches[-1]:
                p.remove()
            del self.patches[-1]
            del self.peaks[-1]
            if len(self.peak_numbers) >= len(self.peaks) + 1:
                del self.peak_numbers[-1]

    def visible(self):
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        x0 = index(self.x, cur_xlim[0])
        x1 = index(self.x, cur_xlim[1])
        y0 = index(self.y, cur_ylim[0])
        y1 = index(self.y, cur_ylim[1])

        values = self.values[x0:x1, y0:y1]
        return values, (x0, x1), (y0, y1)

    def pick(self):
        values, (x0, x1), (y0, y1) = self.visible()
        (xm, ym), (x, y) = centroid(values)
        # Shift the points back
        xm = int(xm + x0)
        ym = int(ym + y0)
        self.add_peak([xm, ym])
        self.add_rectangles(x, y)
        scat = self.ax.scatter(self.x[xm], self.y[ym], color="r", marker="x")
        self.patches[-1].append(scat)
        self.ax.figure.canvas.draw_idle()

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

    def click_factory(self):
        def handle_click(event):
            if (ax := event.inaxes) is None:
                return

            if event.button == 1:  # LEFT
                self.pick()
            elif event.button == 3:  # RIGHT
                pass
        return handle_click

    def press_factory(self):
        press_wait_time = 2  # seconds

        def handle_keypress(event):
            if event.key == 'r':
                self.pick()
            if event.key == 'x':
                self.remove_peak()
                self.ax.figure.canvas.draw_idle()
            if event.key.isdigit():
                digit = int(event.key)
                # Handle composites
                if handle_keypress.do_wait:
                    if time.time() - handle_keypress.last_press > press_wait_time:
                        print("Disregarding composite")
                        handle_keypress.do_wait = False
                        return

                    handle_keypress.digits += event.key
                    if len(handle_keypress.digits) == 2:
                        handle_keypress.do_wait = False
                        digit = int(handle_keypress.digits)
                    else:
                        return
                self.add_peak_number(digit)
            if event.key == '/':
                print("Waiting for composite")
                handle_keypress.digits = ''
                handle_keypress.do_wait = True
                handle_keypress.last_press = time.time()
            if event.key == 'n':
                plt.close('all')
                return

        handle_keypress.last_press = 0
        handle_keypress.do_wait = False
        handle_keypress.digits = ''
        return handle_keypress


    def zoom_factory(self):
        def handle_zoom(event_click, event_release):
            xmin, xmax, ymin, ymax = self.rs.extents
            self.ax.set_xlim([xmin, xmax])
            self.ax.set_ylim([ymin, ymax])
            self.ax.figure.canvas.draw_idle()

        return handle_zoom

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
            xm, ym = centerofmass(values, x, y)
            xm = int(xm + x0)
            ym = int(ym + y0)
            self.add_rectangles(x, y)
            scat = self.ax.scatter(self.x[xm], self.y[ym], color="r", marker="x")
            self.patches[-1].append(scat)
            self.ax.figure.canvas.draw_idle()
            self.add_peak([xm, ym])
            # ax.pcolormesh(self.x[x0:x1], self.y[y0:y1], m.T, norm=LogNorm(), cmap='turbo')
            # plt.show()

        return handle_onselect


    @staticmethod
    def handle_scroll(event, base_scale=1.5):
        if ax := event.inaxes:
            # get the current x and y limits
            cur_xlim = ax.get_xlim()
            cur_ylim = ax.get_ylim()
            cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
            cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
            xdata = event.xdata # get event x location
            ydata = event.ydata # get event y location
            if event.button == 'up':
                # deal with zoom in
                scale_factor = 1/base_scale
            elif event.button == 'down':
                # deal with zoom out
                scale_factor = base_scale
            else:
                # deal with something that should never happen
                scale_factor = 1
            # set new limits
            ax.set_xlim([xdata - cur_xrange*scale_factor,
                         xdata + cur_xrange*scale_factor])
            ax.set_ylim([ydata - cur_yrange*scale_factor,
                        ydata + cur_yrange*scale_factor])
            ax.figure.canvas.draw_idle()

def index(arr, e) -> int:
    return np.argmin(abs(arr-e))

def centroid(matrix, numbins=0.1, mincount=100):
    flat = matrix.flatten()
    # Use numbins% of the highest bins
    numbins = int(len(flat)*numbins)
    # Use at least mincount number or bins, all if necessary
    mincount = len(flat) if len(flat) < mincount else mincount
    numbins = mincount if numbins < mincount else numbins
    indices = flat.argpartition(-numbins)[-numbins:]

    x = indices//matrix.shape[1]
    y = indices%matrix.shape[1]
    xm, ym = centerofmass(matrix, x, y)
    return (xm, ym), (x, y)

def centerofmass(matrix, x, y):
    total = matrix[x, y].sum()
    xm = x@matrix[x, y]/total
    ym = y@matrix[x, y]/total
    return xm, ym


def read(fname: Path):
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
