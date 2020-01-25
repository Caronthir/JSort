from clicky2d import Clicky2D
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from qkinz import Qkinz
from itertools import product
import threading
from queue import Queue
from sklearn.linear_model import LinearRegression


class Handler2D:
    def __init__(self, path, qkinz_path, pattern=r"edef*.bin",
                qkinz_pattern=r"28Si_strip*.txt"):
        self.filename_pairs = self.load(path, qkinz_path, pattern, qkinz_pattern)

    def load(self, path, qkinz_path, pattern, qkinz_pattern):
        data = Path(path).glob(pattern)
        qkinz = Path(qkinz_path).glob(qkinz_pattern)

        # Pair up the data with Qkinz calculations
        pairs = []
        for d, q in product(data, qkinz):
            f = int(d.stem[-1])-1
            fq = int(q.stem[-1])
            if f == fq:
                pairs.append((f, (d, q)))
        pairs = sorted(pairs, key=lambda x: x[0])
        return [(d, q) for _, (d, q) in pairs]

    def doall(self):
        for front in range(8):
            self.pick(front)

    def pick(self, front: int):
        data_fname, qkinz_fname = self.filename_pairs[front]
        qkinz = Qkinz(qkinz_fname)
        data = read(data_fname)
        hist, *edges = np.histogram2d(data[:, 0], data[:, 1], bins=1000)
        clicky = Clicky2D(hist, *edges)

        # Quality plot
        fig, ax = plt.subplots(nrows=2)

        fit_thread = threading.Thread(target=Regressor, args=(clicky, qkinz, ax))
        fit_thread.start()
        clicky.ax.set_title(data_fname.stem)
        qkinz.plot(ax=clicky.ax)

        plt.show()
        clicky.queue.put(("STOP", ))


class Regressor:
    def __init__(self, clicky: Clicky2D, qkinz: Qkinz, ax):
        self.plot = clicky
        self.qkinz = qkinz
        self.queue = clicky.queue
        self.ax = ax
        self.run()

    def run(self):
        while True:
            peaks = self.queue.get()
            if peaks[0] == "STOP":
                return

            coefficients, (calib_e, calib_Δe) = self.regress(peaks, True)
            self.plot_quality(peaks, calib_e, calib_Δe)
            coefficients, (calib_e, calib_Δe) = self.regress(peaks, False)
            self.qkinz.calibrate(calib_e, calib_Δe)

    def regress(self, peaks, etoq=True):
        e, Δe, qe, qΔe = self.peaks_to_e(peaks)

        coeff_e = []
        coeff_Δe = []

        def fitshift(X, y):
            return y[0] - X[0], lambda x: x + (y[0] - X[0])
        def fitlinear(X, y):
            X = X.reshape(-1, 1)
            reg = LinearRegression().fit(X, y)
            return [reg.intercept_, *reg.coef_], lambda x: reg.predict(x.reshape(-1, 1))

        # Handle a shift separately
        fit = fitshift if len(e) == 1 else fitlinear
        if not etoq:
            e, qe = qe, e
            qΔe, Δe = Δe, qΔe

        coeff_e, predict_e  = fit(e, qe)
        coeff_Δe, predict_Δe = fit(Δe, qΔe)

        return (coeff_e, coeff_Δe), (predict_e, predict_Δe)

    def peaks_to_e(self, peaks):
        peaks = list(zip(*peaks))
        peaks = sorted(peaks, key=lambda x: x[0])
        peaks = list(zip(*peaks))
        qkinznum = peaks[0]
        peaks = peaks[1]

        e, Δe = [], []
        qe, qΔe = [], []

        for i, n in enumerate(qkinznum):
            row = self.qkinz.data_.iloc[n]
            qe.append(row.e)
            qΔe.append(row.Δe)
            # Convert from indices to energy
            i, j = peaks[i]
            x = self.plot.x[i]
            y = self.plot.y[j]
            e.append(x)
            Δe.append(y)

        e, Δe = np.array(e), np.array(Δe)
        qe, qΔe = np.array(qe), np.array(qΔe)

        return e, Δe, qe, qΔe

    def plot_quality(self, peaks, calib_e, calib_Δe):
        e, Δe, qe, qΔe = self.peaks_to_e(peaks)
        self.ax[0].cla()
        self.ax[1].cla()

        self.ax[0].scatter(e, qe)
        self.ax[1].scatter(Δe, qΔe)

        if len(e) > 1:
            E = np.linspace(e[0], e[-1])
            dE = np.linspace(Δe[0], Δe[-1])
            self.ax[0].plot(E, calib_e(E))
            self.ax[1].plot(dE, calib_Δe(dE))

        self.ax[0].figure.canvas.draw_idle()



def read(fname: Path):
    if fname.suffix == ".bin":
        data = np.fromfile(fname, dtype="float32")
        return data[1:].reshape((int(data[0]), 2))
    elif fname.suffix == ".csv":
        return np.loadtxt(fname, delimiter=",", skiprows=1)


if __name__ == '__main__':
    handler = Handler2D("/home/erdos/master/sortering/sirius/",
                        "/home/erdos/master/sortering/qkinz/")
    handler.doall()
