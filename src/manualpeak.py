from clicky2d import Clicky2D
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from qkinz import Qkinz
from itertools import zip_longest
import threading
from sklearn.linear_model import LinearRegression
import argparse
import yaml
from typing import Tuple, List, Callable


Peaklist = Tuple[List[int], List[Tuple[int, int]]]
array = np.ndarray
Calibrator = Callable[[array], array]


class Handler2D:
    def __init__(self, path, qkinz_path, isotope, pattern=r"edef*.bin"):
        path = Path(path)
        self.filename_pairs = self.load(path, qkinz_path, pattern, isotope)
        self.peaks_e_handle = (path/"peaks_e_m.csv").open("w")
        self.peaks_Δe_handle = (path/"peaks_de_m.csv").open("w")

    def load(self, path, qkinz_path, pattern, isotope):
        data = sorted(Path(path).glob(pattern))
        protons = read_qkinz(qkinz_path, isotope+"_strip*.txt")
        deutrons = read_qkinz(qkinz_path, isotope+"_d_strip*.txt")
        tritons = read_qkinz(qkinz_path, isotope+"_t_strip*.txt")

        # Pair up the data with Qkinz calculations
        return [paired for paired in
                zip_longest(data, protons, deutrons, tritons, fillvalue=None)]

    def doall(self):
        for front in range(8):
            peaks = self.pick(front)
            self.save_peaks(peaks, front)
        self.peaks_e_handle.close()
        self.peaks_Δe_handle.close()

    def pick(self, front: int):
        data_fname, *qkinz_fnames = self.filename_pairs[front]
        qkinz = [Qkinz(path) for path in qkinz_fnames if path is not None]
        data = read(data_fname)
        hist, *edges = np.histogram2d(data[:, 0], data[:, 1], bins=1000)
        clicky = Clicky2D(hist, *edges)

        # Quality plot
        #fig, ax = plt.subplots(nrows=2)
        ax = None

        fit_thread = threading.Thread(target=Regressor, args=(clicky, qkinz, ax))
        fit_thread.start()
        clicky.ax.set_title(data_fname.stem)
        for q in qkinz:
            q.plot(ax=clicky.ax)

        clicky.ax.set_xlim([0, np.max(edges[0])])
        clicky.ax.set_ylim([0, np.max(edges[1])])
        plt.show()

        clicky.queue.put(("STOP", ))
        fit_thread.join()
        peaks = clicky.queue.get()
        return peaks

    def save_peaks(self, peaks, front):
        def serialize(x, y):
            return ','.join([str(x), str(y)])
        for e, Δe, qe, qΔe in zip(*peaks):
            print(f"Saving {e, Δe, qe, qΔe}")
            self.peaks_e_handle.write(f"{front},{serialize(e, qe)}\n")
            self.peaks_Δe_handle.write(f"{front},{serialize(Δe, qΔe)}\n")
            self.peaks_e_handle.flush()
            self.peaks_Δe_handle.flush()

def read_qkinz(path: Path, pattern: str):
    if not path.exists():
        return []
    qkinz = Path(path).glob(pattern)
    qkinz = list(sorted(qkinz))
    return qkinz


class Regressor:
    def __init__(self, clicky: Clicky2D, qkinz: Qkinz, ax):
        self.plot = clicky
        self.qkinz = qkinz
        self.queue = clicky.queue
        self.ax = ax
        self.peaks_e = []
        self.run()

    def run(self):
        while True:
            peaks: Peaklist = self.queue.get()
            if peaks[0] == "STOP":
                self.queue.put(self.peaks_e)
                return
            self.peaks_e = self.peaks_to_e(peaks)

            coefficients, (calib_e, calib_Δe) = self.regress(peaks, True)
            #self.plot_quality(peaks, calib_e, calib_Δe)
            coefficients, (calib_e, calib_Δe) = self.regress(peaks, False)
            for qkinz in self.qkinz:
                qkinz.calibrate(calib_e, calib_Δe)

    def regress(self, peaks: Peaklist, etoq=True) ->\
        Tuple[Tuple[List[float], List[float]],
              Tuple[Calibrator, Calibrator]]:
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

        coeff_e, predict_e = fit(e, qe)
        coeff_Δe, predict_Δe = fit(Δe, qΔe)

        return (coeff_e, coeff_Δe), (predict_e, predict_Δe)

    def peaks_to_e(self, peaks: Peaklist) ->\
        Tuple[array, array, array, array]:
        peaks = list(zip(*peaks))
        peaks = sorted(peaks, key=lambda x: x[0])
        peaks = list(zip(*peaks))
        qkinznum = peaks[0]
        peaks = peaks[1]

        e, Δe = [], []
        qe, qΔe = [], []

        for i, n in enumerate(qkinznum):
            row = self.qkinz[0].data_.iloc[n]
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

    def plot_quality(self, peaks: Peaklist,
                     calib_e: Calibrator, calib_Δe: Calibrator):
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
        return data.reshape((-1, 2))
    elif fname.suffix == ".csv":
        return np.loadtxt(fname, delimiter=",", skiprows=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("configfile")
    args = parser.parse_args()
    config = yaml.load(open(args.configfile, "r"), Loader=yaml.Loader)
    path = Path(config["read path"])
    # handler = Handler2D("/home/erdos/master/sortering/sirius/",
    #                     "/home/erdos/master/sortering/qkinz/")
    qkinzpath = path/"qkinz"
    assert qkinzpath.exists()
    handler = Handler2D(path, qkinzpath, isotope=config['isotope'])
    handler.doall()
