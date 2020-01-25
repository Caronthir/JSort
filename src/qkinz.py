import numpy as np
from pathlib import Path
from re import split
import pandas as pd
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)


class Qkinz:
    def __init__(self, path: Path):
        path = Path(path)
        self.data_ = self.read(path)
        self.data = self.data_.copy()
        self.artists = []
        self.ax = None

    @staticmethod
    def read(path: Path):
        with open(path) as infile:
            lines = infile.readlines()

        a0, a1, a2 = 0, 0, 0
        rows = []
        for line in lines:
            tokens = split(r",?\s+", line.strip().rstrip())
            if len(tokens) != 6:
                if tokens and tokens[0].startswith("chi"):
                    a0 = tokens[4][:-3]
                    a1 = tokens[7]
                    a2 = tokens[10][:-8]
            try:
                rows.append([float(token) for token in tokens])
            except ValueError:
                continue
        return pd.DataFrame(rows, columns=['ex', 'Δe', 'dΔe', 'e', 'de', 'etot'])

    def plot(self, ax):
        self.ax = ax
        line = ax.errorbar(self.data.e, self.data.Δe,
                    xerr=self.data.de, yerr=self.data.dΔe,
                    fmt='k-o')
        self.artists.append(line)
        for i, row in self.data.iterrows():
            box = TextArea(str(i))
            ab = AnnotationBbox(box, (row.e, row.Δe-150))
            artist = ax.add_artist(ab)
            self.artists.append(artist)

    def calibrate(self, calib_e, calib_Δe):
        self.data.e = calib_e(self.data_.e.to_numpy())
        self.data.Δe = calib_Δe(self.data_.Δe.to_numpy())
        for artist in self.artists:
            artist.remove()
        self.artists = []
        self.plot(self.ax)
        self.ax.figure.canvas.draw_idle()



