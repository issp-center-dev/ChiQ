import numpy as np
import pathlib
import sys

from .eigen_path import EigenPath



class ChiQWEigenPath(EigenPath):
    """
    Generate list of (x,y) coodinates to be plotted

    How to use:
        E = ChiQEigenPath(file_qpath)
        xarray = E.get_x()
        chiqw_re, chiqw_im = E.get_chiqw_on_path(dir_chiqw=".", file_prefix="chiq_w", ielem=0)
        xticks, xlabels = E.get_xticks()
    """

    def __init__(self, _file_qpath):
        super().__init__(_file_qpath)

    def get_chiqw_on_path(self, dir_chiqw=".", file_prefix="chiq_w", ielem=0, wmin=None, wmax=None):
        """
        Get array of y values

        Return
            ws: np.ndarray(n_w)
            chiqw_re: np.ndarray(n_q, n_w)
            chiqw_im: np.ndarray(n_q, n_w)
        """

        if ielem < 0:
            raise ValueError(f"ielem={ielem} must be non-negative")

        path_chiqw = pathlib.Path(dir_chiqw)
        chiqw_re = np.zeros((0,0), dtype=float)
        chiqw_im = np.zeros((0,0), dtype=float)

        for iq,q in enumerate(self.generate_q()):
            file_chiqw = path_chiqw / f"{file_prefix}_{q}.dat"
            if not file_chiqw.exists():
                raise FileNotFoundError(f"File {file_chiqw} not found")
            with open(file_chiqw, "r") as f:
                chiw_re = []
                chiw_im = []
                ws_local = []
                for line in f:
                    line = line.strip()
                    if line.startswith("#"):
                        continue
                    array = line.split()
                    w = float(array[0])
                    if wmin is not None and w < wmin:
                        continue
                    if wmax is not None and w > wmax:
                        continue
                    nelem = (len(array) - 1) // 2
                    if ielem >= nelem:
                        raise ValueError(f"ielem={ielem} is out of range for {file_chiqw} (nelem={nelem})")
                    ws_local.append(w)
                    chiw_re.append(float(array[ielem*2+1]))
                    chiw_im.append(float(array[ielem*2+2]))
                wnum = len(chiw_re)
                if iq == 0:
                    chiqw_re = np.zeros((self.len_q(), wnum), dtype=float)
                    chiqw_im = np.zeros((self.len_q(), wnum), dtype=float)
                    ws = np.array(ws_local, dtype=float)
                if wnum != chiqw_re.shape[1]:
                    raise ValueError(f"wnum={wnum} at q={q} is not consistent with wnum={chiqw_re.shape[0]} at first q")
                chiqw_re[iq,:] = chiw_re
                chiqw_im[iq,:] = chiw_im
        return ws, chiqw_re, chiqw_im
