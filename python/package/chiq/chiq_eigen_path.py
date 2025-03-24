import numpy as np
from collections import defaultdict
import sys


def _str2latex(str):
    greek = ['Gamma', 'Sigma', 'Delta', 'Lambda']
    if str in greek:
        # print '$\\' + str + '$'
        return rf'$\{str}$'
    else:
        return str


# Read chi(q) data (duplicate of q not allowed)
# return dict(np.ndarray)
def _load_eigen(_file_eigen):
    suscep = {}  # suscep[(w,q)] = [suscep1, suscep2, ...]
    with open(_file_eigen, "r") as f:
        for line in f:
            if line[0] == '#':  # skip comment
                continue
            array = line.split()
            if len(array) == 0:  # skip empty line
                continue
            wq = (array[0], array[1])
            if wq in suscep:
                raise Exception(f"ERROR: duplicate key (w,q)={wq}")
            suscep[wq] = np.array(array[2:], dtype=float)
    return suscep


# Read chi(q) data (duplicate of q allowed)
# return dict(list(np.ndarray))
def _load_eigen_dup(_file_eigen):
    suscep_dup = defaultdict(list)
    with open(_file_eigen, "r") as f:
        for line in f:
            if line[0] == '#':  # skip comment
                continue
            array = line.split()
            if len(array) == 0:  # skip empty line
                continue
            suscep_dup[(array[0], array[1])].append(np.array(array[2:], dtype=float))
    return suscep_dup


class ChiQEigenPath(object):
    """
    Generate list of (x,y) coodinates to be plotted

    How to use:
        E = ChiQEigenPath(file_qpath, file_eigen)
        xarray = E.get_x()
        yarray = E.get_y()
        xticks, xlabels = E.get_xticks()
    """

    def __init__(self, _file_qpath):
        self.__set_qlabels(_file_qpath)

    def __set_qlabels(self, _file_qpath):
        """Read q-points on the path"""

        self.__qlabels = []  # __qlabels[i] = (x, name)
        self.__wq_list = []  # __wq_list[i] = (w,q)
        self.__x_list = []
        with open(_file_qpath, "r") as f:
            for line in f:
                array = line.split()
                # print array
                self.__wq_list.append((array[0], array[1]))
                self.__x_list.append(float(array[2]))
                if len(array) >= 4:  # pick up q-point having label
                    self.__qlabels.append((float(array[2]), array[3]))
        print("qlabels =", self.__qlabels)

    def get_x(self):
        """
        Get array of x coordinate for plot

        Return
            x: np.ndarray(n_q)
        """
        return np.array(self.__x_list, dtype=float)

    def get_y_on_path(self, _file_eigen):
        """
        Get array of y values

        Return
            y: np.ndarray(n_q, n_col)
        """
        eigen_wq = _load_eigen(_file_eigen)
        assert isinstance(eigen_wq, dict)
        if not eigen_wq:  # empty
            sys.exit(f"ERROR: File '{_file_eigen}' is empty.")

        for wq in self.__wq_list:
            if not wq in eigen_wq:
                eigen_wq[wq] = [np.nan] * len(list(eigen_wq.values())[0])  # put NaN if data not exist

        eigen_array = np.array([eigen_wq[wq] for wq in self.__wq_list])
        return eigen_array

    def get_points(self, _file_eigen):
        """
        Get array of x-y values

        Return
            points: np.ndarray(n_points, 1+n_y): 2d array of (x, y1, y2, ...)
        """
        eigen_list_wq = _load_eigen_dup(_file_eigen)
        assert isinstance(eigen_list_wq, dict)

        eigen_list = []
        x_list = []
        for x, wq in zip(self.__x_list, self.__wq_list):
            if wq in eigen_list_wq:
                assert isinstance(eigen_list_wq[wq], list)
                eigen_list += eigen_list_wq[wq]
                x_list += [x,] * len(eigen_list_wq[wq])

        eigen_array = np.array(eigen_list, dtype=float)
        x_array = np.array(x_list, dtype=float)
        return np.hstack([x_array[:, None], eigen_array])

    def get_xticks(self, latex=False):
        """
        Get lists x positions and names

        Return
            x_ticks: list(float)
            x_labels: list(str)
        """

        # list of position/name of high-symmetry points
        x_ticks = [x for x, _ in self.__qlabels]
        x_labels = [_str2latex(label) if latex else label for _, label in self.__qlabels]

        return x_ticks, x_labels
