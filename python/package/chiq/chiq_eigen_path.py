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
def _load_eigen(_file_eigen, wb=0):
    suscep = {}  # suscep[q] = [suscep1, suscep2, ...]
    with open(_file_eigen, "r") as f:
        for line in f:
            if line[0].startswith('#'):  # skip comment
                continue
            array = line.split()
            if len(array) == 0:  # skip empty line
                continue
            w = int(array[0])
            if w != wb:
                continue
            q = array[1]
            if q in suscep:
                raise Exception(f"ERROR: duplicate key q={q}")
            suscep[q] = np.array(array[2:], dtype=float)
    return suscep


# Read chi(q) data (duplicate of q allowed)
# return dict(list(np.ndarray))
def _load_eigen_dup(_file_eigen, wb=0):
    suscep_dup = defaultdict(list)
    with open(_file_eigen, "r") as f:
        for line in f:
            if line[0].startswith('#'):  # skip comment
                continue
            array = line.split()
            if len(array) == 0:  # skip empty line
                continue
            w = int(array[0])
            if w != wb:
                continue
            q = array[1]
            suscep_dup[q].append(np.array(array[2:], dtype=float))
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

    def __init__(self, _file_qpath, wb=0):
        self.__set_qlabels(_file_qpath)
        self.wb = wb

    def __set_qlabels(self, _file_qpath):
        """Read q-points on the path"""

        self.__qlabels = []  # __qlabels[i] = (x, name)
        self.__q_list = []  # __q_list[i] = q
        self.__x_list = []
        with open(_file_qpath, "r") as f:
            for line in f:
                array = line.split()

                q = array[0]
                ql = q.split(".")
                if len(ql) != 3:
                    raise Exception(f"ERROR: q={q} is not a valid q-point.")
                if not all(map(lambda x: x.isdigit(), ql)):
                    raise Exception(f"ERROR: q={q} is not a valid q-point.")

                x = float(array[1])
                label = array[2] if len(array) >= 3 else ""
                # print array
                self.__q_list.append(q)
                self.__x_list.append(x)
                if label:  # pick up q-point having label
                    self.__qlabels.append((float(x), label))
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
        eigen_q = _load_eigen(_file_eigen, wb=self.wb)
        assert isinstance(eigen_q, dict)
        if not eigen_q:  # empty
            sys.exit(f"ERROR: No data for wb={self.wb} in the File '{_file_eigen}'.")

        for q in self.__q_list:
            if not q in eigen_q:
                eigen_q[q] = [np.nan] * len(list(eigen_q.values())[0])  # put NaN if data not exist

        eigen_array = np.array([eigen_q[q] for q in self.__q_list])
        return eigen_array

    def get_points(self, _file_eigen):
        """
        Get array of x-y values

        Return
            points: np.ndarray(n_points, 1+n_y): 2d array of (x, y1, y2, ...)
        """
        eigen_list_q = _load_eigen_dup(_file_eigen, wb=self.wb)
        assert isinstance(eigen_list_q, dict)

        eigen_list = []
        x_list = []
        for x, q in zip(self.__x_list, self.__q_list):
            if q in eigen_list_q:
                assert isinstance(eigen_list_q[q], list)
                eigen_list += eigen_list_q[q]
                x_list += [x,] * len(eigen_list_q[q])

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
