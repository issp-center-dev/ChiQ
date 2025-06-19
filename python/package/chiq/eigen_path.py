import numpy as np


def _str2latex(str):
    greek = ['Gamma', 'Sigma', 'Delta', 'Lambda']
    if str in greek:
        # print '$\\' + str + '$'
        return rf'$\{str}$'
    else:
        return str


class EigenPath(object):
    """
    Generate x coodinates along the path

    How to use:
        E = EigenPath(file_qpath)
        xarray = E.get_x()
        xticks, xlabels = E.get_xticks()
    """

    def __init__(self, _file_qpath):
        self.__set_qlabels(_file_qpath)

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

    def generate_q(self):
        for q in self.__q_list:
            yield q

    def len_q(self):
        return len(self.__q_list)


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
