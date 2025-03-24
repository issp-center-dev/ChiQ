import numpy as np
import os
import argparse
from scipy import interpolate, optimize

kind_list = ['linear', 'quadratic', 'cubic']

def find_tc(x, y):
    assert x.shape == y.shape

    y_inv = 1./y

    # interp1 = interpolate.interp1d(x, y_inv, kind='linear')
    interp = { kind : interpolate.interp1d(x, y_inv, kind=kind) for kind in kind_list }

    def get_tc(_interp, xmin, xmax):
        # find root in the interval [xmin, xmax]
        return optimize.brentq(lambda z: _interp(z), xmin, xmax)

    tc_list = []
    n = y.shape[0]
    for i in range(1,n):
        if np.sign(y[i-1]) != np.sign(y[i]):
            # get tc
            tc = { kind : get_tc(interp[kind], x[i-1], x[i]) for kind in kind_list }
            tc_list.append( tc )
    return tc_list


if __name__ == '__main__':
    P = argparse.ArgumentParser()
    P.add_argument('file', help="file name")
    P.add_argument('--format', default='detail', choices=['detail', 'tc_only'], help="In 'detail' format (default), all infomation is displayed. In 'tc_only' format, only values of Tc are printed.")
    P.add_argument('--algo', default='quadratic', choices=kind_list, help="default='quadratic'")

    args = P.parse_args()

    flag_print = True if args.format == 'detail' else False
    if flag_print:
        print(args)

    if not os.path.exists(args.file):
        raise Exception("File '%s' not exists" %args.file)

    # read data
    data = np.loadtxt(args.file, ndmin=2)

    if flag_print:
        print(" data shape", data.shape)
        print("mode: {'interpolation method': Tc, ...}")

    for j in range(1, data.shape[1]):
        # print j
        tc_list = find_tc(data[:,0], data[:,j])
        # print tc_list

        if args.format == 'detail':
            if len(tc_list):
                print("%3d :" %j, end=' ')
                for k, tc in enumerate(tc_list):
                    if k:
                        print("     ", end=' ')
                    print(tc)

        elif args.format == 'tc_only':
            if len(tc_list):
                print("", tc_list[0][args.algo], end=' ')
            else:
                print(" ?", end=' ')
