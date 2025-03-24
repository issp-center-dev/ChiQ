import numpy as np
import os
import sys
import math
from itertools import product

# add path to BSE_ROOT/python
# print os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append( os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) )
from bse_tools.h5bse import h5BSE


filename='dmft_bse.h5'
basename='Pq'  # output filename like basename_unif.dat
Qs={ 'unif' : '00.00.00', 'stag' : '16.16.00' }
block=(0,0)


if __name__ == '__main__':
    print("=================================")
    print("START get_pq")

    if not os.path.exists(filename):
        raise Exception("file '%s' not found" %filename)

    BS = h5BSE(filename)
    beta = BS.get(key=('beta'))
    print("beta =", beta)

    print("\ngetting X0_loc")
    X0_loc = BS.get(key=('X0_loc',0))
    # print type(X0_loc)
    print(list(X0_loc.keys()))
    print(X0_loc[block].shape)
    assert X0_loc[block].shape[0:2] == (1,1)  # single-orbital
    x0_loc = X0_loc[block][0,0,:]
    # print x0_loc.shape

    nw=x0_loc.shape[0]
    print("NW =", nw)
    # print range(-nw/2,nw/2)
    freqs = np.array([ (2*n+1)*math.pi/beta for n in range(-nw/2,nw/2) ]).reshape((nw,1))  # column vector
    print(freqs.shape)


    print("\ngetting X0_q")
    for label, q in list(Qs.items()):
        X0_q = BS.get(key=('X0_q',0,q))
        # print X0_q.keys()

        x0_q = X0_q[block][0,0,:]

        pq = 1./x0_loc - 1./x0_q  # element-wise
        # print pq

        fileout = basename + '_' + label + '.dat'
        print(fileout)
        # Re P_q  Im P_q
        # np.savetxt(fileout, pq.view(float).reshape(-1,2))
        # w_n  Re P_q  Im P_q
        # np.savetxt( fileout, np.hstack((freqs, pq.view(float).reshape(-1,2))) )

        data = freqs.copy()
        data = np.hstack((data, pq.view(float).reshape(-1,2)))
        data = np.hstack((data, x0_loc.view(float).reshape(-1,2)))
        data = np.hstack((data, x0_q.view(float).reshape(-1,2)))
        print(data.shape)
        # w_n  Re P_q  Im P_q  Re X0_loc  Im X0_loc  Re X0_q  Im X0_q
        np.savetxt(fileout, data)


    print("END get_pq")
    print("=================================")
