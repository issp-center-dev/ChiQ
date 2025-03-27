#!/usr/bin/env python3

# coding: utf-8
import numpy as np
import sys
import math
import subprocess
import csv
import argparse
import itertools
import logging
from collections import defaultdict
from more_itertools import divide
import os

from chiq import h5bse
from chiq import bse_toml
import chiq_main # in the same directory


# ============================================================================

class bse_post_base(object):
    def __init__(self, h5bse_obj):
        self.h5 = h5bse_obj
        self.beta = self.h5.get("beta")
        self.output_file_list_for_bsetool = {}
        self.outputfile = {}
        self._output_header()
        # self._set_data_info()
        self.datalist = self._set_datalist()  # dict
        # self.header = "chi_q"
        self.header_type = None  # defined in derived class

    def _output_header(self):
        # Make output file
        for key, name in list(self.output_file_list_for_bsetool.items()):
            self.outputfile[key] = open(name, "w")
            self.outputfile[key].write("#Temperature: {} \n".format(str(1.0 / self.beta)))
            print(f"    Save to '{name}'")

    # def _output_eigen_value(self, outputfile, omegalabel, qlabel, outputdata):
    #     strOutputLine = "{} {}".format(str(omegalabel), str(qlabel)) if qlabel is not None else "{}".format(str(omegalabel))
    #     outputdata_sort = np.sort(outputdata)[::-1]
    #     for _data in outputdata_sort:
    #         strOutputLine += " {}".format(str(np.real(_data)))
    #     outputfile.write("{}\n".format(strOutputLine))

    def _set_datalist(self):
        # Get Data List
        datalist = defaultdict(list)
        for key in self.h5.get_keylist_data(input_output='output'):
            # key = (type, w, q)
            assert isinstance(key, tuple)
            datatype = key[0]
            datalist[datatype].append(key)
        return datalist

    @property
    def datalist_q(self):
        return self.datalist[self.header_type]

    def close_files(self):
        for output in list(self.outputfile.values()):
            output.close()


class bse_post(bse_post_base):
    def __init__(self, bse, mode):
        super(bse_post, self).__init__(bse)

        self.output_file_list_for_bsetool.update({
            "chi_q": f"chi_q_{mode}.dat",
            "chi_loc": f"chi_loc_{mode}.dat",
            "I_q": f"I_q_{mode}.dat",
            "I_r": f"I_r_{mode}.dat",
        })
        self._output_header()
        self.header_type = "chi_q"  # for chi_post.save_vec


class rpa_post(bse_post_base):
    def __init__(self, bse, mode):
        super(rpa_post, self).__init__(bse)

        self.output_file_list_for_bsetool.update({
            "chi_q_rpa": f"chi_q_rpa_{mode}.dat",
        })
        self._output_header()
        self.header_type = "chi_q_rpa"


class scl_post(bse_post_base):
    def __init__(self, bse, mode):
        super(scl_post, self).__init__(bse)

        self.output_file_list_for_bsetool.update({
            "chi_q_scl": f"chi_q_scl_{mode}.dat",
            "chi_loc": f"chi_loc_{mode}.dat",
            "I_q_scl": f"I_q_scl_{mode}.dat",
            "I_r_scl": f"I_r_scl_{mode}.dat",
        })
        self._output_header()
        self.header_type = "chi_q_scl"


class chi0_post(bse_post_base):
    def __init__(self, bse, mode):
        super(chi0_post, self).__init__(bse)
        self.output_file_list_for_bsetool.update({
            "chi0_q": f"chi0_q_{mode}.dat",
        })
        self._output_header()
        # self.header = "chi0_q"
        self.header_type = "chi0_q"


class rrpa_post(bse_post_base):
    def __init__(self, bse, mode):
        super(rrpa_post, self).__init__(bse)
        self.output_file_list_for_bsetool.update({
            "chi_q_rrpa": f"chi_q_rrpa_{mode}.dat",
            "chi_loc": f"chi_loc_{mode}.dat",
        })
        self._output_header()
        self.header_type = "chi_q_rrpa"

# ============================================================================

class chipost_base(object):
    def __init__(self, file_out, datatype, open_h5=False, mode="matrix_element"):
        self.arrayflag = False
        self.h5 = h5bse.h5BSE(file_out)
        if open_h5:
            self.h5.open('r')

        _bse_post = {
            "bse": bse_post,
            "rpa": rpa_post,
            "rrpa": rrpa_post,
            "chi0": chi0_post,
            "scl": scl_post,
        }.get(datatype, bse_post_base)

        self.solver_post = _bse_post(self.h5, mode)

        self.block_dim = self._get_block_size("block_name")
        self.inner_dim = self._get_block_size("inner_name")
        # self.h5.close()

    # Basic functions
    def get_matrix(self, key):
        return (self.h5.get(key))  # dictionary {(block1, block2): ndarray.complex}

    def _reformat_data(self, data):
        _inner_dim = self.inner_dim
        _block_dim = self.block_dim
        _reformat_data = np.zeros((_block_dim * _inner_dim, _block_dim * _inner_dim), dtype=complex)
        for _block in list(data.keys()):
            block_number1 = _block[0]
            block_number2 = _block[1]
            _data = data[_block]
            for i, j in itertools.product(list(range(_inner_dim)), list(range(_inner_dim))):
                try:
                    _reformat_data[block_number1 * _inner_dim + i, block_number2 * _inner_dim + j] = _data[i][j]
                except IndexError:
                    pass
        return (_reformat_data)

    def _get_block_size(self, keyword):
        return (len(self.h5.get(keyword)))

    def post(self, omega, q, data, header=None):
        self.data = self._reformat_data(data)
        if header == None:
            pass
        else:
            # q_str = "" if q is None else f".q{q}"
            # # filename = f"{header}_matrix_element.w{omega}.q{q}.dat"
            # filename = f"{header}_matrix_element.w{omega}" + q_str + ".dat"
            # # np.savetxt(filename, self.data)
            # np.savetxt(filename, self.data.view(float))
            # print(f"    Save to '{filename}'")

            SuscepList = self.data.reshape(-1)
            self._output_suscep(self.solver_post.outputfile[header], omega, q, SuscepList)

    def save_vec(self, header, q):
        pass

    def _load_coefs_file(self, filename):
        print(f"\n    Load {filename!r}")
        if not os.path.exists(filename):
            sys.exit(f"ERROR: File not found: {filename!r}")
        rowVector = np.loadtxt(filename, ndmin=2).view(complex)
        print(f"      shape={rowVector.shape}")
        _, dim = rowVector.shape
        if dim != self.block_dim * self.inner_dim:
            sys.exit(f"ERROR: dimension={dim} is inconsistent with the data size, block*inner={self.block_dim * self.inner_dim}")
        return rowVector.T.copy()  # return column vectors

    def _output_suscep(self, outputfile, omegalabel, qlabel, outputdata):
        assert isinstance(outputdata, np.ndarray)
        assert outputdata.ndim == 1
        strOutputLine = str(omegalabel)
        if qlabel is not None:
            strOutputLine += " " + str(qlabel)

        if outputdata.dtype == float:
            for suscep in outputdata:
                strOutputLine += " " + str(suscep)
        elif outputdata.dtype == complex:
            for suscep in outputdata:
                strOutputLine += " " + str(suscep.real) + " " + str(suscep.imag)
        else:
            raise TypeError(f"outputdata.dtype={outputdata.dtype} is not supported. Only float or complex is supported.")

        outputfile.write(strOutputLine + "\n")

class chipost_eigen(chipost_base):
    def __init__(self, file_out, datatype, open_h5, eigen_order, order_file):
        mode = "eigen"
        super(chipost_eigen, self).__init__(file_out, datatype, open_h5, mode)

        if eigen_order == "file":
            self.arrayflag = True
            self.EigenVecOld = self._load_coefs_file(order_file)
            dim, _ = self.EigenVecOld.shape
            self.IndexListOld = [i for i in range(dim)]
            # check if matrix loaded is unitary
            print("      check if the matrix loaded is unitary", end="")
            if np.allclose(self.EigenVecOld @ self.EigenVecOld.T.conj(), np.identity(dim)):
                print(" --> OK")
            else:
                print(" --> not unitary")
                sys.exit(f"ERROR: Unitary matrix is expected: order_file={order_file!r}")

        self.eigen_order = eigen_order

    def post(self, omega, q, data, header):
        super(chipost_eigen, self).post(omega, q, data)
        self.EigenValue, self.EigenVector = np.linalg.eigh(self.data)
        IndexList = self._calcIndexList()
        self._output_eigen_value_with_label(self.solver_post.outputfile[header], omega, q,
                                                self.EigenValue, IndexList)

    def save_vec(self, header, q):
        rowVector = np.array([self.EigenVecOld[:, i] for i in self.IndexListOld])
        q_str = "" if q is None else f".{q}"
        # filename = f"{header}_eigenvec.{q}.dat"
        filename = f"{header}_eigenvec{q_str}.dat"
        # np.savetxt("{}_eigenvec.{}.dat".format(header, q), rowVector.view(float))
        np.savetxt(filename, rowVector.view(float))
        print(f"    Save to '{filename}'")

    def _output_eigen_value_with_label(self, outputfile, omegalabel, qlabel, outputdata, indexList):
        strOutputLine = "{} {}".format(str(omegalabel), str(qlabel)) if qlabel is not None else "{}".format(
            str(omegalabel))
        for index in indexList:
            strOutputLine += " " + str(np.real(outputdata[index]))
        outputfile.write(strOutputLine + "\n")

    def _calcIndexList(self):
        if self.arrayflag is False:
            IndexList = list(reversed([i for i in range(self.EigenVector.shape[0])]))
            # IndexList = [i for i, eigenvec in enumerate(self.EigenVector)]
            self.arrayflag = True
        else:
            # if self.eigen_order == "overlap":
            if self.eigen_order == "overlap" or self.eigen_order == "file":
                n = len(self.IndexListOld)
                if (n == 0):
                    IndexList = [i for i in range(self.EigenVector.shape[0])]
                else:
                    IndexList = []
                    # find a vector which has the largest overlap with each old vector
                    for j in range(n):
                        # j-th column vector(old)
                        eigenvec_old = self.EigenVecOld[:, self.IndexListOld[j]]
                        # compute overlap with k-th vector(new)
                        # dotList[k] = | < vec_old[j] | vec[k] > |
                        dotList = [np.absolute(np.vdot(eigenvec_old, self.EigenVector[:, k])) for k in range(n)]
                        # sort indices in descending order
                        index_sorted = np.argsort(dotList)[::-1]
                        # get the first index that does not included in the IndexList
                        for i in index_sorted:
                            if i not in IndexList:
                                IndexList.append(i)
                                break
            else:
                IndexList = self.IndexListOld

        if self.eigen_order != "file":
            self.EigenVecOld = self.EigenVector
            self.IndexListOld = IndexList
        return IndexList

class chipost_linear_combination(chipost_base):
    def __init__(self, file_out, datatype, open_h5, coefs_file):
        mode = "linear_combination"
        super(chipost_linear_combination, self).__init__(file_out, datatype, open_h5, mode)

        self.WeightVector = self._load_coefs_file(coefs_file)
        dim, ncol = self.WeightVector.shape
        print(f"      {ncol} kinds of susceptibilities are evaluated.")

    def post(self, omega, q, data, header):
        super(chipost_linear_combination, self).post(omega, q, data)

        # compute susceptibilities and get a real 1d-array
        SuscepList = self._calc_suscep()

        # save the susceptibilities
        self._output_suscep(self.solver_post.outputfile[header],
                            omega, q, SuscepList)

    def _calc_suscep(self):
        SuscepMatrix = self.WeightVector.conjugate().T @ self.data @ self.WeightVector
        SuscepList = np.diagonal(SuscepMatrix)

        # check if susceptibilities are real
        def isreal(array):
            return np.allclose(np.imag(array), np.zeros(array.shape))

        if not isreal(SuscepList):
            # TODO: warning or error
            print(f"Warning: Susceptibilities are not real. Imagnary part (maximum absolute value {np.abs(np.imag(SuscepList)).max():.1e}) is neglected.", file=sys.stderr)

        return np.real(SuscepList)

# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Calculate chi0_q and chi_q.',
        add_help=True,
    )

    parser.add_argument('toml', type=str, help='Parameter file in toml format')

    args = parser.parse_args()

    # Load parameters from toml file
    dict_common, _, dict_post = bse_toml.load_params_from_toml(args.toml)
    file_in = dict_common["input"]
    file_out = dict_common["output"]
    type_list = dict_common["type"]
    path_to_target_list = dict_common["omega_q"]
    OutputVectorFlag = dict_post["vector"]
    order = dict_post["order"]
    order_file = dict_post["order_file"]
    mode_list = dict_post["mode"]
    coefs_file = dict_post["coefs_file"]
    output_dir = dict_post["output_dir"]

    # Convert filename to absolute path
    file_in = os.path.abspath(file_in)
    file_out = os.path.abspath(file_out)
    path_to_target_list = os.path.abspath(path_to_target_list)
    order_file = os.path.abspath(order_file)
    coefs_file = os.path.abspath(coefs_file)

    # Move into output directory
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        os.chdir(output_dir)

    print(f"\nNumerical data are saved in directory {os.getcwd()!r}.")

    target_lists_all = []
    if (path_to_target_list is None):
        flag_calc_all = True
    else:
        flag_calc_all = False
        with open(path_to_target_list, "r") as file:
            lines = file.readlines()
            for _line in lines:
                target_lists_all.append(_line.split())

    # logger = logging.getLogger("bse")
    fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
    # logging.basicConfig(level=logging.DEBUG, format=fmt)
    # logging.basicConfig(level=logging.INFO, filename="logger_post.log", format=fmt)
    done_chi_loc = False
    for calc_type in type_list:
        print(f"\ntype={calc_type!r}")

        _chi_bse = {
            "bse": chiq_main.chi_bse,
            "rpa": chiq_main.chi_rpa,
            "rrpa": chiq_main.chi_rrpa,
            "chi0": chiq_main.chi_chi0,
            "scl": chiq_main.chi_scl,
        }[calc_type]

        for mode in mode_list:
            print(f"\n  mode={mode!r}")
            chi_post_args = dict(
                file_out = file_out,
                datatype = calc_type,
                open_h5 = True,  # Keep the HDF5 file open. Close manurally.
            )
            if mode == "eigen":
                chi_post = chipost_eigen(**chi_post_args, eigen_order=order, order_file=order_file)
            elif mode == "matrix_element":
                chi_post = chipost_base(**chi_post_args)
            elif mode == "linear_combination":
                chi_post = chipost_linear_combination(**chi_post_args, coefs_file=coefs_file)
            else:
                sys.exit("Error: mode = {} is wrong.".format(mode))

            # Common Output chi_loc
            # FIXME: When done_chi_loc==True, chi_loc file is initialized in bse_post_base._output_header() but no data are filled in.
            # if not done_chi_loc:
            if True:
                # logger.info("Output chi_loc")
                if "chi_loc" in chi_post.solver_post.output_file_list_for_bsetool:
                    for _X0LocInfo in chi_post.solver_post.datalist["chi_loc"]:
                        omega = int(_X0LocInfo[1])
                        if flag_calc_all == False:
                            if chiq_main.get_calc_flg(omega, target_lists_all) == False:
                                continue

                        tmp_chi_loc = chi_post.get_matrix(("chi_loc", omega))
                        if tmp_chi_loc != False:
                            chi_post.post(omega, None, tmp_chi_loc, "chi_loc")
                            done_chi_loc = True
                            if OutputVectorFlag:
                                chi_post.save_vec("chi_loc", None)

            # for _X0qInfo in chi_post.solver_post.datalist["chi0_q"]:
            for _X0qInfo in chi_post.solver_post.datalist_q:
                # Calculate chi_q_eigen
                omega = _X0qInfo[1]
                q = _X0qInfo[2]
                if flag_calc_all is False:
                    if chiq_main.get_calc_flg(omega, target_lists_all) == False:
                        continue
                    if chiq_main.get_calc_flg(q, target_lists_all, target="q") == False:
                        continue

                header_type = chi_post.solver_post.header_type
                tmp_chi_q = chi_post.get_matrix((header_type, omega, q))
                if tmp_chi_q != False:
                    # header = chi_post.solver_post.header
                    chi_post.post(omega, q, tmp_chi_q, header_type)
                    if OutputVectorFlag:
                        chi_post.save_vec(header_type, q)

                # Diagonalize Iq: Post process
                for header in ["I_q", "I_r", "I_q_scl", "I_r_scl"]:
                    if header in chi_post.solver_post.output_file_list_for_bsetool:
                        tmp_I = chi_post.get_matrix((header, omega, q))
                        if tmp_I != False:
                            chi_post.post(omega, q, tmp_I, header)
                            if OutputVectorFlag:
                                chi_post.save_vec(header, q)

            chi_post.solver_post.close_files()
            chi_post.h5.close()


if __name__ == "__main__":
    main()
