#!/usr/bin/env python3

from bse import h5bse
import os
import sys
import numpy as np
import math
import subprocess
import csv
import argparse
import itertools
import logging
from bse.mpi import COMM_WORLD as comm
from more_itertools import divide
from bse import bse_toml
from bse import __version__ as version


def get_calc_flg(w_or_q, target_lists, target="w"):
    """Return True if w_or_q is included in target_lists

    Args:
        w_or_q (int or str): The value of omega (e.g., 0 or '0') or the label for q (e.g., '00.00.00')
        target_lists (list): list of ('0', '00.00.00')
        target (str, optional): 'w' or 'q'. Defaults to "w".

    Returns:
        bool:
    """
    index = {"w": 0, "q": 1}[target]
    for _targetomega in target_lists:
        if _targetomega[index] == str(w_or_q):
            return True
    return False


def launch_subprocess(cmd, rank=None):
    isinstance(cmd, (str, list))

    logfile = f"bse_solver_{rank}.log"
    errfile = f"bse_solver_{rank}.err"
    with open(logfile, "a") as fout, open(errfile, "a") as ferr:
        print(f"\nStart subprocess: {cmd}", file=fout, flush=True)
        proc = subprocess.run(cmd, stdout=fout, stderr=ferr, shell=True)
        print(f"End subprocess", file=fout, flush=True)

    if proc.returncode:
        print(f"ERROR at rank {rank!r}: Subprocess call {proc.args!r} aborted with returncode {proc.returncode}. See {os.path.abspath(errfile)}", file=sys.stderr)
        comm.Abort(proc.returncode)

class chi_base(object):
    #For bse
    def __init__(self, fileName, rank=0, h5_open_mode=None):
        self.rank = rank
        self.h5 = h5bse.h5BSE(fileName)
        if h5_open_mode is not None:
            self.h5.open(h5_open_mode)

        self.input_file_list_for_bse_solver_all = {}
        self.data_diagonal_type_from_hdf5 = {}
        self.data_type_for_bse_solver_at_omega = []
        self.dataset_to_obtain_omega_sum = {"sum":[], "direct":[]}
        self.outputfile = {}
        self.output_chiq = ""
        self.calcflag = True
        self.sumflagChi0Loc = False
        self.sumflagChiLoc = False
        self.datasetname = None
        self.type = "test"
        self.beta = self.h5.get("beta")
        self.block_dim = self._get_block_size("block_name")
        self.inner_dim = self._get_block_size("inner_name")
        self.file_name_init = []
        # self._set_data_info()

    #Basic functions
    def get_matrix(self, key):
        return (self.h5.get(key)) #dictionary {(block1, block2): ndarray.complex}

    def get_matrix_omega(self, omega):
        data_omega = {}
        for key in self.dataset_to_obtain_omega_sum["direct"]:
            data_omega[key] = self.get_matrix((key, omega))
        for key_pair in self.dataset_to_obtain_omega_sum["sum"]:
            data_omega[key_pair[0]] = self._get_sum_matrix((key_pair[0], omega), (key_pair[1], omega))
        return data_omega

    def get_results(self, omega, q):
        data = self._get_chi0q(omega, q)
        return data

    def print_init(self, initFileName, omega):
        with open(initFileName, "w") as initFile:
            initFile.write("suffix {}\n".format(self.rank))
            initFile.write("calctype " + self.type + "\n")
            initFile.write("Matrix_info matrix_info.tmp\n")
            for key, filename in list(self.input_file_list_for_bse_solver_all.items()):
                initFile.write("{} {} \n".format(key, filename).format(omega))
            initFile.write("{} {} \n".format("Beta", str(self.beta)))

    def print_file_w(self, data_list, omega, logger=None):
        assert isinstance(data_list, dict)

        for key in self.data_type_for_bse_solver_at_omega:
            assert key in data_list
            if data_list[key] == False:  # Check if data to write has been loaded from HDF5 file
                print(f"ERROR: Data not found: {key!r} is required.", file=sys.stderr)
                comm.Abort(1)
            self.print_matrix(key, data_list[key], omega=omega, logger=logger)

    def print_matrix_info(self):
        with open("matrix_info.tmp", "w") as file:
            # dimension of block matrix (number of blocks on the column(=row))
            file.write(str(self.block_dim) + "\n")
            # dimensions of inner matrices
            # This is written for each block for future extension to systems like d-f orbital system
            for b in range(self.block_dim):
                file.write(str(self.inner_dim) + "\n")

    def print_matrix(self, datatype, data, omega=None, outputFileName="", logger=None):
        assert isinstance(data, dict)

        mat_info = self.h5.get_mat_info(datatype)  # no access to HDF5 file
        assert mat_info in ['dictA', 'dictB', 'dictC']

        # mat_shape will be used for checking data structure
        if mat_info in ['dictA', 'dictB']:
            niw = data[list(data.keys())[0]].shape[2]
        else:
            niw = 1
        mat_shape_ABC = {
            "dictA" : (self.inner_dim, self.inner_dim, niw, niw),
            "dictB" : (self.inner_dim, self.inner_dim, niw),
            "dictC" : (self.inner_dim, self.inner_dim),
        }
        mat_shape = mat_shape_ABC[mat_info]

        # Convert the matrix structure
        if mat_info == 'dictA':
            # convert [inner, inner, iw, iw]
            #      to [inner, iw, inner, iw]
            data_save = {}
            for _block, _data in data.items():
                assert _data.shape == mat_shape
                data_save[_block] = _data.transpose((0, 2, 1, 3))
            mat_shape = (self.inner_dim, niw, self.inner_dim, niw)
        elif mat_info == 'dictB':
            # convert (block, block) [inner, inner, iw]
            #      to (iw-block, iw-block)] [inner, inner]
            data_save = {}
            for _block, _data in data.items():
                assert _data.shape == mat_shape
                for iw in range(niw):
                    _block_iw = (_block[0] + self.block_dim*iw, _block[1] + self.block_dim*iw)
                    _data_iw = _data[:, :, iw]
                    data_save[_block_iw] = _data_iw
            mat_shape = (self.inner_dim, self.inner_dim)
        else:
            data_save = data

        if outputFileName =="":
            outputFileName = self.input_file_list_for_bse_solver_all[datatype].format(omega)

        # Save data
        with open(outputFileName, "w") as file:
            # Print header
            file.write(str(niw) + "\n")
            # Print data for each block
            for _block, _data in data_save.items():
                assert _data.shape == mat_shape
                file.write("block " + str(_block[0]) + ", " + str(_block[1]) + "\n")
                # for d in _data.ravel():  # flatten
                #     file.write(str(d.real) + ", " + str(d.imag) + "\n")
                _data_1d = np.ascontiguousarray(_data).view(float).reshape(-1, 2)
                np.savetxt(file, _data_1d, delimiter=',')

        # if logger is not None:
        #     logger.info("Filesize: {:,} Bytes {}".format(os.path.getsize(outputFileName), outputFileName))
        # print("Filesize: {:,} Bytes {}".format(os.path.getsize(outputFileName), outputFileName))

    def delete_tmp(self):
        subprocess.call("rm *.tmp", shell=True)

        # subprocess.call("rm init_{}.tmp".format(self.rank), shell=True)
        # subprocess.call("rm matrix_info.tmp".format(self.rank), shell=True)
        # for filename in ["X0_loc", "X0_q"]:
        #     if os.path.isfile(self.input_file_list_for_bse_solver_all[filename]):
        #         subprocess.call("rm {}".format(self.input_file_list_for_bse_solver_all[filename]), shell=True)

    def close_files(self):
        for output in list(self.outputfile.values()):
            output.close()

    # Internal functions
    def _get_block_size(self, keyword):
        return (len(self.h5.get(keyword)))

    def _get_sum(self, key, factor):
        data = self.h5.get(key)
        mat_info = self.h5.get_mat_info(key[0])  # no access to HDF5 file
        if data is False:
            return False
        data_sum = {}
        for _block in list(data.keys()):
            _data = data[_block]
            if mat_info == "dictA":
                data_sum[_block] = np.sum(_data, axis = (2, 3)) * factor
            elif mat_info == "dictB":
                data_sum[_block] = np.sum(_data, axis = (2)) * factor
            else: #dictC
                data_sum[_block] = _data * factor
        return data_sum

    def _calc_susceptibility(self, data, data_sum = None):
        factor = 1.0/self.beta if data_sum is not None else 1.0
        dim = self.inner_dim * self.block_dim
        chi_q = np.zeros((dim, dim), dtype=complex)
        iret = True
        for _data in data:
            real_value = float(_data[2])
            imag_value = float(_data[3])
            # if not a number, not calculate
            iret = not (math.isnan(real_value) or math.isnan(imag_value))
            real_value = 0 if math.fabs(real_value) < math.pow(10.0, -12) else real_value
            imag_value = 0 if math.fabs(imag_value) < math.pow(10.0, -12) else imag_value
            chi_q[int(_data[0])][int(_data[1])] = complex(real_value, imag_value) * factor
            if data_sum is not None:
                chi_q[int(_data[0])][int(_data[1])] += data_sum[int(_data[0])][int(_data[1])]
        return (iret, chi_q)

    def _get_sum_matrix(self, keysum, keyorg):
        data_sum = self.get_matrix(keysum)
        if data_sum is False:
            data_sum = self._get_sum(keyorg, 1.0/self.beta)
        return data_sum

    def output2hdf5(self, key, outputlist, outputdata):
        # outputdata: dim1, dim2, re_value, im_value
        block_dim = self.block_dim
        inner_dim = self.inner_dim
        data_array = np.zeros([block_dim, inner_dim, block_dim, inner_dim], dtype=complex)
        for data in outputlist:
            dim1_block = int(int(data[0]) / inner_dim)
            dim1_inner = int(int(data[0]) % inner_dim)
            dim2_block = int(int(data[1]) / inner_dim)
            dim2_inner = int(int(data[1]) % inner_dim)
            data_array[dim1_block][dim1_inner][dim2_block][dim2_inner] = outputdata[int(data[0])][int(data[1])]
        data = {}
        for dim1, dim2 in itertools.product(list(range(block_dim)), list(range(block_dim))):
            data_array_child = np.zeros([inner_dim, inner_dim], dtype=complex)
            idx = 0
            for inner1, inner2 in itertools.product(list(range(inner_dim)), list(range(inner_dim))):
                if (abs(data_array[dim1][inner1][dim2][inner2]) > pow(10.0, -12)):
                    data_array_child[inner1][inner2] = data_array[dim1][inner1][dim2][inner2]
                    idx += 1
            if (idx != 0):
                data[(dim1, dim2)] = data_array_child
        self.h5.save(key, data)

    def _get_chi0q(self, omega, q):
        outputdata_chi0q = self._read_results("chi0_q_{}.tmp".format(self.rank))
        subprocess.call("rm chi0_q_{}.tmp".format(self.rank), shell=True)
        # 2-3 Calculate chi0q
        # iret, chi0_q = self._calc_susceptibility(outputdata_chi0q, self._reformat_data(data_omega["chi0_loc"]))
        iret, chi0_q = self._calc_susceptibility(outputdata_chi0q)
        if iret is False:
            return None
        return {("chi0_q", str(omega), str(q)): [outputdata_chi0q, chi0_q]}

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

    def _read_results(self, outputFileName):
        outputdata = []
        with open(outputFileName, "r") as file:
            dataReader = csv.reader(file)
            for row in dataReader:
                outputdata.append(row)
        return outputdata

    def _set_data_info(self):
        # Get Data List
        for data_list in self.h5.get_keylist_data():
            dataType = data_list[0]
            if self.data_diagonal_type_from_hdf5.get(dataType) != None:
                self.data_diagonal_type_from_hdf5[dataType].append(data_list)

class chi_bse(chi_base):

    def __init__(self, fileName, rank=0, h5_open_mode=None):
        super(chi_bse, self).__init__(fileName, rank, h5_open_mode)
        self.input_file_list_for_bse_solver_all = {"X0_loc":"X0_loc_w{}.tmp",
                         "X0_q":"X0_q_{}.tmp".format(rank),
                         "X_loc":"X_loc_w{}.tmp",
                         "chi0_loc":"chi0_loc_w{}.tmp",
                         "chi_loc":"chi_loc_w{}.tmp",}
        # "_w{}" will be replaced with "_w0" etc
        self.data_diagonal_type_from_hdf5 = {"X0_loc": [],
                             "X0_q":[],
                             "X_loc":[]}
        self.data_type_for_bse_solver_at_omega = ["X0_loc", "X_loc", "chi_loc", "chi0_loc"]
        self.type = "bse"
        self.readtmpdata = "chi_q_{}.tmp".format(rank)
        self.datasetname = "chi_q"
        self.dataset_to_obtain_omega_sum = {"sum":[("chi0_loc", "X0_loc"), ("chi_loc", "X_loc")],
                                     "direct":["X0_loc", "X_loc" ]}
        self._set_data_info()

    def get_results(self, omega, q):
        data = super(chi_bse, self).get_results(omega, q)
        outputdata = self._read_results(self.readtmpdata)
        subprocess.call("rm " + self.readtmpdata, shell=True)
        # iret, chi_q = self._calc_susceptibility(outputdata, self._reformat_data(data_omega["chi_loc"]))
        iret, chi_q = self._calc_susceptibility(outputdata)
        data[(self.datasetname, str(omega), str(q))] = [outputdata, chi_q]

        # I_q
        outputdata = self._read_results("Iq_{}.tmp".format(self.rank))
        iret, Iq = self._calc_susceptibility(outputdata)
        data[("I_q", str(omega), str(q))] = [outputdata, Iq]

        return data

    # def delete_tmp(self):
    #     super(chi_bse, self).delete_tmp()
    #     if os.path.isfile(self.input_file_list_for_bse_solver_all["X_loc"]):
    #         subprocess.call("rm {}".format(self.input_file_list_for_bse_solver_all["X_loc"]), shell=True)

class chi_rpa(chi_base):
    def __init__(self, fileName, rank=0, h5_open_mode=None):
        super(chi_rpa, self).__init__(fileName, rank, h5_open_mode)
        self.input_file_list_for_bse_solver_all = {"X0_loc":"X0_loc_w{}.tmp",
                         "X0_q":"X0_q_{}.tmp".format(rank),
                         "chi0_loc":"chi0_loc_w{}.tmp",
                         "gamma0":"gamma0.tmp",
                         }
        self.data_diagonal_type_from_hdf5 = {"X0_loc": [],
                             "X0_q":[],
                             "gamma0": []
                             }
        self.data_type_for_bse_solver_at_omega = ["X0_loc", "chi0_loc"]
        self.dataset_to_obtain_omega_sum = {"sum":[("chi0_loc", "X0_loc")],
                                     "direct":["X0_loc"]}
        self.readtmpdata = "chi_q_rpa_{}.tmp".format(rank)
        self.datasetname = "chi_q_rpa"
        self.type = "rpa"
        # Output gamma0
        # self.gamma0 = self.get_matrix("gamma0")
        # self.print_matrix("gamma0", self.gamma0)
        self._set_data_info()

    def get_results(self, omega, q):
        data = super(chi_rpa, self).get_results(omega, q)
        outputdata = self._read_results(self.readtmpdata)
        subprocess.call("rm " + self.readtmpdata, shell=True)
        iret, chi_q = self._calc_susceptibility(outputdata)
        data[(self.datasetname, str(omega), str(q))] = [outputdata, chi_q]
        return data

    def close_files(self):
        for filename in ["gamma0", "chi0_loc"]:
            if os.path.isfile(self.input_file_list_for_bse_solver_all[filename]):
                subprocess.call("rm {}".format(self.input_file_list_for_bse_solver_all[filename]), shell=True)
        super(chi_rpa, self).close_files()


class chi_chi0(chi_base):
    def __init__(self, fileName, rank=0, h5_open_mode=None):
        super(chi_chi0, self).__init__(fileName, rank, h5_open_mode)
        self.input_file_list_for_bse_solver_all = {"X0_loc":"X0_loc_w{}.tmp",
                         "X0_q":"X0_q_{}.tmp".format(rank),
                         "chi0_loc":"chi0_loc_w{}.tmp",
                         }
        self.data_diagonal_type_from_hdf5 = {"X0_q":[],"X0_loc": []}
        self.data_type_for_bse_solver_at_omega = ["X0_loc", "chi0_loc"]
        self.dataset_to_obtain_omega_sum = {"sum":[("chi0_loc", "X0_loc")], "direct":["X0_loc"]}
        #self.datasetname = "chi0_q"
        self.type = "chi0"
        self._set_data_info()


class chi_scl(chi_base):
    def __init__(self, fileName, rank=0, h5_open_mode=None):
        super(chi_scl, self).__init__(fileName, rank, h5_open_mode)
        self.input_file_list_for_bse_solver_all = {"X0_loc":"X0_loc_w{}.tmp",
                         "X0_q":"X0_q_{}.tmp".format(rank),
                         "chi0_loc":"chi0_loc_w{}.tmp",
                         "chi_loc":"chi_loc_w{}.tmp",
                         "Phi":"Phi_w{}.tmp",
                         "Phi_sum":"Phi_sum_w{}.tmp",}
        self.data_diagonal_type_from_hdf5 = {"X0_loc": [],
                             "X0_q":[],
                             "Phi":[],
                             "X_loc":[]}
        self.readtmpdata = "chi_q_scl_{}.tmp".format(rank)
        self.datasetname = "chi_q_scl"
        self.data_type_for_bse_solver_at_omega = ["X0_loc", "Phi", "Phi_sum", "chi0_loc"]
        self.dataset_to_obtain_omega_sum = {"sum":[("chi0_loc", "X0_loc"), ("chi_loc", "X_loc")],
                                     "direct":["X0_loc", "X_loc" ]}
        self.type = "scl"
        self._set_data_info()

    def get_matrix_omega(self, omega):
        data_omega = super(chi_scl, self).get_matrix_omega(omega)
        data_omega["Phi"], data_omega["Phi_sum"] = self.calc_scl(data_omega["chi_loc"], omega)
        return data_omega

    def calc_scl(self, chi_loc, omega):
        data_Phi = self.get_matrix(("Phi", omega))
        data_PhiSum = {}
        # Calculate Phi_sum
        l, P = np.linalg.eigh(self._reformat_data(chi_loc))
        L_sqrt = np.sqrt(l.astype("complex")) * np.identity(P.shape[0])
        data_phisum = np.dot(P, np.dot(L_sqrt, np.conjugate(P.T)))
        data_phisum = data_phisum.reshape(self.block_dim, self.inner_dim, self.block_dim, self.inner_dim).swapaxes(1, 2)
        for _block1, _block2 in itertools.product(list(range(data_phisum.shape[0])), list(range(data_phisum.shape[1]))):
            data_PhiSum[(_block1, _block2)] = data_phisum[_block1][_block2]
        return data_Phi, data_PhiSum

    def get_results(self, omega, q):
        data = super(chi_scl, self).get_results(omega, q)
        outputdata = self._read_results(self.readtmpdata)
        subprocess.call("rm " + self.readtmpdata, shell=True)
        iret, chi_q = self._calc_susceptibility(outputdata)
        data[(self.datasetname, str(omega), str(q))] = [outputdata, chi_q]
        if omega == 0:
            outputdata = self._read_results("Iq_scl_{}.tmp".format(self.rank))
            iret, Iq = self._calc_susceptibility(outputdata)
            data[("I_q_scl", str(omega), str(q))] = [outputdata, Iq]
        return data

    # def delete_tmp(self):
    #     super(chi_scl, self).delete_tmp()
    #     for filename in ["Phi", "Phi_sum"]:
    #         if os.path.isfile(self.input_file_list_for_bse_solver_all[filename]):
    #             subprocess.call("rm {}".format(self.input_file_list_for_bse_solver_all[filename]), shell=True)
    #     subprocess.call("rm {}".format("Iq_scl_{}.tmp".format(self.rank)), shell=True)


class chi_rrpa(chi_base):
    def __init__(self, fileName, rank=0, h5_open_mode=None):
        super(chi_rrpa, self).__init__(fileName, rank, h5_open_mode)
        self.input_file_list_for_bse_solver_all = {"X0_loc":"X0_loc_w{}.tmp",
                         "X0_q":"X0_q_{}.tmp".format(rank),
                         "chi0_loc":"chi0_loc_w{}.tmp",
                         "chi_loc":"chi_loc_w{}.tmp",}
        self.data_type_for_bse_solver_at_omega = ["X0_loc", "chi0_loc", "chi_loc"]
        self.data_diagonal_type_from_hdf5 = {"X0_loc": [],
                             "X0_q":[]
                             }
        self.readtmpdata = "chi_q_rrpa_{}.tmp".format(rank)
        self.datasetname = "chi_q_rrpa"
        self.dataset_to_obtain_omega_sum = {"sum":[("chi0_loc", "X0_loc"), ("chi_loc", "X_loc")],
                                     "direct":["X0_loc", "X_loc"]}
        self.type = "rrpa"
        self._set_data_info()

    def get_results(self, omega, q):
        data = super(chi_rrpa, self).get_results(omega, q)
        outputdata = self._read_results(self.readtmpdata)
        subprocess.call("rm " + self.readtmpdata, shell=True)
        iret, chi_q = self._calc_susceptibility(outputdata)
        data[(self.datasetname, str(omega), str(q))] = [outputdata, chi_q]
        return data

    def close_files(self):
        for filename in ["chi_loc", "chi0_loc"]:
            if os.path.isfile(self.input_file_list_for_bse_solver_all[filename]):
                subprocess.call("rm {}".format(self.input_file_list_for_bse_solver_all[filename]), shell=True)
        super(chi_rrpa, self).close_files()


_version_message = f'BSE version {version}'


def main():
    rank = comm.Get_rank()
    size = comm.Get_size()

    parser = argparse.ArgumentParser(
        prog='bse_tool.py',
        description='Calculate chi0_q and chi_q.',
        add_help=True,
    )

    parser.add_argument('toml', type=str, help='Parameter file in toml format')
    parser.add_argument('--version', action='version', version=_version_message)

    args = parser.parse_args()

    if rank==0:
        print(f"{_version_message}\n")

    # Load parameters from toml file
    dict_common, dict_tool, _ = bse_toml.load_params_from_toml(args.toml, print_summary=(rank==0))
    input_file = dict_common["input"]
    type_list = dict_common["type"]
    path_to_target_list  = dict_common["omega_q"]
    path_to_bse_solver = dict_tool["solver"]
    work_dir = dict_tool["work_dir"]

    # Convert filename to absolute path
    input_file = os.path.abspath(input_file)
    path_to_target_list = os.path.abspath(path_to_target_list)

    # Move into working directory
    if work_dir:
        if rank == 0:
            os.makedirs(work_dir, exist_ok=True)
        comm.barrier()
        os.chdir(work_dir)

    if rank==0:
        print(f"\nStart BSE calculations. See log files in the working directory {os.getcwd()!r}.")

    # Read (w, q) to compute
    #   w : bosonic frequency (0, 1, ...)
    #   q : label for q (e.g., 00.00.00)
    #
    # The file contains space-sparated data like
    #   0 00.00.00
    #   0 00.00.01 (more columns are neglected)
    #   # any comment
    #   ...
    target_lists_all = []  # list of tuple (w, q)
    if (path_to_target_list is None):
        flag_calc_all = True
    else:
        flag_calc_all = False
        with open(path_to_target_list, "r") as file:
            lines = file.readlines()
            for _line in lines:
                # target_lists_all.append(_line.split())
                str_list = _line.split()
                if len(str_list) < 2:
                    continue
                if str_list[0] == '#':
                    continue
                target_lists_all.append(tuple(str_list[:2]))
        # remove duplicate target
        target_lists_all = sorted(set(target_lists_all))

    # Slice target_lists into processes (target_lists is rank dependent)
    target_lists = [tuple(i) for i in divide(size, target_lists_all)][rank]

    logger = logging.getLogger("bse")
    fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
    # logging.basicConfig(level=logging.DEBUG, format=fmt)
    logging.basicConfig(level=logging.INFO, filename="bse_tool_{}.log".format(rank), format=fmt, filemode='w')

    # omega list
    omega_list = sorted(set([w for w, _ in target_lists_all]))

    # debug
    # if rank == 0:
    #     print(f"target_lists_all = {target_lists_all}")
    # comm.barrier()
    # for i in range(size):
    #     if rank == i:
    #         print(f"# rank{rank}: target_lists = {target_lists}")
    #     comm.barrier()
    # print("omega_list =", omega_list)
    # return

    for calc_type in type_list:
        logger.info("Start type = {}".format(calc_type))
        _chi_bse = {
            "bse": chi_bse,
            "rpa": chi_rpa,
            "rrpa": chi_rrpa,
            "chi0": chi_chi0,
            "scl": chi_scl,
        }[calc_type]

        bse_solver = _chi_bse(input_file, rank, h5_open_mode='r')
        # HDF5 file is kept open with file mode 'r' (specified by h5_open_mode)
        # The HDF5 file should be closed manually.

        logger.info("Get Dimension.")
        logger.info("  (block_dim, inner_dim)=(%d, %d)" % (bse_solver.block_dim, bse_solver.inner_dim))
        # logger.info("Get Data Info form bse/input/data_list .")

        # TODO: Specify the number of Matsubara frequencies to be used
        #   default: None (use all data in HDF5 file)
        # n_iw = None

        # Load data (q-independent quantities) from HDF5 and output them into text file. Only rank0 do it. All ranks share the files.
        if rank == 0:
            logger.info("  Load data from HDF5 file.")
            for omega in omega_list:
                # Load data from HDF5 file
                data_omega = bse_solver.get_matrix_omega(omega)
                # Output data into text file
                bse_solver.print_file_w(data_omega, omega, logger=logger)
                # release memory
                del data_omega
            bse_solver.print_matrix_info()

            # Output gamma0
            # FIXME: do this in class method (but not in __init__ because bse_post also instantiates this class)
            if bse_solver.type == "rpa":
                gamma0 = bse_solver.get_matrix("gamma0")
                bse_solver.print_matrix("gamma0", gamma0)
        comm.barrier()

        # Get Matrix for chi0Loc
        # Loop(1): Matsubara
        data = {}
        logger.info("omega-loop start")
        for _X0LocInfo in bse_solver.data_diagonal_type_from_hdf5["X0_loc"]:
            flag_calc = False
            omega = int(_X0LocInfo[1])
            logger.info("  omega = {}".format(str(omega)))

            if flag_calc_all == False :
                if get_calc_flg(omega,target_lists) == False:
                    continue

            # Loop(2): wave_number
            for _X0qInfo in bse_solver.data_diagonal_type_from_hdf5["X0_q"]:
                q = _X0qInfo[2]
                flag_calc = False
                #Skip if omega, q is not selected
                if flag_calc_all == False :
                    omega = _X0qInfo[1]
                    if get_calc_flg(omega, target_lists) == False:
                        continue
                    if get_calc_flg(q, target_lists, target = "q") == False:
                        continue
                logger.info("      q = {}".format(q))
                bse_solver.print_matrix("X0_q", bse_solver.get_matrix(("X0_q", omega, q)), logger=logger)
                # 2-0 Make input file for BSE solver (init.tmp)
                initFileName = "init_{}.tmp".format(rank)
                bse_solver.print_init(initFileName, omega)
                # 2-1 Call BSE solver
                logger.info("        Start: BSE calculation.")
                launch_subprocess("{} {}".format(path_to_bse_solver, initFileName), rank=rank)
                logger.info("          End: BSE calculation.\n")
                # 2-2 Output results (calculate susceptibility)
                data.update(bse_solver.get_results(omega, q))

            # delete tmp files
            # bse_solver.delete_tmp()

        bse_solver.h5.close()

        logger.info("omega-loop end")
        comm.barrier()

        # delete tmp files
        if rank == 0:
            bse_solver.delete_tmp()

        # logger.info("Gather Data")
        # tmp_data_all = comm.gather(data, root=0)

        # if rank == 0:
        #     logger.info("Output Data to HDF5 file.")
        #     data_all = {}
        #     for _data in tmp_data_all:
        #         data_all.update(_data)
        #     for key, value in list(data_all.items()):
        #         bse_solver.output2hdf5(key, value[0], value[1])

        logger.info("Output Data to HDF5 file.")
        # Avoid simultaneous write access to HDF5 file
        for p in range(size):
            if rank == p:
                logger.info("  Start: output2hdf5")
                bse_solver.h5.open('a')
                for key, value in list(data.items()):
                    # logger.info(f"         key = {key}")
                    bse_solver.output2hdf5(key, value[0], value[1])
                bse_solver.h5.close()
                logger.info("    End: output2hdf5")
            comm.barrier()

        bse_solver.close_files()
        logger.info("End type = {}".format(calc_type))
        comm.barrier()

    if rank==0:
        print(f"\nEnd BSE calculations.")


if __name__ == "__main__":
    main()
