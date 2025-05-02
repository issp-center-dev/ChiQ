#!/usr/bin/env python3

import os
import sys
import numpy as np
import argparse
import itertools
import logging
import re
from more_itertools import divide

from chiq import __version__ as version
from chiq import bse_toml
from chiq import h5bse
from chiq.mpi import COMM_WORLD as comm

import bse_solver  # shared library implemented with C++


def get_calc_flg(q, target_lists):
    """Return True if w_or_q is included in target_lists

    Args:
        w_or_q (int or str): The value of omega (e.g., 0 or '0') or the label for q (e.g., '00.00.00')
        target_lists (list): list of ('0', '00.00.00')
        target (str, optional): 'w' or 'q'. Defaults to "w".

    Returns:
        bool:
    """
    for _targetomega in target_lists:
        if _targetomega[0] == q:
            return True
    return False


# key = 'X0_q'          -> return 'X0_q',   False
# key = 'output/chi0_q' -> return 'chi0_q', True
def _parse_input(key: str):
    # print("_parse_input:", input)
    if re.match(r'^output/', key):
        # print(input[7:], True)
        return key[7:], True
    else:
        # print(input, False)
        return key, False


class chi_base(object):
    #For bse
    def __init__(self, file_in, file_out, rank=0, open_h5in=False, n_iw=None):
        self.rank = rank

        # attributes defined in derived class
        self.input_q = None  # 'X0_q' or 'chi0_q'
        self.input_local = []
        self.input_chi_sum = []
        self.read_h5out_for_input = False  # True for RPA and RRPA
        self.output_local = []  # output computed within this class
        self.output_q = []  # output from solver

        # overwrite the above attributes
        self._init_mode_dep()

        # set h5_in
        self.h5_in = h5bse.h5BSE(file_in)
        if open_h5in:
            # open h5_in
            try:
                self.h5_in.open('r')
            except IOError as e:
                print(f"\nError in opening HDF5 input file.", file=sys.stderr)
                print(e, file=sys.stderr)
                sys.exit(1)

        # set h5_out
        if file_in == file_out:
            self.h5_out = self.h5_in
        else:
            self.h5_out = h5bse.h5BSE(file_out)
            # open h5_out
            if open_h5in and self.read_h5out_for_input:
                try:
                    self.h5_out.open('r')
                except IOError as e:
                    print(f"\nError in opening HDF5 output file. RPA and RRPA calculations require chi0(q) data contained in HDF5 output file. Run type='chi0' calculation before RPA and RRPA calculations.", file=sys.stderr)
                    print(e, file=sys.stderr)
                    sys.exit(1)

        # get info
        self.beta = self.h5_in.get("beta")
        self._block_names = self.h5_in.get("block_name")
        self._inner_names = self.h5_in.get("inner_name")
        self.block_dim = len(self._block_names)
        self.inner_dim = len(self._inner_names)
        # self.block_dim = self._get_block_size("block_name")
        # self.inner_dim = self._get_block_size("inner_name")

        self.datalist = self._set_datalist()

        if n_iw is None:
            self.n_iw = self._get_niw()
        else:
            self.n_iw = n_iw

    def _init_mode_dep(self):
        raise NotImplementedError

    def h5_close(self):
        self.h5_in.close()
        if self.read_h5out_for_input and self.h5_in is not self.h5_out:
            self.h5_out.close()

    def check_input(self):
        datalist = self.datalist

        # input_q
        try:
            key, _ = _parse_input(self.input_q)
            assert key in datalist, f"Subgroup '{self.input_q}' not found."

            # input_local
            for input_local in self.input_local:
                key, _ = _parse_input(input_local)
                assert key in datalist, f"Subgroup '{input_local}' not found."

            # input_chi_sum
            for key1, key2, key3 in self.input_chi_sum:
                assert key2 in datalist or key3 in datalist, f"Subgroup '{key2}' or '{key3}' (alternative) not found. Either is required."
        except AssertionError as e:
            print("Error in HDF5 file(s):", e, file=sys.stderr)
            sys.exit(1)


    #Basic functions
    def get_matrix(self, key, from_output=False):
        """Get block matrix from HDF5 file

        Args:
            key: (tuple)
                key = (type,)
                      (type, w,)
                      (type, w, q)

        Returns:
            bm: dict(ndarray)
                Block matrix in 'bse_solver' format
        """
        assert isinstance(key, tuple)
        # return (self.h5.get(key)) #dictionary {(block1, block2): ndarray.complex}

        h5 = self.h5_in if from_output==False else self.h5_out

        bm_h5bse = h5.get(key)
        if bm_h5bse is False:
            return False
        assert isinstance(bm_h5bse, dict)
        # {(block1, block2): ndarray(complex)}

        # dictA, dictB, dictC
        mat_info = h5.get_mat_info(key[0])  # no access to HDF5 file

        # Convert block matrix format from 'h5bse' to 'bse_solver'
        return self._convert_block_matrix(bm_h5bse, mat_info)

    def get_matrix_local(self, omega):
        bms_local = {}
        # for key in self.dataset_to_obtain_omega_sum["direct"]:
        for key in self.input_local:
            _key, from_output = _parse_input(key)
            bms_local[_key] = self.get_matrix((_key, omega), from_output=from_output)
        # for key1, key2, key3 in self.dataset_to_obtain_omega_sum["sum"]:
        for key1, key2, key3 in self.input_chi_sum:
            bms_local[key1] = self._get_sum_matrix((key2, omega), (key3, omega))
        return bms_local

    # Internal functions
    def _get_block_size(self, keyword):
        return (len(self.h5_in.get(keyword)))

    def _get_niw(self):
        # "dictB" : (self.inner_dim, self.inner_dim, niw),
        key = ('X0_loc', 0)
        bm = self.h5_in.get(key)
        if bm == False:
            print(f"\nERROR: Cannot determine num_wf. Input HDF5 file does not have 'X0_loc'. You may specify num_wf manually.", file=sys.stderr)
            exit(1)
        return bm[list(bm.keys())[0]].shape[2]

    def _convert_block_matrix(self, bm_h5bse, mat_info):
        """
        Convert block matrix format and dtype to complex

        from 'h5bse'
        "dictA" : (self.inner_dim, self.inner_dim, niw, niw),
        "dictB" : (self.inner_dim, self.inner_dim, niw),
        "dictC" : (self.inner_dim, self.inner_dim),

        to 'bse_solver' (C++)
        // matrix-A
        //   block: [block], inner: [in, iw]
        // matrix-B
        //   block: [iw, block], inner: [in]
        // matrix-C
        //   block: [block], inner: [in]
        """

        n_block = self.block_dim
        n_inner = self.inner_dim
        n_iw = self.n_iw
        if mat_info == "dictA":
            bm_solver = {}
            for block, array in bm_h5bse.items():
                # cut iw
                array = array[:, :, :n_iw, :n_iw]
                assert array.shape == (n_inner, n_inner, n_iw, n_iw)

                # [in1, in2, iw1, iw2] -> [(in1, iw1), (in2, iw2)]
                # bm_solver[block] = array.transpose([0, 2, 1, 3]).reshape([n_inner*n_iw, n_inner*n_iw]).copy()
                bm_solver[block] = np.array(array.transpose([0, 2, 1, 3]).reshape([n_inner*n_iw, n_inner*n_iw]), dtype=complex)
            del bm_h5bse

        elif mat_info == "dictB":
            bm_solver = {}
            for block, array in bm_h5bse.items():
                # cut iw
                array = array[:, :, :n_iw]
                assert array.shape == (n_inner, n_inner, n_iw)

                # [i, j][in1, in2, iw] -> [(iw, i), (iw, j)][in1, in2]
                i, j = block
                for iw in range(n_iw):
                    # bm_solver[(iw * n_block + i, iw * n_block + j)] = array[:, :, iw].copy()
                    bm_solver[(iw * n_block + i, iw * n_block + j)] = np.array(array[:, :, iw], dtype=complex)
            del bm_h5bse

        elif mat_info == "dictC":
            bm_solver = {}
            for block, array in bm_h5bse.items():
                assert array.shape == (n_inner, n_inner)
                bm_solver[block] = np.array(array, dtype=complex)
            # bm_solver = bm_h5bse

        return bm_solver

    def _get_sum(self, key, factor):
        bm = self.h5_in.get(key)
        if bm is False:
            return False
        mat_info = self.h5_in.get_mat_info(key[0])  # no access to HDF5 file
        bm_sum = {}
        for _block in list(bm.keys()):
            _data = bm[_block]
            if mat_info == "dictA":
                bm_sum[_block] = np.sum(_data, axis = (2, 3)) * factor
            elif mat_info == "dictB":
                bm_sum[_block] = np.sum(_data, axis = (2)) * factor
            else: #dictC
                bm_sum[_block] = _data * factor
        return bm_sum

    def _get_sum_matrix(self, keysum, keyorg):
        bm = self.get_matrix(keysum)
        if bm is False:
            bm = self._get_sum(keyorg, 1.0/self.beta)
        return bm

    def output2hdf5(self, key, bm):
        # key = (type, omega, q)

        # Check NaN, inf
        for mat in bm.values():
            if not np.all(np.isfinite(mat)):
                print(f"Warning: Result {key} contains NaN or infinity", file=sys.stderr)
                return

        # Save only non-zero blocks
        bm_save = {block: mat for block, mat in bm.items() if np.any(np.abs(mat) > 1e-12)}
        self.h5_out.save(key, bm_save)

    def output_info(self):
        self.h5_out.save("beta", self.beta)
        self.h5_out.save("block_name", self._block_names)
        self.h5_out.save("inner_name", self._inner_names)

    def _reformat_data(self, bm):
        _inner_dim = self.inner_dim
        _block_dim = self.block_dim
        _reformat_data = np.zeros((_block_dim * _inner_dim, _block_dim * _inner_dim), dtype=complex)
        for _block in list(bm.keys()):
            block_number1 = _block[0]
            block_number2 = _block[1]
            _data = bm[_block]
            for i, j in itertools.product(list(range(_inner_dim)), list(range(_inner_dim))):
                try:
                    _reformat_data[block_number1 * _inner_dim + i, block_number2 * _inner_dim + j] = _data[i][j]
                except IndexError:
                    pass
        return (_reformat_data)

    # def _set_datalist(self):
    #     # Get Data List
    #     datalist = {'X0_loc': [], 'X0_q': [],}
    #     for key in self.h5_in.get_keylist_data(input_output='input'):
    #         # key = (type, w, q)
    #         assert isinstance(key, tuple)
    #         datatype = key[0]
    #         # print(datatype)
    #         if datatype in datalist:
    #             datalist[datatype].append(key)
    #     return datalist

        # Get Data List
    def _set_datalist(self):
        datalist = {}
        keys = self.h5_in.get_keylist_data(input_output='input')
        if self.read_h5out_for_input:
            keys += self.h5_out.get_keylist_data(input_output='output')
        for key in keys:
            # key = (type, w, q)
            assert isinstance(key, tuple)
            datatype = key[0]
            if datatype not in datalist:
                datalist[datatype] = []
            # print(datatype)
            datalist[datatype].append(key)
        return datalist


class chi_bse(chi_base):
    def _init_mode_dep(self):
        self.input_q = 'X0_q'
        self.input_local = ["X0_loc", "X_loc"]
        self.input_chi_sum = [
            # (key1, key2, key3)
            #   key1: key to address
            #   key2: input
            #   key3: alternative input when key2 does not exists
            ("chi0_loc", "chi0_loc_in", "X0_loc"),
            ("chi_loc", "chi_loc_in", "X_loc"),
        ]
        self.output_local = ['chi0_loc', 'chi_loc']
        self.output_q = ['chi0_q','chi_q', 'I_q']


class chi_chi0(chi_base):
    def _init_mode_dep(self):
        self.input_q = 'X0_q'
        self.input_local = ["X0_loc",]
        self.input_chi_sum = [
            ("chi0_loc", "chi0_loc_in", "X0_loc"),
        ]
        self.output_local = ['chi0_loc',]
        self.output_q = ['chi0_q',]


class chi_rpa(chi_base):
    def _init_mode_dep(self):
        self.input_q = 'output/chi0_q'
        self.input_local = ["gamma0",]
        self.read_h5out_for_input = True
        self.output_q = ['chi_q_rpa',]


class chi_rrpa(chi_base):
    def _init_mode_dep(self):
        self.input_q = 'output/chi0_q'
        self.input_local = ["output/chi0_loc",]
        self.input_chi_sum = [
            ("chi_loc", "chi_loc_in", "X_loc"),
        ]
        self.read_h5out_for_input = True
        self.output_local = ['chi_loc',]
        self.output_q = ['chi_q_rrpa',]

class chi_scl(chi_base):
    def _init_mode_dep(self):
        self.input_q = 'X0_q'
        self.input_local = ["X0_loc", "X_loc"]
        self.input_chi_sum = [
            ("chi0_loc", "chi0_loc_in", "X0_loc"),
            ("chi_loc", "chi_loc_in", "X_loc"),
        ]
        self.output_local = ['chi0_loc', 'chi_loc']
        self.output_q = ['chi0_q', 'chi_q_scl', 'I_q_scl']

    def get_matrix_local(self, omega):
        bms_local = super().get_matrix_local(omega)
        bms_local["Phi"], bms_local["Phi_sum"] = self.calc_scl(bms_local["chi_loc"], omega)
        return bms_local

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


_version_message = f'ChiQ version {version}'


def main():
    rank = comm.Get_rank()
    size = comm.Get_size()

    parser = argparse.ArgumentParser(
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
    file_in = dict_common["input"]
    file_out = dict_common["output"]
    type_list = dict_common["type"]
    path_to_target_list  = dict_common["q_points"]
    num_wb = dict_common["num_wb"]
    if num_wb < 0:
        num_wb = sys.maxsize
    work_dir = dict_tool["work_dir"]
    n_iw = dict_tool["num_wf"]

    # Convert filename to absolute path
    file_in = os.path.abspath(file_in)
    file_out = os.path.abspath(file_out)
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
                # target_lists_all.append(tuple(str_list[:2]))
                target_lists_all.append(str_list[0])
        # remove duplicate target
        target_lists_all = sorted(set(target_lists_all))

    # Slice target_lists into processes (target_lists is rank dependent)
    target_lists = [tuple(i) for i in divide(size, target_lists_all)][rank]

    logger = logging.getLogger("bse")
    fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
    # logging.basicConfig(level=logging.DEBUG, format=fmt)
    logging.basicConfig(level=logging.INFO, filename="chiq_main_{}.log".format(rank), format=fmt, filemode='w')

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

        # Invoke worker class
        chi_worker = _chi_bse(file_in, file_out, rank, open_h5in=True, n_iw=n_iw)
        # HDF5 file 'file_in' is kept open with file mode 'r'
        # The HDF5 file should be closed manually.

        # Check if HDF5 file(s) contains necessary data
        chi_worker.check_input()

        logger.info("Get Dimension.")
        logger.info("  (block_dim, inner_dim)=(%d, %d)" % (chi_worker.block_dim, chi_worker.inner_dim))
        # logger.info("Get Data Info form bse/input/data_list .")

        n_block = chi_worker.block_dim
        n_inner = chi_worker.inner_dim
        n_iw = chi_worker.n_iw
        matinfo_A = np.array([n_iw*n_inner,] * n_block)
        matinfo_B = np.array([n_inner,] * (n_block * n_iw))
        matinfo_C = np.array([n_inner,] * n_block)
        # print("matinfo_A", matinfo_A)
        # print("matinfo_B", matinfo_B)
        # print("matinfo_C", matinfo_C)

        # Make omega list
        # type_q = chi_worker.input_q  # "X0_q" or "chi0_q"
        type_q, from_output_q = _parse_input(chi_worker.input_q)  # "X0_q" or "chi0_q"
        omega_set = {_X0qInfo[1] for _X0qInfo in chi_worker.datalist[type_q]}

        # Get Matrix for chi0Loc
        # Loop(1): Matsubara
        results = {}
        logger.info("omega-loop start")
        # for _X0LocInfo in chi_worker.datalist["X0_loc"]:
        #     # _X0LocInfo = ('X0_loc', omega)
        #     omega = int(_X0LocInfo[1])
        for omega in sorted(omega_set):
            if omega >= num_wb:
                continue
            logger.info("  omega = {}".format(str(omega)))

            # BSE solver
            solver = bse_solver.BSESolver(chi_worker.beta, matinfo_A, matinfo_B, matinfo_C)

            # Set local quantities (X_loc, X0_loc, etc)
            bms_local = chi_worker.get_matrix_local(omega)
            for datatype, bm in bms_local.items():
                solver.set(bm, datatype)
            logger.info(f"    Loaded {list(bms_local.keys())}")

            # local quantities to be saved
            for datatype in chi_worker.output_local:
                results[(datatype, omega)] = bms_local[datatype]

            # Loop(2): wave_number
            for _X0qInfo in chi_worker.datalist[type_q]:
                # _X0LocInfo = ('X0_q', omega, q)
                if _X0qInfo[1] != omega:
                    continue
                q = _X0qInfo[2]
                #Skip if q is not selected
                if q not in target_lists:
                    continue
                logger.info("      q = {}".format(q))

                # Set X0_q or chi0_q
                x0_q = chi_worker.get_matrix((type_q, omega, q), from_output=from_output_q)
                solver.set(x0_q, type_q)
                logger.info(f"        Loaded '{type_q}'")

                # 2-1 Call BSE solver
                logger.info("        Start: BSE calculation.")

                solver.calc(calc_type)
                logger.info("          End: BSE calculation.")

                # Get results
                # 2-2 Output results (calculate susceptibility)
                # results.update(chi_worker.get_results(omega, q))
                for datatype in chi_worker.output_q:
                    bm = solver.get(datatype)
                    # dictC
                    for mat in bm.values():
                        assert mat.shape == (chi_worker.inner_dim, chi_worker.inner_dim)
                    results[(datatype, omega, q)] = bm
                logger.info(f"        Get '{chi_worker.output_q}'")

        chi_worker.h5_close()

        logger.info("omega-loop end")
        comm.barrier()

        # logger.info("Gather Data")
        # tmp_data_all = comm.gather(data, root=0)

        # if rank == 0:
        #     logger.info("Output Data to HDF5 file.")
        #     data_all = {}
        #     for _data in tmp_data_all:
        #         data_all.update(_data)
        #     for key, value in list(data_all.items()):
        #         chi_worker.output2hdf5(key, value[0], value[1])

        logger.info("Output Data to HDF5 file.")
        # Avoid simultaneous write access to HDF5 file
        for p in range(size):
            if rank == p:
                logger.info("  Start: output2hdf5")
                chi_worker.h5_out.open('a')

                # output info
                if rank == 0:
                    logger.info("      info")
                    chi_worker.output_info()

                # output results
                for key, bm in list(results.items()):
                    # print("key =", key)
                    logger.info(f"      key = {key}")
                    chi_worker.output2hdf5(key, bm)
                chi_worker.h5_out.close()
                logger.info("    End: output2hdf5")
            comm.barrier()

        logger.info("End type = {}".format(calc_type))
        comm.barrier()

    if rank==0:
        print(f"\nEnd BSE calculations.")


if __name__ == "__main__":
    main()
