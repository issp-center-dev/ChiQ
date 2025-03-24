import os
import sys
import numpy as np
import h5py as h5
import platform

is_python3 = int(platform.python_version_tuple()[0]) == 3

def _from_bytes_to_utf8(s):
    """
    from bytes to string
    :param s:
    :return:
    """
    if is_python3 and isinstance(s, bytes):
        return s.decode('utf-8')
    else:
        return s


class _h5py_File(object):
    """
    An extension of h5py.File. Two different usages:

    (i)

      with _h5py_File(filename, mode) as f:

    is the same as h5py.File.

    (ii)

      with _h5py_File(filename, mode, h5file) as f:

    h5file is an instance of h5py.File. In this case, h5file enters into f. h5file is not closed when the program gets out of the 'with' block.
    """
    def __init__(self, filename, mode, h5file=None):
        self._filename = filename
        self._mode = mode
        self._h5file = h5file

    def __enter__(self):
        if self._h5file is None:
            self._h5file_local = h5.File(self._filename, self._mode)
            return self._h5file_local
        else:
            return self._h5file

    def __exit__(self, exc_type, exc_value, traceback):
        if self._h5file is None:
            self._h5file_local.close()
        else:
            pass


class h5BSE(object):
    def __init__(self, filename, groupname="", compression="gzip", compression_opts=None):
        self._filename = filename
        self._groupname = groupname

        if compression == "gzip":
            self._compression = "gzip"
            if compression_opts is None:
                self._compression_opts = 4
            elif isinstance(compression_opts, int):
                self._compression_opts = compression_opts
            else:
                self._compression_opts = compression_opts["compression_opts"]
        elif compression == "lzf":
            self._compression = "lzf"
        else:
            try:
                import hdf5plugin
                f_cls = hdf5plugin.get_filters(compression)[0]
                if compression_opts is None:
                    self._compression = f_cls()
                else:
                    self._compression = f_cls(**compression_opts)
            except ImportError:
                print(f"Error: h5py does not support the requested compression method '{compression}'. Use 'gzip' or 'lzf', or install hdf5plugin.")
                sys.exit(1)
            except ValueError:
                print(f"Error: Neither h5py nor hdf5plugin support the requested compression method '{compression}'.")
                print(f"Available compression methods:")
                print(f"  h5py:")
                for cm in ("gzip", "lzf"):
                    print(f"    {cm}")
                print(f"  hdf5plugin:")
                for cm in hdf5plugin.get_filters():
                    print(f"    {cm.filter_name}")
                sys.exit(1)
            except TypeError as e:
                print(f"TypeError: {e}")
                print(f"INFO: There is an invalid compression option.")
                sys.exit(1)

        # self._h5file = h5.File(filename, status)
        self._DataType = ["dictA", "dictB", "dictC"]
        self._keyList = {
            "beta": ("info", "float"),
            "block_name": ("info", "list"),
            "inner_name": ("info", "list"),
            "X_loc": ("input", "dictA", "w"),
            "X0_loc": ("input", "dictB", "w"),
            "X0_q": ("input", "dictB", "wq"),
            "chi_loc_in": ("input", "dictC", "w"),
            "chi0_loc_in": ("input", "dictC", "w"),
            "gamma0": ("input", "dictC"),
            "Phi" : ("input", "dictB", "w"),
            "Phi_sum": ("input", "dictC", "w"),
            "chi_q": ("output", "dictC", "wq"),
            "chi0_q": ("output", "dictC", "wq"),
            "chi_loc": ("output", "dictC", "w"),
            "chi0_loc": ("output", "dictC", "w"),
            "chi_q_rpa": ("output", "dictC", "wq"),
            "chi_q_rrpa": ("output", "dictC", "wq"),
            "chi_q_scl": ("output", "dictC", "wq"),
            "I_q":("output", "dictC", "wq"),
            "I_r":("output", "dictC", "wq"),
            "I_q_scl":("output", "dictC", "wq"),
            "I_r_scl":("output", "dictC", "wq"),
        }

        # HDF5 file object; Used by open/close method (optional)
        self._h5file = None

    def open(self, mode='r'):
        """
        To improve peformance, HDF5 file can be opened manually. If this method is called, close() must be called to close the file manually. If this method is NOT called, the file is opened/closed every time.
        """
        # check if file exists
        if mode=='r':
            if not os.path.isfile(self._filename):
                raise IOError(f"File not found: '{self._filename}'.")

        self._h5file = h5.File(self._filename, mode)

    def close(self):
        self._h5file.close()
        self._h5file = None

    def save(self, key, data):
        if type(key) == type(""):
            key = tuple([key])
        try:
            keyinfo = self._get_keyinfo(key)
        except KeyError:
            # print(f"key={key!r} is not defined.")
            # return False
            raise Exception(f"key={key!r} is not defined.")

        if self._check_key(key, keyinfo) is False:
            return False

        if keyinfo[0] == "info":
            if self._save_info(key, data) is False:
                return False
        elif keyinfo[0] == "input" or keyinfo[0] == "output":
            if self._save_data(key, data) is False:
                return False
        else:
            return False
        return True

    def _get_keyinfo(self, key):
        if type(key) == type(()):
            key = key[0]
        keyinfo = self._keyList[_from_bytes_to_utf8(key)]
        return keyinfo

    def get_mat_info(self, key):
        # if type(key) == type(()):
        if isinstance(key, tuple):
            key = key[0]
        keyinfo = self._keyList[_from_bytes_to_utf8(key)]
        if keyinfo[0] == "info":
            return False
        else:  # 'input' or 'output'
            return keyinfo[1]  # 'dictA', 'dictB', or 'dictC'

    def _check_key(self, key, keyinfo):
        # Check key
        if keyinfo[0] == "info":
            if len(key) != 1:
                return False
            return True
        elif keyinfo[0] == "input" or keyinfo[0] == "output":
            if (len(keyinfo) == 2 and len(key) == 1) \
                    or (keyinfo[2] == "w" and len(key) == 2) \
                    or (keyinfo[2] == "wq" and len(key) == 3):
                return True
            else:
                return False
        else:
            # print(("key: " + str(key) + " is invalid."))
            # return False
            raise Exception("key: " + str(key) + " is invalid.")

    def _save_info(self, key, data):
        keyinfo = self._get_keyinfo(key)
        # make dirname
        if self._groupname != "":
            dirname = "bse/" + self._groupname + "/" + keyinfo[0] + "/" + key[0]
        else:
            dirname = "bse/" + keyinfo[0] + "/" + key[0]
        # save info
        if keyinfo[1] == "list":
            for i, _data in enumerate(data):
                childdir = dirname + "/" + str(i)
                self._write_data(childdir, _data)
        else:
            self._write_data(dirname, data)

    # def _save_data(self, key, data):
    #     # key = (type, w, q)
    #     # check type
    #     keyinfo = self._get_keyinfo(key)
    #     if self._check_type(keyinfo[1], data) is False:
    #         # print("Key type is invalid.")
    #         # return False
    #         raise Exception(f"Key type is invalid: keyinfo={keyinfo}, data={data}")
    #     # make dirname
    #     dirname = self._make_dirname(key, keyinfo)
    #     # save info (key, bosonic_freq, q)
    #     self._write_data(dirname + "/info/key", _from_bytes_to_utf8(key[0]))
    #     if len(keyinfo) != 2:
    #         self._write_data(dirname + "/info/bosonic_freq", int(key[1]))
    #         if keyinfo[2] == "wq":
    #             self._write_data(dirname + "/info/q", _from_bytes_to_utf8(str(key[2])))
    #     # save data
    #     for _data in list(data.keys()):
    #         groupname = str(_data[0]) + "_" + str(_data[1])
    #         self._write_data(dirname + "/data/" + groupname, data[_data])

    def _save_data(self, key, data):
        # key = (type, w, q)
        keyinfo = self._get_keyinfo(key)
        if self._check_type(keyinfo[1], data) is False:
            raise Exception(f"Key type is invalid: keyinfo={keyinfo}, data={data}")
        dirname = self._make_dirname(key, keyinfo)

        # save data
        with _h5py_File(self._filename, "a", self._h5file) as _h5file:
            if dirname not in _h5file:
                grp = _h5file.create_group(dirname)
            else:
                grp = _h5file[dirname]
            for (i, j), _data in data.items():
                dsetname = f"{i}_{j}"  # (0, 1) -> '0_1'
                if dsetname in grp:
                    del grp[dsetname]
                if self._compression == "gzip":
                    grp.create_dataset(dsetname, data=_data, compression=self._compression, compression_opts=self._compression_opts)
                else:
                    grp.create_dataset(dsetname, data=_data, compression=self._compression)

    def _check_type(self, keytype, data):
        # Check type is dictionary
        if type(data) != type({}):
            return False
        else:
            # Check type is tupple
            for _data in list(data.keys()):
                if type(_data) != type(()):
                    return False
            # Check Data Type
            for _data in list(data.values()):
                if keytype == "dictA":
                    if (len(_data.shape) != 4):
                        return False
                elif keytype == "dictB":
                    if (len(_data.shape) != 3):
                        return False
                elif keytype == "dictC":
                    if (len(_data.shape) != 2):
                        return False
        return True

    def _make_dirname(self, key, keyinfo):
        keytype = key[0]
        if self._groupname != "":
            dirname = "bse/" + self._groupname + "/" + keyinfo[0]
        else:
            dirname = "bse/" + keyinfo[0]

        if len(keyinfo) == 2:
            subdir = "/" + key[0]
        elif keyinfo[2] == "w":
            # subdir = "/" + key[0] + "_w" + str(key[1])
            subdir = "/" + key[0] + "/w" + str(key[1])
        elif keyinfo[2] == "wq":
            # subdir = "/" + str(key[0]) + "_w" + str(key[1]) + "_q_" + str(key[2])
            subdir = "/" + str(key[0]) + "/w" + str(key[1]) + "/q_" + str(key[2])
        else:
            return False

        dirname = dirname + subdir
        return (dirname)

    def _write_data(self, path, data):
        with _h5py_File(self._filename, "a", self._h5file) as _h5file:
            if path in _h5file:
                del _h5file[path]
            _h5file[path] = data

    def get(self, key):
        if type(key) == type(""):
            key = tuple([key])
        try:
            keyinfo = self._get_keyinfo(key)
        except KeyError:
            # print(f"key={key!r} is not defined.")
            # return False
            raise Exception(f"key={key!r} is not defined.")
        if keyinfo[0] == "info":
            data = self._get_info(key)
        elif keyinfo[0] == "input" or keyinfo[0] == "output":
            data = self._get_data(key)
        else:
            return False
        return data

    # def get_keylist_data(self):
    #     dirname = "bse/" + self._groupname + "/input"
    #     with _h5py_File(self._filename, "r", self._h5file) as _h5file:
    #         datalist = list(_h5file[dirname].keys())
    #         key_list = []
    #         for _data in datalist:
    #             childdir = dirname + "/" + str(_data) + "/info/key"
    #             _key = _from_bytes_to_utf8(_h5file[childdir][()])
    #             if(_key != "gamma0"):
    #                 childdir = dirname + "/" + str(_data) + "/info/bosonic_freq"
    #                 _bosonic_freq = int(_h5file[childdir][()])
    #                 keyinfo = self._get_keyinfo(_key)
    #                 if keyinfo[2] == "wq":
    #                     childdir = dirname + "/" + str(_data) + "/info/q"
    #                     _q = _from_bytes_to_utf8(_h5file[childdir][()])
    #                     key_list.append((_key, _bosonic_freq, _q))
    #                 else:
    #                     key_list.append((_key, _bosonic_freq))
    #             else:
    #                key_list.append((_key))
    #     print(key_list)
    #     return(key_list)

    def get_keylist_data(self, input_output='input'):
        assert input_output in ['input', 'output']
        key_list = []
        with _h5py_File(self._filename, "r", self._h5file) as _h5file:
            # dirname = "bse/" + self._groupname + "/input"
            dirname = "bse/" + self._groupname + "/" + input_output
            if dirname in _h5file:
                grp = _h5file[dirname]
                # list(~.keys()) to get all keys once
                for _key in list(grp.keys()):  # 'X0_q', 'X_loc', etc
                    subgrp_key = grp[_key]
                    if(_key != "gamma0"):
                        for dir_w in list(subgrp_key.keys()):  # 'w0', 'w1', ...
                            _bosonic_freq = int(dir_w.lstrip('w'))  # 'w0' -> 0
                            subgrp_w = subgrp_key[dir_w]
                            keyinfo = self._get_keyinfo(_key)
                            if keyinfo[2] == "wq":
                                for dir_q in list(subgrp_w.keys()):  # 'q_00.00.00', ...
                                    _q = dir_q.lstrip('q_')  # 'q_00.00.00' -> '00.00.00'
                                    key_list.append((_key, _bosonic_freq, _q))
                            else:
                                key_list.append((_key, _bosonic_freq))
                    else:
                        key_list.append((_key,))
        return(key_list)

    def _get_info(self, key):
        keyinfo = self._get_keyinfo(key)
        dirname = "bse/" + self._groupname + "/" + keyinfo[0] + "/" + key[0]
        data = []
        with _h5py_File(self._filename, "r", self._h5file) as _h5file:
            if keyinfo[1] == "list":
                keys = list(_h5file[dirname].keys())
                for i in range(len(keys)):
                    childdir = dirname + "/" + str(i)
                    _data = _h5file[childdir][()]
                    data.append(_data)
            else:
                data = _h5file[dirname][()]
        return data

    # def _get_data(self, key):
    #     keyinfo = self._get_keyinfo(key)
    #     dirname = self._make_dirname(key, keyinfo)
    #     try:
    #         assert key[0] == _from_bytes_to_utf8(self._load_data(dirname + "/info/key")), "key is wrong."
    #         if len(keyinfo) != 2:
    #             if keyinfo[2] == "wq" or keyinfo[2] == "w":
    #                 assert key[1] == self._load_data(dirname + "/info/bosonic_freq"), "bosonic_freq is wrong."
    #             if keyinfo[2] == "wq":
    #                 assert key[2] == _from_bytes_to_utf8(self._load_data(dirname + "/info/q")), "q is wrong."
    #     except AssertionError as err:
    #         # print(f"Warning : {err}")
    #         return False
    #     data = {}
    #     with _h5py_File(self._filename, "r", self._h5file) as _h5file:
    #         if "data" not in _h5file[dirname]:
    #             raise RuntimeError(f"Invalid data structure: Group {dirname!r} in HDF5 file {self._filename!r} contains 'info' but no 'data'.")
    #             # return False
    #         keys = list(_h5file[dirname + "/data"].keys())
    #     for _key in keys:
    #         _innerlist = _key.split("_")
    #         inner = tuple([int(_innerlist[0]), int(_innerlist[1])])
    #         data[inner] = self._load_data(dirname + "/data/" + _key)
    #     return data

    def _get_data(self, key):
        keyinfo = self._get_keyinfo(key)
        dirname = self._make_dirname(key, keyinfo)
        data = {}

        # load data
        with _h5py_File(self._filename, "r", self._h5file) as _h5file:
            if dirname not in _h5file:
                return False
            grp = _h5file[dirname]
            for _key in list(grp.keys()):  # list() to get all keys once
                i, j = list(map(int, _key.split("_")))  # '0_1' -> ['0', '1']
                data[(i, j)] = grp[_key][()]  # load
        return data

    def make(self, key, path):
        pass

    def delete(self, key):
        pass

    def get_path(self, key):
        pass

    def _load_data(self, path):
        with _h5py_File(self._filename, "r", self._h5file) as _h5file:
            if path in _h5file:
                data = _h5file[path][()]
            else:
                # print(path)
                # print("data does not exist.")
                return False
            return data
