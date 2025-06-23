# coding: utf-8
import toml
import sys
import os
from collections import OrderedDict


def _check_if_dict_empty(dict_obj, block):
    if dict_obj:
        sys.exit(f"ERROR: Unrecognized parameter(s) in [{block}]: {tuple(dict_obj.keys())}")


def _obsolete_param(dict_obj, block, name):
    param = dict_obj[block].pop(name, None)
    if param is not None:
        print(f"Warning: Parameter [{block}]{name} is obsolete.", file=sys.stderr)


def load_params_from_toml(file_name, print_summary=True):
    if print_summary:
        print(f"Loading file {file_name!r}")

    try:
        dict_toml = toml.load(open(file_name))
    except Exception as e:
        print(f"ERROR in loading TOML format file: {file_name!r}", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)

    dict_common = OrderedDict()
    dict_main = OrderedDict()
    dict_post = OrderedDict()
    dict_anacont = OrderedDict()

    #set common parameters
    params = dict_toml["chiq_common"]
    dict_common["input"] = params.pop("input", "dmft_bse.h5")
    dict_common["output"] = params.pop("output", "dmft_bse.out.h5")
    dict_common["type"] = params.pop("type", ["bse"])
    if "omega_q" in params:
        print("Error: chiq_common.omega_q is removed. Use chiq_common.q_points instead.")
        sys.exit(1)
    dict_common["q_points"] = params.pop("q_points", None)
    dict_common["num_wb"] = params.pop("num_wb", -1)
    _check_if_dict_empty(params, block="chiq_common")

    #set main parameters
    params = dict_toml["chiq_main"]
    dict_main["work_dir"] = params.pop("work_dir", "")
    dict_main["num_wf"] = params.pop("num_wf", None)
    _obsolete_param(dict_toml, "chiq_main", "solver")
    _check_if_dict_empty(params, block="chiq_main")

    #set post parameters
    if "chiq_post" not in dict_toml:
        dict_toml["chiq_post"] = {}
    params = dict_toml["chiq_post"]
    dict_post["output_dir"] = params.pop("output_dir", "")
    dict_post["mode"] = params.pop("mode", "eigen")
    dict_post["vector"] = params.pop("vector", False)
    dict_post["order"] = params.pop("order", "descend")
    dict_post["order_file"] = params.pop("order_file", "eigenvec.in")
    dict_post["coefs_file"] = params.pop("coefs_file", "coefs.in")
    _check_if_dict_empty(params, block="chiq_post")

    #set anacont parameters
    params = dict_toml["chiq_anacont"]
    dict_anacont["input_file"] = params.pop("input_file", "chi_q_eigen.dat")
    dict_anacont["output_dir"] = params.pop("output_dir", os.path.join(dict_post["output_dir"], "anacont"))
    dict_anacont["output_prefix"] = params.pop("output_prefix", "chi_q_w")
    dict_anacont["output_prefix_chi_q_iw"] = params.pop("output_prefix_chi_q_iw", "chi_q_iw")
    dict_anacont["wmax"] = params.pop("wmax", 10.0)
    dict_anacont["wmin"] = params.pop("wmin", 0.0)
    dict_anacont["wnum"] = params.pop("wnum", 101)
    dict_anacont["method"] = params.pop("method", "pade")
    dict_anacont["print_chi_q_iw"] = params.pop("print_chi_q_iw", False)
    # for pade
    dict_anacont["eta"] = params.pop("eta", 1e-5)
    # for SpM
    dict_anacont["loglambda"] = params.pop("loglambda", 0.0)
    dict_anacont["maxiter"] = params.pop("maxiter", 1000)
    dict_anacont["initial_mu"] = params.pop("initial_mu", 1.0)
    dict_anacont["wmax_factor"] = params.pop("wmax_factor", 1.2)
    _check_if_dict_empty(params, block="chiq_anacont")

    # Print summary
    if print_summary:
        print("=========================================")
        print("Summary of parameters")
        for block, params in [("chiq_common", dict_common), ("chiq_main", dict_main), ("chiq_post", dict_post), ("chiq_anacont", dict_anacont)]:
            print(f"\n[{block}]")
            for key, val in params.items():
                print(f"{key} = {val!r}")
        print("=========================================")

    result_dict = OrderedDict()
    result_dict["common"] = dict_common
    result_dict["main"] = dict_main
    result_dict["post"] = dict_post
    result_dict["anacont"] = dict_anacont

    return result_dict


if __name__ == "__main__":
    with open("test.toml", "w") as fw:
        fw.write("[chiq_common]\n")
        fw.write("input = \"test.in\"\n")
        fw.write("[chiq_main]\n")
        fw.write("num_wf = 10\n")
        fw.write("[chiq_post]\n")
        fw.write("output_dir = \"output/\"\n")
    dict_objs = load_params_from_toml("test.toml")
    print(dict_objs)
