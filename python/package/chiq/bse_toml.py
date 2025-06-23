# coding: utf-8
import toml
import sys
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

    #set tool parameters
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

    # Print summary
    if print_summary:
        print("=========================================")
        print("Summary of parameters")
        for block, params in [("chiq_common", dict_common), ("chiq_main", dict_main), ("chiq_post", dict_post)]:
            print(f"\n[{block}]")
            for key, val in params.items():
                print(f"{key} = {val!r}")
        print("=========================================")

    result_dict = OrderedDict()
    result_dict["common"] = dict_common
    result_dict["main"] = dict_main
    result_dict["post"] = dict_post

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
