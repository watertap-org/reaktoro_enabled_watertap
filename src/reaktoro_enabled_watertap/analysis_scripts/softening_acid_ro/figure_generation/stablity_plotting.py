#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://https://github.com/watertap-org/reaktoro_enabled_watertap"
#################################################################################

from psPlotKit.data_manager.ps_data_manager import psDataManager
from psPlotKit.data_plotter.fig_generator import figureGenerator

from reaktoro_enabled_watertap.utils.report_util import get_lib_path

import numpy as np

__author__ = "Alexander V. Dudchenko"


def merge_data(manager, step, hessian, water_sources, key):
    data = None
    for w in water_sources:
        d = manager[
            ("rkt_numerical_step", step),
            ("rkt_hessian_type", hessian),
            ("sim_cases", w),
            key,
        ]
        if data is None:
            data = d.data
        else:
            data = np.hstack((data, d.data))
    return data


def get_data(manager, step, hessian, water_sources, key):
    data = merge_data(manager, step, hessian, water_sources, key)
    # print(data)
    return data


if __name__ == "__main__":
    work_path = get_lib_path()
    save_location = (
        work_path
        / "analysis_scripts/softening_acid_ro/figure_generation/stability_figures/"
    )
    data_manager = psDataManager(
        str(
            work_path
            / "analysis_scripts/softening_acid_ro/data_generation/output/stability_sweep_analysisType_stability_sweep.h5"
        ),
    )

    import_keys = [
        # {
        #     "filekey": "fs.costing.LCOW",
        #     "return_key": "LCOW",
        # },
        {
            "filekey": f"fs.ipopt_iterations[Objective function evaluations]",
            "return_key": "iterations",
            "units": "dimensionless",
        },
        {
            "filekey": f"fs.unscaled_ipopt_result[dual_infeasibility]",
            "return_key": "Dual infeasibility",
            "units": "dimensionless",
        },
        {
            "filekey": f"fs.unscaled_ipopt_result[overall_nlp_error]",
            "return_key": "Overall NLP error",
            "units": "dimensionless",
        },
    ]
    data_manager.load_data(import_keys)
    data_manager.display()
    # assert False
    hessians = [
        "zero_hs",
        "gauss_newton",
        "lmt_sc1",
        "bfgs",
        "cbfgs",
        "bfgs_ipopt",
        "bfgs_damp",
        "lbfgs",
    ]
    activation = {"bfgs": ["gn", "sc1"]}  # , "lmt": ["sc1"]}
    scalar_names = {
        "gn": "Gauss-Newton",
        "sc1": "scalar1",
        "sc2": "Scaled 2",
        "cn": "Diagonal",
    }
    sc_idx = {sc: i for i, sc in enumerate(scalar_names.values())}
    keys = {}
    names = {
        "zero_hs": "No hessian",
        "gauss_newton": "Gauss-Newton",
        "lmt_sc1": "Ipopt w/ limited-memory",
        "bfgs": "BFGS",
        "cbfgs": " BFGS-cautious ",
        "bfgs_ipopt": "BFGS-ipopt",
        "bfgs_damp": "BFGS-damped",
        "lbfgs": "L-BFGS",
    }

    for h in hessians:
        keys[names[h]] = {}
        act_found = False
        for act in activation:
            if act in h:
                for a in activation[act]:
                    key = f"{h}_{a}"
                    keys[names[h]][a] = key
                    act_found = True
        if not act_found:
            keys[names[h]] = h
    print(keys, sc_idx)
    water_sources = [
        "BGW",
        "BGW_1500",
        "BGW_500",
        "SW_RO",
        "SW_HPRO",
    ]  # ["BGW1", "sample_500_hardness", "sample_1500_hardness"]

    def agg_water_case_data(dtype, hess_type):
        data = []
        for w in water_sources:
            d = data_manager[
                ("hessian_sim_cases", hess_type), ("water_sim_cases", w), dtype
            ].data
            data += list(d)
        return np.array(data)

    data_results = {}
    for keyt in ["iterations", "Overall NLP error", "Dual infeasibility"]:
        for key, info in keys.items():
            if isinstance(info, dict):
                for k, v in info.items():
                    data = agg_water_case_data(keyt, v)
                    if v not in data_results:
                        data_results[v] = {}
                    data_results[v][keyt] = data
            else:
                data = agg_water_case_data(keyt, info)
                if info not in data_results:
                    data_results[info] = {}
                data_results[info][keyt] = data
    # print(data_results)
    colors = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99"]
    figures = {
        "iterations": {
            "ylabel": "Objective function evaluations (#)",
            "yticks": [0, 100, 200, 300],
            "yscale": "linear",
        },
        "Overall NLP error": {
            "ylabel": "log$_{10}$(Overall NLP error)",
            "yticks": [-8, -6, -4, -2, 0],
            "yscale": "log",
        },
        "Dual infeasibility": {
            "ylabel": "log$_{10}$(Dual infeasibility)",
            "yticks": [-8, -6, -4, -2, 0],
            "yscale": "log",
        },
        "Solve rate": {
            "ylabel": "Successful solves (%)",
            "yticks": [0, 25, 50, 75, 100],
            "yscale": "linear",
        },
        # "Objective": {
        #     "ylabel": "Objective difference from 10$^{-4}$",
        #     "yticks": [0, 10, 20],
        #     "yscale": "linear",
        # },
    }

    def process_data(
        dtype,
        data,
    ):
        if dtype == "Solve rate":
            data = data["iterations"]
            print(
                "Solve rate",
                len(data[data == data]) / (len(data) - 8) * 100,
                "num samples",
                (len(data) - 8),
            )
            return len(data[data == data]) / (len(data) - 8) * 100
        elif dtype in ["iterations", "Overall NLP error", "Dual infeasibility"]:
            data = data[dtype]
            # print("data", data)
            return data[data == data]

    def plot_result(dtype, pos, data, color, label=None, w=0.4):
        if dtype == "Solve rate":
            if data is not None:
                fig.plot_bar(
                    pos,
                    data,
                    width=w,
                    color=color,
                    label=label,
                )
        else:
            if len(data[data == data]) > 0:
                fig.plot_box(
                    pos,
                    data[data == data],
                    width=w,
                    color=color,
                    vertical=True,
                    label=label,
                )

    for fig_type, fig_settings in figures.items():
        fig = figureGenerator()
        fig.init_figure()
        created_label = []
        for i, hess in enumerate(keys):
            if isinstance(keys[hess], dict):
                for j, (key, v) in enumerate(keys[hess].items()):
                    pos = i + j * 0.4 + 0.2 - 0.4
                    data = process_data(fig_type, data_results[v])
                    if fig_settings["yscale"] == "log":
                        data = np.log10(data)
                        scale = "linear"
                    elif fig_settings["yscale"] == "mlog":
                        scale = "log"
                    else:
                        scale = fig_settings["yscale"]
                    label = scalar_names[key]
                    if label not in created_label:
                        created_label.append(label)
                    else:
                        label = None
                    plot_result(
                        fig_type,
                        pos,
                        data,
                        colors[sc_idx[scalar_names[key]]],
                        label=label,
                    )

            else:
                pos = i  # + 0.2
                data = process_data(fig_type, data_results[keys[hess]])
                # print("data", fig_type, data)
                if fig_settings["yscale"] == "log":
                    data = np.log10(data)
                color = "gray"
                if hess == "Ipopt w/ limited-memory":
                    color = "#1f78b4"
                plot_result(
                    fig_type,
                    pos,
                    data,
                    color,
                    w=0.6,
                    # label=scalar_names[key],
                )
        fig.set_axis_ticklabels(xticklabels=[i for i, n in keys.items()], rotate=True)
        fig.set_axis(
            yticks=fig_settings["yticks"],
            ylabel=fig_settings["ylabel"],
            yscale=scale,
            # yscale=fig_settings["yscale"],
        )
        fig.add_legend()
        fig.save_fig(str(save_location) + "\\" + fig_type)
    fig.show()
