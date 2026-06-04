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

from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
import reaktoro_enabled_watertap.flowsheets.softening_acid_ro.softening_acid_ro as sar
import time

from reaktoro_enabled_watertap.utils.report_util import get_lib_path

import yaml
import os

__author__ = "Alexander V. Dudchenko"


def solve_with_ma27(m, tee=False, **kwargs):
    result = sar.solve_model(m, tee=tee, linear_solver="ma27")
    return result


def initialize_ma27(m, **kwargs):
    sar.initialize(m, linear_solver="ma27", tee=False)


def update_config(config_location, config_name, new_config_name, num_samples=11):

    with open(os.path.join(config_location, config_name), "r") as f:
        config = yaml.safe_load(f)

    config["stability_sweep"]["build_loop"]["build_loop"]["sweep_param_loop"][
        "water_recovery"
    ]["num_samples"] = num_samples
    with open(os.path.join(config_location, new_config_name), "w") as f:
        yaml.dump(config, f, sort_keys=False)


def main(
    save_location=None,
    config_location=None,
    num_loop_workers=3,
    use_ma27=True,
    num_samples=11,
):
    ts = time.time()
    work_path = get_lib_path()
    work_path = str(work_path) + "/analysis_scripts/softening_acid_ro/data_generation"
    if save_location is None:
        save_location = work_path
    if config_location is None:
        config_location = work_path
    solver = solve_with_ma27 if use_ma27 else sar.solve_model
    initializer = initialize_ma27 if use_ma27 else sar.initialize

    update_config(
        config_location,
        config_location + "/stability_sweep.yaml",
        config_location + "/_temp_stability_sweep.yaml",
        num_samples=num_samples,
    )

    loopTool(
        config_location + "/_temp_stability_sweep.yaml",
        build_function=sar.build_model,
        initialize_function=initializer,
        optimize_function=solver,
        save_name="stability_sweep",
        probe_function=sar.test_func,
        saving_dir=save_location,
        number_of_subprocesses=1,
        num_loop_workers=num_loop_workers,
    )

    print("Total time: ", time.time() - ts)


if __name__ == "__main__":
    main()
