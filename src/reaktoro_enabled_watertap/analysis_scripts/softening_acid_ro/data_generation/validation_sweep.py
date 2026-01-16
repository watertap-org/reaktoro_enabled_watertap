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

__author__ = "Alexander V. Dudchenko"


def solve_with_ma27(m, tee=False, **kwargs):
    result = sar.solve_model(m, tee=tee, linear_solver="mumps")
    return result


def initialize_ma27(m, **kwargs):
    sar.initialize(m, linear_solver="mumps", tee=False)


def main():

    ts = time.time()
    cwd = get_working_dir()
    loopTool(
        cwd + "/validation_soda_ash_hcl_h2so4_sweep.yaml",
        build_function=sar.build_model,
        initialize_function=initialize_ma27,
        optimize_function=solve_with_ma27,
        save_name="validation_soda_ash_hcl_h2so4_sweep",
        probe_function=sar.test_func,
        saving_dir=cwd,
        number_of_subprocesses=1,
        num_loop_workers=1,
    )

    print("Total time: ", time.time() - ts)


if __name__ == "__main__":
    main()
