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

__author__ = "Alexander V. Dudchenko"


def solve_with_ma27(m, tee=True, **kwargs):
    result = sar.solve_model(m, tee=tee, linear_solver="ma27")
    return result


def main():
    cwd = get_working_dir()
    loopTool(
        cwd + "/stability_sweep.yaml",
        build_function=sar.build_model,
        initialize_function=sar.initialize,
        optimize_function=solve_with_ma27,
        save_name="stability_sweep",
        probe_function=sar.test_func,
        saving_dir=cwd,
        number_of_subprocesses=1,
        num_loop_workers=3,
    )


if __name__ == "__main__":
    main()
