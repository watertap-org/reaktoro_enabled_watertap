import pytest

from reaktoro_enabled_watertap.analysis_scripts.softening_acid_ro.data_generation import (
    treatment_sweep,
    validation_sweep,
    stability_sweep,
)
from reaktoro_enabled_watertap.utils.report_util import get_lib_path

from psPlotKit.data_manager.ps_data_manager import psDataManager

import os
import numpy


@pytest.mark.component
def test_stability_sweep():
    work_path = get_lib_path()
    filename = str(
        work_path
        / "analysis_scripts/softening_acid_ro/data_generation/output/stability_sweep_analysisType_stability_sweep.h5"
    )
    if os.path.exists(filename):
        os.remove(filename)

    stability_sweep.main()

    data_manager = psDataManager(
        filename,
    )
    data_manager.register_data_key(
        file_key="fs.water_recovery", return_key="Water recovery", units="%"
    )
    data_manager.register_data_key(file_key="fs.costing.LCOW", return_key="LCOW")
    data_manager.load_data()
    data_manager.display()
    test_data = {
        (
            ("hessian_sim_cases", "bfgs_damp_gn"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "bfgs_damp_gn"), ("water_sim_cases", "BGW"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_gn"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_gn"),
            ("water_sim_cases", "BGW_1500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_gn"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_gn"),
            ("water_sim_cases", "BGW_500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_gn"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_gn"),
            ("water_sim_cases", "SW_HPRO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_gn"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_gn"),
            ("water_sim_cases", "SW_RO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_sc1"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_sc1"),
            ("water_sim_cases", "BGW"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_sc1"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_sc1"),
            ("water_sim_cases", "BGW_500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_sc1"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_sc1"),
            ("water_sim_cases", "SW_HPRO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_sc1"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_damp_sc1"),
            ("water_sim_cases", "SW_RO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_gn"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "bfgs_gn"), ("water_sim_cases", "BGW"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "bfgs_gn"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "bfgs_gn"), ("water_sim_cases", "BGW_1500"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "bfgs_gn"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "bfgs_gn"), ("water_sim_cases", "BGW_500"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "bfgs_gn"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "bfgs_gn"), ("water_sim_cases", "SW_HPRO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "bfgs_gn"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "bfgs_gn"), ("water_sim_cases", "SW_RO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_gn"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_gn"),
            ("water_sim_cases", "BGW"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_gn"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_gn"),
            ("water_sim_cases", "BGW_1500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_gn"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_gn"),
            ("water_sim_cases", "BGW_500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_gn"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_gn"),
            ("water_sim_cases", "SW_HPRO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_gn"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_gn"),
            ("water_sim_cases", "SW_RO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_sc1"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_sc1"),
            ("water_sim_cases", "BGW"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_sc1"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_sc1"),
            ("water_sim_cases", "BGW_500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_sc1"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_sc1"),
            ("water_sim_cases", "SW_HPRO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_sc1"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_ipopt_sc1"),
            ("water_sim_cases", "SW_RO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_sc1"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "bfgs_sc1"), ("water_sim_cases", "BGW"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "bfgs_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "bfgs_sc1"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "bfgs_sc1"), ("water_sim_cases", "BGW_500"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "bfgs_sc1"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "bfgs_sc1"), ("water_sim_cases", "SW_HPRO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "bfgs_sc1"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "bfgs_sc1"), ("water_sim_cases", "SW_RO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "cbfgs_gn"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "cbfgs_gn"), ("water_sim_cases", "BGW"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "cbfgs_gn"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "cbfgs_gn"),
            ("water_sim_cases", "BGW_1500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "cbfgs_gn"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "cbfgs_gn"), ("water_sim_cases", "BGW_500"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "cbfgs_gn"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "cbfgs_gn"), ("water_sim_cases", "SW_HPRO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "cbfgs_gn"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "cbfgs_gn"), ("water_sim_cases", "SW_RO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "cbfgs_sc1"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "cbfgs_sc1"), ("water_sim_cases", "BGW"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "cbfgs_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "cbfgs_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "cbfgs_sc1"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "cbfgs_sc1"),
            ("water_sim_cases", "BGW_500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "cbfgs_sc1"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "cbfgs_sc1"),
            ("water_sim_cases", "SW_HPRO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "cbfgs_sc1"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "cbfgs_sc1"), ("water_sim_cases", "SW_RO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "gauss_newton"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "gauss_newton"), ("water_sim_cases", "BGW"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "gauss_newton"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "gauss_newton"),
            ("water_sim_cases", "BGW_1500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "gauss_newton"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "gauss_newton"),
            ("water_sim_cases", "BGW_500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "gauss_newton"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "gauss_newton"),
            ("water_sim_cases", "SW_HPRO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "gauss_newton"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "gauss_newton"),
            ("water_sim_cases", "SW_RO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "lbfgs_gn"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lbfgs_gn"), ("water_sim_cases", "BGW"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "lbfgs_gn"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "lbfgs_gn"),
            ("water_sim_cases", "BGW_1500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "lbfgs_gn"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lbfgs_gn"), ("water_sim_cases", "BGW_500"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "lbfgs_gn"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lbfgs_gn"), ("water_sim_cases", "SW_HPRO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "lbfgs_gn"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lbfgs_gn"), ("water_sim_cases", "SW_RO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "lbfgs_sc1"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lbfgs_sc1"), ("water_sim_cases", "BGW"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "lbfgs_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "lbfgs_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "lbfgs_sc1"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "lbfgs_sc1"),
            ("water_sim_cases", "BGW_500"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "lbfgs_sc1"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (
            ("hessian_sim_cases", "lbfgs_sc1"),
            ("water_sim_cases", "SW_HPRO"),
            "LCOW",
        ): 11,
        (
            ("hessian_sim_cases", "lbfgs_sc1"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lbfgs_sc1"), ("water_sim_cases", "SW_RO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "lmt_sc1"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lmt_sc1"), ("water_sim_cases", "BGW"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "lmt_sc1"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lmt_sc1"), ("water_sim_cases", "BGW_1500"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "lmt_sc1"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lmt_sc1"), ("water_sim_cases", "BGW_500"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "lmt_sc1"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lmt_sc1"), ("water_sim_cases", "SW_HPRO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "lmt_sc1"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "lmt_sc1"), ("water_sim_cases", "SW_RO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "zero_hs"),
            ("water_sim_cases", "BGW"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "zero_hs"), ("water_sim_cases", "BGW"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "zero_hs"),
            ("water_sim_cases", "BGW_1500"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "zero_hs"), ("water_sim_cases", "BGW_1500"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "zero_hs"),
            ("water_sim_cases", "BGW_500"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "zero_hs"), ("water_sim_cases", "BGW_500"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "zero_hs"),
            ("water_sim_cases", "SW_HPRO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "zero_hs"), ("water_sim_cases", "SW_HPRO"), "LCOW"): 11,
        (
            ("hessian_sim_cases", "zero_hs"),
            ("water_sim_cases", "SW_RO"),
            "Water recovery",
        ): 11,
        (("hessian_sim_cases", "zero_hs"), ("water_sim_cases", "SW_RO"), "LCOW"): 11,
    }
    for key in data_manager:
        # test_data[key] = len(data_manager[key].data)
        # solved_data = data_manager[key].data
        assert (
            len(data_manager[key].data) == test_data[key]
        )  # just test all data is there
        # for i, val in enumerate(test_data[key]):
        #     if val is None:
        #         assert solved_data[i] != solved_data[i]
        #     else:
        #         assert val == pytest.approx(solved_data[i], rel=1e-5)
    # print(test_data)
