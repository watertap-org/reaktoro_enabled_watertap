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
from psPlotKit.data_plotter.ps_break_down_plotter import breakDownPlotter
from reaktoro_enabled_watertap.utils.report_util import get_lib_path

__author__ = "Alexander V. Dudchenko"


if __name__ == "__main__":
    work_path = get_lib_path()
    save_location = (
        work_path / "analysis_scripts/softening_acid_ro/figure_generation/cost_figures/"
    )
    costing_data = psDataManager(
        str(
            work_path
            / "analysis_scripts/softening_acid_ro/data_generation/output/treatment_lime_soda_ash_hcl_h2so4_sweep_analysisType_treatment_sweep.h5"
        ),
    )
    device_groups = {
        "ERD": {
            "units": {
                "ERD": "erd_unit",
            },
        },
        "RO": {
            "units": {
                "ro_unit": "ro_unit",
                "pump": "pump_unit",
            },
        },
        "HPRO": {
            "units": {
                "ro_unit": "hpro_unit",
                "pump": "hp_pump_unit",
            },
        },
        "Softening": {
            "CAPEX": {
                "units": "precipitation_reactor",
            },
            "OPEX": {
                "units": {
                    "fs.costing.aggregate_flow_costs[fs_softening_unit_reagent_CaO]",
                    "fs.costing.aggregate_flow_costs[fs_softening_unit_reagent_Na2CO3]",
                }
            },
        },
        "Acid addition": {
            "CAPEX": {
                "units": "chemical_reactor",
            },
            "OPEX": {
                "units": {
                    "fs.costing.aggregate_flow_costs[fs_acidification_unit_reagent_H2SO4]",
                    "fs.costing.aggregate_flow_costs[fs_acidification_unit_reagent_HCl]",
                }
            },
        },
    }

    costing_data.load_data(
        [
            {
                "filekey": "fs.water_recovery",
                "return_key": "Water recovery",
                "units": "%",
            },
        ],
    )
    costing_data.get_costing(
        device_groups,
        default_flow="fs.product.product.properties[0.0].flow_vol_phase[Liq]",
    )

    regions = {
        "seawater": {
            "Region 1": {"range": [50, 74], "color": "#fef0d9"},
            "Region 2": {"range": [74, 75], "color": "#fdd49e"},
            "Region 3": {"range": [75, 76], "color": "#fdbb84"},
            "Region 4": {"range": [76, 86], "color": "#fc8d59"},
        },
        ("water_sim_cases", "BGW"): {
            "Region 1": {"range": [50, 56], "color": "#fef0d9"},
            "Region 2": {"range": [56, 64], "color": "#fdd49e"},
            "Region 3": {"range": [64, 74], "color": "#fdbb84"},
            "Region 4": {"range": [74, 90], "color": "#fc8d59"},
        },
        ("water_sim_cases", f"BGW_500"): {
            "Region 1": {"range": [50, 59], "color": "#fef0d9"},
            "Region 2": {"range": [59, 64], "color": "#fdd49e"},
            "Region 3": {"range": [64, 68], "color": "#fdbb84"},
            "Region 4": {"range": [68, 90], "color": "#fc8d59"},
        },
        ("water_sim_cases", f"BGW_1500"): {
            "Region 4": {"range": [50, 90], "color": "#fc8d59"},
        },
    }
    for case in [
        ("water_sim_cases", "BGW"),
        ("water_sim_cases", f"BGW_500"),
        ("water_sim_cases", f"BGW_1500"),
    ]:
        costing_data.select_data(case, exact_keys=True)
        wr = costing_data.get_selected_data()
        cost_plotter = breakDownPlotter(
            wr,
            save_folder=save_location,
            save_name="{}".format(case[1]),
        )
        cost_plotter.define_area_groups(
            [
                {"ERD": {"label": None, "color": "#d9f0d3"}},
                {"RO": {"label": None, "color": "#a6cee3"}},
                {"HPRO": {"label": None, "color": "#1f78b4"}},
                {"Acid addition": {"label": None, "color": "#fb9a99"}},
                {"Softening": {"label": None, "color": "#33a02c"}},
            ]
        )
        cost_plotter.define_hatch_groups(
            {"TOTAL": {}}
        )  # {"CAPEX": {}, "OPEX": {"hatch": "//"}})
        cost_plotter.plotbreakdown(
            xdata="Water recovery",
            ydata=["cost_breakdown", "levelized"],
            axis_options={
                "yticks": [0, 0.5, 1, 1.5],
                "xticks": [50, 60, 70, 80, 90],
            },
            legend_loc="upper left",
            generate_figure=False,
            fig_options={"width": 1.8, "height": 1.6},
        )
        for r_name, r in regions[case].items():
            print(r_name, r["range"], r["color"], cost_plotter.fig)
            cost_plotter.fig.plot_area(
                r["range"],
                [0, 0],
                y2data=[1.5, 1.5],
                label=r_name,
                zorder=-10,
                color=r["color"],
                edgecolor=None,
            )
        cost_plotter.generate_figure()
