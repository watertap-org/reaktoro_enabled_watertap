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

from psPlotKit.data_manager.ps_data_manager import PsDataManager
from psPlotKit.data_plotter.ps_break_down_plotter import BreakDownPlotter
from reaktoro_enabled_watertap.utils.report_util import get_lib_path
from psPlotKit.data_manager.costing_packages.watertap_costing import (
    WaterTapCostingPackage,
)
from psPlotKit.data_manager.ps_costing import (
    PsCostingPackage,
    PsCostingGroup,
    PsCostingManager,
)

__author__ = "Alexander V. Dudchenko"


def main(show_figs=True):
    work_path = get_lib_path()
    save_location = (
        work_path / "analysis_scripts/softening_acid_ro/figure_generation/cost_figures/"
    )
    costing_data = PsDataManager(
        str(
            work_path
            / "analysis_scripts/softening_acid_ro/data_generation/output/treatment_lime_soda_ash_hcl_h2so4_sweep_analysisType_treatment_sweep.h5"
        ),
    )
    RO = PsCostingGroup("RO")
    RO.add_unit(
        unit_name="ro_unit.ro_unit",
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    RO.add_unit(
        unit_name="pump_unit.pump",
        capex_keys="capital_cost",
        flow_keys={"electricity": "control_volume.work"},
    )
    HPRO = PsCostingGroup("HPRO")
    HPRO.add_unit(
        unit_name="hpro_unit.ro_unit",
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    HPRO.add_unit(
        unit_name="hp_pump_unit.pump",
        capex_keys="capital_cost",
        flow_keys={"electricity": "control_volume.work"},
    )
    erd = PsCostingGroup("ERD")
    erd.add_unit(
        unit_name="erd_unit.ERD",
        capex_keys="capital_cost",
        flow_keys={"electricity": "control_volume.work"},
    )
    softening = PsCostingGroup("Softening")
    softening.add_unit(
        unit_name="precipitation_reactor",
        capex_keys="capital_cost",
        flow_keys={
            "Lime_cost": "flow_mass_reagent[CaO]",
            "Sodaash_cost": "flow_mass_reagent[Na2CO3]",
        },
    )
    acid_addition = PsCostingGroup("Acid addition")
    acid_addition.add_unit(
        unit_name="chemical_reactor",
        capex_keys="capital_cost",
        flow_keys={
            "Sulfuric_acid_cost": "flow_mass_reagent[H2SO4]",
            "Hydrochloric_acid_cost": "flow_mass_reagent[HCl]",
        },
    )

    costing_data.register_data_key(
        file_key="fs.water_recovery",
        return_key="Water recovery",
        units="%",
    )

    pkg = WaterTapCostingPackage(
        costing_block="fs.costing", validation_key="fs.costing.LCOW"
    )
    pkg.add_flow_cost(
        "Lime_cost",
        cost_parameter_key="fs_softening_unit_reagent_CaO_cost",
        parameter_assign_units="USD/kg",
    )
    pkg.add_flow_cost(
        "Sodaash_cost",
        cost_parameter_key="fs_softening_unit_reagent_Na2CO3_cost",
        parameter_assign_units="USD/kg",
    )
    pkg.add_flow_cost(
        "Sulfuric_acid_cost",
        cost_parameter_key="fs_acidification_unit_reagent_H2SO4_cost",
        parameter_assign_units="USD/kg",
    )
    pkg.add_flow_cost(
        "Hydrochloric_acid_cost",
        cost_parameter_key="fs_acidification_unit_reagent_HCl_cost",
        parameter_assign_units="USD/kg",
    )
    pkg.register_product_flow(
        file_key="fs.product.product.properties[0.0].flow_vol_phase[Liq]"
    )
    cm = PsCostingManager(
        costing_data, pkg, [RO, HPRO, erd, softening, acid_addition]
    )  # , acid_addition])
    cm.build(error_on_validation_failure=True)

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
        cost_plotter = BreakDownPlotter(
            wr,
            save_folder=save_location,
            save_name="{}".format(case[1]),
            show_fig=show_figs,
        )
        cost_plotter.define_area_groups(
            [
                {"ERD": {"label": None, "color": "#d9f0d3"}},
                {"RO": {"label": None, "color": "#a6cee3"}},
                {"HPRO": {"label": None, "color": "#1f78b4"}},
                {"Acid addition": {"label": None, "color": "#a7a7a7"}},
                {"Softening": {"label": None, "color": "#33a02c"}},
            ]
        )
        # cost_plotter.define_hatch_groups(
        #     {"LCOW": {}}
        # )  # {"CAPEX": {}, "OPEX": {"hatch": "//"}})
        cost_plotter.plotbreakdown(
            xdata="Water recovery",
            ydata="LCOW",
            axis_options={
                "yticks": [0, 0.5, 1, 1.5],
                "xticks": [50, 60, 70, 80, 90],
                "ylabel": "LCOW ($\$/m^3$)",
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


if __name__ == "__main__":
    main()
