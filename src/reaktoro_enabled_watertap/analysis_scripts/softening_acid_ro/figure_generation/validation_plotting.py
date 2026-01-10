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
import pandas as pd

import numpy as np


__author__ = "Alexander V. Dudchenko"


def get_csv_data_from_pd(file_path):
    df = pd.read_csv(file_path)
    return df


def get_specific_data(df, key, acid_type):
    data = df[key].to_numpy()
    idx = np.where(df["Acid choice"].to_numpy() == acid_type)
    return data[idx]


def add_legend_and_save(
    fig,
    save_name,
    acid_legend,
    work_legend,
    ylabel,
    yticks,
    xlabel="Water recovery (%)",
    xticks=[50, 60, 70, 80, 90],
):
    for label, marker in acid_legend.items():
        fig.plot_line(
            [],
            [],
            marker=marker,
            markersize=4,
            color="black",
            label=label,
        )
    for label, color in work_legend.items():
        fig.plot_line(
            [],
            [],
            marker="",
            markersize=4,
            color=color,
            label=label,
        )
    fig.set_axis(
        xticks=xticks,
        xlabel=xlabel,
        yticks=yticks,
        ylabel=ylabel,
    )
    fig.add_legend()
    fig.save_fig(f"{save_location}/{save_name}")
    # fig.show()
    return fig


def get_different(data_in, data_out):

    return (data_in.data - data_out.data) / data_in.data * 100


def get_val_data(dm, case, acid, data_key):
    if case == "stacked_data":
        x = dm[
            case,
            ("acidification_reagents", acid),
            "Water recovery",
        ].data
        y = dm[case, ("acidification_reagents", acid), data_key].data
    else:
        x = dm[
            ("acidification_reagents", acid),
            case,
            "Water recovery",
        ].data
        y = dm[
            ("acidification_reagents", acid),
            case,
            data_key,
        ].data
    return x, y


if __name__ == "__main__":
    work_path = get_lib_path()
    save_location = (
        work_path
        / "analysis_scripts/softening_acid_ro/figure_generation/validation_figures/"
    )
    data_manager = psDataManager(
        str(
            work_path
            / "analysis_scripts/softening_acid_ro/data_generation/output/validation_soda_ash_hcl_h2so4_sweep_analysisType_validation_sweep.h5"
        ),
    )

    seawater_df = get_csv_data_from_pd(
        str(
            work_path
            / "analysis_scripts/softening_acid_ro/figure_generation/validation_data/sw_onestage_medium_capacity_rightdensity.csv"
        )
    )

    bgw_df = get_csv_data_from_pd(
        str(
            work_path
            / "analysis_scripts/softening_acid_ro/figure_generation/validation_data/bw_onestage_medium_capacity_rightdensity.csv"
        )
    )

    import_keys = [
        {
            "filekey": "fs.water_recovery",
            "return_key": "Water recovery",
            "units": "%",
        },
        {
            "filekey": "fs.costing.LCOW",
            "return_key": "LCOW",
        },
        {
            "filekey": f"fs.softening_unit.precipitation_reactor.reagent_dose[Na2CO3]",
            "return_key": ("softening chemical dose", "Soda ash"),
            "units": "PPM",
        },
        {
            "filekey": f"fs.acidification_unit.chemical_reactor.reagent_dose[HCl]",
            "return_key": ("acid addition dose", "HCl"),
            "units": "PPM",
        },
        {
            "filekey": f"fs.acidification_unit.chemical_reactor.reagent_dose[H2SO4]",
            "return_key": ("acid addition dose", "H2SO4"),
            "units": "PPM",
        },
        {
            "filekey": f"fs.ro_unit.ro_unit.scaling_tendency[Calcite]",
            "return_key": ("Scaling tendency", "Calcite"),
            "units": "dimensionless",
        },
        {
            "filekey": f"fs.ro_unit.ro_unit.scaling_tendency[Gypsum]",
            "return_key": ("Scaling tendency", "Gypsum"),
            "units": "dimensionless",
        },
        {
            "filekey": f"fs.hpro_unit.ro_unit.scaling_tendency[Calcite]",
            "return_key": ("HP Scaling tendency", "Calcite"),
            "units": "dimensionless",
        },
        {
            "filekey": f"fs.hpro_unit.ro_unit.scaling_tendency[Gypsum]",
            "return_key": ("HP Scaling tendency", "Gypsum"),
            "units": "dimensionless",
        },
    ]
    data_manager.load_data(import_keys)
    data_manager.select_data(("water_sim_cases", "SW_RO"))
    data_manager.select_data(("water_sim_cases", "SW_HPRO"), add_to_existing=True)
    data_manager.display()
    sea_water = data_manager.get_selected_data()
    sea_water.reduce_data(
        stack_keys="water_sim_cases", data_key="LCOW", reduction_type="min"
    )
    plot_cases = {
        "Seawater": {
            "dm": sea_water,
            "case": "stacked_data",
            "val_data": seawater_df,
        },
        "BGW": {
            "dm": data_manager,
            "case": ("water_sim_cases", "BGW"),
            "val_data": bgw_df,
        },
    }
    acids = {"HCl": "o", "H2SO4": "s"}
    study = {"This work": "red", "Amusat et al. (2024)": "black"}
    for case, case_details in plot_cases.items():
        fig = figureGenerator()
        fig.init_figure(width=3, height=2)
        for acid, marker in acids.items():
            dm = case_details["dm"]
            vd = case_details["val_data"]
            x, y = get_val_data(dm, case_details["case"], acid, "LCOW")
            val_x = get_specific_data(vd, "System recovery", acid.upper()) * 100
            val_y = get_specific_data(vd, "LCOW", acid.upper())

            if case == "Seawater":
                yticks = [0, 0.5, 1.0, 1.5, 2.0]
            else:
                yticks = [0, 0.2, 0.4, 0.6, 0.8]
            fig.plot_line(
                x,
                y,
                marker=marker,
                markersize=4,
                color=study["This work"],
                zorder=10,
            )
            fig.plot_line(
                val_x,
                val_y,
                marker=marker,
                markersize=4,
                color=study["Amusat et al. (2024)"],
            )

        add_legend_and_save(
            fig,
            f"case_{case}_lcow_comparison",
            acid_legend=acids,
            work_legend=study,
            ylabel="LCOW ($\$$/$m^3$)",
            yticks=yticks,
        )
        fig = figureGenerator()
        fig.init_figure(width=3, height=2)
        for acid, marker in acids.items():
            dm = case_details["dm"]
            vd = case_details["val_data"]

            x, y = get_val_data(
                dm,
                case_details["case"],
                acid,
                ("softening chemical dose", "Soda ash"),
            )
            val_x = get_specific_data(vd, "System recovery", acid.upper()) * 100
            val_y = get_specific_data(vd, "Optimal soda ash dose", acid.upper())
            if case == "Seawater":
                yticks = [0, 200, 400, 600, 800, 1000]
            else:
                yticks = [0, 200, 400, 600]
            fig.plot_line(
                x,
                y,
                marker=marker,
                markersize=4,
                color=study["This work"],
                zorder=10,
            )
            fig.plot_line(
                val_x,
                val_y,
                marker=marker,
                markersize=4,
                color=study["Amusat et al. (2024)"],
            )

        add_legend_and_save(
            fig,
            f"case_{case}_soda_ash_comparison",
            acid_legend=acids,
            work_legend=study,
            ylabel="Soda ash dose (mg/L)",
            yticks=yticks,
        )
        fig = figureGenerator()
        fig.init_figure(width=3, height=2)
        for acid, marker in acids.items():
            dm = case_details["dm"]
            vd = case_details["val_data"]
            x, y = get_val_data(
                dm,
                case_details["case"],
                acid,
                ("acid addition dose", acid),
            )
            val_x = get_specific_data(vd, "System recovery", acid.upper()) * 100
            val_y = get_specific_data(vd, "Optimal acid dose", acid.upper())
            if acid == "H2SO4":
                val_y = (
                    val_y / 0.93
                )  # Amusat reports pure acid dose, while we report dose based actual purity
            if acid == "HCl":
                val_y = (
                    val_y / 0.3
                )  # Amusat reports pure acid dose, while we report dose based actual purity
            if case == "Seawater":
                yticks = [0, 20, 40, 60, 80, 100]
            else:
                yticks = [0, 50, 100, 150, 200, 250, 300]
            fig.plot_line(
                x,
                y,
                marker=marker,
                markersize=4,
                color=study["This work"],
                zorder=10,
            )
            fig.plot_line(
                val_x,
                val_y,
                marker=marker,
                markersize=4,
                color=study["Amusat et al. (2024)"],
            )

        add_legend_and_save(
            fig,
            f"case_{case}_acid_dose_comparison",
            acid_legend=acids,
            work_legend=study,
            ylabel="Acid dose (mg/L)",
            yticks=yticks,
        )
        fig = figureGenerator()
        fig.init_figure(width=3, height=2)
        for acid, marker in acids.items():
            dm = case_details["dm"]
            vd = case_details["val_data"]
            if case == "Seawater":

                x, y = get_val_data(
                    dm,
                    case_details["case"],
                    acid,
                    ("HP Scaling tendency", "Calcite"),
                )

                x2, y2 = get_val_data(
                    dm,
                    case_details["case"],
                    acid,
                    ("Scaling tendency", "Calcite"),
                )
                # HPRO only has 1 stage upto recovery of ~65 % or so, so we merge HPRO and RO data together to get SI
                y[y == 0] = y2[y == 0]
            else:
                x, y = get_val_data(
                    dm,
                    case_details["case"],
                    acid,
                    ("Scaling tendency", "Calcite"),
                )
            val_x = get_specific_data(vd, "System recovery", acid.upper()) * 100
            val_y = get_specific_data(vd, "SI CaCO3", acid.upper())
            yticks = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
            fig.plot_line(
                x,
                y,
                marker=marker,
                markersize=4,
                color=study["This work"],
                zorder=10,
            )
            fig.plot_line(
                val_x,
                val_y,
                marker=marker,
                markersize=4,
                color=study["Amusat et al. (2024)"],
            )

        add_legend_and_save(
            fig,
            f"case_{case}_scaling_tendency_calcite_comparison",
            acid_legend=acids,
            work_legend=study,
            ylabel="Scaling Tendency (Calcite)",
            yticks=yticks,
        )
        fig = figureGenerator()
        fig.init_figure(width=3, height=2)
        for acid, marker in acids.items():
            dm = case_details["dm"]
            vd = case_details["val_data"]
            if case == "Seawater":

                x, y = get_val_data(
                    dm,
                    case_details["case"],
                    acid,
                    ("HP Scaling tendency", "Gypsum"),
                )

                x2, y2 = get_val_data(
                    dm,
                    case_details["case"],
                    acid,
                    ("Scaling tendency", "Gypsum"),
                )
                # HPRO only has 1 stage upto recovery of ~65 % or so, so we merge HPRO and RO data together to get SI
                y[y == 0] = y2[y == 0]
            else:
                x, y = get_val_data(
                    dm,
                    case_details["case"],
                    acid,
                    ("Scaling tendency", "Gypsum"),
                )
            val_x = get_specific_data(vd, "System recovery", acid.upper()) * 100
            val_y = get_specific_data(vd, "SI Gypsum", acid.upper())
            yticks = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
            fig.plot_line(
                x,
                y,
                marker=marker,
                markersize=4,
                color=study["This work"],
                zorder=10,
            )
            fig.plot_line(
                val_x,
                val_y,
                marker=marker,
                markersize=4,
                color=study["Amusat et al. (2024)"],
            )

        add_legend_and_save(
            fig,
            f"case_{case}_scaling_tendency_gypsum_comparison",
            acid_legend=acids,
            work_legend=study,
            ylabel="Scaling Tendency (Gypsum)",
            yticks=yticks,
        )
