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
        / "analysis_scripts/property_comparison/figure_generation/validation_figures/"
    )
    data_manager = psDataManager(
        str(
            work_path
            / "analysis_scripts/property_comparison/data_generation/output/prop_sweep_analysisType_prop_sweep.h5"
        ),
    )

    import_keys = [
        {
            "filekey": "fs.water_recovery",
            "return_key": "Water recovery",
            "units": "%",
        },
        {
            "filekey": "fs.feed_tds",
            "return_key": "TDS",
        },
        {
            "filekey": f"fs.modified_properties[osmoticPressure,H2O]",
            "return_key": ("Reaktoro", "osmotic pressure"),
            "assign_units": "Pa",
            "conversion_factor": 1e-5,
        },
        {
            "filekey": f"fs.seawater_feed.properties[0.0].pressure_osm_phase[Liq]",
            "return_key": ("Watertap Seawater property package", "osmotic pressure"),
            "units": "bar",
        },
        {
            "filekey": f"fs.nacl_feed.properties[0.0].pressure_osm_phase[Liq]",
            "return_key": ("Watertap NaCl property package", "osmotic pressure"),
            "units": "bar",
        },
        {
            "filekey": f"fs.multicomp_feed.properties[0.0].dens_mass_phase[Liq]",
            "return_key": ("Reaktoro", "density"),
        },
        {
            "filekey": f"fs.seawater_feed.properties[0.0].dens_mass_phase[Liq]",
            "return_key": ("Watertap Seawater property package", "density"),
        },
        {
            "filekey": f"fs.nacl_feed.properties[0.0].dens_mass_phase[Liq]",
            "return_key": ("Watertap NaCl property package", "density"),
        },
    ]
    data_manager.load_data(import_keys)
    data_manager.display()
    water_reference = {
        "SW": "seawater",
        "BGW": "USBR BGW",
        "BGW_500": "USGS case 1",
        "BGW_1500": "USGS case 2",
    }
    packages = {
        "Reaktoro": "#a6cee3",
        "Watertap Seawater property package": "#b2df8a",
        "Watertap NaCl property package": "#fb9a99",
    }
    for water in water_reference.keys():
        fig = figureGenerator()
        fig.init_figure()
        for case, color in packages.items():
            fig.plot_line(
                data_manager[("water_sim_cases", water), "Water recovery"].data,
                data_manager[
                    ("water_sim_cases", water), (case, "osmotic pressure")
                ].data,
                marker="o",
                markersize=4,
                color=color,
                label=case,
            )
        yticks = [0, 5, 10, 15, 20, 25, 30]
        if water == "SW":
            yticks = [0, 50, 100, 150, 200, 250, 300, 350, 400]
        fig.set_axis(
            xticks=[50, 60, 70, 80, 90],
            xlabel="Water recovery (%)",
            yticks=yticks,
            ylabel="Osmotic pressure (bar)",
        )
        fig.add_legend()
        fig.save_fig(f"figures/{water_reference[water]}_osmotic_pressure_comparison")
    for water in water_reference.keys():
        fig = figureGenerator()
        fig.init_figure()
        for case, color in packages.items():
            fig.plot_line(
                data_manager[("water_sim_cases", water), "Water recovery"].data,
                data_manager[("water_sim_cases", water), (case, "density")].data,
                marker="o",
                markersize=4,
                color=color,
                label=case,
            )
        yticks = [1000, 1005, 1010, 1015, 1020, 1025, 1030]
        if water == "SW":
            yticks = [
                1000,
                1050,
                1100,
                1150,
                1200,
                1250,
            ]
        fig.set_axis(
            xticks=[50, 60, 70, 80, 90],
            xlabel="Water recovery (%)",
            yticks=yticks,
            ylabel="Density (kg/m³)",
        )
        fig.add_legend()
        fig.save_fig(f"figures/{water_reference[water]}_density_comparison")
    fig.show()
