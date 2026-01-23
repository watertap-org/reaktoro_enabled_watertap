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

__author__ = "Alexander V. Dudchenko"


def plot_fig(
    xdata,
    ydata,
    ylabel,
    yticks,
    xlabel="Water recovery (%)",
    xticks=[50, 60, 70, 80, 90],
    save_name="blank",
    fig=None,
    save_fig=True,
    data_label=None,
    width=3.6,
    height=3.6,
    regions=None,
    add_region_lbl=False,
    marker="o",
    color="black",
):
    if fig is None:
        fig = figureGenerator()
        fig.init_figure(width=width, height=height)
    if hasattr(ydata, "data"):
        ydata = ydata.data
    fig.plot_line(
        xdata.data,
        ydata,
        marker=marker,
        markersize=4,
        color=color,
        label=data_label,
    )
    if regions is not None:
        for r_name, r in regions.items():
            if add_region_lbl:
                lb = r_name
            else:
                lb = None
            fig.plot_area(
                r["range"],
                [yticks[0], yticks[0]],
                y2data=[yticks[-1], yticks[-1]],
                label=lb,
                zorder=-10,
                color=r["color"],
                edgecolor=None,
            )
    if save_fig:
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


def main(show_figs=True):
    water_reference = {
        "seawater": "seawater",
        "BGW": "USBR BGW",
        "BGW_500": "USGS case 1",
        "BGW_1500": "USGS case 2",
    }
    work_path = get_lib_path()
    global save_location
    save_location = (
        work_path
        / "analysis_scripts/softening_acid_ro/figure_generation/treatment_figures/"
    )
    data_manager = psDataManager(
        str(
            work_path
            / "analysis_scripts/softening_acid_ro/data_generation/output/treatment_lime_soda_ash_hcl_h2so4_sweep_analysisType_treatment_sweep.h5"
        ),
    )

    data_manager.register_data_key(
        file_key="fs.water_recovery",
        return_key="Water recovery",
        units="%",
    )
    data_manager.register_data_key(
        file_key="fs.costing.LCOW",
        return_key="LCOW",
    )
    data_manager.register_data_key(
        file_key="fs.hpro_unit.ro_unit.flux_mass_phase_comp_avg[0.0,Liq,H2O]",
        return_key=("HPRO1", "flux"),
    )
    data_manager.register_data_key(
        file_key="fs.hpro_unit.ro_unit.flux_mass_phase_comp[0.0,0.1,Liq,H2O]",
        return_key=("HPRO1", "flux inlet"),
    )
    data_manager.register_data_key(
        file_key="fs.hpro_unit.ro_unit.flux_mass_phase_comp[0.0,1.0,Liq,H2O]",
        return_key=("HPRO1", "flux outlet"),
    )
    data_manager.register_data_key(
        file_key="fs.hpro_unit.ro_unit.area",
        return_key=("HPRO1", "area"),
        units="m**2",
    )
    data_manager.register_data_key(
        file_key="fs.ro_unit.ro_retentate.pH",
        return_key=("RO1", "pH"),
        units="dimensionless",
    )
    data_manager.register_data_key(
        file_key="fs.hpro_unit.ro_retentate.pH",
        return_key=("HPRO1", "pH"),
        units="dimensionless",
    )
    data_manager.register_data_key(
        file_key="fs.softening_unit.precipitation_reactor.pH[outlet]",
        return_key=("softening", "pH"),
        units="dimensionless",
    )
    data_manager.register_data_key(
        file_key="fs.hpro_unit.ro_unit.feed_side.properties[0.0,0.0].pressure",
        return_key=("HPRO1", "pressure"),
        units="bar",
    )
    data_manager.register_data_key(
        file_key="fs.hpro_unit.ro_unit.feed_side.velocity[0.0,0.0]",
        return_key=("HPRO1", "velocity"),
        units="m/s",
    )
    data_manager.register_data_key(
        file_key="fs.ro_unit.ro_unit.flux_mass_phase_comp_avg[0.0,Liq,H2O]",
        return_key=("RO1", "flux"),
    )
    data_manager.register_data_key(
        file_key="fs.ro_unit.ro_unit.flux_mass_phase_comp[0.0,0.1,Liq,H2O]",
        return_key=("RO1", "flux inlet"),
    )
    data_manager.register_data_key(
        file_key="fs.ro_unit.ro_unit.flux_mass_phase_comp[0.0,1.0,Liq,H2O]",
        return_key=("RO1", "flux outlet"),
    )
    data_manager.register_data_key(
        file_key="fs.ro_unit.ro_unit.area",
        return_key=("RO1", "area"),
        units="m**2",
    )
    data_manager.register_data_key(
        file_key="fs.ro_unit.ro_unit.feed_side.properties[0.0,0.0].pressure",
        return_key=("RO1", "pressure"),
        units="bar",
    )
    data_manager.register_data_key(
        file_key="fs.ro_unit.ro_unit.feed_side.velocity[0.0,0.0]",
        return_key=("RO1", "velocity"),
        units="m/s",
    )
    data_manager.register_data_key(
        file_key="fs.softening_unit.precipitation_reactor.reagent_dose[CaO]",
        return_key=("Softening_RKT_1", "Chemical dose", "Lime"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.softening_unit.precipitation_reactor.reagent_dose[Na2CO3]",
        return_key=("Softening_RKT_1", "Chemical dose", "Soda ash"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.softening_unit.precipitation_reactor.dissolution_reactor.properties_in[0.0].conc_mass_phase_comp[Liq,Ca_2+]",
        return_key=("Softening_RKT_1", "Effluent in", "Ca"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.softening_unit.precipitation_reactor.dissolution_reactor.properties_in[0.0].conc_mass_phase_comp[Liq,HCO3_-]",
        return_key=("Softening_RKT_1", "Effluent in", "HCO3"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.softening_unit.precipitation_reactor.dissolution_reactor.properties_in[0.0].conc_mass_phase_comp[Liq,Mg_2+]",
        return_key=("Softening_RKT_1", "Effluent in", "Mg"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.softening_unit.precipitation_reactor.precipitation_reactor.properties_out[0.0].conc_mass_phase_comp[Liq,Ca_2+]",
        return_key=("Softening_RKT_1", "Effluent out", "Ca"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.softening_unit.precipitation_reactor.precipitation_reactor.properties_out[0.0].conc_mass_phase_comp[Liq,HCO3_-]",
        return_key=("Softening_RKT_1", "Effluent out", "HCO3"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.softening_unit.precipitation_reactor.precipitation_reactor.properties_out[0.0].conc_mass_phase_comp[Liq,Mg_2+]",
        return_key=("Softening_RKT_1", "Effluent out", "Mg"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.softening_unit.precipitation_reactor.alkalinity",
        return_key=("Softening_RKT_1", "Alkalinity"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.acidification_unit.chemical_reactor.reagent_dose[HCl]",
        return_key=("acid_addition", "Chemical dose", "HCl"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.acidification_unit.chemical_reactor.reagent_dose[H2SO4]",
        return_key=("acid_addition", "Chemical dose", "H2SO4"),
        units="PPM",
    )
    data_manager.register_data_key(
        file_key="fs.ro_unit.ro_unit.scaling_tendency[Calcite]",
        return_key=("Scaling tendency", "Calcite"),
        units="dimensionless",
    )
    data_manager.register_data_key(
        file_key="fs.ro_unit.ro_unit.scaling_tendency[Gypsum]",
        return_key=("Scaling tendency", "Gypsum"),
        units="dimensionless",
    )
    data_manager.register_data_key(
        file_key="fs.hpro_unit.ro_unit.scaling_tendency[Calcite]",
        return_key=("HP Scaling tendency", "Calcite"),
        units="dimensionless",
    )
    data_manager.register_data_key(
        file_key="fs.hpro_unit.ro_unit.scaling_tendency[Gypsum]",
        return_key=("HP Scaling tendency", "Gypsum"),
        units="dimensionless",
    )
    data_manager.load_data()
    data_manager.display()
    data_manager.select_data(("water_sim_cases", "SW_RO"))
    data_manager.select_data(("water_sim_cases", "SW_HPRO"), add_to_existing=True)
    sea_water = data_manager.get_selected_data()
    sea_water.display()
    sea_water.reduce_data(
        stack_keys="water_sim_cases", data_key="LCOW", reduction_type="min"
    )
    sea_water.display()
    # need to merge stage 1 and stage 2 data to single data set for HPRO
    for t in [
        "area",
        "pressure",
        "velocity",
        "flux",
        "pH",
        "flux inlet",
        "flux outlet",
    ]:

        hpro_area = sea_water[("water_sim_cases", "SW_HPRO"), ("HPRO1", t)].data
        ro_area = sea_water[("water_sim_cases", "SW_RO"), ("RO1", t)].data
        ssro_area = sea_water[("water_sim_cases", "SW_HPRO"), ("RO1", t)].data
        hpro_area[ro_area == ro_area] = 0
        sea_water[("water_sim_cases", "SW_HPRO"), ("HPRO1", t)].set_data(hpro_area)
        ro_area[ro_area != ro_area] = ssro_area[ro_area != ro_area]
        # print(ro_area, sea_water[("water_sim_cases", "SW_RO"), ("RO1", t)])
        sea_water[("water_sim_cases", "SW_RO"), ("RO1", t)].set_data(ro_area)
    for scalant in ["Gypsum", "Calcite"]:
        ro_data = sea_water[
            ("water_sim_cases", "SW_RO"), ("Scaling tendency", scalant)
        ].data
        hpro_data = sea_water[
            ("water_sim_cases", "SW_HPRO"), ("HP Scaling tendency", scalant)
        ].data
        ro_data[ro_data != ro_data] = hpro_data[ro_data != ro_data]

        sea_water["stacked_data", ("HP Scaling tendency", scalant)].set_data(ro_data)
    bgw_cases = [
        ("water_sim_cases", "BGW"),
        ("water_sim_cases", f"BGW_500"),
        ("water_sim_cases", f"BGW_1500"),
    ]

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
            "Region 3": {"range": [64, 69], "color": "#fdbb84"},
            "Region 4": {"range": [69, 90], "color": "#fc8d59"},
        },
        ("water_sim_cases", f"BGW_1500"): {
            "Region 4": {"range": [50, 90], "color": "#fc8d59"},
        },
    }

    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water[("water_sim_cases", "SW_RO"), ("RO1", "flux")].data * 3600,
        ylabel="Flux (LMH)",
        yticks=[0, 10, 20, 30],
        data_label="RO",
        save_fig=False,
        width=1.8,
        height=1.6,
        save_name="seawater/hpro_avg_flux",
    )
    plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water[("water_sim_cases", "SW_HPRO"), ("HPRO1", "flux")].data * 3600,
        ylabel="Flux (LMH)",
        yticks=[0, 10, 20, 30],
        data_label="HPRO",
        save_fig=True,
        fig=fig,
        marker="s",
        color="red",
        save_name="seawater/hpro_avg_flux",
        regions=regions["seawater"],
    )
    for point in ["inlet", "outlet"]:
        fig = plot_fig(
            sea_water["stacked_data", "Water recovery"],
            sea_water[("water_sim_cases", "SW_RO"), ("RO1", f"flux {point}")].data
            * 3600,
            ylabel="Flux (LMH)",
            yticks=[0, 10, 20, 30, 40, 50, 60],
            data_label="RO",
            save_fig=False,
            width=1.8,
            height=1.6,
            save_name=f"seawater/hpro_{point}_flux",
        )
        plot_fig(
            sea_water["stacked_data", "Water recovery"],
            sea_water[("water_sim_cases", "SW_HPRO"), ("HPRO1", f"flux {point}")].data
            * 3600,
            ylabel="Flux (LMH)",
            yticks=[0, 10, 20, 30, 40, 50, 60],
            data_label="HPRO",
            save_fig=True,
            fig=fig,
            marker="s",
            color="red",
            save_name=f"seawater/hpro_{point}_flux",
            regions=regions["seawater"],
        )
    # fig = plot_fig(
    #     sea_water["stacked_data", "Water recovery"],
    #     sea_water[("water_sim_cases", "SW_RO"), ("RO1", "area")],
    #     ylabel="Area (m$^2$)",
    #     yticks=[0, 50, 100, 150, 200],
    #     data_label="RO",
    #     save_fig=False,
    #     width=1.8,
    #     height=1.6,
    #     save_name="hpro_area",
    # )
    # plot_fig(
    #     sea_water["stacked_data", "Water recovery"],
    #     sea_water[("water_sim_cases", "SW_HPRO"), ("HPRO1", "area")],
    #     ylabel="Area (m$^2$)",
    #     yticks=[0, 50, 100, 150, 200],
    #     data_label="HPRO",
    #     save_fig=True,
    #     color="red",
    #     fig=fig,
    #     marker="s",
    #     save_name="hpro_area",
    #     regions=regions["seawater"],
    # )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water[("water_sim_cases", "SW_RO"), ("RO1", "velocity")],
        ylabel="Velocity (m/s)",
        yticks=[0, 0.1, 0.2, 0.3],
        data_label="RO",
        save_fig=False,
        width=1.8,
        height=1.6,
        save_name="hpro_Velocity",
    )
    plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water[("water_sim_cases", "SW_HPRO"), ("HPRO1", "velocity")],
        ylabel="Velocity (m/s)",
        yticks=[0, 0.1, 0.2, 0.3],
        data_label="HPRO",
        save_fig=True,
        fig=fig,
        color="red",
        marker="s",
        save_name="seawater/hpro_Velocity",
        regions=regions["seawater"],
    )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water[("water_sim_cases", "SW_RO"), ("RO1", "pressure")],
        ylabel="Pressure (bar)",
        yticks=[0, 100, 200, 300, 400],
        data_label="RO",
        save_fig=False,
        width=1.8,
        height=1.6,
        save_name="hpro_Pressure",
    )
    plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water[("water_sim_cases", "SW_HPRO"), ("HPRO1", "pressure")],
        ylabel="Pressure (bar)",
        yticks=[0, 100, 200, 300, 400],
        data_label="HPRO",
        save_fig=True,
        fig=fig,
        color="red",
        marker="s",
        save_name="seawater/hpro_Pressure",
        regions=regions["seawater"],
    )

    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water["stacked_data", ("Softening_RKT_1", "Chemical dose", "Soda ash")],
        ylabel="Soda ash (PPM)",
        yticks=[0, 400, 800, 1200],
        save_fig=True,
        width=1.8,
        height=1.6,
        save_name=f"seawater/hpro_soda ash",
        regions=regions["seawater"],
        # data_label="Soda ash",
    )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water["stacked_data", ("Softening_RKT_1", "Chemical dose", "Lime")],
        ylabel="Lime (PPM)",
        yticks=[0, 50, 100, 150, 200],
        # save_fig=True,
        width=1.8,
        height=1.6,
        data_label="Lime",
        marker="s",
        save_name=f"seawater/hpro_lime",
        regions=regions["seawater"],
        # fig=fig,
    )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water["stacked_data", ("acid_addition", "Chemical dose", "HCl")].data
        * 0.35,
        ylabel="Acid dose (PPM)",
        yticks=[0, 20, 40, 60, 80],
        save_fig=False,
        width=1.8,
        height=1.6,
        data_label="HCl",
        save_name=f"seawater/hpro_Velocity",
    )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water["stacked_data", ("acid_addition", "Chemical dose", "H2SO4")].data
        * 0.92,
        ylabel="Acid dose (PPM)",
        yticks=[0, 20, 40, 60, 80],
        save_fig=True,
        width=1.8,
        height=1.6,
        data_label="H$_2$SO$_4$",
        color="red",
        marker="s",
        regions=regions["seawater"],
        save_name=f"hpro_acid",
        fig=fig,
    )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water["stacked_data", ("Softening_RKT_1", "Alkalinity")],
        ylabel="Alkalinity (PPM)",
        yticks=[0, 200, 400, 600],
        save_fig=True,
        width=1.8,
        height=1.6,
        marker="s",
        regions=regions["seawater"],
        save_name=f"seawater/hpro_alkalinity",
    )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        get_different(
            sea_water["stacked_data", ("Softening_RKT_1", "Effluent in", "Ca")],
            sea_water["stacked_data", ("Softening_RKT_1", "Effluent out", "Ca")],
        ),
        ylabel="Ion removal (%)",
        yticks=[0, 25, 50, 75, 100],
        save_fig=False,
        width=1.8,
        height=1.6,
        data_label="Ca",
        save_name=f"seawater/hpro_velocity",
    )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        get_different(
            sea_water["stacked_data", ("Softening_RKT_1", "Effluent in", "HCO3")],
            sea_water["stacked_data", ("Softening_RKT_1", "Effluent out", "HCO3")],
        ),
        ylabel="Ion removal (%)",
        yticks=[0, 25, 50, 75, 100],
        save_fig=True,
        width=1.8,
        height=1.6,
        data_label="HCO$_3$",
        color="red",
        marker="s",
        regions=regions["seawater"],
        save_name=f"seawater/hpro_ion_removal",
        fig=fig,
    )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        get_different(
            sea_water["stacked_data", ("Softening_RKT_1", "Effluent in", "Mg")],
            sea_water["stacked_data", ("Softening_RKT_1", "Effluent out", "Mg")],
        ),
        ylabel="Ion removal (%)",
        yticks=[0, 25, 50, 75, 100],
        save_fig=True,
        width=1.8,
        height=1.6,
        data_label="Mg",
        color="blue",
        marker="d",
        regions=regions["seawater"],
        save_name=f"seawater/hpro_ion_removal",
        fig=fig,
    )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water["stacked_data", ("HP Scaling tendency", "Calcite")],
        ylabel="Scaling tendency",
        yticks=[0.6, 0.8, 1.0, 1.2],
        save_fig=False,
        width=1.8,
        height=1.6,
        data_label="Calcite",
        save_name=f"seawater/hpro_scaling",
    )
    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water["stacked_data", ("HP Scaling tendency", "Gypsum")],
        ylabel="Scaling tendency",
        yticks=[0.6, 0.8, 1.0, 1.2],
        save_fig=True,
        width=1.8,
        height=1.6,
        # data_label="H$_2$SO$_4$",
        marker="s",
        regions=regions["seawater"],
        save_name=f"seawater/hpro_scaling",
        fig=fig,
        color="red",
        data_label="Gypsum",
    )

    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water["stacked_data", ("softening", "pH")],
        ylabel="Softened pH",
        yticks=[6, 7, 8, 9, 10, 11],
        save_fig=True,
        width=1.8,
        height=1.6,
        save_name="seawater/hpro_soft_ph",
        regions=regions["seawater"],
    )

    fig = plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water[("water_sim_cases", "SW_RO"), ("RO1", "pH")],
        ylabel="Retentate pH",
        yticks=[6, 7, 8, 9, 10, 11],
        data_label="RO",
        save_fig=False,
        width=1.8,
        height=1.6,
        save_name="seawater/hpro_ph",
    )
    plot_fig(
        sea_water["stacked_data", "Water recovery"],
        sea_water[("water_sim_cases", "SW_HPRO"), ("HPRO1", "pH")],
        ylabel="Retentate pH",
        yticks=[6, 7, 8, 9, 10, 11],
        data_label="HPRO",
        color="red",
        save_fig=True,
        fig=fig,
        marker="s",
        save_name="seawater/hpro_ph",
        regions=regions["seawater"],
    )
    for i, hd in enumerate(bgw_cases):
        water_case = water_reference[hd[1]]
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("RO1", "flux")].data * 3600,
            ylabel="Flux (LMH)",
            yticks=[0, 10, 20, 30, 40],
            save_fig=True,
            width=1.8,
            height=1.6,
            save_name=f"{water_case}/{water_case}_avg_flux",
            regions=regions[hd],
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("RO1", "area")],
            ylabel="Area (m$^2$)",
            yticks=[0, 50, 100, 150, 200],
            save_fig=True,
            width=1.8,
            height=1.6,
            save_name=f"{water_case}/{water_case}_area",
            regions=regions[hd],
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("RO1", "pressure")],
            ylabel="Pressure (bar)",
            yticks=[0, 10, 20, 30],
            save_fig=True,
            width=1.8,
            height=1.6,
            save_name=f"{water_case}/{water_case}_Pressure",
            regions=regions[hd],
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("RO1", "velocity")],
            ylabel="Velocity (m/s)",
            yticks=[0, 0.1, 0.2, 0.3],
            save_fig=True,
            width=1.8,
            height=1.6,
            save_name=f"{water_case}/{water_case}_Velocity",
            regions=regions[hd],
        )

        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("Softening_RKT_1", "Chemical dose", "Soda ash")],
            ylabel="Soda ash (PPM)",
            yticks=[0, 400, 800, 1200],
            save_fig=True,
            width=1.8,
            height=1.6,
            save_name=f"{water_case}/{water_case}_soda ash",
            regions=regions[hd],
            # data_label="Soda ash",
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("Softening_RKT_1", "Chemical dose", "Lime")],
            ylabel="Lime (PPM)",
            yticks=[0, 50, 100, 150, 200],
            # save_fig=True,
            width=1.8,
            height=1.6,
            data_label="Lime",
            marker="s",
            save_name=f"{water_case}/{water_case}_lime",
            regions=regions[hd],
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("acid_addition", "Chemical dose", "HCl")].data * 0.35,
            ylabel="Acid dose (PPM)",
            yticks=[0, 20, 40, 60, 80],
            save_fig=False,
            width=1.8,
            height=1.6,
            data_label="HCl",
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("acid_addition", "Chemical dose", "H2SO4")].data * 0.92,
            ylabel="Acid dose (PPM)",
            yticks=[0, 20, 40, 60, 80],
            save_fig=True,
            width=1.8,
            height=1.6,
            data_label="H$_2$SO$_4$",
            color="red",
            marker="s",
            regions=regions[hd],
            save_name=f"{water_case}/{water_case}_acid",
            fig=fig,
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            get_different(
                data_manager[hd, ("Softening_RKT_1", "Effluent in", "Ca")],
                data_manager[hd, ("Softening_RKT_1", "Effluent out", "Ca")],
            ),
            ylabel="Ion removal (%)",
            yticks=[0, 50, 100],
            save_fig=False,
            width=1.8,
            height=1.6,
            data_label="Ca",
            save_name=f"{water_case}/{water_case}_ion_removal",
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            get_different(
                data_manager[hd, ("Softening_RKT_1", "Effluent in", "HCO3")],
                data_manager[hd, ("Softening_RKT_1", "Effluent out", "HCO3")],
            ),
            ylabel="Ion removal (%)",
            yticks=[0, 25, 50, 75, 100],
            save_fig=True,
            width=1.8,
            height=1.6,
            data_label="HCO$_3$",
            color="red",
            marker="s",
            regions=regions[hd],
            save_name=f"{water_case}/{water_case}_ion_removal",
            fig=fig,
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            get_different(
                data_manager[hd, ("Softening_RKT_1", "Effluent in", "Mg")],
                data_manager[hd, ("Softening_RKT_1", "Effluent out", "Mg")],
            ),
            ylabel="Ion removal (%)",
            yticks=[0, 25, 50, 75, 100],
            save_fig=True,
            width=1.8,
            height=1.6,
            data_label="Mg",
            color="blue",
            marker="d",
            regions=regions[hd],
            save_name=f"{water_case}/{water_case}_ion_removal",
            fig=fig,
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("Scaling tendency", "Calcite")],
            ylabel="Scaling tendency",
            yticks=[0.6, 0.8, 1.0, 1.2],
            save_fig=False,
            width=1.8,
            height=1.6,
            data_label="Calcite",
            regions=regions[hd],
            save_name=f"{water_case}/{water_case}_scaling",
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("Scaling tendency", "Gypsum")],
            ylabel="Scaling tendency",
            yticks=[0.6, 0.8, 1.0, 1.2],
            save_fig=True,
            width=1.8,
            height=1.6,
            # data_label="H$_2$SO$_4$",
            marker="s",
            regions=regions[hd],
            save_name=f"{water_case}/{water_case}_scaling",
            fig=fig,
            data_label="Gypsum",
            color="red",
        )

        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("RO1", "pH")],
            ylabel="Retentate pH",
            yticks=[6, 7, 8, 9, 10, 11],
            save_fig=True,
            width=1.8,
            regions=regions[hd],
            height=1.6,
            save_name=f"{water_case}/{water_case}_ro_ph",
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("softening", "pH")],
            ylabel="Softened pH",
            yticks=[6, 7, 8, 9, 10, 11],
            save_fig=True,
            regions=regions[hd],
            width=1.8,
            height=1.6,
            save_name=f"{water_case}/{water_case}_soft_ph",
        )
        fig = plot_fig(
            data_manager[hd, "Water recovery"],
            data_manager[hd, ("Softening_RKT_1", "Alkalinity")],
            ylabel="Alkalinity (PPM)",
            yticks=[0, 200, 400, 600],
            save_fig=True,
            width=1.8,
            height=1.6,
            marker="s",
            regions=regions[hd],
            save_name=f"{water_case}/{water_case}_alkalinity",
        )
    if show_figs:
        fig.show()


if __name__ == "__main__":
    main()
