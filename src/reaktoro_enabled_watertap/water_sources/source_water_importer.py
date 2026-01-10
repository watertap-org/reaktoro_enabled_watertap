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

import yaml
from pyomo.environ import (
    units as pyunits,
)
from reaktoro_enabled_watertap.utils.report_util import get_lib_path

__author__ = "Alexander V. Dudchenko"


def get_source_water_data(water_source, file_location=None):
    """simple function to load feed water compostion from yaml file"""
    if file_location is None:
        file_location = get_lib_path() / "water_sources" / water_source
    with open(file_location, "r") as ymlfile:
        data_dict = yaml.safe_load(ymlfile)
    # Converts yaml structure to dict structure for use with MCAS
    mcas_param_dict = {}
    mcas_param_dict["solute_list"] = get_solute_dict(data_dict)
    mcas_param_dict["diffusivity_data"] = gen_diffusivity_dict(data_dict)
    mcas_param_dict["mw_data"] = gen_mw_dict(data_dict)
    mcas_param_dict["stokes_radius_data"] = gen_stoke_dict(data_dict)
    mcas_param_dict["charge"] = gen_charge_dict(data_dict)

    # Creats dict with feed properties to pass into multi_comp_feed
    mass_comp_dict = get_feed_comp(data_dict)
    pH = float(data_dict["pH"])
    feed_temperature = data_dict.get("temperature", 293.15)
    alkalinity = data_dict.get("alkalinity_as_CaCO3", None)

    if alkalinity != None:
        alkalinity = float(alkalinity) * pyunits.mg / pyunits.L
    feed_spec_dict = {
        "ion_concentrations": mass_comp_dict,
        "pH": pH,
        "temperature": feed_temperature,
        "alkalinity_as_CaCO3": alkalinity,
    }
    if data_dict.get("flow_mass", None) is not None:
        feed_spec_dict["mass_flowrate"] = (
            data_dict.get("flow_mass", None) * pyunits.kg / pyunits.s
        )
    if data_dict.get("volumetric_flowrate", None) is not None:
        feed_spec_dict["volumetric_flowrate"] = (
            data_dict.get("volumetric_flowrate") * pyunits.L / pyunits.s
        )
    return mcas_param_dict, feed_spec_dict


def get_solute_dict(data_dict):
    solute_list = list(data_dict["solute_list"].keys())
    return solute_list


def gen_diffusivity_dict(data_dict):
    diff_dict = {}
    for solute in data_dict["solute_list"].keys():
        diff_dict[("Liq", solute)] = float(
            data_dict["solute_list"][solute].get("diffusivity", 0)
        )
    return diff_dict


def gen_mw_dict(data_dict):
    mw_dict = {}
    for solute in data_dict["solvent_list"].keys():
        mw_dict[solute] = float(
            data_dict["solvent_list"][solute].get("molecular_weight (kg/mol)", 0)
        )
    for solute in data_dict["solute_list"].keys():
        mw_dict[solute] = float(
            data_dict["solute_list"][solute].get("molecular_weight (kg/mol)", 0)
        )
    return mw_dict


def gen_stoke_dict(data_dict):
    stokes_dict = {}
    for solute in data_dict["solute_list"].keys():
        stokes_dict[solute] = float(
            data_dict["solute_list"][solute].get("stokes_radius (m)", 0)
        )
    return stokes_dict


def gen_charge_dict(data_dict):
    charge_dict = {}
    for solute in data_dict["solute_list"].keys():
        charge_dict[solute] = float(
            data_dict["solute_list"][solute].get("elemental charge", 0)
        )
    return charge_dict


def get_feed_comp(data_dict):
    mass_loading_dict = {}
    for solute in data_dict["solute_list"].keys():
        value = float(
            data_dict["solute_list"][solute].get("concentration (mg/L)", None)
        )
        if value is None:
            raise ValueError(
                f"Concentration for {solute} not found in the data dictionary."
            )
        mass_loading_dict[solute] = value * pyunits.mg / pyunits.L

    return mass_loading_dict
