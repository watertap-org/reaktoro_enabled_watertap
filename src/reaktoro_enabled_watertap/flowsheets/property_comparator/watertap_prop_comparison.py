from reaktoro_pse.reaktoro_block import ReaktoroBlock
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
from pyomo.environ import (
    ConcreteModel,
    Var,
    value,
    Constraint,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import idaes.core.util.scaling as iscale

import reaktoro as rkt

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
)
from idaes.models.unit_models import (
    Product,
    Feed,
    StateJunction,
)
import watertap.property_models.seawater_prop_pack as props_seawater

import watertap.property_models.NaCl_T_dep_prop_pack as props_nacl
from idaes.core import (
    FlowsheetBlock,
)

import numpy as np
from idaes.core.util.model_statistics import degrees_of_freedom
import csv

from reaktoro_enabled_watertap.utils.report_util import get_lib_path
from reaktoro_enabled_watertap.water_sources.source_water_importer import (
    get_source_water_data,
)


def main():
    m = build_model("sample_1500_hardness.yaml")
    m.fs.water_recovery.fix(0.5)
    result = solve_model(m)
    assert_optimal_termination(result)
    print_comparison(m)


def print_comparison(m):
    header = [
        "Temperature (K)",
        "TDS",
        "Reaktoro_osmoticPressure (Pa)",
        "Seawater_osmoticPressure (Pa)",
        "NaCl_osmoticPressure( Pa)",
        "Reaktoro_vaporPressure (Pa)",
        "Seawater_vaporPressure (Pa)",
        "NaCl_vaporPressure (Pa)",
        "Reaktoro_density (kg/m3)",
        "Seawater_density (kg/m3)",
        "NaCl_density (kg/m3)",
    ]
    for phase, ion in m.fs.multicomp_feed.properties[0].conc_mass_phase_comp:
        header.append(f"{ion} (mg/L)")

    data_row = [
        m.fs.multicomp_feed.properties[0].temperature.value,
        m.fs.feed_tds.value,
        m.fs.modified_properties["osmoticPressure", "H2O"].value,
        m.fs.seawater_feed.properties[0].pressure_osm_phase["Liq"].value,
        m.fs.nacl_feed.properties[0].pressure_osm_phase["Liq"].value,
        m.fs.modified_properties["vaporPressure", "H2O(g)"].value,
        m.fs.seawater_feed.properties[0].pressure_sat.value,
        m.fs.nacl_feed.properties[0].pressure_sat.value,
        m.fs.multicomp_feed.properties[0].dens_mass_phase["Liq"].value,
        m.fs.seawater_feed.properties[0].dens_mass_phase["Liq"].value,
        m.fs.nacl_feed.properties[0].dens_mass_phase["Liq"].value,
    ]
    for phase, ion in m.fs.multicomp_feed.properties[0].conc_mass_phase_comp:
        data_row.append(
            value(
                pyunits.convert(
                    m.fs.multicomp_feed.properties[0].conc_mass_phase_comp[phase, ion],
                    to_units=pyunits.mg / pyunits.L,
                )
            )
        )
    result = {}
    for header, val in zip(header, data_row):
        print(f"{header}: {float(val)}")
        result[header] = float(val)
    return result


def build_model(water_case):

    mcas_props, feed_specs = get_source_water_data(water_case)
    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.mcas_props = MCASParameterBlock(**mcas_props)
    m.fs.seawater_props = props_seawater.SeawaterParameterBlock()
    m.fs.nacl_props = props_nacl.NaClParameterBlock()
    m.fs.multicomp_feed = Feed(property_package=m.fs.mcas_props)
    m.fs.multicomp_feed.properties[0].total_dissolved_solids[...]
    m.fs.multicomp_feed.properties[0].flow_vol_phase[...]
    m.fs.seawater_feed = Feed(property_package=m.fs.seawater_props)
    m.fs.seawater_feed.properties[0].pressure_osm_phase[...]
    m.fs.seawater_feed.properties[0].pressure_sat[...]
    m.fs.seawater_feed.properties[0].flow_vol_phase[...]
    m.fs.nacl_feed = Feed(property_package=m.fs.nacl_props)
    m.fs.nacl_feed.properties[0].pressure_osm_phase[...]
    m.fs.nacl_feed.properties[0].pressure_sat[...]
    m.fs.nacl_feed.properties[0].flow_vol_phase[...]
    m.fs.feed_pH = Var(
        initialize=feed_specs["pH"], bounds=(4, 12), units=pyunits.dimensionless
    )
    m.fs.feed_pH.fix()

    set_feed_composition(m, feed_specs)
    m.fs.water_recovery = Var(
        initialize=1e-8,
        units=pyunits.dimensionless,
    )
    m.fs.eq_water_recovery = Constraint(
        expr=(1 - m.fs.water_recovery)
        == m.fs.multicomp_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    # add global TDS constraint
    m.fs.feed_tds = Var(
        initialize=1e-8,
        units=pyunits.g / pyunits.L,
    )
    m.fs.feed_tds.fix(m.fs.multicomp_feed.properties[0].total_dissolved_solids)

    # make sure TDS is equal across all property packages
    m.fs.multicomp_feed.eq_TDS = Constraint(
        expr=m.fs.feed_tds
        == pyunits.convert(
            m.fs.multicomp_feed.properties[0].total_dissolved_solids,
            to_units=pyunits.g / pyunits.L,
        )
    )

    m.fs.seawater_feed.eq_TDS = Constraint(
        expr=m.fs.feed_tds
        == pyunits.convert(
            m.fs.seawater_feed.properties[0].conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    m.fs.nacl_feed.eq_TDS = Constraint(
        expr=m.fs.feed_tds
        == pyunits.convert(
            m.fs.nacl_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
            to_units=pyunits.g / pyunits.L,
        )
    )

    # fix the mass flows of each ion and unfix water mass flow, it wil be adjusted to match TDS,
    # similarly we unfix Cl in multicomp to allow charge balance
    for idx in m.fs.multicomp_feed.properties[0].conc_mass_phase_comp:
        m.fs.multicomp_feed.properties[0].conc_mass_phase_comp[idx].unfix()
        m.fs.multicomp_feed.properties[0].flow_mass_phase_comp[idx].fix()

    m.fs.multicomp_feed.properties[0].flow_mass_phase_comp["Liq", "Cl_-"].unfix()
    m.fs.multicomp_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].unfix()

    m.fs.nacl_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].unfix()
    m.fs.nacl_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix()

    m.fs.multicomp_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.nacl_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].unfix()

    m.fs.seawater_feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix()
    m.fs.seawater_feed.properties[0].conc_mass_phase_comp["Liq", "TDS"].unfix()
    m.fs.seawater_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    print(degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    result = solve_model(m)
    assert_optimal_termination(result)
    m.fs.feed_tds.unfix()
    return m


def set_feed_composition(m, feed_specs):
    temperature = feed_specs["temperature"]
    m.fs.multicomp_feed.properties[0].temperature.fix(temperature * pyunits.K)
    m.fs.multicomp_feed.properties[0].pressure.fix(101325 * pyunits.Pa)
    m.fs.nacl_feed.properties[0].temperature.fix(temperature * pyunits.K)
    m.fs.nacl_feed.properties[0].pressure.fix(101325 * pyunits.Pa)
    m.fs.seawater_feed.properties[0].temperature.fix(temperature * pyunits.K)
    m.fs.seawater_feed.properties[0].pressure.fix(101325 * pyunits.Pa)

    for ion in feed_specs["ion_concentrations"]:
        m.fs.multicomp_feed.properties[0].conc_mass_phase_comp["Liq", ion].fix(
            feed_specs["ion_concentrations"][ion]
        )
        m.fs.multicomp_feed.properties[0].flow_mol_phase_comp["Liq", ion].unfix()
    m.fs.multicomp_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        1 * pyunits.kg / pyunits.s
    )
    m.fs.multicomp_feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].unfix()
    assert degrees_of_freedom(m.fs.multicomp_feed) == 0
    result = solve_model(m.fs.multicomp_feed)
    scale_model(m)
    assert degrees_of_freedom(m.fs.multicomp_feed) == 0
    result = solve_model(m.fs.multicomp_feed)
    assert_optimal_termination(result)
    m.fs.nacl_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        1 * pyunits.kg / pyunits.s
    )
    m.fs.seawater_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        1 * pyunits.kg / pyunits.s
    )

    m.fs.nacl_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
    m.fs.seawater_feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].unfix()
    m.fs.nacl_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
        m.fs.multicomp_feed.properties[0].total_dissolved_solids
    )
    m.fs.seawater_feed.properties[0].conc_mass_phase_comp["Liq", "TDS"].fix(
        m.fs.multicomp_feed.properties[0].total_dissolved_solids
    )
    assert degrees_of_freedom(m.fs.seawater_feed) == 0
    result = solve_model(m.fs.seawater_feed)
    assert_optimal_termination(result)
    assert degrees_of_freedom(m.fs.nacl_feed) == 0
    result = solve_model(m.fs.nacl_feed)
    assert_optimal_termination(result)

    add_standard_properties(m)

    m.fs.eq_modified_properties.initialize()
    m.fs.modified_properties[("charge", None)].fix(0)
    assert degrees_of_freedom(m) == 0
    result = solve_model(m)
    assert_optimal_termination(result)
    print(
        "reconciled Cl concentration from:",
        m.fs.initial_neutral_ion,
        "to",
        m.fs.multicomp_feed.properties[0].conc_mass_phase_comp["Liq", "Cl_-"].value,
    )


def add_standard_properties(
    m,
    database="pitzer.dat",
    activity_model="ActivityModelPitzer",
):
    m.fs.modified_properties = Var(
        [
            ("vaporPressure", "H2O(g)"),
            ("osmoticPressure", "H2O"),
            ("charge", None),
        ],
        initialize=1,
    )
    prop_dict = {
        ("density", None): m.fs.multicomp_feed.properties[0].dens_mass_phase["Liq"],
    }
    for key in m.fs.modified_properties:
        prop_dict[key] = m.fs.modified_properties[key]
    print(prop_dict)
    m.fs.multicomp_feed.properties[0].dens_mass_phase["Liq"].unfix()
    # this will charge nutralize the solution and remove water as vapor, allowing assessment of how pressure (osmotic/or vapor changes as function of tempearuter and recovery)
    m.fs.eq_modified_properties = ReaktoroBlock(
        aqueous_phase={
            "composition": m.fs.multicomp_feed.properties[0].flow_mol_phase_comp,
            "convert_to_rkt_species": True,
            "activity_model": activity_model,
        },
        database_file=database,
        system_state={
            "temperature": m.fs.multicomp_feed.properties[0].temperature,
            "pressure": m.fs.multicomp_feed.properties[0].pressure,
            "pH": m.fs.feed_pH,
        },
        gas_phase={
            "phase_components": ["H2O(g)", "Ntg(g)"],
            "activity_model": rkt.ActivityModelPengRobinsonPhreeqc(),
        },
        outputs=prop_dict,
        assert_charge_neutrality=False,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=False,  # direct calculations here
    )
    m.fs.multicomp_feed.properties[0].conc_mass_phase_comp["Liq", "Cl_-"].unfix()
    m.fs.initial_neutral_ion = (
        m.fs.multicomp_feed.properties[0].conc_mass_phase_comp["Liq", "Cl_-"].value
    )

    m.fs.modified_properties[("charge", None)].fix(0)
    m.fs.multicomp_feed.properties[0].eq_dens_mass_phase["Liq"].deactivate()
    iscale.set_scaling_factor(m.fs.modified_properties[("charge", None)], 1e8)


def scale_model(m):

    tds_scale = 0
    for idx in m.fs.multicomp_feed.properties[0].flow_mol_phase_comp:
        scale = 1 / m.fs.multicomp_feed.properties[0].flow_mol_phase_comp[idx].value
        if "H2O" not in idx[1]:
            tds_scale += (
                m.fs.multicomp_feed.properties[0].flow_mass_phase_comp[idx].value
            )
        m.fs.mcas_props.set_default_scaling("flow_mol_phase_comp", scale, index=idx)
    m.fs.nacl_props.set_default_scaling(
        "flow_mass_phase_comp", (1 / (tds_scale / 1000)), index=("Liq", "NaCl")
    )
    m.fs.nacl_props.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.seawater_props.set_default_scaling(
        "flow_mass_phase_comp", (1 / (tds_scale / 1000)), index=("Liq", "NaCl")
    )
    m.fs.seawater_props.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    iscale.calculate_scaling_factors(m)


def solve_model(m, **kwargs):
    cy_solver = get_cyipopt_watertap_solver()  # get_solver(solver="cyipopt-watertap")
    result = cy_solver.solve(m, tee=True)
    return result


if __name__ == "__main__":
    main()
