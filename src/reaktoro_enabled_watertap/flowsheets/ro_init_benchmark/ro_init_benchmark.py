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
"""
Benchmark flowsheet (feed -> pump -> RO -> ERD, no reaktoro) for testing how the
initial pump pressure guess and initial RO area affect RO initialization and the
first flowsheet optimization solves.

Initialization sequence mirrors softening_acid_ro.initialize:
    1. sequential unit initialization (pump pressure fixed at guess, RO init objective)
    2. square solve
    3. set optimization operation, fix system recovery at 0.5 (or current if < 0.5), solve
    4. step to final recovery (0.8 brackish / 0.6 seawater), solve

Run:
    python ro_init_benchmark.py
"""

import csv
import itertools
import os
import time

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Objective,
    Reals,
    TransformationFactory,
    Var,
    check_optimal_termination,
    value,
    units as pyunits,
)
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.costing import WaterTAPCosting
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
)
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)

from reaktoro_enabled_watertap.water_sources.source_water_importer import (
    get_source_water_data,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_feed_unit import MultiCompFeed
from reaktoro_enabled_watertap.unit_models.multi_comp_pump_unit import (
    MultiCompPumpUnit,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_ro_unit import MultiCompROUnit
from reaktoro_enabled_watertap.unit_models.multi_comp_erd_unit import (
    MultiCompERDUnit,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_product_unit import (
    MultiCompProduct,
)
from reaktoro_enabled_watertap.utils import scale_utils as scu

__author__ = "Alexander V. Dudchenko"

WATERS = [
    "USDA_brackish.yaml",
    "sample_500_hardness.yaml",
    "sample_1500_hardness.yaml",
    "Seawater.yaml",
]


def build_model(
    water_case,
    osmotic_over_pressure=2,
    osmotic_pressure_offset=2e5 * pyunits.Pa,
    initialization_pressure="osmotic_pressure",
    initial_area=None,
    ro_target_recovery=0.5,
    feed_flow_rate=5000 * pyunits.m**3 / pyunits.day,
    maximum_pressure=85 * pyunits.bar,
):
    """Builds feed -> pump -> RO -> ERD flowsheet without reaktoro.

    Args:
        water_case: water source yaml file
        osmotic_over_pressure: pump guess = osmotic pressure * this + offset
        osmotic_pressure_offset: pump guess offset (pyomo units)
        initialization_pressure: 'osmotic_pressure' or fixed pressure with pyomo units
        initial_area: initial RO area guess in m2 (None keeps RO unit default of 100 m2)
        ro_target_recovery: RO init objective target recovery
        feed_flow_rate: volumetric feed flow
        maximum_pressure: pump maximum pressure
    """
    mcas_props, feed_specs = get_source_water_data(water_case)
    feed_specs["volumetric_flowrate"] = feed_flow_rate
    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m = ConcreteModel()
    m.water_case = water_case
    m.fs = FlowsheetBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.ro_properties = SeawaterParameterBlock()

    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        reconcile_using_reaktoro=False,
        **feed_specs,
    )
    m.fs.pump_unit = MultiCompPumpUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        initialization_pressure=initialization_pressure,
        osmotic_over_pressure=osmotic_over_pressure,
        osmotic_pressure_offset=osmotic_pressure_offset,
        maximum_pressure=maximum_pressure,
    )
    m.fs.ro_unit = MultiCompROUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        ro_property_package=m.fs.ro_properties,
        add_reaktoro_chemistry=False,
        target_recovery=ro_target_recovery,
    )
    m.fs.erd_unit = MultiCompERDUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
    )
    m.fs.product = MultiCompProduct(default_property_package=m.fs.properties)
    m.fs.brine = MultiCompProduct(default_property_package=m.fs.properties)

    m.flowsheet_unit_order = [
        m.fs.feed,
        m.fs.pump_unit,
        m.fs.ro_unit,
        m.fs.erd_unit,
        m.fs.product,
        m.fs.brine,
    ]

    m.fs.feed.outlet.connect_to(m.fs.pump_unit.inlet)
    m.fs.pump_unit.outlet.connect_to(m.fs.ro_unit.feed)
    m.fs.ro_unit.retentate.connect_to(m.fs.erd_unit.inlet)
    m.fs.ro_unit.product.connect_to(m.fs.product.inlet)
    m.fs.erd_unit.outlet.connect_to(m.fs.brine.inlet)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.product.product.properties[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.product.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.product.properties[0].flow_vol
    )
    m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0.0, 1),
        domain=Reals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )
    m.fs.eq_water_recovery = Constraint(
        expr=sum(m.fs.feed.feed.properties[0].flow_mass_phase_comp["Liq", :])
        * m.fs.water_recovery
        == m.fs.product.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    iscale.set_scaling_factor(m.fs.water_recovery, 1)
    iscale.constraint_scaling_transform(m.fs.eq_water_recovery, 1)

    for unit in m.flowsheet_unit_order:
        unit.fix_and_scale()
    iscale.calculate_scaling_factors(m)
    scu.scale_costing_block(m.fs.costing)
    # override initial area guess (set_fixed_operation sets 100 m2)
    if initial_area is not None:
        m.fs.ro_unit.ro_unit.area.fix(initial_area)
    assert degrees_of_freedom(m) == 0
    return m


def solve_model(m, tee=False):
    solver = get_cyipopt_watertap_solver(
        linear_solver="mumps", max_iter=1000, pivtol=1e-3, pivtolmax=1e0
    )
    result = solver.solve(m, tee=tee)
    return check_optimal_termination(result)


def ro_state(m):
    ro = m.fs.ro_unit.ro_unit
    to_lmh = 3600  # kg/m2/s -> kg/m2/hr (~LMH)
    # flux is not defined at x=0 for BACKWARD scheme, use first interior node
    n0 = ro.length_domain.next(ro.length_domain.first())
    n1 = ro.length_domain.last()
    return {
        "pump_pressure_bar": value(m.fs.pump_unit.pump.outlet.pressure[0]) / 1e5,
        "feed_osm_bar": value(
            m.fs.pump_unit.pump.control_volume.properties_in[0].pressure_osm_phase[
                "Liq"
            ]
        )
        / 1e5,
        "ro_area_m2": value(ro.area),
        "ro_recovery": value(ro.recovery_vol_phase[0, "Liq"]),
        "flux_avg_lmh": to_lmh * value(ro.flux_mass_phase_comp_avg[0, "Liq", "H2O"]),
        "flux_in_lmh": to_lmh * value(ro.flux_mass_phase_comp[0, n0, "Liq", "H2O"]),
        "flux_out_lmh": to_lmh * value(ro.flux_mass_phase_comp[0, n1, "Liq", "H2O"]),
        "inlet_velocity": value(ro.feed_side.velocity[0, 0]),
        "system_recovery": value(m.fs.water_recovery),
        "lcow": value(m.fs.costing.LCOW),
    }


def run_case(water, final_recovery=None, tee=False, **build_kwargs):
    """Builds and initializes a case, returning a dict with stage status and state.

    stages: build -> unit_init -> square_solve -> opt_init_solve -> final_solve
    """
    if final_recovery is None:
        final_recovery = 0.6 if "Seawater" in water else 0.8
    build_kwargs.setdefault("osmotic_over_pressure", 2)
    build_kwargs.setdefault("osmotic_pressure_offset", 2e5 * pyunits.Pa)
    result = {"water": water, **{k: _fmt(v) for k, v in build_kwargs.items()}}
    result["failed_stage"] = None
    t0 = time.time()
    stage = "build"
    try:
        m = build_model(water, **build_kwargs)
        result["init_area_m2"] = value(m.fs.ro_unit.ro_unit.area)

        stage = "unit_init"
        for unit in m.flowsheet_unit_order:
            unit.initialize()
        m.fs.costing.initialize()
        result.update({f"init_{k}": v for k, v in ro_state(m).items()})

        stage = "square_solve"
        if not solve_model(m, tee=tee):
            raise RuntimeError("not optimal")

        stage = "opt_init_solve"
        for unit in m.flowsheet_unit_order:
            unit.set_optimization_operation()
        if m.fs.water_recovery.value < 0.5:
            m.fs.water_recovery.fix()
            if not solve_model(m, tee=tee):
                raise RuntimeError("not optimal")
        m.fs.water_recovery.fix(0.5)
        if not solve_model(m, tee=tee):
            raise RuntimeError("not optimal")
        result.update({f"opt50_{k}": v for k, v in ro_state(m).items()})

        stage = "final_solve"
        if "Seawater" not in water:
            m.fs.water_recovery.fix(0.7)
            if not solve_model(m, tee=tee):
                raise RuntimeError("not optimal")
        m.fs.water_recovery.fix(final_recovery)
        if not solve_model(m, tee=tee):
            raise RuntimeError("not optimal")
        result.update({f"final_{k}": v for k, v in ro_state(m).items()})
    except Exception as e:  # record failure and move on
        result["failed_stage"] = stage
        result["error"] = f"{type(e).__name__}: {str(e).strip()[:200]}"
    result["time_s"] = round(time.time() - t0, 1)
    return result


def _fmt(v):
    try:
        return value(v)
    except Exception:
        return v


def run_sweep(
    waters=WATERS,
    over_pressures=(1.5, 2, 3, 4, 6.5),
    offsets_bar=(2, 5, 10),
    initial_areas=(None,),
    output_file="ro_init_benchmark_results.csv",
    tee=False,
):
    idaeslog.getLogger("idaes").setLevel(idaeslog.ERROR)
    results = []
    for water, ovp, off, area in itertools.product(
        waters, over_pressures, offsets_bar, initial_areas
    ):
        r = run_case(
            water,
            tee=tee,
            osmotic_over_pressure=ovp,
            osmotic_pressure_offset=off * 1e5 * pyunits.Pa,
            initial_area=area,
        )
        results.append(r)
        print(_summary_line(r), flush=True)
    keys = []
    for r in results:
        keys += [k for k in r if k not in keys]
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(results)
    return results


def _summary_line(r):
    status = "OK" if r["failed_stage"] is None else f"FAIL@{r['failed_stage']}"
    g = lambda k: r.get(k, float("nan"))
    return (
        f"{r['water']:<26} ovp={r['osmotic_over_pressure']:<4} "
        f"off={r['osmotic_pressure_offset'] / 1e5:<5.1f}bar "
        f"area0={g('init_area_m2'):<8.0f} | "
        f"P={g('init_pump_pressure_bar'):6.1f}bar "
        f"rec={g('init_ro_recovery'):.3f} "
        f"Jin={g('init_flux_in_lmh'):5.1f} Javg={g('init_flux_avg_lmh'):5.1f} "
        f"A={g('init_ro_area_m2'):7.0f} | "
        f"LCOW={g('final_lcow'):.4f} {status} ({r['time_s']}s)"
    )


if __name__ == "__main__":
    run_sweep(initial_areas=(None, 1000, 5000))
