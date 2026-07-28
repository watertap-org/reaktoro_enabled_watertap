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

import pytest
from reaktoro_enabled_watertap.unit_models.nanofiltration_unit import (
    Nanofiltration,
)

from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
from pyomo.environ import (
    assert_optimal_termination,
)
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.core.util.scaling as iscale
from pyomo.environ import (
    TransformationFactory,
)

from reaktoro_enabled_watertap.unit_models.tests.test_multi_comp_feed_product import (
    build_case,
)
from idaes.core.util.model_statistics import large_residuals_set
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

__author__ = "Chenyu Wang"


DEFAULT_REJECTIONS = {
    "Na_+": 0.3,
    "Ca_2+": 0.95,
    "Mg_2+": 0.90,
    "SO4_2-": 0.98,
    "Cl_-": 0.25,
    "K_+": 0.35,
    "HCO3_-": 0.50,
}


@pytest.mark.core
@pytest.mark.component
def test_nanofiltration_no_reaktoro():
    m = build_case("USDA_brackish", True)

    m.fs.NF = Nanofiltration(
        default_property_package=m.fs.properties,
        add_reaktoro_chemistry=False,
        component_rejections=DEFAULT_REJECTIONS,
    )

    m.fs.NF.fix_and_scale()
    m.fs.feed.outlet.connect_to(m.fs.NF.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.NF.nanofiltration.pH["feed"].fix()
    m.fs.NF.initialize()
    m.fs.NF.report()

    dg = DiagnosticsToolbox(m)
    dg.report_structural_issues()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_watertap_solver(linear_solver="mumps", max_iter=50)
    result = solver.solve(m, tee=True)

    m.fs.NF.report()
    assert_optimal_termination(result)

    assert (
        pytest.approx(m.fs.NF.nanofiltration.flux_vol_solvent[0, "H2O"].value, abs=1e-6)
        == 9.9768e-6
    )

    assert (
        pytest.approx(m.fs.NF.nanofiltration.pH["retentate"].value, abs=1e-6)
        == m.fs.NF.nanofiltration.pH["feed"].value
    )

    assert (
        pytest.approx(
            m.fs.NF.nanofiltration.recovery_vol_phase[0, "Liq"].value, rel=1e-4
        )
        == 0.5
    )


@pytest.mark.core
@pytest.mark.component
def test_nanofiltration_with_reaktoro():
    m = build_case("USDA_brackish", True)

    m.fs.NF = Nanofiltration(
        default_property_package=m.fs.properties,
        add_reaktoro_chemistry=True,
        selected_scalants={"Calcite": 1, "Gypsum": 1},
        component_rejections=DEFAULT_REJECTIONS,
    )

    m.fs.NF.fix_and_scale()
    m.fs.feed.outlet.connect_to(m.fs.NF.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    for props in [
        m.fs.NF.nanofiltration.feed_side.properties_in[0],
        m.fs.NF.nanofiltration.feed_side.properties_out[0],
    ]:
        iscale.set_scaling_factor(props.flow_mol_phase_comp["Liq", "H2O"], 0.1)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()

    m.fs.NF.nanofiltration.pH["feed"].fix(7)

    m.fs.NF.initialize()
    m.fs.NF.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_watertap_solver(linear_solver="mumps", max_iter=50)
    result = solver.solve(m, tee=True)

    m.fs.NF.report()

    # Diagnostics
    for lr in large_residuals_set(m, tol=1e-8, return_residual_values=True):
        print(lr)
    for i in iscale.list_unscaled_variables(m):
        print(i)
    for i in iscale.list_unscaled_constraints(m):
        print(i)
    for var, scale in iscale.badly_scaled_var_generator(m):
        print(
            "Badly scaled variable:",
            var.name,
            var.value,
            iscale.get_scaling_factor(var),
        )

    assert_optimal_termination(result)

    assert pytest.approx(m.fs.NF.nanofiltration.pH["feed"].value, rel=1e-4) == 7
    osm_feed = (
        m.fs.NF.nanofiltration.feed_side.properties_in[0]
        .pressure_osm_phase["Liq"]
        .value
    )
    assert pytest.approx(osm_feed, rel=1e-4) == 184453.4
    assert (
        pytest.approx(
            m.fs.NF.nanofiltration.scaling_tendency["feed", "Calcite"].value, rel=1e-4
        )
        == 0.9474
    )
    assert (
        pytest.approx(
            m.fs.NF.nanofiltration.scaling_tendency["feed", "Gypsum"].value, rel=1e-4
        )
        == 0.29900
    )

    assert pytest.approx(m.fs.NF.nanofiltration.pH["retentate"].value, rel=1e-4) == 7
    osm_retentate = (
        m.fs.NF.nanofiltration.feed_side.properties_out[0]
        .pressure_osm_phase["Liq"]
        .value
    )
    assert pytest.approx(osm_retentate, rel=1e-4) == 254546.7
    assert (
        pytest.approx(
            m.fs.NF.nanofiltration.scaling_tendency["retentate", "Calcite"].value,
            rel=1e-4,
        )
        == 2.1503
    )
    assert (
        pytest.approx(
            m.fs.NF.nanofiltration.scaling_tendency["retentate", "Gypsum"].value,
            rel=1e-4,
        )
        == 0.74324
    )
