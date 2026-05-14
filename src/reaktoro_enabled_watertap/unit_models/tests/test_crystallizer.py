import pytest
from reaktoro_enabled_watertap.unit_models.crystallizer import (
    CrystallizationUnit,
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
    units as pyunits,
)
from watertap.costing import WaterTAPCosting

# from brackish_valorization_reaktoro.property_models.tests.test_mcas_with_enthalpy import (
#     build_case,
# )
from reaktoro_enabled_watertap.unit_models.tests.test_multi_comp_feed_product import (
    build_case,
)
from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap.core.util.model_diagnostics.infeasible import *
from idaes.core.util.model_statistics import large_residuals_set

from reaktoro_enabled_watertap.utils import scale_utils as scu
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

__author__ = "Alexander V. Dudchenko"


@pytest.mark.core
@pytest.mark.component
def test_crystallizer_no_reaktoro():
    m = build_case("USDA_brackish", True)

    m.fs.water_properties_vapor = WaterParameterBlock()
    # m.fs.costing = WaterTAPCosting()
    m.fs.crystallizer = CrystallizationUnit(
        default_property_package=m.fs.properties,
        # default_costing_package=m.fs.costing,
        vapor_property_package=m.fs.water_properties_vapor,
        add_reaktoro_chemistry=False,
    )

    m.fs.crystallizer.fix_and_scale()
    m.fs.feed.outlet.connect_to(m.fs.crystallizer.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.crystallizer.initialize()
    m.fs.crystallizer.report()
    dg = DiagnosticsToolbox(m)
    dg.report_structural_issues()
    dg.display_underconstrained_set()
    dg.display_overconstrained_set()
    dg.display_potential_evaluation_errors()
    dg.display_variables_fixed_to_zero()
    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_watertap_solver(linear_solver="mumps", max_iter=50)
    result = solver.solve(m, tee=True)

    m.fs.crystallizer.report()
    m.fs.crystallizer.crystallizer.temperature_operating.fix(315)
    m.fs.crystallizer.crystallizer.evaporation_percent.fix(90)

    result = solver.solve(m, tee=True)
    m.fs.crystallizer.report()
    assert_optimal_termination(result)
    # m.fs.precipitation.report()
    # assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.crystallizer.crystallizer.flow_mass_sol_comp_apparent["Calcite"].value,
            1e-6,
        )
        == 6.4432e-08
    )
    assert (
        pytest.approx(
            m.fs.crystallizer.crystallizer.pH["outlet"].value,
            1e-5,
        )
        == 7.07
    )

    assert (
        pytest.approx(
            m.fs.crystallizer.crystallizer.pressure_operating.value,
            1e-2,
        )
        == 8143.4887616
    )


@pytest.mark.core
@pytest.mark.component
def test_reaktoro_crystallizer():
    m = build_case("USDA_brackish", True)

    m.fs.water_properties_vapor = WaterParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.crystallizer = CrystallizationUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        vapor_property_package=m.fs.water_properties_vapor,
        selected_precipitants=[
            "Calcite",
            # "Anhydrite",
            "Gypsum",
            "Halite",
            "Brucite",
            "Mirabilite",
            "Sylvite",
            "Bischofite",
            "Epsomite",
        ],
    )

    m.fs.costing.cost_process()

    m.fs.feed.outlet.connect_to(m.fs.crystallizer.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.crystallizer.fix_and_scale()

    iscale.calculate_scaling_factors(m)
    scu.scale_costing_block(m.fs.costing)
    m.fs.feed.initialize()
    m.fs.crystallizer.initialize()
    m.fs.crystallizer.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_watertap_solver(linear_solver="mumps", max_iter=50)
    result = solver.solve(m, tee=True)

    m.fs.crystallizer.report()
    m.fs.crystallizer.crystallizer.temperature_operating.fix(315)
    m.fs.crystallizer.crystallizer.evaporation_percent.fix(90)

    result = solver.solve(m, tee=True)

    m.fs.crystallizer.report()
    print("------vars_close_to_bound-tests---------")
    print_variables_close_to_bounds(m)
    print("------constraints_close_to_bound-tests---------")
    print_constraints_close_to_bounds(m)
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
    print(large_residuals_set(m, tol=1e-8, return_residual_values=True))
    for lr in large_residuals_set(m):
        print(lr)

    assert_optimal_termination(result)
    assert (
        pytest.approx(
            m.fs.crystallizer.crystallizer.flow_mass_sol_comp_apparent["Calcite"].value,
            abs=1e-5,
        )
        == 0.00015684482
    )
    assert (
        pytest.approx(
            m.fs.crystallizer.crystallizer.pH["outlet"].value,
            abs=1e-2,
        )
        == 6.0309
    )
    for r in [95, 98]:
        m.fs.crystallizer.crystallizer.evaporation_percent.fix(r)

        result = solver.solve(m, tee=True)
        m.fs.crystallizer.report()

        assert_optimal_termination(result)
    assert (
        pytest.approx(
            m.fs.crystallizer.crystallizer.flow_mass_sol_comp_apparent["Calcite"].value,
            abs=1e-5,
        )
        == 0.00017666
    )
    assert (
        pytest.approx(
            m.fs.crystallizer.crystallizer.pH["outlet"].value,
            abs=1e-2,
        )
        == 5.6588
    )
