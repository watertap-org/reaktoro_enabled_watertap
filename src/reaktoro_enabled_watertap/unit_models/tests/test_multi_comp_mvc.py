import pytest
from reaktoro_enabled_watertap.unit_models.multi_comp_mvc import (
    MultiCompMVC,
)

from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
from pyomo.environ import ConcreteModel, assert_optimal_termination
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import (
    FlowsheetBlock,
)
import idaes.core.util.scaling as iscale
from pyomo.environ import (
    Objective,
    check_optimal_termination,
    TransformationFactory,
    units as pyunits,
)
from watertap.costing import WaterTAPCosting

from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap.core.util.model_diagnostics.infeasible import *

from reaktoro_enabled_watertap.utils import scale_utils as scu

from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    ActivityCoefficientModel,
    DensityCalculation,
)
import watertap.property_models.seawater_prop_pack as props_sw
from reaktoro_enabled_watertap.water_sources.source_water_importer import (
    get_source_water_data,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_feed_unit import (
    MultiCompFeed,
)

__author__ = "Alexander V. Dudchenko"


def build_case(water, reconcile_using_reaktoro=False):
    mcas_props, feed_specs = get_source_water_data(f"{water}.yaml")
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.seawater

    m.fs.properties = MCASParameterBlock(**mcas_props)
    # m.fs.properties = MCASWEParameterBlock(**mcas_props)
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        reconcile_using_reaktoro=reconcile_using_reaktoro,
        **feed_specs,
    )
    m.fs.feed.fix_and_scale()
    m.fs.feed.report()
    return m


def build_nacl_case():
    mcas_props = {
        "solute_list": ["Na_+", "Cl_-"],
        "mw_data": {
            "H2O": 18e-3,
            "Na_+": 23e-3,
            "Cl_-": 35.5e-3,
        },
        "elec_mobility_data": {
            ("Liq", "Na_+"): 5.19e-8,
            ("Liq", "Cl_-"): 7.92e-8,
        },
        "charge": {"Na_+": 1, "Cl_-": -1},
        "diffusivity_data": {
            ("Liq", "Na_+"): 1.33e-9,
            ("Liq", "Cl_-"): 2.03e-9,
        },
    }
    feed_specs = {
        "temperature": 298.15,
        "pressure": 101325,
        "volumetric_flowrate": 1 * pyunits.L / pyunits.s,
        "ion_concentrations": {
            "Na_+": 15 * pyunits.g / pyunits.L,
            "Cl_-": 15 * pyunits.g / pyunits.L,
        },  # "H_+": 1e-7, "OH_-": 1e-7},
    }
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        reconcile_using_reaktoro=False,
        **feed_specs,
    )
    m.fs.feed.fix_and_scale()
    m.fs.feed.report()
    return m


@pytest.mark.core
@pytest.mark.component
def test_enthalpy_term():
    m = build_case("USDA_brackish", True)
    m.fs.feed.feed.properties[0].get_enthalpy_flow_terms(["Liq"])
    solver = get_cyipopt_watertap_solver()
    r = solver.solve(m.fs.feed, tee=True)
    assert_optimal_termination(r)
    m.fs.feed.feed.properties[0].display()
    assert m.fs.feed.feed.properties[0].enth_mass_phase["Liq"].value == pytest.approx(
        83554.43, rel=1e-1
    )


@pytest.mark.core
@pytest.mark.component
def test_mvc_no_reaktoro_sea_water_prop_with_translators():
    m = build_nacl_case()

    m.fs.water_properties_vapor = WaterParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.sea_water_props = props_sw.SeawaterParameterBlock()
    m.fs.MVC = MultiCompMVC(
        default_property_package=m.fs.properties,
        mvc_property_package=m.fs.sea_water_props,
        default_costing_package=m.fs.costing,
        vapor_property_package=m.fs.water_properties_vapor,
        add_reaktoro_chemistry=False,
    )

    m.fs.MVC.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.MVC.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.MVC.tb_distillate.properties_out[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.MVC.tb_distillate.properties_out[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.MVC.tb_distillate.properties_out[0].flow_vol
    )

    iscale.calculate_scaling_factors(m)

    scu.scale_costing_block(m.fs.costing)
    print("Degrees of freedom before initialization: ", degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    m.fs.feed.initialize()
    m.fs.MVC.initialize()
    assert degrees_of_freedom(m) == 0

    m.fs.MVC.set_optimization_operation()

    print("Degrees of freedom optimization: ", degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 5
    m.fs.MVC.recovery.fix()
    assert degrees_of_freedom(m) == 4
    m.fs.cost_objective = Objective(expr=m.fs.costing.LCOW)
    solver = get_cyipopt_watertap_solver(linear_solver="mumps", max_iter=3000)
    for r in [0.5, 0.6, 0.7, 0.8]:

        m.fs.MVC.recovery.fix(r)
        result = solver.solve(m, tee=True)
        print(f"Water removal: {r}")
        m.fs.MVC.report()

        assert check_optimal_termination(result)


@pytest.mark.core
@pytest.mark.component
def test_mvc_no_reaktoro_mcas():
    m = build_nacl_case()

    m.fs.water_properties_vapor = WaterParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.MVC = MultiCompMVC(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        vapor_property_package=m.fs.water_properties_vapor,
        add_reaktoro_chemistry=False,
    )

    m.fs.MVC.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.MVC.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.MVC.tb_distillate.properties_out[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.MVC.tb_distillate.properties_out[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.MVC.tb_distillate.properties_out[0].flow_vol
    )

    iscale.calculate_scaling_factors(m)

    scu.scale_costing_block(m.fs.costing)

    print("Degrees of freedom before initialization: ", degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    m.fs.feed.initialize()
    m.fs.MVC.initialize()

    m.fs.MVC.set_optimization_operation()
    print("Degrees of freedom optimization: ", degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 5
    m.fs.MVC.recovery.fix()
    assert degrees_of_freedom(m) == 4
    m.fs.cost_objective = Objective(expr=m.fs.costing.LCOW)
    solver = get_cyipopt_watertap_solver(linear_solver="mumps", max_iter=3000)
    for r in [0.5, 0.6, 0.7, 0.8]:

        m.fs.MVC.recovery.fix(r)
        result = solver.solve(m, tee=True)
        print(f"Water removal: {r}")
        m.fs.MVC.report()
        assert check_optimal_termination(result)


@pytest.mark.core
@pytest.mark.component
def test_mvc_with_reaktoro_sea_water_prop_with_translators():
    m = build_case("Seawater", True)

    m.fs.water_properties_vapor = WaterParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.sea_water_props = props_sw.SeawaterParameterBlock()
    m.fs.MVC = MultiCompMVC(
        default_property_package=m.fs.properties,
        mvc_property_package=m.fs.sea_water_props,
        default_costing_package=m.fs.costing,
        vapor_property_package=m.fs.water_properties_vapor,
        add_reaktoro_chemistry=True,
    )

    m.fs.MVC.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.MVC.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.MVC.tb_distillate.properties_out[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.MVC.tb_distillate.properties_out[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.MVC.tb_distillate.properties_out[0].flow_vol
    )

    iscale.calculate_scaling_factors(m)

    scu.scale_costing_block(m.fs.costing)
    print("Degrees of freedom before initialization: ", degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    m.fs.feed.initialize()
    m.fs.MVC.initialize()
    m.fs.MVC.set_optimization_operation()
    m.fs.MVC.deactivate_scaling_constraints()
    assert degrees_of_freedom(m) == 5
    m.fs.MVC.recovery.fix()
    assert degrees_of_freedom(m) == 4
    m.fs.cost_objective = Objective(expr=m.fs.costing.LCOW)
    for r in [0.5, 0.8]:

        m.fs.MVC.recovery.fix(r)
        solver = get_cyipopt_watertap_solver(linear_solver="mumps", max_iter=3000)
        result = solver.solve(m, tee=True)
        print(f"Water removal: {r}")
        m.fs.MVC.report()
        assert check_optimal_termination(result)


@pytest.mark.core
@pytest.mark.component
@pytest.mark.parametrize("r", [0.5, 0.6, 0.7, 0.8])
def test_mvc_with_reaktoro_mcas(r):
    m = build_case("Seawater", True)

    m.fs.water_properties_vapor = WaterParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.MVC = MultiCompMVC(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        vapor_property_package=m.fs.water_properties_vapor,
        add_reaktoro_chemistry=True,
    )

    m.fs.MVC.fix_and_scale()
    m.fs.feed.outlet.connect_to(m.fs.MVC.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.MVC.tb_distillate.properties_out[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.MVC.tb_distillate.properties_out[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.MVC.tb_distillate.properties_out[0].flow_vol
    )

    iscale.calculate_scaling_factors(m)
    scu.scale_costing_block(m.fs.costing)

    assert degrees_of_freedom(m) == 0
    m.fs.feed.initialize()
    m.fs.MVC.initialize()
    m.fs.MVC.set_optimization_operation()
    m.fs.MVC.deactivate_scaling_constraints()

    assert degrees_of_freedom(m) == 5
    m.fs.MVC.recovery.fix(r)
    assert degrees_of_freedom(m) == 4

    m.fs.cost_objective = Objective(expr=m.fs.costing.LCOW)

    solver = get_cyipopt_watertap_solver(linear_solver="mumps", max_iter=3000)
    result = solver.solve(m, tee=True)

    print(f"Water removal: {r}")
    m.fs.MVC.report()

    assert check_optimal_termination(result)
