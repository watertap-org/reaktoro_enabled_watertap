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
from reaktoro_enabled_watertap.unit_models.multi_comp_feed_unit import (
    MultiCompFeed,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_ph_splitter import (
    SplitterPhUnit,
)
from reaktoro_enabled_watertap.water_sources.source_water_importer import (
    get_source_water_data,
)
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
from pyomo.environ import ConcreteModel
from idaes.core import (
    FlowsheetBlock,
)

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
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

__author__ = "Alexander V. Dudchenko"


@pytest.mark.core
@pytest.mark.component
def test_splitter():
    mcas_props, USDA_feed_specs = get_source_water_data(f"USDA_brackish.yaml")

    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        reconcile_using_reaktoro=True,
        **USDA_feed_specs,
    )
    m.fs.feed.fix_and_scale()

    m.fs.splitter = SplitterPhUnit(
        default_property_package=m.fs.properties,
        outlet_ports=["feed_a", "feed_b"],
        splitter_initialization_guess=0.4,
    )
    m.fs.splitter.fix_and_scale()
    m.fs.feed.outlet.connect_to(m.fs.splitter.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    assert degrees_of_freedom(m) == 0
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()

    m.fs.splitter.initialize()
    m.fs.splitter.report()
    assert degrees_of_freedom(m) == 0
    solver = get_cyipopt_watertap_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.splitter.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.splitter.splitter.pH.value,
            1e-5,
        )
        == 7.07
    )
    assert (
        pytest.approx(
            m.fs.splitter.splitter.pH.value,
            1e-5,
        )
        == 7.07
    )
    assert (
        pytest.approx(
            m.fs.splitter.splitter.feed_a_state[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .value,
            1e-5,
        )
        == 0.39863
    )
    assert (
        pytest.approx(
            m.fs.splitter.splitter.feed_b_state[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .value,
            1e-5,
        )
        == 0.59795
    )
    m.fs.splitter.set_optimization_operation()
    for port in m.fs.splitter.config.outlet_ports:
        assert m.fs.splitter.splitter.split_fraction[0, port].fixed == False


@pytest.mark.core
@pytest.mark.component
def test_splitter_3_ports():
    mcas_props, USDA_feed_specs = get_source_water_data(f"USDA_brackish.yaml")

    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        reconcile_using_reaktoro=True,
        **USDA_feed_specs,
    )
    m.fs.feed.fix_and_scale()

    m.fs.splitter = SplitterPhUnit(
        default_property_package=m.fs.properties,
        outlet_ports=["feed_a", "feed_b", "port_c"],
        splitter_initialization_guess={"feed_a": 0.4, "feed_b": 0.4},
    )
    m.fs.splitter.fix_and_scale()
    m.fs.feed.outlet.connect_to(m.fs.splitter.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    assert degrees_of_freedom(m) == 0
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()

    m.fs.splitter.initialize()
    m.fs.splitter.report()
    assert degrees_of_freedom(m) == 0
    solver = get_cyipopt_watertap_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.splitter.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.splitter.splitter.pH.value,
            1e-3,
        )
        == 7.07
    )
    assert (
        pytest.approx(
            m.fs.splitter.splitter.pH.value,
            1e-3,
        )
        == 7.07
    )
    assert (
        pytest.approx(
            m.fs.splitter.splitter.feed_a_state[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .value,
            1e-3,
        )
        == 0.39863
    )
    assert (
        pytest.approx(
            m.fs.splitter.splitter.feed_b_state[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .value,
            1e-3,
        )
        == 0.39863
    )
    assert (
        pytest.approx(
            m.fs.splitter.splitter.port_c_state[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .value,
            1e-3,
        )
        == 0.19932
    )
    m.fs.splitter.set_optimization_operation()
    for port in m.fs.splitter.config.outlet_ports:
        assert m.fs.splitter.splitter.split_fraction[0, port].fixed == False
