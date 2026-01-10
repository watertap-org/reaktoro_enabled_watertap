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
from reaktoro_enabled_watertap.unit_models.multi_comp_erd_unit import (
    MultiCompERDUnit,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_pump_unit import (
    MultiCompPumpUnit,
)
from reaktoro_enabled_watertap.unit_models.tests.test_multi_comp_feed_product import (
    build_case,
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

from watertap.costing import WaterTAPCosting

__author__ = "Alexander V. Dudchenko"


@pytest.mark.component
def test_erd_pressure():
    m = build_case("USDA_brackish", True)
    m.fs.costing = WaterTAPCosting()
    m.fs.pump_unit = MultiCompPumpUnit(
        default_property_package=m.fs.properties,
        initialization_pressure="osmotic_pressure",
        default_costing_package=m.fs.costing,
    )
    m.fs.pump_unit.fix_and_scale()
    m.fs.ERD_unit = MultiCompERDUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
    )
    m.fs.ERD_unit.fix_and_scale()
    m.fs.feed.outlet.connect_to(m.fs.pump_unit.inlet)
    m.fs.pump_unit.outlet.connect_to(m.fs.ERD_unit.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.costing.cost_process()
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.pump_unit.initialize()
    m.fs.ERD_unit.initialize()
    m.fs.ERD_unit.report()
    m.fs.ERD_unit.report(use_default_units=True)
    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_watertap_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.ERD_unit.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.ERD_unit.inlet.port.pressure[0].value,
            1e-5,
        )
        == 6.1710e05
    )
    assert (
        pytest.approx(
            m.fs.ERD_unit.outlet.port.pressure[0].value,
            1e-5,
        )
        == 101325
    )
    assert (
        pytest.approx(
            m.fs.ERD_unit.ERD.costing.capital_cost.value,
            1e-1,
        )
        == 3852.0
    )
