from reaktoro_enabled_watertap.property_models.mcas_with_enthalpy import (
    MCASWEParameterBlock,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    ActivityCoefficientModel,
    DensityCalculation,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_feed_unit import (
    MultiCompFeed,
)
from pyomo.environ import ConcreteModel
from idaes.core import (
    FlowsheetBlock,
)

from pyomo.environ import (
    assert_optimal_termination,
)


from reaktoro_enabled_watertap.water_sources.source_water_importer import (
    get_source_water_data,
)
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
import pytest


def build_case(water, reconcile_using_reaktoro=False):
    mcas_props, feed_specs = get_source_water_data(f"{water}.yaml")
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.seawater

    m.fs.properties = MCASWEParameterBlock(**mcas_props)
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        reconcile_using_reaktoro=reconcile_using_reaktoro,
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
