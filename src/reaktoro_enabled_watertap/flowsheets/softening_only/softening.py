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

from reaktoro_enabled_watertap.unit_models.multi_comp_feed_unit import (
    MultiCompFeed,
)
from reaktoro_enabled_watertap.water_sources.source_water_importer import (
    get_source_water_data,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_product_unit import (
    MultiCompProduct,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_ph_mixer_unit import (
    MixerPhUnit,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_erd_unit import (
    MultiCompERDUnit,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_pump_unit import (
    MultiCompPumpUnit,
)
from reaktoro_enabled_watertap.unit_models.multi_comp_ro_unit import (
    MultiCompROUnit,
)
from reaktoro_enabled_watertap.unit_models.precipitation_unit import (
    PrecipitationUnit,
)
from reaktoro_enabled_watertap.unit_models.chemical_addition_unit import (
    ChemicalAdditionUnit,
)

from reaktoro_enabled_watertap.utils import ipopt_performance_utils as ipopt_perf_utils
from pyomo.environ import (
    TransformationFactory,
    units as pyunits,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
)
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.costing import WaterTAPCosting
from pyomo.environ import ConcreteModel, Var, Reals, Constraint, Objective
from idaes.core import (
    FlowsheetBlock,
)
import idaes.core.util.scaling as iscale

from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
from pyomo.environ import (
    assert_optimal_termination,
)
from reaktoro_enabled_watertap.utils.report_util import (
    build_report_table,
)
from idaes.core.util.model_statistics import degrees_of_freedom

from reaktoro_enabled_watertap.utils import scale_utils as scu
from watertap.core.util.model_diagnostics.infeasible import *

from reaktoro_enabled_watertap.utils.report_util import get_lib_path
from reaktoro_enabled_watertap.costing import (
    amusat_2024_costing as ams,
)
import numpy as np

__author__ = "Alexander V. Dudchenko"


def main():
    m = build_model(water_case="USDA_brackish.yaml", softening_reagents=["NaOH"])
    initialize_model(m)
    for ph in np.arange(9.0, 13.5, 0.5):
        m.fs.softening_unit.precipitation_reactor.pH["outlet"].fix(ph)
        solve_model(m, limited_memory=False)
        m.fs.softening_unit.report()


def build_model(
    water_case,
    rkt_hessian_type="LBFGS",
    bfgs_initialization_type="GaussNewton",
    softening_reagents=["Na2CO3", "CaO"],
):
    mcas_props, feed_specs = get_source_water_data(water_case)
    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m = ConcreteModel()
    m.water_case = water_case
    rkt_options = {
        "hessian_options": {
            "hessian_type": rkt_hessian_type,
            "bfgs_initialization_type": bfgs_initialization_type,
        }
    }
    m.fs = FlowsheetBlock()
    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.ro_properties = SeawaterParameterBlock()
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        reconcile_using_reaktoro=True,
        **feed_specs,
    )
    m.fs.softening_unit = PrecipitationUnit(
        default_property_package=m.fs.properties,
        selected_precipitants=[
            "Calcite",
            "Brucite",
        ],
        selected_reagents=softening_reagents,
        add_alkalinity=True,
        reaktoro_options=rkt_options,
    )

    m.fs.feed.outlet.connect_to(m.fs.softening_unit.inlet)

    m.fs.feed.fix_and_scale()
    m.fs.softening_unit.fix_and_scale()
    TransformationFactory("network.expand_arcs").apply_to(m)

    iscale.calculate_scaling_factors(m)
    return m


def initialize_model(m):
    m.fs.feed.initialize()
    m.fs.softening_unit.initialize()
    m.fs.softening_unit.set_optimization_operation()


def solve_model(m, tee=True, linear_solver="mumps", **kwargs):
    if linear_solver == "mumps":
        pivtol = 1e-3
        maxpivtol = 1e0
    else:
        pivtol = None
        maxpivtol = None
    solver = get_cyipopt_watertap_solver(
        linear_solver=linear_solver,
        max_iter=1000,
        limited_memory=False,
        pivtol=pivtol,
        pivtolmax=maxpivtol,
    )

    result = solver.solve(m, tee=tee)
    if tee:
        print("------vars_close_to_bound-tests---------")
        print_variables_close_to_bounds(m)
        print("------constraints_close_to_bound-tests---------")
        print_constraints_close_to_bounds(m)

    assert_optimal_termination(result)
    return result


if __name__ == "__main__":
    main()
