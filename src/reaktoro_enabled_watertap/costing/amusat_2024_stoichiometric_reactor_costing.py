#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

# Import Pyomo libraries
from pyomo.environ import (
    Param,
    Var,
    Constraint,
    units as pyunits,
)
from watertap.costing.util import (
    make_capital_cost_var,
)


def build_stoichiometric_reactor_cost_param_block(blk):
    blk.capital_cost_softening = Param(
        ["const", "exp"],
        initialize=1,
        units=pyunits.dimensionless,
        mutable=True,
        doc="Cost for softening from Amusat et al. (2024)",
    )
    blk.capital_cost_softening["const"] = 12985
    blk.capital_cost_softening["exp"] = 0.5901
    blk.capital_cost_acid_addition = Param(
        ["a", "b", "c"],
        initialize=127.8,
        units=pyunits.dimensionless,
        doc="Cost for acid mixer from Amusat et al. (2024)",
    )
    blk.capital_cost_acid_addition["a"] = -0.0029
    blk.capital_cost_acid_addition["b"] = 48.434
    blk.capital_cost_acid_addition["c"] = 22648


def cost_stoichiometric_reactor(blk):
    build_stoichiometric_reactor_cost_param_block(blk)
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    # if softening reactor
    if (
        blk.unit_model.has_precipitation_reaction
        and blk.unit_model.has_dissolution_reaction
    ):
        mass_softening_reagents = sum(
            pyunits.convert(
                obj,
                to_units=pyunits.lb / pyunits.day,
            )
            for reagent, obj in blk.unit_model.flow_mass_reagent.items()
        )
        blk.capital_cost_constraint = Constraint(
            expr=blk.capital_cost
            == blk.cost_factor
            * blk.capital_cost_softening["const"]
            * mass_softening_reagents ** blk.capital_cost_softening["exp"],
        )
    # if acid addition reactor
    elif blk.unit_model.has_dissolution_reaction:
        volume_acid = sum(
            pyunits.convert(
                obj,
                to_units=pyunits.gallon / pyunits.day,
            )
            for reagent, obj in blk.unit_model.flow_vol_reagent.items()
        )
        blk.capital_cost_constraint = Constraint(
            expr=blk.capital_cost
            == blk.cost_factor * blk.capital_cost_acid_addition["a"] * volume_acid**2
            + blk.capital_cost_acid_addition["b"] * volume_acid
            + blk.capital_cost_acid_addition["c"],
        )
