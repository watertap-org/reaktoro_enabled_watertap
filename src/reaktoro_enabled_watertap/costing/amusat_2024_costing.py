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

import idaes.core.util.scaling as iscale

##################################################################
# Costing from Amusat et al. (2024) "Cost-optimal selection of pH control for mineral scaling
# prevention in high recovery reverse osmosis desalination"
##################################################################


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
        )  # 1e1 stable for validation
        # iscale.set_scaling_factor(blk.capital_cost, 1 / 1e1)
        # iscale.constraint_scaling_transform(blk.capital_cost_constraint, 1 / 1e1)
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
            == blk.cost_factor
            * (
                blk.capital_cost_acid_addition["a"] * volume_acid**2
                + blk.capital_cost_acid_addition["b"] * volume_acid
                + blk.capital_cost_acid_addition["c"]
            ),
        )
        # # 1e4 stable for validation
        # iscale.set_scaling_factor(blk.capital_cost, 1 / 1e4)
        # iscale.constraint_scaling_transform(blk.capital_cost_constraint, 1 / 1e4)


def build_high_pressure_pump_cost_param_block(blk):

    blk.pump_cost = Param(
        initialize=700,
        mutable=True,
        doc="High pressure pump cost",
        units=pyunits.USD_2018 / pyunits.kW,
    )


def cost_high_pressure_pump(blk, cost_electricity_flow=True):
    """
    High pressure pump costing method

    `TODO: describe equations`

    Args:
        cost_electricity_flow (bool): if True, the Pump's work_mechanical will
            be converted to kW and costed as an electricity. Defaults to True.
    """
    t0 = blk.flowsheet().time.first()
    make_capital_cost_var(blk)
    build_high_pressure_pump_cost_param_block(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyunits.convert(
            blk.pump_cost
            * pyunits.convert(blk.unit_model.work_mechanical[t0], pyunits.W),
            to_units=blk.costing_package.base_currency,
        )
    )
    if cost_electricity_flow:
        # grab lower bound of mechanical work
        lb = blk.unit_model.work_mechanical[t0].lb
        # set lower bound to 0 to avoid negative defined flow warning when lb is not >= 0
        blk.unit_model.work_mechanical.setlb(0)
        blk.costing_package.cost_flow(
            pyunits.convert(blk.unit_model.work_mechanical[t0], to_units=pyunits.kW),
            "electricity",
        )
        # set lower bound back to its original value that was assigned to lb
        blk.unit_model.work_mechanical.setlb(lb)
