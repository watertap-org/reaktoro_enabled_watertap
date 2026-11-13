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

from pyomo.environ import Constraint, assert_optimal_termination
from pyomo.core.expr.visitor import identify_variables
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.core.solvers import get_solver

__author__ = "Alexander V. Dudchenko"


def init_translator_block(block, additional_vars_to_fix=None):
    """Initializes a translator block.

    Translator equations are explicit mappings from inlet to outlet, so each
    outlet variable is computed directly from its translator constraint before
    solving, giving the solver a near-converged starting point.

    Args:
        block: Translator block to initialize
        additional_vars_to_fix: variable (or list of variables) outside the
            translator block that appear in its constraints (e.g. upstream
            composition or pH). These are fixed during initialization and
            unfixed afterwards.
    """
    if additional_vars_to_fix is not None:
        if not isinstance(additional_vars_to_fix, list):
            additional_vars_to_fix = [additional_vars_to_fix]
        for var in additional_vars_to_fix:
            var.fix()
    flags = fix_state_vars(block.properties_in)
    block.properties_in.initialize()
    for con in block.component_data_objects(
        Constraint, active=True, descend_into=False
    ):
        free_vars = [v for v in identify_variables(con.body) if not v.fixed]
        if len(free_vars) == 1:
            calculate_variable_from_constraint(free_vars[0], con)
    # Only initialize outlet properties if it is square with its state vars fixed,
    # outlet can be over-specified if translator fixes non-state outlet vars
    # (e.g. MVC distillate solute flows), translator solve below still computes
    # outlet properties in that case
    out_flags = fix_state_vars(block.properties_out)
    if degrees_of_freedom(block.properties_out) == 0:
        block.properties_out.initialize()
    revert_state_vars(block.properties_out, out_flags)
    solver = get_solver()
    results = solver.solve(block)
    assert_optimal_termination(results)
    revert_state_vars(block.properties_in, flags)
    if additional_vars_to_fix is not None:
        for var in additional_vars_to_fix:
            var.unfix()
