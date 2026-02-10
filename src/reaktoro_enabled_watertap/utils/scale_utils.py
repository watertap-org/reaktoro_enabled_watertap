from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale
from pyomo.environ import value


def calculate_scale_from_dependent_vars(var, constraint, dependent_vars):
    """
    Calculate the scale factor for a variable based on the scales of other variables
    in a constraint.

    Args:
        var: The variable to calculate the scale for.
        constraint: The constraint that relates the variable to other variables.
        dependent_vars: A list of other variables in the constraint that have known scales.

    Returns:
        The calculated scale factor for the variable.
    """
    if isinstance(dependent_vars, list) == False:
        dependent_vars = [dependent_vars]
    initial_values = [v.value for v in dependent_vars]
    fixed_states = [v.fixed for v in dependent_vars]
    # Get the current scale factors for the dependent variables
    scales = [iscale.get_scaling_factor(v) for v in dependent_vars]
    # Calculate the value of the variable using the constraint
    for v, scale in zip(dependent_vars, scales):
        if scale is None:
            raise ValueError(
                f"Cannot calculate scale for variable {v.name} with unknown scale factor"
            )
        v.fix(1 / scale)
        # print("Fixed ", v.name, 1 / scale)
    calculate_variable_from_constraint(var, constraint)
    scale_var = var.value
    if scale_var is None or scale_var == 0:
        raise ValueError(
            f"Cannot calculate scale for variable {var.name} with value {scale_var}"
        )
    # Use the average of the scales of the dependent variables as a basis
    # print("Applying scale to ", var.name, constraint.name, "of", 1 / scale_var)
    scale = 1 / scale_var
    if iscale.get_scaling_factor(var) is None:
        iscale.set_scaling_factor(var, scale)
    else:
        print("Skiping, as scaling far already existsf for ", var.name)
    if iscale.get_constraint_transform_applied_scaling_factor(constraint) is None:
        iscale.constraint_scaling_transform(constraint, scale)
    else:
        print(
            "Skiping, as scaling far already existsf for ",
            var.name,
            iscale.get_constraint_transform_applied_scaling_factor(constraint),
        )
    for v, (initial_values, fixed_states) in zip(
        dependent_vars, zip(initial_values, fixed_states)
    ):
        # print("Setting back ", v.name, initial_values, fixed_states)
        v.fix(initial_values)
        if not fixed_states:
            v.unfix()


def get_vars_from_expr(var_list, expr):
    if isinstance(expr, (int, float)) == False and expr.is_expression_type():
        for arg in expr.args:
            get_arg = get_vars_from_expr(var_list, arg)
            if get_arg is not None:
                var_list.append(get_arg)
    else:
        return expr


def get_scale_from_expr(expr):
    var_list = []
    get_vars_from_expr(var_list, expr)
    found_vars = []
    preset_values = []
    preset_states = []
    for v in var_list:
        if isinstance(v, (int, float)) == False and v.is_variable_type():
            found_vars.append(v)
            preset_values.append(v.value)
            preset_states.append(v.fixed)
            scale = iscale.get_scaling_factor(v)
            if scale is None or scale == 0:
                scale = 1
            v.fix(1 / scale)
    expr_scale = value(expr)
    for v, val, state in zip(found_vars, preset_values, preset_states):
        v.fix(val)
        if not state:
            v.unfix()
    return 1 / expr_scale


def scale_costing_block(costing_block):
    """
    Scale the costing block based on the registered unit costing blocks
    This requires that all capital costs and operating costs are scaled, including
    flow costs.
    Args:
        costing_block: The costing block to scale.
    """
    capital_costing_factor = 0
    operating_costing_factor = 0
    variable_costing_factor = 0
    flow_cost_types = {}
    flow_costs = {}
    for unit in costing_block._registered_unit_costing:
        if hasattr(unit, "capital_cost"):
            capital_costing_factor += iscale.get_scaling_factor(unit.capital_cost)
        if hasattr(unit, "fixed_operating_cost"):
            operating_costing_factor += iscale.get_scaling_factor(
                unit.fixed_operating_cost
            )
        if hasattr(unit, "variable_operating_cost"):
            variable_costing_factor += iscale.get_scaling_factor(
                unit.variable_operating_cost
            )
    for ftype in costing_block.used_flows:
        flow_cost_types[ftype] = 0
        for flow in costing_block._registered_flows[ftype]:
            if flow.is_variable_type():
                scale = iscale.get_scaling_factor(flow)
            else:
                scale = get_scale_from_expr(flow)
            flow_cost_types[ftype] += scale
        cost_value = getattr(costing_block, f"{ftype}_cost")
        cost_scale = iscale.get_scaling_factor(cost_value)
        if cost_scale is None:
            cost_scale = 1 / value(cost_value)
            iscale.set_scaling_factor(cost_value, cost_scale)
        flow_costs[ftype] = cost_scale

    if capital_costing_factor == 0:
        capital_costing_factor = 1
    if operating_costing_factor == 0:
        operating_costing_factor = 1
    if variable_costing_factor == 0:
        variable_costing_factor = 1
    operating_costing_factor = 1 / (
        1 / operating_costing_factor + 1 / capital_costing_factor * 0.03
    )  # account for small fixed cost due to capital cost
    print("Capital cost scaling factor: ", capital_costing_factor)
    print("Operating cost scaling factor: ", operating_costing_factor)
    print("Variable cost scaling factor: ", variable_costing_factor)
    print("Flow cost types: ", flow_cost_types)
    # iscale.set_scaling_factor(costing_block.electricity_cost, 100)

    iscale.set_scaling_factor(costing_block.utilization_factor, 1)
    iscale.set_scaling_factor(costing_block.maintenance_labor_chemical_factor, 100)

    iscale.set_scaling_factor(costing_block.total_investment_factor, 1)
    iscale.set_scaling_factor(costing_block.plant_lifetime, 1)
    iscale.set_scaling_factor(costing_block.wacc, 10)

    iscale.set_scaling_factor(costing_block.capital_recovery_factor, 10)
    iscale.set_scaling_factor(costing_block.TPEC, 1)

    iscale.set_scaling_factor(costing_block.TIC, 1)
    # total capital costs
    iscale.set_scaling_factor(
        costing_block.aggregate_capital_cost, capital_costing_factor
    )
    iscale.constraint_scaling_transform(
        costing_block.aggregate_capital_cost_constraint, capital_costing_factor
    )

    iscale.set_scaling_factor(costing_block.total_capital_cost, capital_costing_factor)
    iscale.constraint_scaling_transform(
        costing_block.total_capital_cost_constraint, capital_costing_factor
    )

    # operating costs
    iscale.set_scaling_factor(
        costing_block.aggregate_fixed_operating_cost, operating_costing_factor
    )
    iscale.constraint_scaling_transform(
        costing_block.aggregate_fixed_operating_cost_constraint,
        operating_costing_factor,
    )
    iscale.set_scaling_factor(
        costing_block.total_operating_cost, operating_costing_factor
    )
    iscale.constraint_scaling_transform(
        costing_block.total_operating_cost_constraint, operating_costing_factor
    )

    iscale.set_scaling_factor(
        costing_block.aggregate_variable_operating_cost, variable_costing_factor
    )

    iscale.constraint_scaling_transform(
        costing_block.aggregate_variable_operating_cost_constraint,
        variable_costing_factor,
    )
    iscale.constraint_scaling_transform(
        costing_block.capital_recovery_factor_constraint,
        1,
    )
    for ftype, scale in flow_cost_types.items():
        if "electricity" in ftype:
            mc = 1
        else:
            # this assume our mass flow is based on seconds! This needs to be refined!
            mc = 1 / 3.15576e7
        iscale.set_scaling_factor(
            costing_block.aggregate_flow_costs[ftype],
            scale * flow_costs[ftype] * mc,
        )

        iscale.constraint_scaling_transform(
            costing_block.aggregate_flow_costs_constraint[ftype],
            scale * flow_costs[ftype] * mc,
        )
        agg_flow_cost = costing_block.find_component(f"aggregate_flow_{ftype}")
        agg_flow_constraint = costing_block.find_component(
            f"aggregate_flow_{ftype}_constraint"
        )
        iscale.set_scaling_factor(agg_flow_cost, scale)
        iscale.constraint_scaling_transform(agg_flow_constraint, scale)
