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
import os
import tempfile

from reaktoro_enabled_watertap.utils.report_util import get_lib_path
from reaktoro_enabled_watertap.costing import (
    amusat_2024_costing as ams,
)

__author__ = "Alexander V. Dudchenko"


def main(
    multi_process_reaktoro=True, result_save_location="default", linear_solver="mumps"
):
    if result_save_location == "default":
        result_save_location = os.path.join(get_lib_path(), "flowsheets/results")

    for water in [
        "USDA_brackish.yaml",
        "sample_500_hardness.yaml",
        "sample_1500_hardness.yaml",
        "Seawater.yaml",
    ]:
        if (
            water == "sample_500_hardness.yaml" or water == "sample_1500_hardness.yaml"
        ) and linear_solver == "mumps":
            import time

            print(f"\n\n\n\n{water}, do not solve well with mumps!!!!!!!!!!!!\n\n\n\n")

            time.sleep(4)
        if "Seawater" in water:
            hpro = True
        else:
            hpro = False

        m = build_model(
            water,
            multi_process_reaktoro=multi_process_reaktoro,
            hpro=hpro,
            softening_reagents=["Na2CO3", "CaO"],
            acidification_reagents=["HCl", "H2SO4"],
            rkt_hessian_type="LBFGS",
            bfgs_initialization_type="GaussNewton",
            system_costing="Amusat_et_al_2024",
        )
        initialize(m, linear_solver=linear_solver, tee=True)
        for r in [60, 70, 80]:
            m.fs.water_recovery.fix(r / 100)
            print(f"\n\n------------Solving for water recovery: {r}%------------")
            solve_model(m, linear_solver=linear_solver, tee=True)
            report_all_units(m)
        if m.find_component("reaktoro_manager") is not None:
            m.reaktoro_manager.terminate_workers()


def enable_multi_process_reaktoro(
    m, rkt_hessian_type="LBFGS", bfgs_initialization_type="GaussNewton"
):
    """Enables use of parallel solves for reaktoro blocks,
    in RO mode there will be 3 reaktoro blocks
        1 for Softening,
        1 for Acidification,
        1 for RO
        requiring 3 logical cores
    in HPRO mode there will be 4 reaktoro blocks
        1 for Softening,
        1 for Acidification,
        1 for RO,
        1 for HPRO
        requiring 4 logical cores
    """
    from reaktoro_pse.parallel_tools.reaktoro_block_manager import (
        ReaktoroBlockManager,
    )

    rkt_options = {}
    opt = {
        "hessian_type": rkt_hessian_type,
        "bfgs_initialization_type": bfgs_initialization_type,
    }
    if bfgs_initialization_type == "constant":
        opt["bfgs_init_const_hessian_value"] = 1e-16

    m.reaktoro_manager = ReaktoroBlockManager(
        hessian_options=opt,
    )

    rkt_options["reaktoro_block_manager"] = m.reaktoro_manager
    return rkt_options


def build_model(
    water_case,
    multi_process_reaktoro=True,
    hpro=False,
    rkt_hessian_type="LBFGS",
    bfgs_initialization_type="GaussNewton",
    softening_reagents=["Na2CO3", "CaO"],
    acidification_reagents=["HCl", "H2SO4"],
    feed_flow_rate=5000 * pyunits.m**3 / pyunits.day,
    system_costing="watertap_default",
):
    """Builds the flowsheet model for the softening-acidification-RO process.
    Args:
        water_case (str): location of water source data file with yaml format.
        multi_process_reaktoro (bool): If True, enables parallel processing for reaktoro blocks.
        hpro (bool): If True, includes high-pressure RO unit in the flowsheet.
        rkt_hessian_type (str): Type of Hessian to use in reaktoro blocks.
            Options are:
                - limited-memory: Uses limited memory BFGS Hessian approximation within ipopt solver.
                - ZeroHessian - no hessian
                - GaussNewton - Naive Gauss-Newton Hessian approximation (Jacobian^T * Jacobian)
                - LBFGS - Limited Memory Broyden-Fletcher-Goldfarb-Shanno
                - BFGS - Broyden-Fletcher-Goldfarb-Shanno
                - CBFGS - conditional BFGS
                - BFGS_mod - modified BFGS
                - BFGS_damp - damped BFGS
                - BFGS_ipopt - BFGS with ipopt update step
        softening_reagents (list): List of reagents to use in the softening unit.
        acidification_reagents (list): List of reagents to use in the acidification unit.
        feed_flow_rate: volumetric flow rate of the feed water to the system.
        system_costing (str): Costing method for softening and acidification units. (watertap_default, or Amusat_et_al_2024)
    """

    mcas_props, feed_specs = get_source_water_data(water_case)
    if feed_flow_rate is not None:
        if isinstance(feed_flow_rate, dict):
            feed_specs["volumetric_flowrate"] = feed_flow_rate["value"] * getattr(
                pyunits, feed_flow_rate["units"]
            )
        else:
            feed_specs["volumetric_flowrate"] = feed_flow_rate
    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m = ConcreteModel()
    m.water_case = water_case
    if rkt_hessian_type == "limited-memory":
        rkt_hessian_type = "ZeroHessian"
        m.solver_limited_memory = True
    else:
        m.solver_limited_memory = False
    m.solver_limited_memory_scalar = bfgs_initialization_type
    rkt_options = {
        "hessian_options": {
            "hessian_type": rkt_hessian_type,
            "bfgs_initialization_type": bfgs_initialization_type,
        }
    }
    if bfgs_initialization_type == "constant":
        rkt_options["hessian_options"]["bfgs_init_const_hessian_value"] = 1e-16
    if multi_process_reaktoro:
        rkt_options = enable_multi_process_reaktoro(
            m, rkt_hessian_type, bfgs_initialization_type
        )

    m.fs = FlowsheetBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.ro_properties = SeawaterParameterBlock()
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        reconcile_using_reaktoro=True,
        **feed_specs,
    )
    if isinstance(softening_reagents, str):
        softening_reagents = [softening_reagents]
    if isinstance(acidification_reagents, str):
        acidification_reagents = [acidification_reagents]
    if system_costing == "watertap_default":
        chemical_costing_type = {}
        pump_costing_type = {}
    elif system_costing == "Amusat_et_al_2024":
        chemical_costing_type = {"costing_method": ams.cost_stoichiometric_reactor}
        pump_costing_type = {"costing_method": ams.cost_high_pressure_pump}
    else:
        raise (
            ValueError(
                f"Unknown costing method {system_costing} for softening and acidification units."
            )
        )
    m.fs.softening_unit = PrecipitationUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        selected_precipitants=[
            "Calcite",
            "Brucite",
        ],
        selected_reagents=softening_reagents,
        add_alkalinity=True,
        reaktoro_options=rkt_options,
        default_costing_package_kwargs=chemical_costing_type,
    )

    m.fs.acidification_unit = ChemicalAdditionUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        selected_reagents=acidification_reagents,
        reaktoro_options=rkt_options,
        default_costing_package_kwargs=chemical_costing_type,
    )
    if "Seawater" in water_case:
        ovp = 1.5
    else:
        ovp = 6.5
    m.fs.pump_unit = MultiCompPumpUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        initialization_pressure="osmotic_pressure",
        osmotic_over_pressure=ovp,
        maximum_pressure=85 * pyunits.bar,
        default_costing_package_kwargs=pump_costing_type,
    )

    m.fs.ro_unit = MultiCompROUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        ro_property_package=m.fs.ro_properties,
        selected_scalants={"Calcite": 1, "Gypsum": 1},
        use_interfacecomp_for_effluent_pH=True,
        reaktoro_options=rkt_options,
        target_recovery=0.5,
    )

    m.fs.erd_unit = MultiCompERDUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
    )

    if hpro:
        m.fs.hp_pump_unit = MultiCompPumpUnit(
            default_property_package=m.fs.properties,
            default_costing_package=m.fs.costing,
            initialization_pressure="osmotic_pressure",
            maximum_pressure=400 * pyunits.bar,
            default_costing_package_kwargs=pump_costing_type,
        )
        m.fs.hpro_unit = MultiCompROUnit(
            default_property_package=m.fs.properties,
            default_costing_package=m.fs.costing,
            ro_property_package=m.fs.ro_properties,
            selected_scalants={"Calcite": 1, "Gypsum": 1},
            reaktoro_options=rkt_options,
            use_interfacecomp_for_effluent_pH=True,
            default_costing_package_kwargs={
                "costing_method_arguments": {"ro_type": "high_pressure"}
            },
            target_recovery=0.3,
        )
        m.fs.product_mixer = MixerPhUnit(
            default_property_package=m.fs.properties,
            default_costing_package=m.fs.costing,
            inlet_ports=["ro_inlet", "hpro_inlet"],
            add_reaktoro_chemistry=False,
        )

    m.fs.product = MultiCompProduct(
        default_property_package=m.fs.properties,
    )
    m.fs.brine = MultiCompProduct(
        default_property_package=m.fs.properties,
    )
    if multi_process_reaktoro:
        m.reaktoro_manager.build_reaktoro_blocks()
    # to simplify initialization and fixing model state
    m.flowsheet_unit_order = []
    m.flowsheet_unit_order.append(m.fs.feed)
    m.flowsheet_unit_order.append(m.fs.softening_unit)
    m.flowsheet_unit_order.append(m.fs.acidification_unit)
    m.flowsheet_unit_order.append(m.fs.pump_unit)
    m.flowsheet_unit_order.append(m.fs.ro_unit)

    if hpro:
        m.flowsheet_unit_order.append(m.fs.hp_pump_unit)
        m.flowsheet_unit_order.append(m.fs.hpro_unit)
        m.flowsheet_unit_order.append(m.fs.product_mixer)

    m.flowsheet_unit_order.append(m.fs.erd_unit)
    m.flowsheet_unit_order.append(m.fs.product)
    m.flowsheet_unit_order.append(m.fs.brine)

    m.fs.reaktoro_blocks = []
    m.fs.reaktoro_blocks.append(m.fs.softening_unit.precipitation_block)
    m.fs.reaktoro_blocks.append(m.fs.acidification_unit.chemistry_block)
    m.fs.reaktoro_blocks.append(m.fs.ro_unit.scaling_block)
    if hpro:
        m.fs.reaktoro_blocks.append(m.fs.hpro_unit.scaling_block)

    # build all connections
    m.fs.feed.outlet.connect_to(m.fs.softening_unit.inlet)
    m.fs.softening_unit.outlet.connect_to(m.fs.acidification_unit.inlet)
    m.fs.acidification_unit.outlet.connect_to(m.fs.pump_unit.inlet)
    m.fs.pump_unit.outlet.connect_to(m.fs.ro_unit.feed)

    if hpro:
        m.fs.ro_unit.retentate.connect_to(m.fs.hp_pump_unit.inlet)
        m.fs.hp_pump_unit.outlet.connect_to(m.fs.hpro_unit.feed)
        m.fs.hpro_unit.retentate.connect_to(m.fs.erd_unit.inlet)

        m.fs.ro_unit.product.connect_to(m.fs.product_mixer.ro_inlet)
        m.fs.hpro_unit.product.connect_to(m.fs.product_mixer.hpro_inlet)
        m.fs.product_mixer.outlet.connect_to(m.fs.product.inlet)
    else:
        m.fs.ro_unit.retentate.connect_to(m.fs.erd_unit.inlet)
        m.fs.ro_unit.product.connect_to(m.fs.product.inlet)
    m.fs.erd_unit.outlet.connect_to(m.fs.brine.inlet)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.product.product.properties[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.product.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.product.properties[0].flow_vol
    )
    m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)
    TransformationFactory("network.expand_arcs").apply_to(m)
    add_global_constraints(m)
    fix_and_scale(m)
    scu.scale_costing_block(m.fs.costing)
    add_perfomance_tracking_vars(m)
    return m


def add_global_constraints(m):
    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0.0, 1),
        domain=Reals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )
    m.fs.eq_water_recovery = Constraint(
        expr=sum(m.fs.feed.feed.properties[0].flow_mass_phase_comp["Liq", :])
        * m.fs.water_recovery
        == m.fs.product.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    iscale.set_scaling_factor(m.fs.water_recovery, 1)
    iscale.constraint_scaling_transform(m.fs.eq_water_recovery, 1)


def add_perfomance_tracking_vars(m):
    m.fs.ipopt_iterations = Var(
        [
            "Number of iterations",
            "Objective function evaluations",
            "Objective gradient evaluations",
            "Equality constraint evaluations",
            "Inequality constraint evaluations",
            "Equality constraint Jacobian evaluations",
            "Inequality constraint Jacobian evaluations",
            "Lagrangian Hessian evaluations",
        ],
        initialize=0,
        domain=Reals,
        doc="Number of ipopt iterations",
    )
    m.fs.ipopt_iterations.fix()
    m.fs.scaled_ipopt_result = Var(
        [
            "dual_infeasibility",
            "constraint_violation",
            "complementarity_error",
            "variable_bound_violation",
            "overall_nlp_error",
        ],
        initialize=0,
        domain=Reals,
        doc="Scaled dual infeasibility from ipopt result",
    )
    m.fs.unscaled_ipopt_result = Var(
        [
            "dual_infeasibility",
            "constraint_violation",
            "complementarity_error",
            "variable_bound_violation",
            "overall_nlp_error",
        ],
        initialize=0,
        domain=Reals,
        doc="unscaled dual infeasibility from ipopt result",
    )
    m.fs.scaled_ipopt_result.fix()


def report_global_state(m):
    data_dict = {"Global results": {}}
    data_dict["Global results"]["DOfs"] = int(degrees_of_freedom(m))
    data_dict["Global results"]["Water recovery"] = m.fs.water_recovery
    data_dict["Global results"]["LCOW"] = m.fs.costing.LCOW

    build_report_table("Global results", data_dict)


def fix_and_scale(m):
    for unit in m.flowsheet_unit_order:
        unit.fix_and_scale()
    # limiiting maximum pH during optimization for
    # stability - in all example waters, pH in softening does not really operate
    # above 11, but ability of softening unit to go above 11 can cuase instability
    # Additonally, when alkalinitry drops bellow <10 ppm we observe poor solvablity and in practice
    # softening does not reduce alkalinity to below 20 ppm, as such
    # we constrain it to 20 ppm
    # (per veolia https://www.watertechnologies.com/handbook/chapter-07-precipitation-softening)
    max_ph = 11
    iscale.calculate_scaling_factors(m)
    m.fs.softening_unit.precipitation_reactor.alkalinity.setlb(20)
    m.fs.ro_unit.ro_feed.pH.setlb(6)
    m.fs.ro_unit.ro_feed.pH.setub(max_ph)
    if m.fs.find_component("hp_pump_unit") is not None:
        m.fs.hpro_unit.ro_feed.pH.setlb(6)
        m.fs.hpro_unit.ro_feed.pH.setub(max_ph)
    m.fs.softening_unit.precipitation_reactor.pH["outlet"].setlb(6)
    m.fs.softening_unit.precipitation_reactor.pH["outlet"].setub(max_ph)
    assert degrees_of_freedom(m) == 0


def initialize(m, linear_solver="mumps", tee=False, **kwargs):
    for unit in m.flowsheet_unit_order:
        unit.initialize()
    m.fs.costing.initialize()
    report_all_units(m)
    solve_model(m, linear_solver=linear_solver, tee=tee)
    set_optimization(m)

    if m.fs.water_recovery.value < 0.5:
        m.fs.water_recovery.fix()
        solve_model(m, linear_solver=linear_solver, tee=tee)
    # else:
    m.fs.water_recovery.fix(0.5)
    solve_model(m, linear_solver=linear_solver, tee=tee)
    print("--------------Initialization complete--------")


def report_all_units(m):
    for unit in m.flowsheet_unit_order:
        unit.report()
    report_global_state(m)


def set_optimization(m):
    for unit in m.flowsheet_unit_order:
        unit.set_optimization_operation()
    report_all_units(m)


def test_func(m, **kwargs):
    """Test function to check if the model is suitable for solving
    For Seawater case, with out HPRO maximum recoveriy is 65% while, with HPRO it is 86%
    """
    print(
        f"Testing func: {m.water_case}, {m.fs.find_component('hpro_unit')}, {m.fs.water_recovery.value}"
    )

    if "Seawater" in m.water_case:
        if (
            m.fs.find_component("hpro_unit") is None
            and m.fs.water_recovery.value > 0.64
        ):
            return False
        elif m.fs.water_recovery.value > 0.86:
            return False
    return True


def solve_model(m, tee=True, linear_solver="mumps", **kwargs):
    if linear_solver == "mumps":
        pivtol = 1e-4
        maxpivtol = 1e0
    else:
        pivtol = None
        maxpivtol = None
    solver = get_cyipopt_watertap_solver(
        linear_solver=linear_solver,
        max_iter=1000,
        limited_memory=m.solver_limited_memory,
        scalar_type=m.solver_limited_memory_scalar,
        pivtol=pivtol,
        pivtolmax=maxpivtol,
    )
    if m.fs.find_component("ipopt_iterations") is not None:
        tmp = tempfile.NamedTemporaryFile(suffix=".txt", delete=False)
        tmp.close()
        solver.options["output_file"] = tmp.name
    else:
        tmp = None
    result = solver.solve(m, tee=tee)
    if tmp is not None:
        matched_keys, parsed_output = ipopt_perf_utils.get_ipopt_performance_data(
            tmp.name
        )
        os.remove(tmp.name)
        iters_keys = list(m.fs.ipopt_iterations.keys())
        for i, k in enumerate(matched_keys.groups()):
            if k is not None:
                k = int(float(k))
            else:
                k = 0

            m.fs.ipopt_iterations[iters_keys[i]] = k
        for key in m.fs.unscaled_ipopt_result:
            m.fs.unscaled_ipopt_result[key] = parsed_output.get(key, 0)
            m.fs.scaled_ipopt_result[key] = parsed_output["final_scaled_results"].get(
                key, 0
            )
        m.fs.ipopt_iterations["Number of iterations"] = int(parsed_output["iters"])
    if tee:
        print("------vars_close_to_bound-tests---------")
        print_variables_close_to_bounds(m)
        print("------constraints_close_to_bound-tests---------")
        print_constraints_close_to_bounds(m)

    assert_optimal_termination(result)
    return result


if __name__ == "__main__":
    main()
