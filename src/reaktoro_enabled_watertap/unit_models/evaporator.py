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
    Block,
    Var,
    check_optimal_termination,
)

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
)
from watertap.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.unit_models.mvc.components.evaporator import EvaporatorData

from idaes.core.util.model_statistics import degrees_of_freedom

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("EvaporatorMCAS")
class EvaporatorMCASData(EvaporatorData):
    def initialize_build(
        self,
        delta_temperature_in=None,
        delta_temperature_out=None,
        recovery=None,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for pressure changer initialization routines

        Keyword Arguments:
            delta_temperature_in : value to fix delta_temperature_in
            delta_temperature_out : value to fix delta_temperature_out
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

        if hasattr(self, "connection_to_condenser"):
            self.connection_to_condenser.deactivate()

        # ---------------------------------------------------------------------
        # Initialize feed side
        self.properties_feed[0].flow_mass_phase_comp[...]  # Ensure its there
        flags_feed = self.properties_feed.initialize(
            solver=solver, optarg=optarg, hold_state=True
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # # ---------------------------------------------------------------------
        # # Initialize brine
        # Set state_args from inlet state
        if state_args is None:
            state_args = {}
            state_dict = self.properties_feed[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        self.properties_brine.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args
        )
        if recovery is not None:
            target_vapor_flow = (
                self.properties_feed[0].flow_mass_phase_comp["Liq", "H2O"].value
                * recovery
            )
        else:
            target_vapor_flow = (
                self.properties_feed[0].flow_mass_phase_comp["Liq", "H2O"].value
            )
        state_args_vapor = {}
        state_args_vapor["pressure"] = 0.5 * state_args["pressure"]
        state_args_vapor["temperature"] = state_args["temperature"]
        state_args_vapor["flow_mass_phase_comp"] = {
            ("Liq", "H2O"): self.properties_vapor[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .lb,
            ("Vap", "H2O"): target_vapor_flow,
        }
        self.properties_vapor.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_vapor,
        )

        init_log.info_high("Initialization Step 2 Complete.")

        # incorporate guessed temperature differences
        has_guessed_delta_temperature_in = False
        if delta_temperature_in is not None:
            if self.delta_temperature_in.is_fixed():
                raise RuntimeError(
                    "A guess was provided for the delta_temperature_in variable in the "
                    "initialization, but it is already fixed. Either do not "
                    "provide a guess for or unfix delta_temperature_in"
                )
            self.delta_temperature_in.fix(delta_temperature_in)
            has_guessed_delta_temperature_in = True

        has_guessed_delta_temperature_out = False
        if delta_temperature_out is not None:
            if self.delta_temperature_out.is_fixed():
                raise RuntimeError(
                    "A guess was provided for the delta_temperature_out variable in the "
                    "initialization, but it is already fixed. Either do not "
                    "provide a guess for or unfix delta_temperature_out"
                )
            self.delta_temperature_out.fix(delta_temperature_out)
            has_guessed_delta_temperature_out = True

        if recovery is not None:
            self.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix(
                target_vapor_flow
            )
        if delta_temperature_in is None and delta_temperature_out is None:
            self.properties_brine[0].temperature.fix()
        # print(
        #     "Degrees of freedom:",
        #     degrees_of_freedom(self),
        # )
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        if recovery is not None:
            self.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
        # ----------------------------------------  -----------------------------
        # Release feed and condenser inlet states and release delta_temperature
        self.properties_feed.release_state(flags_feed, outlvl=outlvl)
        if has_guessed_delta_temperature_in:
            self.delta_temperature_in.unfix()
        if has_guessed_delta_temperature_out:
            self.delta_temperature_out.unfix()
        if hasattr(self, "connection_to_condenser"):
            self.connection_to_condenser.activate()

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        sf = iscale.get_scaling_factor(self.properties_vapor[0].enth_flow_phase["Vap"])
        iscale.constraint_scaling_transform(self.eq_energy_balance[0], sf)
        iscale.constraint_scaling_transform(self.eq_brine_pressure[0], 1e-5)
        iscale.constraint_scaling_transform(self.eq_vapor_pressure[0], 1e-5)
        iscale.constraint_scaling_transform(self.eq_vapor_temperature[0], 1e-2)

        iscale.constraint_scaling_transform(self.eq_lmtd[0], 1)
        sf = iscale.get_scaling_factor(self.heat_transfer)
        iscale.constraint_scaling_transform(self.eq_evaporator_heat[0], sf)

        iscale.constraint_scaling_transform(
            self.connection_to_condenser.eq_delta_temperature_in[0], 1e-2
        )
        iscale.constraint_scaling_transform(
            self.connection_to_condenser.eq_delta_temperature_out[0], 1e-2
        )

        iscale.constraint_scaling_transform(
            self.connection_to_condenser.eq_heat_balance[0], sf
        )
        for phase, ion in self.properties_feed[0].flow_mass_phase_comp.keys():
            sf = iscale.get_scaling_factor(
                self.properties_feed[0].flow_mass_phase_comp[phase, ion]
            )
            iscale.constraint_scaling_transform(self.eq_mass_balance[0, ion], sf)
