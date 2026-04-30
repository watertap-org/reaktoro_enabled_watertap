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

from reaktoro_enabled_watertap.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from idaes.core.util.initialization import propagate_state
from reaktoro_enabled_watertap.utils.reaktoro_utils import (
    ReaktoroOptionsContainer,
)
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import (
    assert_optimal_termination,
)
from pyomo.environ import (
    Var,
    value,
    Constraint,
    Objective,
    units as pyunits,
)
from pyomo.common.config import ConfigValue
from idaes.models.unit_models import Separator, Mixer
from idaes.models.unit_models.translator import Translator
from idaes.models.unit_models.separator import SplittingType
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
)
from watertap.unit_models.mvc.components import Compressor, Condenser

from reaktoro_enabled_watertap.unit_models.evaporator import EvaporatorMCAS
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from watertap.unit_models.pressure_changer import Pump
from idaes.core import (
    declare_process_block_class,
)
from idaes.core import UnitModelCostingBlock

import idaes.core.util.scaling as iscale
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
from reaktoro_enabled_watertap.utils import scale_utils as scu

from reaktoro_pse.reaktoro_block import ReaktoroBlock

from pyomo.network import Arc
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)

__author__ = "Alexander V. Dudchenko"


@declare_process_block_class("MultiCompMVC")
class MultiCompMVCData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "mvc_property_package",
        ConfigValue(
            default=None,
            description="Property package to use with MVC componenetns instead of default property package",
            doc="""
                This will use supplied property package, but also build translator blocks between default property pacakge
                and MVC property pacakge.
            """,
        ),
    )
    CONFIG.declare(
        "vapor_property_package",
        ConfigValue(
            default=None,
            description="Vapor property package that defines all possible reagents",
            doc="""
                Should be a ViableReagents class that contains all of the selected reagents
            """,
        ),
    )
    CONFIG.declare(
        "selected_scalants",
        ConfigValue(
            default={"Calcite": 1, "Gypsum": 1},
            description="Dict of scalants to track",
            doc="""
            Defines which scalants to track, and their maximum potential during optimization. 
            Provide a dictionary with following structure:
            {'Scalant': maximum scaling potential}
            if add_reaktoro_chemistry == True:
                Will add maximum_scaling_potential variable, and use int eq_maximum_scaling_potential constraint to define scaling limits
                The eq_maximum_scaling_potential is deactivated until set_optimization_operation is called.
            """,
        ),
    )
    CONFIG.declare(
        "add_reaktoro_chemistry",
        ConfigValue(
            default=True,
            description="To use Reaktoro-PSE for estimating scaling potential",
            doc="""
            If True, builds a reaktoro block and uses it to calculate scaling potential and pH, which are then used in constraints to limit scaling and track pH. If False, does not build reaktoro block and does not calculate scaling, but still tracks pH and limits it to be the same in feed and brine.
            """,
        ),
    )

    CONFIG.declare(
        "reaktoro_options",
        ConfigValue(
            default=None,
            description="User options for configuring Reaktoro-PSE provided as a dict",
            doc="""
            User can provide additional reaktoro options, or override defaults provided by ReaktoroOptionsContainer class
            """,
        ),
    )

    CONFIG.declare(
        "add_alkalinity",
        ConfigValue(
            default=False,
            description="Defines if to add alkalinity to Reaktoro output",
            doc="""
                Defines if to add alkalinity to Reaktoro output, this will use Reaktoro to calculate alkalinity
            """,
        ),
    )
    CONFIG.declare(
        "track_pH",
        ConfigValue(
            default=True,
            description="if pH should be tracked in the model",
            doc="""
                    Providing True will add pH variable to the model and track it
            """,
        ),
    )
    CONFIG.declare(
        "track_pE",
        ConfigValue(
            default=False,
            description="if pE should be tracked in the model",
            doc="""
                    Providing True will add pE variable to the model and track it
            """,
        ),
    )
    CONFIG.declare(
        "target_recovery",
        ConfigValue(
            default=0.5,
            description="Target recovery for MVC during initialization",
            doc="""
            Sets a target recovery for MVC during initialization
            """,
        ),
    )
    CONFIG.declare(
        "target_temperature",
        ConfigValue(
            default=50,
            description="Target temperature for MVC during initialization",
            doc="""
            Sets a target temperature for MVC during initialization, this is the temperature at which evaporator outlet brine will be fixed during initialization
            """,
        ),
    )
    CONFIG.declare(
        "mvc_material_factor",
        ConfigValue(
            default=9,
            description="Material factor for evaporator used in costing",
            doc="""
            Material factor for evaporator used in costing, this is a multiplier on the evaporator cost to account for more expensive materials needed for high scaling potential and corrosive streams.
            """,
        ),
    )

    def build(self):
        super().build()
        feed_vars = None
        brine_vars = None
        mvc_prop_pack = self.config.default_property_package
        if self.config.mvc_property_package is not None:
            mvc_prop_pack = self.config.mvc_property_package
            if len(self.config.mvc_property_package.solute_set) > 1:
                raise TypeError("MVC prop pack should only be single comp type for now")
            self.mvc_solute_type = list(self.config.mvc_property_package.solute_set)[0]
            self.feed_translator = Translator(
                inlet_property_package=self.config.default_property_package,
                outlet_property_package=self.config.mvc_property_package,
            )
            self.feed_translator.properties_in[0].conc_mass_phase_comp[...]
            self.brine_translator = Translator(
                inlet_property_package=self.config.mvc_property_package,
                outlet_property_package=self.config.default_property_package,
            )
            self.brine_translator.properties_out[0].conc_mass_phase_comp[...]
            self.distillate_translator = Translator(
                inlet_property_package=self.config.mvc_property_package,
                outlet_property_package=self.config.default_property_package,
            )
            self.distillate_translator.properties_out[0].conc_mass_phase_comp[...]
            self.setup_inlet_translator_block(self.feed_translator)

            # we get 100% rejection, so no solutes in outlet
            self.setup_outlet_translator_block(
                self.distillate_translator,
                self.feed_translator.inlet.flow_mol_phase_comp,
                no_solutes=True,
            )
            # we get 100% rejection so all solutes pass
            self.setup_outlet_translator_block(
                self.brine_translator,
                self.feed_translator.inlet.flow_mol_phase_comp,
                all_solutes=True,
            )

        self.pump_feed = Pump(property_package=mvc_prop_pack)
        self.pump_feed.control_volume.properties_in[0].conc_mass_phase_comp[...]
        self.pump_brine = Pump(property_package=mvc_prop_pack)

        self.pump_distillate = Pump(property_package=mvc_prop_pack)
        self.separator_feed = Separator(
            property_package=mvc_prop_pack,
            outlet_list=["hx_distillate_cold", "hx_brine_cold"],
            split_basis=SplittingType.totalFlow,
        )

        self.hx_distillate = HeatExchanger(
            hot_side_name="hot",
            cold_side_name="cold",
            hot={
                "property_package": mvc_prop_pack,
                "has_pressure_change": True,
            },
            cold={
                "property_package": mvc_prop_pack,
                "has_pressure_change": True,
            },
            delta_temperature_callback=delta_temperature_chen_callback,
            flow_pattern=HeatExchangerFlowPattern.countercurrent,
        )
        self.hx_brine = HeatExchanger(
            hot_side_name="hot",
            cold_side_name="cold",
            hot={
                "property_package": mvc_prop_pack,
                "has_pressure_change": True,
            },
            cold={
                "property_package": mvc_prop_pack,
                "has_pressure_change": True,
            },
            delta_temperature_callback=delta_temperature_chen_callback,
            flow_pattern=HeatExchangerFlowPattern.countercurrent,
        )

        self.mixer_feed = Mixer(
            property_package=mvc_prop_pack,
            momentum_mixing_type=MomentumMixingType.none,
            inlet_list=["hx_distillate_cold", "hx_brine_cold"],
        )
        # self.mixer_feed.pressure_equality_constraints[0, 2].deactivate()

        self.mixer_feed.mixer_pressure_eq = Constraint(
            expr=self.mixer_feed.hx_distillate_cold_state[0].pressure
            == self.mixer_feed.mixed_state[0].pressure
        )

        self.evaporator = EvaporatorMCAS(
            property_package_feed=mvc_prop_pack,
            property_package_vapor=self.config.vapor_property_package,
        )
        self.compressor = Compressor(
            property_package=self.config.vapor_property_package
        )

        self.condenser = Condenser(property_package=self.config.vapor_property_package)

        self.tb_distillate = Translator(
            inlet_property_package=self.config.vapor_property_package,
            outlet_property_package=mvc_prop_pack,
        )

        @self.tb_distillate.Constraint(
            list(self.tb_distillate.properties_out[0].flow_mass_phase_comp.keys())
        )
        def eq_flow_mass_comp(blk, key, ion):
            if "Liq" in key:
                if "H2O" in ion:
                    return (
                        blk.properties_in[0].flow_mass_phase_comp["Liq", ion]
                        == blk.properties_out[0].flow_mass_phase_comp[key, ion]
                    )
                else:
                    blk.properties_out[0].flow_mass_phase_comp[key, ion].fix(1e-12)
                    return Constraint.Skip
                    # return  blk.properties_in[0].flow_mass_phase_comp["Liq", ion]
                    #   == blk.properties_out[0].flow_mass_phase_comp[key, ion]
            else:
                return Constraint.Skip

        @self.tb_distillate.Constraint()
        def eq_temperature(blk):
            return blk.properties_in[0].temperature == blk.properties_out[0].temperature

        @self.tb_distillate.Constraint()
        def eq_pressure(blk):
            return blk.properties_in[0].pressure == blk.properties_out[0].pressure

        if self.config.track_pH:
            self.pH = Var(
                ["feed", "brine"],
                initialize=7,
                bounds=(0, 13),
                units=pyunits.dimensionless,
            )

            feed_vars = {"pH": self.pH["feed"]}
            brine_vars = {"pH": self.pH["brine"]}

            if self.config.track_pE:
                self.pE = Var(
                    ["feed", "brine"],
                    initialize=0,
                    units=pyunits.dimensionless,
                    bounds=(None, None),
                )
                feed_vars["pE"] = self.pE["feed"]
                brine_vars["pE"] = self.pE["brine"]

        self.pump_to_splitter = Arc(
            source=self.pump_feed.outlet, destination=self.separator_feed.inlet
        )
        self.sep_to_distillate_hx = Arc(
            source=self.separator_feed.hx_distillate_cold,
            destination=self.hx_distillate.cold_inlet,
        )
        self.sep_to_brine_hx = Arc(
            source=self.separator_feed.hx_brine_cold,
            destination=self.hx_brine.cold_inlet,
        )
        self.hx_distillate_to_mixer = Arc(
            source=self.hx_distillate.cold_outlet,
            destination=self.mixer_feed.hx_distillate_cold,
        )
        self.hx_brine_to_mixer = Arc(
            source=self.hx_brine.cold_outlet, destination=self.mixer_feed.hx_brine_cold
        )
        self.mixer_to_evaporator = Arc(
            source=self.mixer_feed.outlet, destination=self.evaporator.inlet_feed
        )
        self.evap_to_compressor = Arc(
            source=self.evaporator.outlet_vapor, destination=self.compressor.inlet
        )
        self.compressor_to_condenser = Arc(
            source=self.compressor.outlet, destination=self.condenser.inlet
        )
        self.evap_to_brine_pump = Arc(
            source=self.evaporator.outlet_brine, destination=self.pump_brine.inlet
        )
        self.brine_pump_to_hx = Arc(
            source=self.pump_brine.outlet, destination=self.hx_brine.hot_inlet
        )
        self.condenser_to_tb_distillate = Arc(
            source=self.condenser.outlet, destination=self.tb_distillate.inlet
        )
        self.tb_distillate_to_pump = Arc(
            source=self.tb_distillate.outlet, destination=self.pump_distillate.inlet
        )
        self.pump_to_hx_distillate = Arc(
            source=self.pump_distillate.outlet, destination=self.hx_distillate.hot_inlet
        )

        self.evaporator.connect_to_condenser(self.condenser)
        if self.config.mvc_property_package is not None:
            self.feed_to_pump = Arc(
                source=self.feed_translator.outlet, destination=self.pump_feed.inlet
            )
            self.distillate_hx_to_distillate = Arc(
                source=self.hx_distillate.hot_outlet,
                destination=self.distillate_translator.inlet,
            )
            self.brine_hx_to_brine = Arc(
                source=self.hx_brine.hot_outlet, destination=self.brine_translator.inlet
            )
            self.register_port("inlet", self.feed_translator.inlet, feed_vars)
            self.register_port("distillate", self.distillate_translator.outlet)
            self.register_port("brine", self.brine_translator.outlet, brine_vars)
        else:
            self.register_port("inlet", self.pump_feed.inlet, feed_vars)
            self.register_port("distillate", self.hx_distillate.hot_outlet)
            self.register_port("brine", self.hx_brine.hot_outlet, brine_vars)

            self.hx_brine.hot_side.properties_out[0].conc_mass_phase_comp[...]
            self.hx_distillate.hot_side.properties_out[0].conc_mass_phase_comp[...]

        self.build_water_removal_constraint()

        if self.config.add_reaktoro_chemistry and self.config.track_pH:
            self.build_scaling_constraints()
            self.add_reaktoro_chemistry()
        if self.config.add_reaktoro_chemistry == False and self.config.track_pH:
            self.add_retentate_ph_pe_constraint()
        if self.config.default_costing_package is not None:
            self.add_costing()
        self.add_recovery_constraint()
        self.initialized = False
        # self.add_Q_ext()

    def get_brine_state(self):
        if self.config.mvc_property_package is not None:
            return self.brine_translator.properties_out[0]
        else:
            return self.hx_brine.hot_side.properties_out[0]

    def setup_inlet_translator_block(
        self,
        translator_block,
    ):
        """defines inlet translator block, will sum up all mcas massflow and translate them to mass of single solute
        in ro property package
        Args:
            translator_block -- Inlet translator block (should be feed)"""
        ions = []
        for index in translator_block.properties_in[0].flow_mol_phase_comp:
            if "H2O" not in index:
                ions.append(
                    translator_block.properties_in[0].flow_mol_phase_comp[index]
                    * translator_block.properties_in[0].mw_comp[index[1]]
                )

        @translator_block.Constraint(["H2O", self.mvc_solute_type])
        def eq_flow_phase_comp(blk, ion):
            if ion == "H2O":
                return (
                    translator_block.properties_in[0].flow_mol_phase_comp["Liq", "H2O"]
                    * translator_block.properties_in[0].mw_comp["H2O"]
                    == translator_block.properties_out[0].flow_mass_phase_comp[
                        "Liq", "H2O"
                    ]
                )
            else:
                return translator_block.properties_out[0].flow_mass_phase_comp[
                    "Liq", self.mvc_solute_type
                ] == sum(ions)

        translator_block.eq_pressure_equality = Constraint(
            expr=translator_block.properties_in[0].pressure
            == translator_block.properties_out[0].pressure
        )
        translator_block.eq_temperature_equality = Constraint(
            expr=translator_block.properties_in[0].temperature
            == translator_block.properties_out[0].temperature
        )

    def setup_outlet_translator_block(
        self, translator_block, inlet_composition, no_solutes=False, all_solutes=False
    ):
        """defines outlet translator block, will convert single outlet solute to multi solutes, assuming they are
        ratiometrically related to changes between inlet and outlet properties. e.g.
         out_ion=in_total_ion_mass/out_total_ion_mass*in_ion_mass
        Args:
            translator_block -- Outlet translator block (should be feed)
            inlet_composition -- Inlet composition used to get original mass flow of ions entering system
            no_solutes -- If True, will assume that all solute are removed and set mass flow of all solutes to 0
            all_solutes -- If True, will assume that all solutes are present in the same ratio as inlet and set mass flow of all solutes to inlet mass flow ratio*total outlet mass flow of solute (e.g. out_ion=in_ion_mass/in_total_ion_mass*out_total_ion_mass)
        """
        tds_in = []
        ions = []
        for index in inlet_composition:
            if "H2O" not in index:
                tds_in.append(
                    inlet_composition[index]
                    * translator_block.properties_out[0].mw_comp[index[-1]]
                )
            ions.append(index[-1])

        @translator_block.Constraint(ions)
        def eq_flow_phase_comp(blk, ion):
            liq = "Liq"
            if "H2O" in ion:
                return (
                    blk.properties_out[0].flow_mol_phase_comp[liq, ion]
                    == blk.properties_in[0].flow_mass_phase_comp[liq, ion]
                    / blk.properties_out[0].mw_comp[ion]
                )
            else:
                if no_solutes:
                    blk.properties_out[0].flow_mol_phase_comp[liq, ion].fix(0)
                    return Constraint.Skip
                elif all_solutes:
                    return (
                        blk.properties_out[0].flow_mol_phase_comp[liq, ion]
                        == inlet_composition[0.0, liq, ion]
                    )
                else:
                    return blk.properties_out[0].flow_mol_phase_comp[
                        liq, ion
                    ] == inlet_composition[0.0, liq, ion] * blk.properties_in[
                        0
                    ].flow_mass_phase_comp[
                        "Liq", self.mvc_solute_type
                    ] / sum(
                        tds_in
                    )

        translator_block.eq_pressure_equality = Constraint(
            expr=translator_block.properties_in[0].pressure
            == translator_block.properties_out[0].pressure
        )

        translator_block.eq_temperature_equality = Constraint(
            expr=translator_block.properties_in[0].temperature
            == translator_block.properties_out[0].temperature
        )

    def add_costing(self):
        self.pump_feed.costing = UnitModelCostingBlock(
            flowsheet_costing_block=self.config.default_costing_package,
            costing_method_arguments={"pump_type": "low_pressure"},
        )
        self.pump_distillate.costing = UnitModelCostingBlock(
            flowsheet_costing_block=self.config.default_costing_package,
            costing_method_arguments={"pump_type": "low_pressure"},
        )
        self.pump_brine.costing = UnitModelCostingBlock(
            flowsheet_costing_block=self.config.default_costing_package,
            costing_method_arguments={"pump_type": "low_pressure"},
        )
        self.hx_distillate.costing = UnitModelCostingBlock(
            flowsheet_costing_block=self.config.default_costing_package
        )
        self.hx_brine.costing = UnitModelCostingBlock(
            flowsheet_costing_block=self.config.default_costing_package
        )
        self.mixer_feed.costing = UnitModelCostingBlock(
            flowsheet_costing_block=self.config.default_costing_package
        )
        self.evaporator.costing = UnitModelCostingBlock(
            flowsheet_costing_block=self.config.default_costing_package
        )
        self.compressor.costing = UnitModelCostingBlock(
            flowsheet_costing_block=self.config.default_costing_package
        )

    def add_elevated_evap_temp(self):
        self.evap_constraint = Constraint(
            expr=self.evaporator.properties_feed[0].temperature
            <= self.evaporator.properties_brine[0].temperature
        )
        iscale.set_scaling_factor(self.evap_constraint, 1e-2)

    def add_recovery_constraint(self):
        """adds recovery constraint for MVC"""
        self.recovery = Var(
            initialize=self.config.target_recovery,
            bounds=(1e-8, 1),
            units=pyunits.dimensionless,
        )
        self.recovery.fix(self.config.target_recovery)
        self.eq_recovery = Constraint(
            expr=self.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
            == self.recovery
            * sum(
                j
                for i, j in self.pump_feed.control_volume.properties_in[
                    0
                ].flow_mass_phase_comp.items()
            )
        )
        self.split_ratio_recovery_equality = Constraint(
            expr=self.separator_feed.split_fraction[0, "hx_distillate_cold"]
            == self.recovery
        )

    def add_retentate_ph_pe_constraint(self):
        """adds evaporator's brine pH constraint"""
        self.eq_ph_equality = Constraint(expr=self.pH["feed"] == self.pH["brine"])
        if self.config.track_pE:
            self.eq_pE_equality = Constraint(expr=self.pE["feed"] == self.pE["brine"])

    def build_scaling_constraints(self):
        """builds scaling constraints"""
        self.scaling_tendency = Var(
            list(self.config.selected_scalants.keys()),
            initialize=1,
            units=pyunits.dimensionless,
        )

        self.maximum_scaling_tendency = Var(
            list(self.config.selected_scalants.keys()),
            initialize=lambda blk, idx: self.config.selected_scalants[idx],
            units=pyunits.dimensionless,
        )
        self.maximum_scaling_tendency.fix()

        @self.Constraint(list(self.config.selected_scalants.keys()))
        def eq_max_scaling_tendency(blk, scalant):
            return (
                blk.scaling_tendency[scalant] <= blk.maximum_scaling_tendency[scalant]
            )

        self.deactivate_scaling_constraints()

    def activate_scaling_constraints(self):
        """activates scaling constraints"""
        if self.config.add_reaktoro_chemistry:
            for scalant in self.config.selected_scalants.keys():
                self.eq_max_scaling_tendency[scalant].activate()
                self.scaling_tendency[scalant].unfix()
                self.maximum_scaling_tendency[scalant].fix()

    def deactivate_scaling_constraints(self):
        """deactivates scaling constraints"""
        if self.config.add_reaktoro_chemistry:
            for scalant in self.config.selected_scalants.keys():
                self.eq_max_scaling_tendency[scalant].deactivate()
                self.maximum_scaling_tendency[scalant].fix()

    def build_water_removal_constraint(self):
        """builds water removal constraint"""
        self.evaporated_water_amount = Var(
            initialize=1, bounds=(1e-8, None), units=pyunits.mol / pyunits.s
        )

        self.eq_evaporated_water = Constraint(
            expr=(
                self.evaporated_water_amount
                == self.evaporator.properties_vapor[0].flow_mass_phase_comp[
                    "Vap", "H2O"
                ]
                / self.config.default_property_package.mw_comp["H2O"]
            )
        )

    def add_reaktoro_chemistry(self):
        """add relevant reaktoro block"""

        outputs = {("pH", None): self.pH["brine"]}
        if self.config.track_pE:
            outputs[("pE", None)] = self.pE["brine"]
        for scalant in self.scaling_tendency:
            outputs[("scalingTendency", scalant)] = self.scaling_tendency[scalant]
        if self.config.add_alkalinity:
            self.alkalinity = Var(
                initialize=1,
                units=pyunits.mg / pyunits.L,
                doc="Alkalinity (mg/L as CaCO3)",
            )

            outputs[("alkalinityAsCaCO3", None)] = self.alkalinity

        self.reaktoro_options = ReaktoroOptionsContainer()
        if self.config.mvc_property_package is not None:
            inlet_stream = self.feed_translator.properties_in[0]
        else:
            inlet_stream = self.pump_feed.control_volume.properties_in[0]
        self.reaktoro_options.system_state_option(
            "temperature",
            inlet_stream.temperature,
        )
        if self.config.track_pE:
            self.reaktoro_options.system_state_option("pE", self.pE["feed"])
        self.reaktoro_options.system_state_option(
            "pressure",
            inlet_stream.pressure,
        )
        self.reaktoro_options.system_state_option("pH", self.pH["feed"])
        self.reaktoro_options.aqueous_phase_option(
            "composition",
            inlet_stream.flow_mol_phase_comp,
        )

        self.reaktoro_options["chemistry_modifier"] = {
            "H2O_evaporation": self.evaporated_water_amount
        }
        self.reaktoro_options.system_state_modifier_option(
            "pressure", self.evaporator.properties_brine[0].pressure
        )
        self.reaktoro_options.system_state_modifier_option(
            "temperature", self.evaporator.properties_brine[0].temperature
        )
        iscale.set_scaling_factor(
            self.evaporator.properties_brine[0].pressure, 1e-5
        )  # ensure its scaled!
        self.reaktoro_options["outputs"] = outputs
        self.reaktoro_options.update_with_user_options(self.config.reaktoro_options)
        self.scaling_block = ReaktoroBlock(**self.reaktoro_options)

    def set_fixed_operation(self):
        """fixes operation point for pump unit model"""
        self.pump_feed.efficiency_pump[0].fix(0.8)
        self.pump_feed.control_volume.deltaP[0].fix(7e3)

        self.separator_feed.split_fraction[0, "hx_distillate_cold"] = (
            self.recovery.value
        )
        # Distillate HX
        self.hx_distillate.overall_heat_transfer_coefficient[0].fix(2e3)
        self.hx_distillate.area.fix(125)
        self.hx_distillate.cold.deltaP[0].fix(-7e3)
        self.hx_distillate.hot.deltaP[0].fix(-7e3)

        # Set lower bound of approach temperatures
        # self.hx_distillate.delta_temperature_in.setlb(0.1)
        # self.hx_distillate.delta_temperature_out.setlb(0.1)
        self.hx_distillate.area.setlb(1)
        self.hx_distillate.hot_side.properties_in[0].pressure.setlb(1)
        self.hx_distillate.hot_side.properties_out[0].pressure.setlb(1)
        # Brine HX
        self.hx_brine.overall_heat_transfer_coefficient[0].fix(2e3)
        self.hx_brine.area.fix(115)
        self.hx_brine.cold.deltaP[0].fix(-7e3)
        self.hx_brine.hot.deltaP[0].fix(-7e3)
        # Set lower bound of approach temperatures
        # self.hx_brine.delta_temperature_in.setlb(0.1)
        # self.hx_brine.delta_temperature_out.setlb(0.1)
        self.hx_brine.area.setlb(1)

        self.hx_brine.hot_side.properties_in[0].pressure.setlb(1)
        self.hx_brine.hot_side.properties_out[0].pressure.setlb(1)
        # Evaporator
        self.evaporator.inlet_feed.temperature[0] = (
            self.config.target_temperature - 1 + 273.15
        )  # provide guess
        self.evaporator.outlet_brine.temperature[0].fix(
            self.config.target_temperature + 273.15
        )
        self.evaporator.U.fix(3e3)  # W/K-m^2
        self.evaporator.area.setub(1e4)  # m^2
        self.evaporator.delta_temperature_out.setlb(0.01)  # K
        # Compressor
        self.compressor.pressure_ratio.fix(1.6)
        self.compressor.efficiency.fix(0.8)

        self.compressor.pressure_ratio.setlb(1.1)
        self.compressor.pressure_ratio.setub(3)

        # Brine pump
        self.pump_brine.efficiency_pump[0].fix(0.8)
        self.pump_brine.inlet.pressure.setlb(1)
        self.pump_brine.outlet.pressure.setlb(1)
        self.pump_brine.outlet.pressure.fix(101325 + 4e4)
        # Distillate pump
        self.pump_distillate.inlet.pressure.setlb(1)
        self.pump_distillate.outlet.pressure.setlb(1)
        self.pump_distillate.outlet.pressure.fix(101325 + 4e4)
        self.pump_distillate.efficiency_pump[0].fix(0.8)

        # Temperature bounds
        self.evaporator.properties_vapor[0].temperature.setub(75 + 273.15)
        self.evaporator.properties_brine[0].pressure.setlb(1)
        self.compressor.control_volume.properties_out[0].temperature.setub(450)

        # pressure on translator block
        self.tb_distillate.properties_out[0].pressure.setlb(1)
        if self.config.default_costing_package is not None:
            self.config.default_costing_package.evaporator.material_factor_cost.fix(
                self.config.mvc_material_factor
            )

            self.config.default_costing_package.heat_exchanger.material_factor_cost.fix(
                self.config.mvc_material_factor
            )
        self.recovery.fix()
        if self.initialized:
            self.set_brine_distillate_outlet_pressure()
        # self.Q_ext.fix(0)

    def set_brine_distillate_outlet_pressure(self, pressure=101325):
        self.pump_distillate.outlet.pressure.unfix()
        self.pump_brine.outlet.pressure.unfix()
        self.hx_brine.hot_side.properties_out[0].pressure.fix(pressure)
        self.hx_distillate.hot_side.properties_out[0].pressure.fix(pressure)
        self.pump_feed.control_volume.deltaP[0].unfix()
        self.evaporator.properties_feed[0].pressure.fix(pressure)

    def set_optimization_operation(self):
        """unfixes operation for optimization"""
        self.evaporator.area.unfix()
        self.evaporator.outlet_brine.temperature[0].unfix()
        self.evaporator.properties_vapor[0].temperature.setlb(273.15 + 40)
        self.compressor.pressure_ratio.unfix()
        self.hx_distillate.area.unfix()
        self.hx_brine.area.unfix()
        self.recovery.unfix()
        self.hx_distillate.delta_temperature_in.setlb(1)
        self.hx_distillate.delta_temperature_out.setlb(1)

        self.hx_brine.delta_temperature_in.setlb(1)
        self.hx_brine.delta_temperature_out.setlb(1)

        self.evaporator.delta_temperature_in.setlb(1)
        self.evaporator.delta_temperature_out.setlb(1)
        if self.find_component("scaling_block") is not None:
            self.activate_scaling_constraints()

    def get_default_scaling_factors(self):
        """returns scale for ro property package and default property package"""
        scale_factors = {"mol": {}, "mass": {}}

        scale_factors["mass"]["H2O"] = (
            self.config.default_property_package._default_scaling_factors[
                "flow_mol_phase_comp", ("Liq", "H2O")
            ]
            / value(self.config.default_property_package.mw_comp["H2O"])
        )
        scale_factors["mol"]["H2O"] = (
            self.config.default_property_package._default_scaling_factors[
                "flow_mol_phase_comp", ("Liq", "H2O")
            ]
        )
        tds_sf = []
        for ion in self.config.default_property_package.solute_set:
            sf = self.config.default_property_package._default_scaling_factors[
                "flow_mol_phase_comp", ("Liq", ion)
            ]

            scale_factors["mol"][ion] = sf
            scale_factors["mass"][ion] = sf * value(
                self.config.default_property_package.mw_comp[ion]
            )
            tds_sf.append(scale_factors["mass"][ion])
        # scale for ro property package
        tds_sf = sum(tds_sf)

        if self.config.mvc_property_package is not None:
            scale_factors["mass"][self.mvc_solute_type] = tds_sf
        return scale_factors

    def scale_before_initialization(self, **kwargs):
        """scales the unit before initialization"""

        # Get the scaling factors for the MVC property package and default property package
        prop_scaling = self.get_default_scaling_factors()

        iscale.constraint_scaling_transform(self.mixer_feed.mixer_pressure_eq, 1e-5)
        iscale.set_scaling_factor(self.recovery, 1)
        iscale.constraint_scaling_transform(self.eq_recovery, 1)

        iscale.constraint_scaling_transform(self.split_ratio_recovery_equality, 1)
        if self.find_component("eq_ph_equality") is not None:
            iscale.constraint_scaling_transform(self.eq_ph_equality, 1)

        if self.config.track_pH:
            iscale.set_scaling_factor(self.pH, 1)
            if self.config.track_pE:
                iscale.set_scaling_factor(self.pE, 1)
        for phase, ion in self.tb_distillate.eq_flow_mass_comp.keys():
            sf = prop_scaling["mass"][ion]
            iscale.constraint_scaling_transform(
                self.tb_distillate.eq_flow_mass_comp["Liq", ion], sf
            )
            iscale.constraint_scaling_transform(
                self.evaporator.eq_mass_balance[0, ion], sf
            )

        if self.config.mvc_property_package is not None:
            if (
                self.config.mvc_property_package._default_scaling_factors.get(
                    ("flow_mass_phase_comp", ("Liq", "H2O"))
                )
                is None
            ):
                self.config.mvc_property_package.set_default_scaling(
                    "flow_mass_phase_comp",
                    prop_scaling["mass"]["H2O"],
                    index=("Liq", "H2O"),
                )
                self.config.mvc_property_package.set_default_scaling(
                    "flow_mass_phase_comp",
                    prop_scaling["mass"][self.mvc_solute_type],
                    index=("Liq", self.mvc_solute_type),
                )
            # Scale the translator blocks
            for block, scale_type in [
                (self.feed_translator, "mass"),
                (self.distillate_translator, "mol"),
                (self.brine_translator, "mol"),
            ]:
                iscale.constraint_scaling_transform(block.eq_pressure_equality, 1e-5)
                iscale.constraint_scaling_transform(block.eq_temperature_equality, 1e-2)

                iscale.constraint_scaling_transform(
                    block.eq_flow_phase_comp["H2O"],
                    prop_scaling[scale_type]["H2O"],
                )
                # Scale the mass flow of the solute in the translator block,
                # for inlet, we use scale from ro_prop_scaling, for outlet we use default scaling
                # for default property package, Water scaling is same for both

                for ion in block.eq_flow_phase_comp:
                    sf = prop_scaling[scale_type][ion]
                    iscale.constraint_scaling_transform(
                        block.eq_flow_phase_comp[ion],
                        sf,
                    )
        iscale.set_scaling_factor(self.tb_distillate.eq_temperature, 1e-2)
        iscale.set_scaling_factor(self.tb_distillate.eq_pressure, 1e-5)

        fluid_flow_scale = prop_scaling["mass"]["H2O"]
        if (
            self.config.vapor_property_package._default_scaling_factors.get(
                ("flow_mass_phase_comp", ("Liq", "H2O"))
            )
            is None
        ):
            self.config.vapor_property_package.set_default_scaling(
                "flow_mass_phase_comp", fluid_flow_scale, index=("Liq", "H2O")
            )
            self.config.vapor_property_package.set_default_scaling(
                "flow_mass_phase_comp", fluid_flow_scale, index=("Vap", "H2O")
            )

            print("Scaled vapor H2O  flow with scale ", fluid_flow_scale)
        self.config.vapor_property_package.set_default_scaling(
            "enth_mass_phase", 1e-5, index=("Vap")
        )
        if self.find_component("evaporated_water_amount") is not None:
            iscale.constraint_scaling_transform(
                self.eq_evaporated_water,
                self.config.default_property_package._default_scaling_factors[
                    "flow_mol_phase_comp", ("Liq", "H2O")
                ],
            )
            iscale.set_scaling_factor(
                self.evaporated_water_amount,
                self.config.default_property_package._default_scaling_factors[
                    "flow_mol_phase_comp", ("Liq", "H2O")
                ],
            )
        if self.find_component("scaling_tendency") is not None:
            for scalant in self.config.selected_scalants.keys():
                iscale.set_scaling_factor(
                    self.scaling_tendency[scalant],
                    1 / self.config.selected_scalants[scalant],
                )
                iscale.set_scaling_factor(
                    self.maximum_scaling_tendency[scalant],
                    1 / self.config.selected_scalants[scalant],
                )
                iscale.constraint_scaling_transform(
                    self.eq_max_scaling_tendency[scalant],
                    1 / self.config.selected_scalants[scalant],
                )

        work_scale = 1e-3 / (1 / fluid_flow_scale)
        heat_scale = 1e-6 / (1 / fluid_flow_scale)
        iscale.set_scaling_factor(
            self.separator_feed.split_fraction[0.0, "hx_distillate_cold"], 1
        )

        iscale.set_scaling_factor(
            self.separator_feed.split_fraction[0.0, "hx_brine_cold"], 1
        )
        iscale.set_scaling_factor(self.pump_feed.control_volume.work, work_scale)
        iscale.set_scaling_factor(self.pump_brine.control_volume.work, work_scale)
        iscale.set_scaling_factor(self.pump_distillate.control_volume.work, work_scale)

        # distillate HX
        iscale.set_scaling_factor(self.hx_distillate.hot.heat, heat_scale)
        iscale.set_scaling_factor(self.hx_distillate.cold.heat, heat_scale)

        iscale.set_scaling_factor(self.hx_distillate.delta_temperature_in[0.0], 1e-2)
        iscale.set_scaling_factor(self.hx_distillate.delta_temperature_out[0.0], 1e-2)
        iscale.set_scaling_factor(
            self.hx_distillate.overall_heat_transfer_coefficient, 1e-3
        )

        iscale.set_scaling_factor(self.hx_distillate.area, 1e-1)
        iscale.constraint_scaling_transform(
            self.hx_distillate.cold_side.pressure_balance[0], 1e-5
        )
        iscale.constraint_scaling_transform(
            self.hx_distillate.hot_side.pressure_balance[0], 1e-5
        )

        # brine HX
        iscale.set_scaling_factor(self.hx_brine.hot.heat, heat_scale)
        iscale.set_scaling_factor(self.hx_brine.cold.heat, heat_scale)
        iscale.set_scaling_factor(self.hx_brine.delta_temperature_in[0.0], 1e-2)
        iscale.set_scaling_factor(self.hx_brine.delta_temperature_out[0.0], 1e-2)
        iscale.set_scaling_factor(self.hx_brine.overall_heat_transfer_coefficient, 1e-3)
        iscale.set_scaling_factor(self.hx_brine.area, 1e-1)
        iscale.constraint_scaling_transform(
            self.hx_brine.cold_side.pressure_balance[0], 1e-5
        )
        iscale.constraint_scaling_transform(
            self.hx_brine.hot_side.pressure_balance[0], 1e-5
        )

        # evaporator
        iscale.set_scaling_factor(self.evaporator.area, 1e-1)
        iscale.set_scaling_factor(self.evaporator.U, 1e-3)
        iscale.set_scaling_factor(self.evaporator.delta_temperature_in, 1e-2)
        iscale.set_scaling_factor(self.evaporator.delta_temperature_out, 1e-2)
        iscale.set_scaling_factor(self.evaporator.lmtd, 1)

        iscale.set_scaling_factor(self.evaporator.heat_transfer, heat_scale)

        # compressor
        iscale.set_scaling_factor(self.compressor.pressure_ratio, 1)

        iscale.set_scaling_factor(self.compressor.efficiency, 1 / 0.8)

        iscale.set_scaling_factor(self.compressor.control_volume.work, work_scale)
        # condenser
        iscale.set_scaling_factor(self.condenser.control_volume.heat, heat_scale)

        iscale.set_scaling_factor(self.condenser.eq_condenser_pressure_sat, 1e-5)
        if self.config.default_costing_package is not None:
            iscale.set_scaling_factor(
                self.compressor.control_volume.properties_in[0].flow_mass_phase_comp[
                    "Vap", "H2O"
                ],
                fluid_flow_scale,
            )
            costing_package = self.config.default_costing_package
            iscale.set_scaling_factor(costing_package.compressor.unit_cost, 1e-3)
            iscale.set_scaling_factor(
                costing_package.evaporator.material_factor_cost, 1
            )
            scu.calculate_scale_from_dependent_vars(
                self.evaporator.costing.capital_cost,
                self.evaporator.costing.capital_cost_constraint,
                [
                    costing_package.compressor.unit_cost,
                    costing_package.evaporator.material_factor_cost,
                    self.evaporator.area,
                ],
            )
            iscale.set_scaling_factor(costing_package.compressor.unit_cost, 1e-3)
            iscale.set_scaling_factor(costing_package.compressor.exponent, 1)
            scu.calculate_scale_from_dependent_vars(
                self.compressor.costing.capital_cost,
                self.compressor.costing.capital_cost_constraint,
                [
                    costing_package.compressor.unit_cost,
                    costing_package.compressor.exponent,
                    self.compressor.control_volume.properties_in[
                        0
                    ].flow_mass_phase_comp["Vap", "H2O"],
                    self.compressor.pressure_ratio,
                    self.compressor.efficiency,
                ],
            )

            iscale.set_scaling_factor(costing_package.mixer.unit_cost, 1)
            mixer_scale = 1 / (
                costing_package.mixer.unit_cost.value * fluid_flow_scale  # / 1000
            )
            iscale.set_scaling_factor(self.mixer_feed.costing.capital_cost, mixer_scale)
            iscale.constraint_scaling_transform(
                self.mixer_feed.costing.capital_cost_constraint, mixer_scale
            )
            for pump in [self.pump_feed, self.pump_brine, self.pump_distillate]:
                scu.calculate_scale_from_dependent_vars(
                    pump.costing.capital_cost,
                    pump.costing.capital_cost_constraint,
                    [pump.control_volume.work[0]],
                )

            iscale.set_scaling_factor(costing_package.heat_exchanger.unit_cost, 1e-2)
            iscale.set_scaling_factor(
                costing_package.heat_exchanger.material_factor_cost, 1
            )
            for hx in [self.hx_distillate, self.hx_brine]:
                scu.calculate_scale_from_dependent_vars(
                    hx.costing.capital_cost,
                    hx.costing.capital_cost_constraint,
                    [
                        costing_package.heat_exchanger.unit_cost,
                        costing_package.heat_exchanger.material_factor_cost,
                        hx.area,
                    ],
                )

    def init_translator_block(self, block, additonal_var_to_fix=None):
        """initializes translator block"""
        if additonal_var_to_fix is not None:
            if not isinstance(additonal_var_to_fix, list):
                additonal_var_to_fix = [additonal_var_to_fix]
            for var in additonal_var_to_fix:
                var.fix()
        flags = fix_state_vars(block.properties_in)
        solver = get_solver()
        results = solver.solve(block, tee=False)
        assert_optimal_termination(results)
        revert_state_vars(block.properties_in, flags)
        if additonal_var_to_fix is not None:
            for var in additonal_var_to_fix:
                var.unfix()

    def guess_hx_operation(self, initialzied_inlet, guessed_inlet, deltaT_guess):
        """provides an initial guess for heat exchanger operation based on inlet conditions and a delta T guess"""

        if self.config.mvc_property_package is not None:
            for i in initialzied_inlet.flow_mass_phase_comp:
                guessed_inlet.flow_mass_phase_comp[i].value = (
                    initialzied_inlet.flow_mass_phase_comp[i].value
                )
        else:
            for i in initialzied_inlet.flow_mol_phase_comp:
                guessed_inlet.flow_mol_phase_comp[i].value = (
                    initialzied_inlet.flow_mol_phase_comp[i].value
                )
        guessed_inlet.temperature.value = (
            initialzied_inlet.temperature.value + deltaT_guess
        )
        guessed_inlet.pressure.value = initialzied_inlet.pressure.value

    def add_Q_ext(self):
        # Allows additional heat to be added to evaporator so that an initial feasible solution can be found as a starting
        # guess for optimization in case physically infeasible simulation is proposed

        self.Q_ext = Var(initialize=0, units=pyunits.J / pyunits.s)
        self.Q_ext.setlb(0)
        self.Q_ext.fix(0)
        self.evaporator.eq_energy_balance.deactivate()
        self.evaporator.eq_energy_balance_with_additional_Q = Constraint(
            expr=self.evaporator.heat_transfer
            + self.Q_ext
            + self.evaporator.properties_feed[0].enth_flow
            == self.evaporator.properties_brine[0].enth_flow
            + self.evaporator.properties_vapor[0].enth_flow_phase["Vap"]
        )
        iscale.set_scaling_factor(self.Q_ext, 1e-6)

    def initialize_unit(self):
        if self.config.mvc_property_package is not None:
            self.init_translator_block(self.feed_translator)
            propagate_state(self.feed_to_pump)
        self.pump_feed.initialize()
        propagate_state(self.pump_to_splitter)
        self.separator_feed.initialize()
        propagate_state(self.sep_to_distillate_hx)
        propagate_state(self.sep_to_brine_hx)
        self.guess_hx_operation(
            self.hx_distillate.cold_side.properties_in[0],
            self.hx_distillate.hot_side.properties_in[0],
            10,
        )
        self.guess_hx_operation(
            self.hx_brine.cold_side.properties_in[0],
            self.hx_brine.hot_side.properties_in[0],
            10,
        )
        for k in range(1):
            self.hx_distillate.initialize()
            self.hx_brine.initialize()
            propagate_state(self.hx_distillate_to_mixer)
            propagate_state(self.hx_brine_to_mixer)
            self.mixer_feed.initialize()
            propagate_state(self.mixer_to_evaporator)
            self.evaporator.initialize(recovery=self.config.target_recovery)
            # Disitilate side
            propagate_state(self.evap_to_compressor)
            self.compressor.initialize()
            propagate_state(self.compressor_to_condenser)
            self.condenser.initialize(heat=-self.evaporator.heat_transfer.value)
            propagate_state(self.condenser_to_tb_distillate)
            self.init_translator_block(self.tb_distillate)
            # self.tb_distillate.initialize()
            propagate_state(self.tb_distillate_to_pump)
            self.pump_distillate.initialize()
            propagate_state(self.pump_to_hx_distillate)
            # Brine side
            propagate_state(self.evap_to_brine_pump)
            self.pump_brine.initialize()
            propagate_state(self.brine_pump_to_hx)
            # self.hx_brine.initialize()
            # self.hx_distillate.initialize()
        if self.config.mvc_property_package is not None:
            inlet_props = [self.feed_translator.properties_in[0].flow_mol_phase_comp]
            if self.config.track_pH:
                inlet_props.append(self.pH["feed"])
            propagate_state(self.distillate_hx_to_distillate)
            propagate_state(self.brine_hx_to_brine)
            self.init_translator_block(self.brine_translator, inlet_props)
            self.init_translator_block(
                self.distillate_translator,
                self.feed_translator.properties_in[0].flow_mol_phase_comp,
            )
        if self.find_component("scaling_block") is not None:
            self.scaling_block.initialize()
            self.deactivate_scaling_constraints()
        self.inlet.fix()
        solver = get_solver()
        # prepare for initial optimal design
        self.set_brine_distillate_outlet_pressure()
        self.hx_distillate.delta_temperature_in.setlb(1)
        self.hx_distillate.delta_temperature_out.setlb(1)

        self.hx_brine.delta_temperature_in.setlb(1)
        self.hx_brine.delta_temperature_out.setlb(1)

        self.evaporator.delta_temperature_in.setlb(None)
        self.evaporator.delta_temperature_out.setlb(None)
        self.recovery.fix()
        assert degrees_of_freedom(self) == 0
        self.evaporator.area.unfix()
        self.evaporator.outlet_brine.temperature[0].unfix()
        self.compressor.pressure_ratio.unfix()
        self.hx_distillate.area.unfix()
        self.hx_brine.area.unfix()
        self.recovery.unfix()

        self._objective = Objective(
            expr=(
                (self.config.target_temperature + 273.15)
                - self.evaporator.outlet_brine.temperature[0]
            )
            ** 2
            + (self.config.target_recovery - self.recovery) ** 2
            + (self.evaporator.area) ** 2
            + (self.hx_distillate.area) ** 2
            + (self.hx_brine.area) ** 2
            + (self.compressor.pressure_ratio) ** 2
        )

        if self.find_component("scaling_block") is not None:
            self.scaling_block.deactivate()
        result = solver.solve(self, tee=True)
        if self.find_component("scaling_block") is not None:
            self.scaling_block.activate()
            self.scaling_block.initialize()
            solver = get_cyipopt_watertap_solver()
            result = solver.solve(self, tee=True)
        self.del_component("_objective")
        self.report()

        # from idaes.core.util.model_diagnostics import DiagnosticsToolbox

        # dg = DiagnosticsToolbox(self)
        # dg.display_constraints_with_large_residuals()
        # dg.display_variables_at_or_outside_bounds()
        # dg.compute_infeasibility_explanation()
        assert_optimal_termination(result)
        self.evaporator.area.unfix()
        self.evaporator.outlet_brine.temperature[0].fix()
        self.compressor.pressure_ratio.fix()
        self.hx_distillate.area.fix()
        self.hx_brine.area.fix()
        self.recovery.fix()
        assert degrees_of_freedom(self) == 0
        self.inlet.unfix()
        self.initialized = True
        # assert False

    def get_model_state_dict(self):
        """Returns a dictionary with the model state"""

        def get_ion_comp(stream, pH=None, pE=None, mass_flows_only=False):
            data_dict = {}
            data_dict["Mass flow of H2O"] = stream.flow_mass_phase_comp["Liq", "H2O"]
            if mass_flows_only:
                for phase, ion in stream.flow_mass_phase_comp:
                    if ion != "H2O":
                        data_dict[ion] = stream.flow_mass_phase_comp[phase, ion]
            else:
                for phase, ion in stream.conc_mass_phase_comp:
                    if ion != "H2O":
                        data_dict[ion] = stream.conc_mass_phase_comp[phase, ion]
            # data_dict["TDS"] = pyunits.convert(
            #     stream.total_disolved_solids, pyunits.g / pyunits.L
            # )
            if pH is not None:
                data_dict["pH"] = pH
            if pE is not None:
                data_dict["pE"] = pE
            data_dict["Temperature"] = stream.temperature
            data_dict["Pressure"] = stream.pressure
            return data_dict

        model_state_dict = {}
        if self.config.mvc_property_package is not None:
            model_state_dict["MVC Feed"] = get_ion_comp(
                self.feed_translator.properties_in[0]
            )
            model_state_dict["MVC Brine"] = get_ion_comp(
                self.brine_translator.properties_out[0]
            )
            model_state_dict["MVC Distillate"] = get_ion_comp(
                self.distillate_translator.properties_out[0]
            )
        else:
            model_state_dict["MVC Feed"] = get_ion_comp(
                self.pump_feed.control_volume.properties_in[0]
            )
            model_state_dict["MVC Brine"] = get_ion_comp(
                self.hx_brine.hot_side.properties_out[0]
            )
            model_state_dict["MVC Distillate"] = get_ion_comp(
                self.hx_distillate.hot_side.properties_out[0]
            )
        model_state_dict["Evaporator Feed"] = get_ion_comp(
            self.evaporator.properties_feed[0], mass_flows_only=True
        )
        model_state_dict["Evaporator Brine"] = get_ion_comp(
            self.evaporator.properties_brine[0], mass_flows_only=True
        )
        model_state_dict["Evaporator Vapor"] = {
            "Vapor flow H2O": self.evaporator.properties_vapor[0].flow_mass_phase_comp[
                "Vap", "H2O"
            ],
            "Temperature": self.evaporator.properties_vapor[0].temperature,
            "Pressure": self.evaporator.properties_vapor[0].pressure,
        }
        model_state_dict["Evaporator Design"] = {}
        model_state_dict["Evaporator Design"]["Area"] = self.evaporator.area

        model_state_dict["Evaporator Design"][
            "Heat transfer"
        ] = self.evaporator.heat_transfer
        model_state_dict["Evaporator Design"][
            "Delta T in"
        ] = self.evaporator.delta_temperature_in
        model_state_dict["Evaporator Design"][
            "Delta T out"
        ] = self.evaporator.delta_temperature_out
        model_state_dict["Evaporator Design"]["LMTD"] = self.evaporator.lmtd
        model_state_dict["Condenser"] = {
            "Temp in": self.condenser.control_volume.properties_in[0].temperature,
            "Temp out": self.condenser.control_volume.properties_out[0].temperature,
            "Pressure out": self.condenser.control_volume.properties_out[0].pressure,
            "Sat pressure": self.condenser.control_volume.properties_out[
                0
            ].pressure_sat,
        }
        model_state_dict["Compressor"] = {
            "Pressure ratio": self.compressor.pressure_ratio,
            "Work": self.compressor.control_volume.work[0],
            "Inlet pressure": self.compressor.control_volume.properties_in[0].pressure,
            "Inlet temperature": self.compressor.control_volume.properties_in[
                0
            ].temperature,
            "Outlet pressure": self.compressor.control_volume.properties_out[
                0
            ].pressure,
            "Outlet temperature": self.compressor.control_volume.properties_out[
                0
            ].temperature,
        }
        model_state_dict["Brine HX"] = {
            "Area": self.hx_brine.area,
            "Heat": self.hx_brine.hot.heat[0],
            "Delta temp in": self.hx_brine.delta_temperature_in[0],
            "Delta temp out": self.hx_brine.delta_temperature_out[0],
            # "LMTD": self.hx_brine.lmtd,
        }
        model_state_dict["Distillate HX"] = {
            "Area": self.hx_distillate.area,
            "Heat": self.hx_distillate.hot.heat[0],
            "Delta temp in": self.hx_distillate.delta_temperature_in[0],
            "Delta temp out": self.hx_distillate.delta_temperature_out[0],
            # "LMTD": self.hx_distillate.lmtd,
        }
        model_state_dict["MVC operation"] = {}
        if self.find_component("evaporated_water_amount") is not None:
            model_state_dict["MVC operation"][
                "Water Evaporated"
            ] = self.evaporated_water_amount
        model_state_dict["MVC operation"]["Water recovery"] = self.recovery
        if self.config.track_pH:
            model_state_dict["pH"] = self.pH
        if self.config.track_pE:
            model_state_dict["pE"] = self.pE
        if self.config.add_reaktoro_chemistry:
            model_state_dict["Scaling potential"] = {}
            model_state_dict["Maximum scaling potential"] = {}
            for scalant in self.scaling_tendency:
                model_state_dict["Scaling potential"][scalant] = self.scaling_tendency[
                    scalant
                ]
                model_state_dict["Maximum scaling potential"][scalant] = (
                    self.maximum_scaling_tendency[scalant]
                )
            model_state_dict["pH"] = self.pH
            if (
                self.config.add_alkalinity != False
                and self.find_component("alkalinity") is not None
            ):
                model_state_dict["Alkalinity (mg/L as CaCO3)"] = self.alkalinity

        return model_state_dict
