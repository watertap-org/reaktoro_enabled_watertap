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

from idaes.core.util.model_statistics import degrees_of_freedom

from pyomo.environ import (
    Var,
    Constraint,
    value,
    Expression,
    Reals,
    units as pyunits,
)
from pyomo.common.config import ConfigValue
from idaes.core import (
    declare_process_block_class,
)
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale

from reaktoro_enabled_watertap.unit_models.surrogate_crystallizer import (
    SurrogateCrystallizer,
)
from reaktoro_enabled_watertap.utils import scale_utils as scu

from reaktoro_pse.reaktoro_block import ReaktoroBlock
from collections import OrderedDict
from reaktoro_enabled_watertap.utils.reaktoro_utils import (
    ViablePrecipitants,
    ReaktoroOptionsContainer,
)

__author__ = "Alexander V. Dudchenko"


def get_default_viable_precipitants():
    """
    Returns the default set of viable precipitants for the crystallization unit.
    """
    viable_precipitants = ViablePrecipitants()
    viable_precipitants.register_solid(
        "Anhydrite",
        172.17 * pyunits.g / pyunits.mol,
        {"Ca_2+": 1, "SO4_2-": 1},
        "Ca_2+",
        reaktoro_modifier={"Ca": -1, "S": -1, "O": -4},
    )
    viable_precipitants.register_solid(
        "Halite",
        58.3197 * pyunits.g / pyunits.mol,
        {"Na_+": 1, "Cl_-": 1},
        "Na_+",
        reaktoro_modifier={"Na": -1, "Cl": -1},
    )
    viable_precipitants.register_solid(
        "Mirabilite",
        322.9 * pyunits.g / pyunits.mol,
        {"Na_+": 2, "SO4_2-": 2, "H2O": 10},
        "Na_+",
        reaktoro_modifier={"Na": -2, "SO4_2": -2, "H2O": -10},
    )
    viable_precipitants.register_solid(
        "Sylvite",
        74.55 * pyunits.g / pyunits.mol,
        {"K_+": 1, "Cl_-": 1},
        "K_+",
        reaktoro_modifier={"K": -1, "Cl": -1},
    )
    viable_precipitants.register_solid(
        "Bischofite",
        203.30 * pyunits.g / pyunits.mol,
        {"Mg_2+": 1, "Cl_-": 2, "H2O": 6},
        "Mg_2+",
        reaktoro_modifier={"Mg": -1, "Cl": -2, "H2O": -6},
    )
    viable_precipitants.register_solid(
        "Epsomite",
        246.48 * pyunits.g / pyunits.mol,
        {"Mg_2+": 1, "SO4_2-": 1, "H2O": 7},
        "Mg_2+",
        reaktoro_modifier={"Mg": -1, "SO4_2": -1, "H2O": -7},
    )
    return viable_precipitants


@declare_process_block_class("CrystallizationUnit")
class CrystallizationUnitData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
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
        "selected_precipitants",
        ConfigValue(
            default=["Calcite"],
            description="List of precipitants that should form in the reactor",
            doc="""
            This selects precipitants to form in reactor from Viableprecipitants class and adds them to reactor. 
            Default available precipitants are Calcite, Gypsum, and Brucite.
            To add additional precipitants, pass initialized Viableprecipitants
            object with registered new precipitants into selected_precipitants config option.
            """,
        ),
    )
    CONFIG.declare(
        "viable_precipitants",
        ConfigValue(
            default=None,
            description="Viableprecipitants class that defines all possible precipitants that can form",
            doc="""
                Should be a Viableprecipitants class that defines information for all formed precipitants
            """,
        ),
    )
    CONFIG.declare(
        "add_reaktoro_chemistry",
        ConfigValue(
            default=True,
            description="To use Reaktoro-PSE for estimating amount of solids formed",
            doc="""
            If True, builds a reaktoro block and uses it to calculate how much solids form based on feed
            composition adn amount of added reagent. 
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
        "reaktoro_inlet_state",
        ConfigValue(
            default=None,
            description="""dictionary that contains composition, pressure and/or temperature t
            "use for reaktoro inlet calcualtions""",
            doc="""
                    This primarly used in TVC/MVR Crystallizers where we want to pass in temperature from
                    the mixer into the crystallizer, as we do not wish to model effect of temperature on pH in the heat exchanger
            """,
        ),
    )

    def build(self):
        super().build()
        if self.config.viable_precipitants is None:
            self.config.viable_precipitants = get_default_viable_precipitants()
        self.selected_precipitants = {
            key: self.config.viable_precipitants[key]
            for key in self.config.selected_precipitants
        }

        self.crystallizer = SurrogateCrystallizer(
            property_package=self.config.default_property_package,
            vapor_property_package=self.config.vapor_property_package,
            solids_ions_dict={
                key: item["precipitation_stoichiometric"]
                for key, item in self.selected_precipitants.items()
            },
        )
        self.crystallizer.properties_in[0].flow_mol_phase_comp[...]
        self.crystallizer.pH = Var(
            ["inlet", "outlet"],
            initialize=7,
            units=pyunits.dimensionless,
            bounds=(0, 13),
        )
        if self.config.track_pE:
            self.crystallizer.pE = Var(
                ["inlet", "outlet"],
                initialize=0,
                units=pyunits.dimensionless,
                bounds=(None, None),
            )
        if self.config.default_costing_package is not None:
            self.crystallizer.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package,
                **self.config.default_costing_package_kwargs,
            )

        self.build_mol_flow_constraints()
        inlet_vars = {"pH": self.crystallizer.pH["inlet"]}
        outlet_vars = {"pH": self.crystallizer.pH["outlet"]}
        if self.config.track_pE:
            inlet_vars["pE"] = self.crystallizer.pE["inlet"]
            outlet_vars["pE"] = self.crystallizer.pE["outlet"]
        self.register_port(
            "inlet",
            self.crystallizer.inlet,
            inlet_vars,
        )
        self.register_port(
            "outlet",
            self.crystallizer.outlet,
            outlet_vars,
        )
        self.register_port(
            "vapor",
            self.crystallizer.vapor,
            outlet_vars,
        )
        if self.config.add_reaktoro_chemistry:
            self.add_reaktoro_chemistry()
        # elif self.config.add_non_eq_reaktoro_chemistry:
        #     self.add_non_eq_reaktoro_chemistry()
        else:
            self.build_equality_ph_pe_constraints()
            self.build_sat_p_pressure()

    def build_mol_flow_constraints(self):
        self.crystallizer.flow_mol_precipitate = Var(
            self.selected_precipitants.keys(),
            initialize=1,
            units=pyunits.mol / pyunits.s,
            bounds=(None, None),
        )
        self.crystallizer.flow_mol_vapor = Var(
            initialize=1,
            units=pyunits.mol / pyunits.s,
            bounds=(0, None),
        )

        @self.crystallizer.Constraint(self.selected_precipitants.keys())
        def mol_flow_eq(fs, phase):
            return self.crystallizer.flow_mol_precipitate[phase] == pyunits.convert(
                self.crystallizer.flow_mass_sol_comp_apparent[phase]
                / self.selected_precipitants[phase]["mw"],
                to_units=pyunits.mol / pyunits.s,
            )

        self.crystallizer.mol_flow_vapor_eq = Constraint(
            expr=self.crystallizer.flow_mol_vapor
            == pyunits.convert(
                self.crystallizer.flow_mass_vap_total
                / (18.015 * pyunits.g / pyunits.mol),
                to_units=pyunits.mol / pyunits.s,
            )
        )

    def build_equality_ph_pe_constraints(self):
        self.eq_outlet_pH = Constraint(
            expr=self.crystallizer.pH["inlet"] == self.crystallizer.pH["outlet"]
        )
        if self.config.track_pE:
            self.eq_outlet_pE = Constraint(
                expr=self.crystallizer.pE["inlet"] == self.crystallizer.pE["outlet"]
            )

    def build_sat_p_pressure(self):
        self.crystallizer.properties_out_vapor[0].pressure_sat()
        self.crystallizer.eq_sat_op_pressure = Constraint(
            expr=self.crystallizer.pressure_operating
            == self.crystallizer.properties_out_vapor[0].pressure_sat
        )
        iscale.set_scaling_factor(self.crystallizer.eq_sat_op_pressure, 1e-5)

    def add_reaktoro_chemistry(self):
        reagents = {}

        self.rkt_block_outputs = {("pH", None): self.crystallizer.pH["outlet"]}
        if self.config.track_pE:
            self.rkt_block_outputs[("pE", None)] = self.crystallizer.pE["outlet"]
        self.rkt_block_outputs[("vaporPressure", "H2O(g)")] = (
            self.crystallizer.vapor_pressure
        )
        for (
            phase,
            obj,
        ) in self.crystallizer.flow_mol_precipitate.items():
            self.rkt_block_outputs[("speciesAmount", phase)] = obj

        if self.config.add_alkalinity:
            self.crystallizer.alkalinity = Var(
                initialize=1,
                units=pyunits.mg / pyunits.L,
                doc="Alkalinity (mg/L as CaCO3)",
            )
            self.rkt_block_outputs[("alkalinityAsCaCO3", None)] = (
                self.crystallizer.alkalinity
            )

        self.reaktoro_options = ReaktoroOptionsContainer()
        inlet_state = {
            "temperature": self.crystallizer.properties_in[0].temperature,
            "pressure": self.crystallizer.properties_in[0].pressure,
            "pH": self.crystallizer.pH["inlet"],
            "composition": self.crystallizer.properties_in[0].flow_mol_phase_comp,
        }
        if self.config.track_pE:
            inlet_state["pE"] = self.crystallizer.pE["inlet"]
        if self.config.reaktoro_inlet_state is not None:
            inlet_state.update(self.config.reaktoro_inlet_state)
            print(inlet_state)
        self.reaktoro_options.system_state_option(
            "temperature",
            inlet_state["temperature"],
        )
        self.reaktoro_options.system_state_option(
            "pressure",
            inlet_state["pressure"],
        )
        self.reaktoro_options.system_state_option("pH", inlet_state["pH"])
        if self.config.track_pE:
            self.reaktoro_options.system_state_option("pE", inlet_state["pE"])
        self.reaktoro_options.aqueous_phase_option(
            "composition",
            inlet_state["composition"],
        )
        self.reaktoro_options["mineral_phase"] = {
            "phase_components": list(self.selected_precipitants.keys())
        }
        self.reaktoro_options["gas_phase"] = {
            "phase_components": ["H2O(g)", "Ntg(g)"],
            "activity_model": "ActivityModelRedlichKwong",
        }

        self.reaktoro_options["system_state_modifier"] = {
            "temperature": self.crystallizer.temperature_operating,
            "pressure": self.crystallizer.pressure_operating,
        }
        self.reaktoro_options
        self.reaktoro_options["chemistry_modifier"] = {
            "H2O_evaporation": self.crystallizer.flow_mol_vapor
        }
        self.reaktoro_options["outputs"] = self.rkt_block_outputs
        self.reaktoro_options.update_with_user_options(self.config.reaktoro_options)

        self.precipitation_block = ReaktoroBlock(**self.reaktoro_options)

    def set_fixed_operation(self):
        self.crystallizer.temperature_operating.fix(315)
        self.crystallizer.evaporation_percent.fix(50)
        self.crystallizer.evaporation_percent.setlb(1)
        self.crystallizer.evaporation_percent.setub(100)
        self.crystallizer.properties_out_vapor[0].pressure.setlb(0)

        self.crystallizer.properties_out_liq[0].pressure.setlb(0)
        self.crystallizer.eq_P_con1.deactivate()
        # self.crystallizer.eq_P_con3.deactivate()
        # self.crystallizer.pressure_operating.fix(101325)
        self.crystallizer.properties_out_liq[0].pressure.fix(101325)
        for precip in self.selected_precipitants.keys():
            self.crystallizer.flow_mol_precipitate[precip].setlb(0)
            self.crystallizer.flow_mass_sol_comp_apparent[precip].setlb(None)
            if self.config.add_reaktoro_chemistry == False:
                self.crystallizer.flow_mass_sol_comp_apparent[precip].fix(1e-5)

    def set_optimization_operation(self):
        """if we have reaktoro chemistry, we need to unfix the precipitate flow and reagent addition"""
        self.crystallizer.evaporation_percent.unfix()
        self.crystallizer.temperature_operating.unfix()
        if self.config.add_reaktoro_chemistry:
            for phase, data in self.selected_precipitants.items():
                self.crystallizer.flow_mol_precipitate[phase].unfix()
                self.crystallizer.flow_mass_sol_comp_apparent[phase].unfix()

    def rescale_reaktoro_jacobian(self, rescale_factor):
        rkt_jac = self.precipitation_block.display_jacobian_scaling()
        print(rkt_jac)
        new_scaling = {}
        for k in rkt_jac["property_block"]:
            new_scaling[k] = rkt_jac["property_block"][k] * rescale_factor
        self.precipitation_block.update_jacobian_scaling(
            user_scaling_dict=new_scaling, set_on_speciation_block=False
        )
        self.precipitation_block.display_jacobian_scaling()

    def scale_before_initialization(self, flow_factor=1, solid_gen_factor=1, **kwargs):
        # TODO: Identify better auto scaling factors for dose/precipitation
        # generally it appears using the smallest does step is a good scaling factor
        # but it might change depending on chemicals being added.
        sc = (
            self.config.default_property_package._default_scaling_factors[
                "flow_mol_phase_comp", ("Liq", "H2O")
            ]
            / flow_factor
        )
        iscale.set_scaling_factor(self.crystallizer.flow_mol_vapor, sc)
        iscale.constraint_scaling_transform(self.crystallizer.mol_flow_vapor_eq, sc)
        solid_mass_scales = []
        for precip in self.selected_precipitants.keys():
            scales = []
            # use scale for ion and add scale from added chemicals

            for ion, stoich in self.selected_precipitants[precip][
                "precipitation_stoichiometric"
            ].items():
                sc = (
                    self.config.default_property_package._default_scaling_factors[
                        "flow_mol_phase_comp", ("Liq", ion)
                    ]
                    / stoich
                ) / solid_gen_factor
                scales.append(sc)
            # want maximum scale for limiting ion

            precip_scale = max(scales)
            mol_precip_scale = precip_scale
            mass_precip_scale = precip_scale / value(
                pyunits.convert(
                    self.config.viable_precipitants[precip]["mw"],
                    pyunits.kg / pyunits.mol,
                )
            )

            solid_mass_scales.append(mass_precip_scale)
            iscale.set_scaling_factor(
                self.crystallizer.flow_mass_sol_comp_apparent[precip],
                mass_precip_scale,
            )

            iscale.set_scaling_factor(
                self.crystallizer.flow_mol_precipitate[precip],
                mol_precip_scale,
            )
            iscale.constraint_scaling_transform(
                self.crystallizer.mol_flow_eq[precip],
                mol_precip_scale,
            )
        iscale.set_scaling_factor(
            self.crystallizer.flow_mass_sol_total, max(solid_mass_scales)
        )
        iscale.set_scaling_factor(
            self.crystallizer.eq_total_solids_constraint, max(solid_mass_scales)
        )
        if self.config.default_costing_package is not None:
            costing_package = self.config.default_costing_package
            iscale.set_scaling_factor(
                costing_package.surrogate_crystallizer.steam_pressure, 1
            )
            iscale.set_scaling_factor(
                costing_package.surrogate_crystallizer.efficiency_pump, 1
            )
            iscale.set_scaling_factor(
                costing_package.surrogate_crystallizer.pump_head_height, 1
            )
            iscale.set_scaling_factor(
                costing_package.surrogate_crystallizer.fob_unit_cost, 1e-6
            )
            iscale.set_scaling_factor(
                costing_package.surrogate_crystallizer.ref_capacity, 1
            )
            iscale.set_scaling_factor(
                costing_package.surrogate_crystallizer.ref_exponent, 1
            )
            iscale.set_scaling_factor(
                costing_package.surrogate_crystallizer.iec_percent, 1
            )
            print(iscale.get_scaling_factor(self.crystallizer.flow_mass_sol_total))
            scu.calculate_scale_from_dependent_vars(
                self.crystallizer.costing.capital_cost,
                self.crystallizer.costing.capital_cost_constraint,
                [self.crystallizer.flow_mass_sol_total],
            )
        for ind, c in self.crystallizer.eq_P_con3.items():
            iscale.constraint_scaling_transform(c, 1e-5)
        iscale.set_scaling_factor(self.crystallizer.pressure_operating, 1e-5)
        iscale.set_scaling_factor(self.crystallizer.pH, 1)
        if self.config.track_pE:
            iscale.set_scaling_factor(self.crystallizer.pE, 1)
        if self.config.add_reaktoro_chemistry:
            if self.config.add_alkalinity:
                iscale.set_scaling_factor(self.crystallizer.alkalinity, 1)
        else:
            iscale.constraint_scaling_transform(self.eq_outlet_pH, 1)
            if self.config.track_pE:
                iscale.constraint_scaling_transform(self.eq_outlet_pE, 1)

    def initialize_unit(self, **kwargs):
        for phase, data in self.selected_precipitants.items():
            # assume that only fraction of ions will actually precipitate
            flow = (
                self.crystallizer.inlet.flow_mol_phase_comp[
                    0.0, "Liq", data["primary_ion"]
                ].value
                * 0.0001
            )
            self.crystallizer.flow_mol_precipitate[phase].fix(flow)
            self.crystallizer.flow_mass_sol_comp_apparent[phase].unfix()
        self.crystallizer.initialize()
        if self.config.add_reaktoro_chemistry:
            self.precipitation_block.initialize()
            self.precipitation_block.display_jacobian_scaling()
            self.crystallizer.initialize()
            for phase, data in self.selected_precipitants.items():
                self.crystallizer.flow_mol_precipitate[phase].unfix()
                self.crystallizer.flow_mass_sol_comp_apparent[phase].unfix()

    def get_model_state_dict(self):
        def get_ion_comp(stream, pH, pE=None):
            data_dict = OrderedDict()
            data_dict["Mass flow of H2O"] = stream.flow_mass_phase_comp["Liq", "H2O"]
            for phase, ion in stream.conc_mass_phase_comp:
                if ion != "H2O":
                    data_dict[ion] = stream.conc_mass_phase_comp[phase, ion]
            data_dict["pH"] = pH
            if pE is not None:
                data_dict["pE"] = pE
            data_dict["Temperature"] = stream.temperature
            data_dict["Pressure"] = stream.pressure
            return data_dict

        def get_pe(port):
            if self.config.track_pE:
                return self.crystallizer.pE[port]
            else:
                return None

        unit_dofs = degrees_of_freedom(self)

        treated_state = get_ion_comp(
            self.crystallizer.properties_out_liq[0],
            self.crystallizer.pH["outlet"],
            get_pe("outlet"),
        )

        if (
            self.config.add_alkalinity != False
            and self.crystallizer.find_component("alkalinity") is not None
        ):
            treated_state["Alkalinity (mg/L as CaCO3)"] = self.crystallizer.alkalinity

        model_state_dict = {
            "Model": {"DOFs": unit_dofs},
            "Inlet state": get_ion_comp(
                self.crystallizer.properties_in[0],
                self.crystallizer.pH["inlet"],
            ),
            "Inlet mol flow": self.crystallizer.properties_in[0].flow_mol_phase_comp,
            "Outlet mol flow": self.crystallizer.properties_out_liq[
                0
            ].flow_mol_phase_comp,
            "Inlet mol flow": self.crystallizer.properties_in[0].flow_mass_phase_comp,
            "Outlet mol flow": self.crystallizer.properties_out_liq[
                0
            ].flow_mass_phase_comp,
            "Solids ion mass flow": self.crystallizer.flow_mass_sol_comp_true,
            "Solids formed:": self.crystallizer.flow_mass_sol_comp_apparent,
            "Solids formed (mol basis):": self.crystallizer.flow_mol_precipitate,
            "Treated state": treated_state,
            "Operation": {
                "Heater work": self.crystallizer.heat_required,
                "Heat duty": pyunits.convert(self.crystallizer.heat_duty, pyunits.kW),
                "Evaporation %": self.crystallizer.evaporation_percent,
                "Operating temperature": self.crystallizer.temperature_operating,
                "Operating pressure": self.crystallizer.pressure_operating,
                "Vapor pressure": self.crystallizer.vapor_pressure,
            },
        }
        if hasattr(self, "rkt_block_outputs"):
            model_state_dict["Reaktoro outputs"] = self.rkt_block_outputs
        if hasattr(self, "eq_rkt_outputs"):
            model_state_dict["Equilibrium rkt outputs"] = self.eq_rkt_outputs
            model_state_dict["Non-equilibrium rkt outputs"] = self.non_eq_rkt_outputs
            model_state_dict["Reaktoro precipitation efficiency"] = (
                self.crystallizer.non_eq_efficacy
            )
        if self.config.default_costing_package is not None:
            model_state_dict["Costs"] = {
                "Capital cost": self.crystallizer.costing.capital_cost
            }
            # for reagent in self.crystallizer.flow_mass_reagent:
            #     model_state_dict[f"Costs"][f"Reagent {reagent} cost"] = (
            #         self.config.default_costing_package.aggregate_flow_costs[
            #             f"{self.name}_reagent_{reagent}".replace(".", "_")
            #         ]
            #     )

        return model_state_dict
