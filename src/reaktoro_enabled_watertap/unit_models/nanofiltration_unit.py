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
    units as pyunits,
)
from pyomo.common.config import ConfigValue
from idaes.core import (
    declare_process_block_class,
)
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale

from watertap.unit_models.nanofiltration_ZO import NanofiltrationZO
from reaktoro_enabled_watertap.unit_models.multi_comp_pump_unit import (
    MultiCompPumpUnit,
    MultiCompPumpUnitData,
)
from pyomo.network import Arc
from idaes.core.util.initialization import propagate_state

from reaktoro_pse.reaktoro_block import ReaktoroBlock
from collections import OrderedDict
from reaktoro_enabled_watertap.utils.reaktoro_utils import (
    ReaktoroOptionsContainer,
)

__author__ = "Chenyu Wang"


@declare_process_block_class("Nanofiltration")
class NanofiltrationUnitData(WaterTapFlowsheetBlockData):
    LOCATIONS = ["feed", "retentate"]

    CONFIG = WaterTapFlowsheetBlockData.CONFIG()

    CONFIG.declare(
        "selected_scalants",
        ConfigValue(
            default={"Calcite": 1, "Gypsum": 1},
            description="Dict of scalants to track at feed and retentate, keyed by scalant name with maximum scaling tendency as value",
            doc="""
                Selects scalants to track in the model.
                Reaktoro will compute the scaling tendency for each selected scalant
                at both feed and retentate locations. The value represents the maximum
                allowable scaling tendency (1 = scaling threshold).
                Default: {"Calcite": 1, "Gypsum": 1}
            """,
        ),
    )
    CONFIG.declare(
        "add_reaktoro_chemistry",
        ConfigValue(
            default=True,
            description="Use Reaktoro-PSE for estimating osmotic pressure and scaling tendencies",
            doc="""
                If True, builds an indexed ReaktoroBlock over ['inlet', 'retentate'] and
                uses it to calculate:
                  - osmoticPressure at feed and retentate (replaces MCAS prop pack)
                  - scalingTendency for each selected scalant at feed and retentate
                  - pH at feed and retentate
                If False, simple equality constraints are used to pass pH through.
            """,
        ),
    )
    CONFIG.declare(
        "reaktoro_options",
        ConfigValue(
            default=None,
            description="User options for configuring Reaktoro-PSE as a dict",
            doc="""
                User can provide additional reaktoro options, or override defaults
                provided by ReaktoroOptionsContainer.
            """,
        ),
    )
    CONFIG.declare(
        "add_alkalinity",
        ConfigValue(
            default=False,
            description="Defines if to add alkalinity to Reaktoro output",
            doc="""
                If True, adds an alkalinity variable to the NF block and uses Reaktoro
                to calculate alkalinity (mg/L as CaCO3) at feed and retentate.
            """,
        ),
    )
    CONFIG.declare(
        "track_pE",
        ConfigValue(
            default=False,
            description="Track pE (redox potential) in the model",
            doc="""
                If True, adds a pE variable at feed and retentate locations and
                includes it in the Reaktoro inputs/outputs.
            """,
        ),
    )
    CONFIG.declare(
        "add_feed_pump",
        ConfigValue(
            default=False,
            description="Whether to include a feed pump for the NF unit",
            doc="""
                If True, a MultiCompPumpUnit is built and connected between the
                inlet and the NF unit.
            """,
        ),
    )
    CONFIG.declare("pump_options", MultiCompPumpUnitData.CONFIG())
    CONFIG.declare(
        "component_rejections",
        ConfigValue(
            default=None,
            description="Ion rejection fractions as a dict by ion name",
            doc="""
                Dict of {ion: rejection_fraction} applied when set_fixed_operation
                is called.
            """,
        ),
    )

    def build(self):
        super().build()

        self.nanofiltration = NanofiltrationZO(
            property_package=self.config.default_property_package,
        )

        self.nanofiltration.pH = Var(
            self.LOCATIONS,
            initialize=7,
            units=pyunits.dimensionless,
            bounds=(0, 13),
            doc="pH at inlet and retentate",
        )

        if self.config.track_pE:
            self.nanofiltration.pE = Var(
                self.LOCATIONS,
                initialize=0,
                units=pyunits.dimensionless,
                bounds=(None, None),
                doc="pE (redox potential) at inlet and retentate",
            )

        if self.config.default_costing_package is not None:
            self.nanofiltration.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package,
                **self.config.default_costing_package_kwargs,
            )

        # Register ports (inlet, retentate, permeate)
        retentate_vars = {"pH": self.nanofiltration.pH["retentate"]}
        if self.config.track_pE:
            retentate_vars["pE"] = self.nanofiltration.pE["retentate"]

        if self.config.add_feed_pump:
            self.pump = MultiCompPumpUnit(
                default_property_package=self.config.default_property_package,
                default_costing_package=self.config.default_costing_package,
                default_costing_package_kwargs=self.config.pump_options.default_costing_package_kwargs,
                track_pH=False,
                track_pE=False,
                initialization_pressure="osmotic_pressure",
                osmotic_over_pressure=self.config.pump_options.osmotic_over_pressure,
                osmotic_pressure_offset=self.config.pump_options.osmotic_pressure_offset,
                maximum_pressure=self.config.pump_options.maximum_pressure,
                pump_efficiency=self.config.pump_options.pump_efficiency,
            )
            self.register_port("inlet", self.pump.inlet.port)
            self.pump_to_nf = Arc(
                source=self.pump.outlet.port,
                destination=self.nanofiltration.inlet,
            )
        else:
            self.register_port("inlet", self.nanofiltration.inlet)
        self.register_port("retentate", self.nanofiltration.retentate, retentate_vars)
        self.register_port("permeate", self.nanofiltration.permeate)

        if self.config.add_reaktoro_chemistry:
            self.add_reaktoro_chemistry()
        else:
            self.build_equality_ph_pe_constraints()

    def build_equality_ph_pe_constraints(self):
        self.eq_retentate_pH = Constraint(
            expr=self.nanofiltration.pH["feed"] == self.nanofiltration.pH["retentate"]
        )
        if self.config.track_pE:
            self.eq_retentate_pE = Constraint(
                expr=self.nanofiltration.pE["feed"]
                == self.nanofiltration.pE["retentate"]
            )

    def add_reaktoro_chemistry(self):
        if self.config.add_alkalinity:
            self.nanofiltration.alkalinity = Var(
                self.LOCATIONS,
                initialize=1,
                units=pyunits.mg / pyunits.L,
                doc="Alkalinity (mg/L as CaCO3) at inlet and retentate",
            )

        self.nanofiltration.scaling_tendency = Var(
            self.LOCATIONS,
            list(self.config.selected_scalants.keys()),
            initialize=1,
            units=pyunits.dimensionless,
            doc="Scaling tendency for each scalant at inlet and retentate",
        )

        self.nanofiltration.maximum_scaling_tendency = Var(
            self.LOCATIONS,
            list(self.config.selected_scalants.keys()),
            initialize=lambda blk, loc, scalant: self.config.selected_scalants[scalant],
            units=pyunits.dimensionless,
            doc="Maximum allowable scaling tendency for each scalant at feed and retentate",
        )
        self.nanofiltration.maximum_scaling_tendency.fix()

        @self.nanofiltration.Constraint(
            self.LOCATIONS, list(self.config.selected_scalants.keys())
        )
        def eq_max_scaling_tendency(blk, loc, scalant):
            return (
                blk.scaling_tendency[loc, scalant]
                <= blk.maximum_scaling_tendency[loc, scalant]
            )

        self.deactivate_scaling_constraints()

        # Build reaktoro_outputs dict
        self.nanofiltration.reaktoro_outputs = {}

        # Osmotic pressure
        self.nanofiltration.reaktoro_outputs[("feed", "osmoticPressure", "H2O")] = (
            self.nanofiltration.feed_side.properties_in[0].pressure_osm_phase["Liq"]
        )
        self.nanofiltration.reaktoro_outputs[
            ("retentate", "osmoticPressure", "H2O")
        ] = self.nanofiltration.feed_side.properties_out[0].pressure_osm_phase["Liq"]

        # Scaling tendency
        for (loc, scalant), obj in self.nanofiltration.scaling_tendency.items():
            self.nanofiltration.reaktoro_outputs[(loc, "scalingTendency", scalant)] = (
                obj
            )

        # pH
        self.nanofiltration.reaktoro_outputs[("retentate", "pH", None)] = (
            self.nanofiltration.pH["retentate"]
        )

        # pE
        if self.config.track_pE:
            for loc in self.LOCATIONS:
                self.nanofiltration.reaktoro_outputs[(loc, "pE", None)] = (
                    self.nanofiltration.pE[loc]
                )

        # Optional alkalinity
        if self.config.add_alkalinity:
            for loc in self.LOCATIONS:
                self.nanofiltration.reaktoro_outputs[
                    (loc, "alkalinityAsCaCO3", None)
                ] = self.nanofiltration.alkalinity[loc]

        # Build ReaktoroOptionsContainer
        self.reaktoro_options = ReaktoroOptionsContainer()
        self.reaktoro_options.system_state_option(
            "temperature",
            self.nanofiltration.feed_side.properties_in[0].temperature,
        )
        self.reaktoro_options.system_state_option(
            "pressure",
            {
                "feed": self.nanofiltration.feed_side.properties_in[0].pressure,
                "retentate": self.nanofiltration.feed_side.properties_out[0].pressure,
            },
        )
        self.reaktoro_options.system_state_option(
            "pH",
            self.nanofiltration.pH["feed"],
        )
        if self.config.track_pE:
            self.reaktoro_options.system_state_option(
                "pE",
                self.nanofiltration.pE["feed"],
            )

        comp_map = {}
        for ion in self.config.default_property_package.component_list:
            comp_map["feed", ion] = self.nanofiltration.feed_side.properties_in[
                0
            ].flow_mol_phase_comp["Liq", ion]
            comp_map["retentate", ion] = self.nanofiltration.feed_side.properties_out[
                0
            ].flow_mol_phase_comp["Liq", ion]
        self.reaktoro_options.aqueous_phase_option("composition", comp_map)

        self.reaktoro_options["outputs"] = self.nanofiltration.reaktoro_outputs
        self.reaktoro_options.update_with_user_options(self.config.reaktoro_options)

        self.scaling_block = ReaktoroBlock(
            self.LOCATIONS,
            **self.reaktoro_options,
        )

        self.nanofiltration.feed_side.properties_in[0].eq_pressure_osm_phase[
            "Liq"
        ].deactivate()
        self.nanofiltration.feed_side.properties_out[0].eq_pressure_osm_phase[
            "Liq"
        ].deactivate()

    def scale_before_initialization(self, flow_factor=1, **kwargs):
        iscale.set_scaling_factor(self.nanofiltration.pH, 1)
        if self.config.track_pE:
            iscale.set_scaling_factor(self.nanofiltration.pE, 1)

        if self.config.add_reaktoro_chemistry:
            iscale.set_scaling_factor(self.nanofiltration.scaling_tendency, 1)
            if self.config.add_alkalinity:
                iscale.set_scaling_factor(self.nanofiltration.alkalinity, 1e-2)
            for loc in self.LOCATIONS:
                if loc == "feed":
                    osm_var = self.nanofiltration.feed_side.properties_in[
                        0
                    ].pressure_osm_phase["Liq"]
                else:
                    osm_var = self.nanofiltration.feed_side.properties_out[
                        0
                    ].pressure_osm_phase["Liq"]
                iscale.set_scaling_factor(osm_var, 1e-5)

        prop_pkg = self.config.default_property_package
        props_in = self.nanofiltration.feed_side.properties_in[0]
        props_out = self.nanofiltration.feed_side.properties_out[0]
        if hasattr(props_in, "flow_mol_phase_comp"):
            for (var, j), sf in prop_pkg._default_scaling_factors.items():
                if var == "flow_mol_phase_comp":
                    iscale.set_scaling_factor(props_in.flow_mol_phase_comp[j], sf)
        if hasattr(props_out, "flow_mol_phase_comp"):
            for (var, j), sf in prop_pkg._default_scaling_factors.items():
                if var == "flow_mol_phase_comp":
                    iscale.set_scaling_factor(props_out.flow_mol_phase_comp[j], sf)

        if self.config.add_feed_pump:
            self.pump.scale_before_initialization()

    def activate_scaling_constraints(self):
        if self.config.add_reaktoro_chemistry:
            for loc in self.LOCATIONS:
                for scalant in self.config.selected_scalants.keys():
                    self.nanofiltration.eq_max_scaling_tendency[loc, scalant].activate()
                    self.nanofiltration.scaling_tendency[loc, scalant].unfix()
                    self.nanofiltration.maximum_scaling_tendency[loc, scalant].fix()

    def deactivate_scaling_constraints(self):
        if self.config.add_reaktoro_chemistry:
            for loc in self.LOCATIONS:
                for scalant in self.config.selected_scalants.keys():
                    self.nanofiltration.eq_max_scaling_tendency[
                        loc, scalant
                    ].deactivate()
                    self.nanofiltration.maximum_scaling_tendency[loc, scalant].fix()

    def set_fixed_operation(self):
        if self.config.add_feed_pump:
            self.pump.set_fixed_operation()
        self.nanofiltration.flux_vol_solvent[0, "H2O"].fix(1e-5)
        self.nanofiltration.recovery_vol_phase[0, "Liq"].fix(0.5)
        if self.config.component_rejections is not None:
            self.nanofiltration.rejection_phase_comp[...].fix(0.5)
            for ion, rej in self.config.component_rejections.items():
                self.nanofiltration.rejection_phase_comp[0, "Liq", ion].fix(rej)
        else:
            self.nanofiltration.rejection_phase_comp[...].fix(0.5)

        if self.config.add_reaktoro_chemistry:
            self.nanofiltration.pH["feed"].fix(7)

    def set_optimization_operation(self):
        if self.config.add_feed_pump:
            self.pump.set_optimization_operation()
        self.nanofiltration.area.unfix()
        self.nanofiltration.area.setlb(1)

        self.nanofiltration.recovery_vol_phase[0.0, "Liq"].unfix()
        self.nanofiltration.recovery_vol_phase[0.0, "Liq"].setlb(0.25)
        self.nanofiltration.recovery_vol_phase[0.0, "Liq"].setub(0.95)

        self.activate_scaling_constraints()

    def initialize_unit(self, **kwargs):
        if self.config.add_feed_pump:
            self.pump.initialize()
            propagate_state(self.pump_to_nf)
        self.nanofiltration.initialize()

        if self.config.add_reaktoro_chemistry:
            for blk_idx, blk_obj in self.scaling_block.items():
                blk_obj.initialize()
                blk_obj.display_jacobian_scaling()

            self.nanofiltration.initialize()

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

        def get_pE(loc):
            if self.config.track_pE:
                return self.nanofiltration.pE[loc]
            return None

        unit_dofs = degrees_of_freedom(self)

        feed_state = get_ion_comp(
            self.nanofiltration.feed_side.properties_in[0],
            self.nanofiltration.pH["feed"],
            get_pE("feed"),
        )
        retentate_state = get_ion_comp(
            self.nanofiltration.feed_side.properties_out[0],
            self.nanofiltration.pH["retentate"],
            get_pE("retentate"),
        )

        model_state_dict = {
            "Model": {"DOFs": unit_dofs},
            "Feed inlet state": feed_state,
            "Retentate state": retentate_state,
            "Feed mol flow": self.nanofiltration.feed_side.properties_in[
                0
            ].flow_mol_phase_comp,
            "Retentate mol flow": self.nanofiltration.feed_side.properties_out[
                0
            ].flow_mol_phase_comp,
            "Permeate mol flow": self.nanofiltration.properties_permeate[
                0
            ].flow_mol_phase_comp,
            "Operation": {
                "Membrane area": self.nanofiltration.area,
                "Recovery (vol)": self.nanofiltration.recovery_vol_phase[0, "Liq"],
                "Water flux": self.nanofiltration.flux_vol_solvent[0, "H2O"],
            },
        }

        if self.config.add_reaktoro_chemistry:
            model_state_dict["Osmotic pressure"] = {
                "Feed (Pa)": self.nanofiltration.feed_side.properties_in[
                    0
                ].pressure_osm_phase["Liq"],
                "Retentate (Pa)": self.nanofiltration.feed_side.properties_out[
                    0
                ].pressure_osm_phase["Liq"],
            }
            scaling_dict = {}
            for (loc, scalant), obj in self.nanofiltration.scaling_tendency.items():
                scaling_dict[f"{scalant} ({loc})"] = obj
            model_state_dict["Scaling tendencies"] = scaling_dict
            if self.config.add_alkalinity:
                model_state_dict["Alkalinity"] = {
                    f"({loc})": obj
                    for loc, obj in self.nanofiltration.alkalinity.items()
                }

        if self.config.default_costing_package is not None:
            model_state_dict["Costs"] = {
                "Capital cost": self.nanofiltration.costing.capital_cost
            }

        if self.config.add_feed_pump:
            self.pump.report()

        return model_state_dict
