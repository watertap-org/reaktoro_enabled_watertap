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
    value,
    units as pyunits,
)
from pyomo.common.config import ConfigValue
from watertap.unit_models.pressure_changer import EnergyRecoveryDevice
from idaes.core import (
    declare_process_block_class,
)
from idaes.core import UnitModelCostingBlock

from reaktoro_enabled_watertap.utils import scale_utils as scu
import idaes.core.util.scaling as iscale


@declare_process_block_class("MultiCompERDUnit")
class MultiCompERDUnitData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "erd_efficiency",
        ConfigValue(
            default=0.9,
            description="default energy recovery device efficiency",
            doc="""
            default energy recovery device efficiency
            """,
        ),
    )
    CONFIG.declare(
        "erd_outlet_pressure",
        ConfigValue(
            default=1 * pyunits.atm,
            description="ERD outlet pressure",
            doc="""
            Energy recovery device pressure

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

    def build(self):
        """Build a multi-component ERD unit model"""
        super().build()
        self.ERD = EnergyRecoveryDevice(
            property_package=self.config.default_property_package
        )
        if self.config.default_costing_package is not None:
            self.ERD.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package,
                **self.config.default_costing_package_kwargs
            )
        # Add ERD flow rate
        self.ERD.control_volume.properties_in[0].flow_vol_phase[...]
        self.ERD.pH = Var(initialize=7, bounds=(0, 13), units=pyunits.dimensionless)
        tracked_vars = {"pH": self.ERD.pH}
        if self.config.track_pE:
            self.ERD.pE = Var(
                initialize=0,
                units=pyunits.dimensionless,
                bounds=(None, None),
            )
            tracked_vars["pE"] = self.ERD.pE
        iscale.set_scaling_factor(self.ERD.pH, 1)
        self.register_port("inlet", self.ERD.inlet, tracked_vars)
        self.register_port("outlet", self.ERD.outlet, tracked_vars)

    def set_fixed_operation(self):
        """fixes operation point for ERD unit model
        Uses osmotic pressure to initialize ERD outlet pressure or user defined pressure
        """
        self.ERD.efficiency_pump[0].fix(self.config.erd_efficiency)
        self.ERD.control_volume.properties_out[0].pressure.fix(
            self.config.erd_outlet_pressure
        )

    def scale_before_initialization(self, **kwargs):
        h2o_scale = self.config.default_property_package._default_scaling_factors[
            "flow_mol_phase_comp", ("Liq", "H2O")
        ] / value(self.config.default_property_package.mw_comp["H2O"])
        work_scale = 1e-5 / (1 / h2o_scale)
        iscale.set_scaling_factor(self.ERD.inlet.pressure, 1e-5)
        iscale.set_scaling_factor(self.ERD.outlet.pressure, 1e-5)
        iscale.set_scaling_factor(self.ERD.control_volume.work, work_scale)
        iscale.set_scaling_factor(self.ERD.pH, 1)
        if self.config.default_costing_package is not None:
            scu.calculate_scale_from_dependent_vars(
                self.ERD.costing.capital_cost,
                self.ERD.costing.capital_cost_constraint,
                [self.ERD.control_volume.work[0]],
            )
        if self.config.track_pE:
            iscale.set_scaling_factor(self.ERD.pE, 1)

    def initialize_unit(self):
        self.ERD.initialize()

    def get_model_state_dict(self):
        """Returns a dictionary with the model state"""
        self.inlet.fix()
        unit_dofs = degrees_of_freedom(self)
        self.inlet.unfix()

        model_state_dict = {
            "Model": {"DOFs": unit_dofs},
            "Overall": {
                "pH": self.ERD.pH,
                "Temperature": self.ERD.inlet.temperature[0],
                "Flow rate": self.ERD.control_volume.properties_in[0].flow_vol_phase[
                    "Liq"
                ],
            },
            "Inlet": {
                "Pressure": self.ERD.inlet.pressure[0],
            },
            "Outlet": {
                "Pressure": self.ERD.outlet.pressure[0],
            },
        }
        if self.config.default_costing_package is not None:
            model_state_dict["Costs"] = {
                "Capital cost": self.ERD.costing.capital_cost,
                # "Operating cost": self.ERD.costing.operating_cost,
            }
        return self.name, model_state_dict
