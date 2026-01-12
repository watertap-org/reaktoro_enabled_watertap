from reaktoro_enabled_watertap.flowsheets.softening_acid_ro import (
    softening_acid_ro as sar,
)


__author__ = "Alexander V. Dudchenko"
from pyomo.environ import value
import pytest


@pytest.mark.component
def test_softening_acid_ro_seawater_with_hpro():
    m = sar.build_model(
        "Seawater.yaml",
        multi_process_reaktoro=True,
        hpro=True,
        rkt_hessian_type="LBFGS",
        bfgs_initialization_type="GaussNewton",
    )
    sar.initialize(m)
    m.fs.water_recovery.fix(0.75)
    sar.solve_model(m, tee=True)
    m.fs.water_recovery.fix(0.78)
    sar.solve_model(m, tee=True)
    m.fs.water_recovery.fix(0.8)
    sar.solve_model(m, tee=True)
    sar.report_all_units(m)
    m.reaktoro_manager.terminate_workers()
    assert (
        pytest.approx(
            value(m.fs.costing.LCOW),
            1e-3,
        )
        == 1.0542
    )
    assert (
        pytest.approx(
            value(m.fs.softening_unit.precipitation_reactor.reagent_dose["CaO"]),
            1e-3,
        )
        == 0.062647
    )
    assert (
        pytest.approx(
            value(m.fs.softening_unit.precipitation_reactor.reagent_dose["Na2CO3"]),
            1e-3,
        )
        == 0.31345
    )

    assert (
        pytest.approx(
            value(m.fs.acidification_unit.chemical_reactor.reagent_dose["H2SO4"]),
            1e-3,
        )
        == 0.0040669
    )
    assert (
        pytest.approx(
            value(m.fs.acidification_unit.chemical_reactor.reagent_dose["HCl"]),
            1e-3,
        )
        == 1.0018e-05
    )


@pytest.mark.component
def test_softening_acid_ro_bgw():
    m = sar.build_model(
        "USDA_brackish.yaml",
        multi_process_reaktoro=True,
        hpro=False,
        rkt_hessian_type="LBFGS",
        bfgs_initialization_type="GaussNewton",
    )
    sar.initialize(m)
    m.fs.water_recovery.fix(0.65)
    sar.solve_model(m, tee=True)
    m.fs.water_recovery.fix(0.725)
    sar.solve_model(m, tee=True)
    m.fs.water_recovery.fix(0.8)
    sar.solve_model(m, tee=True)

    sar.report_all_units(m)
    m.reaktoro_manager.terminate_workers()
    assert (
        pytest.approx(
            value(m.fs.costing.LCOW),
            1e-3,
        )
        == 0.45611
    )
    assert (
        pytest.approx(
            value(m.fs.softening_unit.precipitation_reactor.reagent_dose["CaO"]),
            1e-3,
        )
        == 0.19795
    )
    assert (
        pytest.approx(
            value(m.fs.softening_unit.precipitation_reactor.reagent_dose["Na2CO3"]),
            1e-3,
        )
        == 0.096213
    )

    assert (
        pytest.approx(
            value(m.fs.acidification_unit.chemical_reactor.reagent_dose["H2SO4"]), 1e-3
        )
        == 0.0013254
    )
    assert (
        pytest.approx(
            value(m.fs.acidification_unit.chemical_reactor.reagent_dose["HCl"]),
            1e-3,
        )
        == 1.0018e-05
    )
