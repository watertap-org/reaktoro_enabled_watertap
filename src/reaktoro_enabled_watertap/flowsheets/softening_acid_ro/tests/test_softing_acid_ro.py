from reaktoro_enabled_watertap.flowsheets.softening_acid_ro import (
    softening_acid_ro as sar,
)
from pyomo.environ import (
    units as pyunits,
)
from pyomo.environ import value
import pytest
from idaes.core.util.model_statistics import degrees_of_freedom


__author__ = "Alexander V. Dudchenko"


@pytest.mark.parametrize("costing", ["watertap_default", "Amusat_et_al_2024"])
def test_softening_acid_ro_seawater_with_hpro(costing):
    solution_results = {
        "watertap_default": {
            "lcow": 1.0472,
            "na2co3_dose": 0.28775,
            "cao_dose": 0.062542,
            "h2so4_dose": 0.0039079,
            "hcl_dose": 1.0018e-05,
        },
        "Amusat_et_al_2024": {
            "lcow": 0.94351,
            "na2co3_dose": 0.29235,
            "cao_dose": 0.062560,
            "h2so4_dose": 0.0039292,
            "hcl_dose": 1.0018e-05,
        },
    }
    m = sar.build_model(
        "Seawater.yaml",
        multi_process_reaktoro=True,
        hpro=True,
        rkt_hessian_type="LBFGS",
        bfgs_initialization_type="GaussNewton",
        feed_flow_rate=5000 * pyunits.m**3 / pyunits.day,
        system_costing=costing,
    )
    sar.initialize(m, tee=True)
    assert degrees_of_freedom(m) == 9
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
        == solution_results[costing]["lcow"]
    )
    assert (
        pytest.approx(
            value(m.fs.softening_unit.precipitation_reactor.reagent_dose["CaO"]),
            abs=1e-4,
        )
        == solution_results[costing]["cao_dose"]
    )
    assert (
        pytest.approx(
            value(m.fs.softening_unit.precipitation_reactor.reagent_dose["Na2CO3"]),
            abs=1e-4,
        )
        == solution_results[costing]["na2co3_dose"]
    )

    assert (
        pytest.approx(
            value(m.fs.acidification_unit.chemical_reactor.reagent_dose["H2SO4"]),
            abs=1e-4,
        )
        == solution_results[costing]["h2so4_dose"]
    )
    assert (
        pytest.approx(
            value(m.fs.acidification_unit.chemical_reactor.reagent_dose["HCl"]),
            abs=1e-4,
        )
        == solution_results[costing]["hcl_dose"]
    )


@pytest.mark.parametrize(
    "water",
    [
        "USDA_brackish.yaml",
        "sample_500_hardness.yaml",
        "sample_1500_hardness.yaml",
        "Seawater.yaml",
    ],
)
@pytest.mark.component
def test_softening_acid_ro_default_costing(water):
    solution_results = {
        "USDA_brackish.yaml": {
            "lcow": 0.47618,
            "na2co3_dose": 0.10193,
            "cao_dose": 0.19797,
            "h2so4_dose": 0.0013513,
            "hcl_dose": 1.0006e-05,
        },
        "sample_500_hardness.yaml": {
            "lcow": 0.37896,
            "na2co3_dose": 0.16254,
            "cao_dose": 0.033480,
            "h2so4_dose": 0.0013017,
            "hcl_dose": 1.0006e-05,
        },
        "sample_1500_hardness.yaml": {
            "lcow": 1.1834,
            "na2co3_dose": 1.0006,
            "cao_dose": 0.057183,
            "h2so4_dose": 0.0018188,
            "hcl_dose": 1.0006e-05,
        },
        "Seawater.yaml": {
            "lcow": 0.58700,
            "na2co3_dose": 1.0018e-05,
            "cao_dose": 1.0006e-05,
            "h2so4_dose": 0.012424,
            "hcl_dose": 1.0006e-05,
        },
    }
    m = sar.build_model(
        water,
        multi_process_reaktoro=True,
        hpro=False,
        rkt_hessian_type="LBFGS",
        bfgs_initialization_type="GaussNewton",
        feed_flow_rate=5000 * pyunits.m**3 / pyunits.day,
        system_costing="watertap_default",
    )
    sar.initialize(m, tee=True)
    assert degrees_of_freedom(m) == 6
    if water == "Seawater.yaml":
        m.fs.water_recovery.fix(0.6)
    else:
        m.fs.water_recovery.fix(0.7)
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
        == solution_results[water]["lcow"]
    )
    assert (
        pytest.approx(
            value(m.fs.softening_unit.precipitation_reactor.reagent_dose["CaO"]),
            abs=1e-4,
        )
        == solution_results[water]["cao_dose"]
    )
    assert (
        pytest.approx(
            value(m.fs.softening_unit.precipitation_reactor.reagent_dose["Na2CO3"]),
            abs=1e-4,
        )
        == solution_results[water]["na2co3_dose"]
    )

    assert (
        pytest.approx(
            value(m.fs.acidification_unit.chemical_reactor.reagent_dose["H2SO4"]),
            abs=1e-4,
        )
        == solution_results[water]["h2so4_dose"]
    )
    assert (
        pytest.approx(
            value(m.fs.acidification_unit.chemical_reactor.reagent_dose["HCl"]),
            abs=1e-4,
        )
        == solution_results[water]["hcl_dose"]
    )


@pytest.mark.parametrize(
    "water",
    [
        "USDA_brackish.yaml",
        "sample_500_hardness.yaml",
        "sample_1500_hardness.yaml",
        "Seawater.yaml",
    ],
)
@pytest.mark.component
def test_softening_acid_ro_amusat_costing(water):
    solution_results = {
        "USDA_brackish.yaml": {
            "lcow": 0.54469,
            "na2co3_dose": 0.10301,
            "cao_dose": 0.19797,
            "h2so4_dose": 0.0013561,
            "hcl_dose": 1.0006e-05,
        },
        "sample_500_hardness.yaml": {
            "lcow": 0.45384,
            "na2co3_dose": 0.16452,
            "cao_dose": 0.033491,
            "h2so4_dose": 0.0013135,
            "hcl_dose": 1.0006e-05,
        },
        "sample_1500_hardness.yaml": {
            "lcow": 1.0618,
            "na2co3_dose": 1.0074,
            "cao_dose": 0.057245,
            "h2so4_dose": 0.0018910,
            "hcl_dose": 1.0006e-05,
        },
        "Seawater.yaml": {
            "lcow": 0.42629,
            "na2co3_dose": 1.0018e-05,
            "cao_dose": 1.0006e-05,
            "h2so4_dose": 0.012710,
            "hcl_dose": 1.0006e-05,
        },
    }
    m = sar.build_model(
        water,
        multi_process_reaktoro=True,
        hpro=False,
        rkt_hessian_type="LBFGS",
        bfgs_initialization_type="GaussNewton",
        feed_flow_rate=5000 * pyunits.m**3 / pyunits.day,
        system_costing="Amusat_et_al_2024",
    )
    sar.initialize(m, tee=True)
    assert degrees_of_freedom(m) == 6
    if water == "Seawater.yaml":
        m.fs.water_recovery.fix(0.6)
    else:
        m.fs.water_recovery.fix(0.7)
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
        == solution_results[water]["lcow"]
    )
    assert (
        pytest.approx(
            value(m.fs.softening_unit.precipitation_reactor.reagent_dose["CaO"]),
            abs=1e-4,
        )
        == solution_results[water]["cao_dose"]
    )
    assert (
        pytest.approx(
            value(m.fs.softening_unit.precipitation_reactor.reagent_dose["Na2CO3"]),
            abs=1e-4,
        )
        == solution_results[water]["na2co3_dose"]
    )

    assert (
        pytest.approx(
            value(m.fs.acidification_unit.chemical_reactor.reagent_dose["H2SO4"]),
            abs=1e-4,
        )
        == solution_results[water]["h2so4_dose"]
    )
    assert (
        pytest.approx(
            value(m.fs.acidification_unit.chemical_reactor.reagent_dose["HCl"]),
            abs=1e-4,
        )
        == solution_results[water]["hcl_dose"]
    )
