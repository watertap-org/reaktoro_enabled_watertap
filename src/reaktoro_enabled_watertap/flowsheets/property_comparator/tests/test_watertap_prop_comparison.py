from reaktoro_enabled_watertap.flowsheets.property_comparator import (
    watertap_prop_comparison as wpc,
)
from pyomo.environ import (
    assert_optimal_termination,
)
import pytest
from idaes.core.util.model_statistics import degrees_of_freedom


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
def test_props(water):
    expected_results = {
        "USDA_brackish.yaml": {
            "Temperature (K)": 293.15,
            "TDS": 6.835665983881908,
            "Reaktoro_osmoticPressure (Pa)": 358518.87111541274,
            "Seawater_osmoticPressure (Pa)": 478462.72467433196,
            "NaCl_osmoticPressure( Pa)": 529506.6166892514,
            "Reaktoro_vaporPressure (Pa)": 2349.294364783399,
            "Seawater_vaporPressure (Pa)": 2331.2925100741427,
            "NaCl_vaporPressure (Pa)": 2039.6936689781523,
            "Reaktoro_density (kg/m3)": 1003.9775231984338,
            "Seawater_density (kg/m3)": 1003.2512267757132,
            "NaCl_density (kg/m3)": 1002.5416769691211,
            "H2O (mg/L)": 997141.8572145521,
            "Na_+ (mg/L)": 1477.2177267548957,
            "K_+ (mg/L)": 17.990472991602246,
            "Ca_2+ (mg/L)": 515.7268924259312,
            "Mg_2+ (mg/L)": 179.90472991602255,
            "Cl_- (mg/L)": 1854.3039066515955,
            "SO4_2- (mg/L)": 2020.9297993899868,
            "HCO3_- (mg/L)": 769.5924557518741,
        },
        "sample_500_hardness.yaml": {
            "Temperature (K)": 293.15,
            "TDS": 4.893645534302204,
            "Reaktoro_osmoticPressure (Pa)": 229676.47488595956,
            "Seawater_osmoticPressure (Pa)": 342409.5745541198,
            "NaCl_osmoticPressure( Pa)": 378868.67227459326,
            "Reaktoro_vaporPressure (Pa)": 2351.536794851882,
            "Seawater_vaporPressure (Pa)": 2333.460981422937,
            "NaCl_vaporPressure (Pa)": 2041.7834136823647,
            "Reaktoro_density (kg/m3)": 1002.6010450981555,
            "Seawater_density (kg/m3)": 1001.7692756897618,
            "NaCl_density (kg/m3)": 1001.1656934057822,
            "H2O (mg/L)": 997707.3995638536,
            "Ca_2+ (mg/L)": 334.7170115020679,
            "Cl_- (mg/L)": 745.0890679060753,
            "HCO3_- (mg/L)": 148.02329965052618,
            "K_+ (mg/L)": 10.99727337671071,
            "Mg_2+ (mg/L)": 39.83012466619589,
            "Na_+ (mg/L)": 1215.5836568273796,
            "SO4_2- (mg/L)": 2399.405100373247,
        },
        "sample_1500_hardness.yaml": {
            "Temperature (K)": 293.15,
            "TDS": 5.355685150322795,
            "Reaktoro_osmoticPressure (Pa)": 162647.1774283672,
            "Seawater_osmoticPressure (Pa)": 374767.7678237279,
            "NaCl_osmoticPressure( Pa)": 414691.5580779625,
            "Reaktoro_vaporPressure (Pa)": 2352.7042486981372,
            "Seawater_vaporPressure (Pa)": 2332.9475218362004,
            "NaCl_vaporPressure (Pa)": 2041.28684392959,
            "Reaktoro_density (kg/m3)": 1003.4896329098508,
            "Seawater_density (kg/m3)": 1002.1222527321229,
            "NaCl_density (kg/m3)": 1001.4933191434558,
            "H2O (mg/L)": 998133.9477595284,
            "Ca_2+ (mg/L)": 979.9665650824727,
            "Cl_- (mg/L)": 156.54988828686618,
            "HCO3_- (mg/L)": 239.99181185693212,
            "K_+ (mg/L)": 31.99890824759094,
            "Mg_2+ (mg/L)": 133.99542828678707,
            "Na_+ (mg/L)": 413.29854725560887,
            "SO4_2- (mg/L)": 3399.8840013065364,
        },
        "Seawater.yaml": {
            "Temperature (K)": 293.15,
            "TDS": 68.18345175258048,
            "Reaktoro_osmoticPressure (Pa)": 4970815.806700132,
            "Seawater_osmoticPressure (Pa)": 4989433.615615857,
            "NaCl_osmoticPressure( Pa)": 5493003.640986841,
            "Reaktoro_vaporPressure (Pa)": 2270.4128943290298,
            "Seawater_vaporPressure (Pa)": 2250.556038088748,
            "NaCl_vaporPressure (Pa)": 1969.128381662586,
            "Reaktoro_density (kg/m3)": 1050.6143611599573,
            "Seawater_density (kg/m3)": 1048.0102465083494,
            "NaCl_density (kg/m3)": 1044.7108412847863,
            "H2O (mg/L)": 982430.9094073769,
            "Na_+ (mg/L)": 20945.2633765292,
            "K_+ (mg/L)": 753.9977342820292,
            "Ca_2+ (mg/L)": 793.6818255600307,
            "Mg_2+ (mg/L)": 2504.066159641897,
            "Cl_- (mg/L)": 37652.49612785001,
            "SO4_2- (mg/L)": 5256.157889771304,
            "HCO3_- (mg/L)": 277.78863894601085,
        },
    }
    m = wpc.build_model(water_case=water)
    m.fs.water_recovery.fix(0.5)
    assert degrees_of_freedom(m) == 0
    result = wpc.solve_model(m)
    assert_optimal_termination(result)
    test_results = wpc.print_comparison(m)
    for key, val in expected_results[water].items():
        assert pytest.approx(test_results[key], rel=1e-3) == val
    print(test_results)
