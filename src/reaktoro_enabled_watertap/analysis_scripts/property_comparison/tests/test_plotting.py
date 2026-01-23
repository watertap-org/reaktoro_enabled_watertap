import pytest

from reaktoro_enabled_watertap.analysis_scripts.property_comparison.figure_generation import (
    prop_comp_plotting,
)


def test_prop_plotting():
    """Test cost breakdown plotting for BGW scenario."""
    try:
        prop_comp_plotting.main(show_figs=False)
    except Exception as e:
        pytest.fail(f"cost_breakdown_bgw.main() raised an exception: {e}")
