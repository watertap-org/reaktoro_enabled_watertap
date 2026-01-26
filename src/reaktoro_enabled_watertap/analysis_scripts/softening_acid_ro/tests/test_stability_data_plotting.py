import pytest

from reaktoro_enabled_watertap.analysis_scripts.softening_acid_ro.figure_generation import (
    stablity_plotting,
)


@pytest.mark.analysis
@pytest.mark.component
def test_validation_plotting():
    """Test cost breakdown plotting for BGW scenario."""
    try:
        stablity_plotting.main(show_figs=False)
    except Exception as e:
        pytest.fail(f"cost_breakdown_bgw.main() raised an exception: {e}")
