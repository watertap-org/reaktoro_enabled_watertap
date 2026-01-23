import pytest

from reaktoro_enabled_watertap.analysis_scripts.softening_acid_ro.figure_generation import (
    cost_breakdown_bgw,
    cost_breakdown_hpro,
    plotting_performance_regions,
    validation_plotting,
)


def test_cost_breakdown_bgw():
    """Test cost breakdown plotting for BGW scenario."""
    try:
        cost_breakdown_bgw.main(show_figs=False)
    except Exception as e:
        pytest.fail(f"cost_breakdown_bgw.main() raised an exception: {e}")


def test_cost_breakdown_hpro():
    """Test cost breakdown plotting for HPro scenario."""
    try:
        cost_breakdown_hpro.main(show_figs=False)
    except Exception as e:
        pytest.fail(f"cost_breakdown_hpro.main() raised an exception: {e}")


def test_plotting_performance_regions():
    """Test plotting performance regions."""
    try:
        plotting_performance_regions.main(show_figs=False)
    except Exception as e:
        pytest.fail(f"plotting_performance_regions.main() raised an exception: {e}")
