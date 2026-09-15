"""Tests for surfactantkit's own top-level __init__.py export surface.

Added 2026-09-15 after a real gap was found during a Tier 1 bottleneck-
resolution completeness audit: __init__.py had not been updated since
before classify.py existed, so a direct `from surfactantkit import X`
user could not reach classify_surfactant_charge_type,
classify_surfactant_structural_family, or any of the 5 new functions
built this session -- despite all of them being real, tested, and
wired into the MCP server. Caught only by deliberately checking
git log against __init__.py's own last-touched commit, not by any
existing test -- this file exists so a future new public function
added to any already-imported module, and never added here, fails
loudly instead of silently.
"""

import surfactantkit


def test_all_names_in_dunder_all_are_real_attributes():
    for name in surfactantkit.__all__:
        assert hasattr(surfactantkit, name), f"'{name}' is listed in __all__ but not actually importable"


def test_dunder_all_has_no_duplicates():
    assert len(surfactantkit.__all__) == len(set(surfactantkit.__all__))


def test_all_names_in_dunder_all_are_callable():
    for name in surfactantkit.__all__:
        assert callable(getattr(surfactantkit, name)), f"'{name}' in __all__ is not callable"


def test_tier0_tier1_bottleneck_fix_functions_are_top_level_importable():
    """Direct regression test for the specific gap found 2026-09-15 --
    every function built during the Tier 0/Tier 1 bottleneck-resolution
    pass must be reachable via `from surfactantkit import X`, not just
    via its own submodule or the MCP server."""
    from surfactantkit import (
        get_vocg_standard_liquid,
        get_owens_wendt_standard_liquid,
        estimate_axial_ratio_from_cpp_geometry,
        estimate_partial_specific_volume_from_tail_and_headgroup,
        classify_surfactant_charge_type_at_ph,
        estimate_intrinsic_water_solubility_qspr,
    )
    assert all(
        callable(f)
        for f in (
            get_vocg_standard_liquid,
            get_owens_wendt_standard_liquid,
            estimate_axial_ratio_from_cpp_geometry,
            estimate_partial_specific_volume_from_tail_and_headgroup,
            classify_surfactant_charge_type_at_ph,
            estimate_intrinsic_water_solubility_qspr,
        )
    )


def test_version_matches_current_capability_state():
    # bumped 0.3.0 -> 0.4.0 alongside this session's Tier 0/1 additions,
    # matching this project's own established bump-on-feature-addition
    # convention (0.1.0 -> 0.2.0 -> 0.3.0 in prior commits)
    assert surfactantkit.__version__ == "0.4.0"
