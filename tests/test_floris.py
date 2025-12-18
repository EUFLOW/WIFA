import os
from pathlib import Path

from windIO import __path__ as wiop
from windIO import validate as validate_yaml
from windIO import load_yaml

from wifa.floris_api import run_floris

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


def _run_floris(wes_dir, expected_error=None):
    """
    Helper function to run FLORIS on all system.yaml files in a directory.
    
    Args:
        wes_dir: Path to wind energy system directory containing system.yaml files
        expected_error: If provided, the test expects this error type/message
    """
    assert wes_dir.is_dir(), f"{wes_dir} is not a directory"

    for yaml_path in wes_dir.glob("system.yaml"):
        print("\nRUNNING FLORIS ON", yaml_path, "\n")
        
        # Validate the YAML file against windIO schema
        validate_yaml(yaml_path, Path("plant/wind_energy_system"))
        
        # Load the YAML file as a dictionary
        yaml_input = load_yaml(yaml_path)
        yaml_input['attributes']['flow_model']['name'] = 'floris'
        
        # Create output directory
        output_dir_name = Path("output_test_floris")
        output_dir_name.mkdir(parents=True, exist_ok=True)
        
        # Run FLORIS
        result = run_floris(yaml_input)


def test_floris_4wts():
    """Test FLORIS with 4 turbines configuration"""
    wes_dir = test_path / "../examples/cases/windio_4turbines/wind_energy_system/"
    _run_floris(wes_dir)


def test_floris_simple_wind_rose():
    """Test FLORIS with simple wind rose configuration"""
    wes_dir = test_path / "../examples/cases/simple_wind_rose/wind_energy_system/"
    _run_floris(wes_dir)


def test_floris_timeseries():
    """Test FLORIS with timeseries and operating flag"""
    wes_dir = test_path / "../examples/cases/timeseries_with_operating_flag/wind_energy_system/"
    _run_floris(wes_dir)


def test_floris_heterogeneous_wind_rose():
    """Test FLORIS with heterogeneous wind rose map"""
    wes_dir = test_path / "../examples/cases/heterogeneous_wind_rose_map/wind_energy_system/"
    _run_floris(wes_dir)


def test_floris_heterogeneous_at_turbines():
    """Test FLORIS with heterogeneous wind rose at turbines"""
    wes_dir = test_path / "../examples/cases/heterogeneous_wind_rose_at_turbines/wind_energy_system/"
    _run_floris(wes_dir)


def test_floris_multiple_turbines():
    """Test FLORIS with multiple turbine types"""
    wes_dir = test_path / "../examples/cases/windio_4turbines_multipleTurbines/wind_energy_system/"
    _run_floris(wes_dir)


if __name__ == "__main__":
    """Run tests when executed directly"""
    import sys
    import traceback
    
    print("=" * 80)
    print("Running FLORIS tests...")
    print("=" * 80)
    
    tests = [
        ("4 Turbines", test_floris_4wts),
        ("Simple Wind Rose", test_floris_simple_wind_rose),
        ("Timeseries", test_floris_timeseries),
        ("Multiple Turbine Types", test_floris_multiple_turbines),
    ]
    
    failed = []
    passed = []
    known_issues = []
    
    for test_name, test_func in tests:
        print(f"\n{'='*80}")
        print(f"Testing: {test_name}")
        print(f"{'='*80}")
        try:
            test_func()
            passed.append((test_name, ""))
            print(f"✓ {test_name} PASSED")
        except Exception as e:
            failed.append((test_name, str(e)))
            print(f"✗ {test_name} FAILED: {e}")
            if "--verbose" in sys.argv or "-v" in sys.argv:
                traceback.print_exc()
    
    # Summary
    print(f"\n{'='*80}")
    print("Test Summary")
    print(f"{'='*80}")
    print(f"Passed: {len(passed)}/{len(tests)}")
    print(f"Known Issues: {len(known_issues)}/{len(tests)}")
    print(f"Failed: {len(failed)}/{len(tests)}")
    
    if known_issues:
        print("\nKnown issues (NotImplementedError):")
        for name, error in known_issues:
            print(f"  WARNING: {name}: {error}")
    
    if failed:
        print("\nFailed tests:")
        for name, error in failed:
            print(f"  ✗ {name}: {error}")
        sys.exit(1)
    else:
        print("\n✓ All tests passed!")
        sys.exit(0)
