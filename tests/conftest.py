"""
Pytest configuration and fixtures for WIFA tests.

Provides:
- Pre-test cleanup of leftover output directories
- Output directory fixtures with conditional cleanup (preserved on failure)
- FOXES engine fixture with proper initialization and teardown
"""

import shutil
from pathlib import Path

import pytest
from foxes import Engine

# Handle different foxes versions - reset_engine location varies
try:
    from foxes import reset_engine
except ImportError:
    from foxes.core import reset_engine

# Store test outcomes for conditional cleanup
_test_outcomes = {}


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Track test outcomes to conditionally preserve output on failure."""
    outcome = yield
    report = outcome.get_result()
    if report.when == "call":
        _test_outcomes[item.nodeid] = report.failed


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )


@pytest.fixture(scope="session", autouse=True)
def cleanup_old_outputs():
    """Remove output directories from previous test runs at session start."""
    patterns = [
        "output_pywake_*",
        "output_test_*",
        "output",
    ]
    for pattern in patterns:
        for path in Path(".").glob(pattern):
            if path.is_dir():
                shutil.rmtree(path)
            elif path.is_file():
                path.unlink()
    yield


@pytest.fixture
def output_dir(request, tmp_path):
    """
    Provide a unique temporary output directory for tests.

    Cleans up automatically on test success, preserves on failure for debugging.
    """
    yield tmp_path

    # Only clean up if test passed
    test_failed = _test_outcomes.get(request.node.nodeid, False)
    if not test_failed:
        if tmp_path.exists():
            shutil.rmtree(tmp_path)


@pytest.fixture
def named_output_dir(request):
    """
    Provide a named output directory based on test name.

    Use this when tests need output in the current working directory
    (e.g., for Code Saturne integration).

    Cleans up automatically on test success, preserves on failure.
    """
    test_name = request.node.name.replace("[", "_").replace("]", "_").rstrip("_")
    output_path = Path(f"output_test_{test_name}")
    output_path.mkdir(parents=True, exist_ok=True)

    yield output_path

    # Only clean up if test passed
    test_failed = _test_outcomes.get(request.node.nodeid, False)
    if not test_failed:
        if output_path.exists():
            shutil.rmtree(output_path)


@pytest.fixture
def cleanup_output_dir(request):
    """
    Cleanup fixture for tests that write to the default 'output/' directory.

    Use this when tests can't control their output location (e.g., when YAML
    specifies output_folder). Cleans up on success, preserves on failure.
    """
    output_path = Path("output")

    yield output_path

    # Only clean up if test passed
    test_failed = _test_outcomes.get(request.node.nodeid, False)
    if not test_failed:
        if output_path.exists():
            shutil.rmtree(output_path)


@pytest.fixture
def foxes_engine():
    """
    Provide a fresh FOXES engine instance per test.

    Properly initializes the engine before the test and resets it after,
    ensuring no state leaks between tests.
    """
    engine = Engine.new("default", verbosity=0)
    engine.initialize()

    yield engine

    # Always reset engine after test
    reset_engine()
