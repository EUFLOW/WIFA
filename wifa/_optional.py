import importlib.util


def require(package_name: str, extra_name: str) -> None:
    """Raise a clear ImportError if an optional dependency is missing."""
    if importlib.util.find_spec(package_name) is None:
        raise ImportError(
            f"'{package_name}' is required for this functionality but is not installed. "
            f"Install it with: pip install wifa[{extra_name}]"
        )
