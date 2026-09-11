"""Track which windIO input keys a flow-model runner actually consumes.

Adapted from the TrackedDict concept in FLORIS's windIO reader
(https://github.com/lejeunemax/floris, floris/read_windio/utils.py): the
input dict is wrapped so that key accesses are recorded, and after the run
every key the runner never read is reported. This gives users a contract:
a parameter they set either influenced the simulation or produced a
warning saying it was ignored.

Unlike the FLORIS implementation (a UserDict), TrackedDict subclasses
dict so that the many isinstance(x, dict) checks in WIFA, windIO, and the
flow models keep working unchanged. Reads are recorded on __getitem__,
get, pop, items, and values (iterating a mapping's items or values counts
as consuming all of them, which is how windIO's dict_to_netcdf reads a
resource). C-level consumers that bypass these methods can only cause
under-reporting, i.e. a spurious "unread" warning, never a wrong
simulation.
"""

import warnings


def _wrap(value, context):
    if isinstance(value, TrackedDict):
        return value
    if isinstance(value, dict):
        return TrackedDict(value, context=context)
    if isinstance(value, list):
        return [_wrap(v, f"{context}[{i}]") for i, v in enumerate(value)]
    return value


class TrackedDict(dict):
    """A dict that records which of its keys have been read."""

    def __init__(self, data, context="wind_energy_system"):
        super().__init__(
            {key: _wrap(value, f"{context}.{key}") for key, value in data.items()}
        )
        self._context = context
        self._read = set()

    def __getitem__(self, key):
        if key in self:
            self._read.add(key)
        return super().__getitem__(key)

    def get(self, key, default=None):
        if key in self:
            self._read.add(key)
        return super().get(key, default)

    def pop(self, key, *default):
        if key in self:
            self._read.add(key)
        return super().pop(key, *default)

    def items(self):
        self._read.update(super().keys())
        return super().items()

    def values(self):
        self._read.update(super().keys())
        return super().values()

    def unread_paths(self):
        """Return the full paths of every key that was never read.

        A subtree whose root key was never read is reported as the root
        path only, matching how a user thinks about their input file.
        """
        paths = []
        for key in super().keys():
            child = super().__getitem__(key)
            if key not in self._read:
                paths.append(f"{self._context}.{key}")
            else:
                paths.extend(_collect_unread(child))
        return paths


def _collect_unread(value):
    if isinstance(value, TrackedDict):
        return value.unread_paths()
    if isinstance(value, list):
        paths = []
        for item in value:
            paths.extend(_collect_unread(item))
        return paths
    return []


def report_unread(tracked, model_name):
    """Warn about input keys the flow-model runner never consumed."""
    unread = tracked.unread_paths()
    if unread:
        warnings.warn(
            f"The following windIO input keys were not used by the "
            f"'{model_name}' runner and did not influence the results: "
            f"{unread}"
        )
    return unread
