import sys
from pathlib import Path

import windIO

windIO.load_yaml(Path(sys.argv[1]))
