import os
from pathlib import Path

from windIO import __path__ as wiop

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])
