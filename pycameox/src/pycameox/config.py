"""
This module provides constants and other package-wide stuff.

"""

from pathlib import Path
from typing import Dict, Counter, NewType, Union

import pandas as pd

# ### Type annotations
Filename = NewType('Filename', Path)
Id = NewType('Id', int)
Seq = NewType('Seq', str)
RunsSet = NewType('RunsSet', Dict[str, pd.DataFrame])
SampleSet = NewType('Sample', Dict[str, pd.DataFrame])







