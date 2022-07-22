#
#     Copyright (C) 2022, Jose Manuel Martí Martínez (LLNL)
#
"""
Python library for CAMEOX: CAMEOs eXtended
"""

__all__ = ['config', 'downstream',
           '__author__', '__date__', '__version__']
__author__ = 'Jose Manuel Martí'
__copyright__ = 'Copyright (C) 2022 Jose Manuel Martí (LLNL)'
__license__ = 'UNLICENSED'
__maintainer__ = 'Jose Manuel Martí'
__status__ = 'Pre-Alpha'
__date__ = 'July 2022'
__version__ = '0.0.5'

import sys

# python
MAJOR, MINOR, *_ = sys.version_info
PYTHON_REL = (MAJOR == 3 and MINOR >= 8)
if not PYTHON_REL:
    raise ImportError('Recentrifuge requires Python 3.8 or later')
