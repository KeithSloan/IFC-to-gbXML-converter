"""
pytest configuration file.

This file helps pytest find the ifcgbxml module and sets up common fixtures.
"""

import sys
import os
from pathlib import Path

# Add the project root to Python path so we can import ifcgbxml
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import the module to make sure it's available
import ifcgbxml
