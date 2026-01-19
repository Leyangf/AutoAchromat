"""
AutoAchromat - Automated achromat doublet lens synthesis.

Modules:
    glass_reader: AGF glass catalog parsing
    synthesis: Core analytical synthesis pipeline
    score: PE scoring and candidate ranking
    io: Input/output utilities
    cli: Command-line interface and config handling
"""

from . import glass_reader
from . import synthesis
from . import score
from . import io
from . import cli

__version__ = "0.1.0"
__all__ = ["glass_reader", "synthesis", "score", "io", "cli"]
