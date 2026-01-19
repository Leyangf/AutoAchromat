"""Test script for src/io.py"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.io import print_summary, ensure_out_dir, extract_rows, save_csv, print_table
from src.score import Candidate


@dataclass
class MockGlass:
    name: str
    catalog: str = "TEST"
    relative_cost: float | None = None


# Create mock glasses
g1 = MockGlass("N-BK7", "SCHOTT", 1.0)
g2 = MockGlass("SF2", "SCHOTT", 2.5)
g3 = MockGlass("FK51", "SCHOTT", None)  # No cost data

# Create candidates (g1/g2 are dynamically attached, as done in synthesis.py)
c1: Any = Candidate("cemented", P2=0.001, W=0.05, R2=50.0, PE=0.01)
c1.g1 = g1
c1.g2 = g2
c1.meta = {
    "phi1": 0.01,
    "phi2": -0.005,
    "n1": 1.5,
    "n2": 1.7,
    "nu1": 60,
    "nu2": 30,
    "delta_nu": 30,
}

c2: Any = Candidate("cemented", P2=0.002, W=0.03, R2=40.0, PE=0.02)
c2.g1 = g1
c2.g2 = g3  # g3 has no cost
c2.meta = {}

c3: Any = Candidate("air_spaced", Pi=[0.1, 0.2, 0.15, 0.12], PE=0.15)
c3.g1 = g3
c3.g2 = g2
c3.meta = {}

# Test ensure_out_dir
print("Testing ensure_out_dir...")
path = ensure_out_dir("test_output")
print(f"  Created/verified: {path}")
print()

# Test extract_rows and save_csv
print("Testing extract_rows and save_csv...")
rows = extract_rows([c1, c2, c3])
save_csv(rows, "test_output/io_test.csv")
print()

# Display CSV contents
print("CSV Output:")
print("-" * 80)
with open("test_output/io_test.csv", "r") as f:
    print(f.read())
print("-" * 80)
print()

# Test print_table
print("Testing print_table...")
print()
print_table(rows, max_rows=2)
print()

# Test print_summary
print("Testing print_summary...")
print()
print_summary(
    {
        "glasses_total": 100,
        "glasses_filtered": 80,
        "pairs_tested": 6400,
        "candidates_generated": 500,
        "candidates_saved": 3,
        "best_PE": 0.01,
        "custom_stat": 42,
    }
)
