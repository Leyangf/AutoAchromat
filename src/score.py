"""
Preliminary Evaluation (PE) scoring module for AUTOACHROMAT.

This module implements the 2022 paper's PE metric for ranking achromat candidates.
At this stage, ranking and Top-N selection depend ONLY on PE.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass
class Candidate:
    """
    Represents an achromat candidate for preliminary evaluation.

    Attributes:
        system_type: "cemented" or "air_spaced"
        P2: Spherical aberration of cemented surface (cemented only)
        W: Coma coefficient (cemented only)
        R2: Radius of cemented surface (cemented only)
        Pi: Spherical aberration coefficients of all surfaces (air-spaced only)
        PE: Preliminary Evaluation score (computed)
    """

    system_type: str
    P2: float | None = None
    W: float | None = None
    R2: float | None = None
    Pi: list[float] = field(default_factory=list)
    PE: float | None = None


def compute_pe(candidate: Candidate) -> float:
    """
    Compute the Preliminary Evaluation (PE) score according to the 2022 paper.

    For cemented doublets:
        PE = |P2| * (3 ^ |W|) / (R2 ^ 2)

    For air-spaced doublets:
        PE = mean(|Pi|)

    Args:
        candidate: The Candidate to evaluate.

    Returns:
        The PE score. Returns +inf if required values are missing or invalid.
        Also stores the result in candidate.PE.
    """
    if candidate.system_type == "cemented":
        # Check for missing or invalid values
        if candidate.P2 is None or candidate.W is None or candidate.R2 is None:
            candidate.PE = math.inf
            return math.inf

        # Check for zero R2 (would cause division by zero)
        if candidate.R2 == 0:
            candidate.PE = math.inf
            return math.inf

        try:
            pe = abs(candidate.P2) * (3 ** abs(candidate.W)) / (candidate.R2**2)
        except (ValueError, OverflowError):
            candidate.PE = math.inf
            return math.inf

        candidate.PE = pe
        return pe

    elif candidate.system_type == "air_spaced":
        # Check for empty or invalid Pi list
        if not candidate.Pi or len(candidate.Pi) == 0:
            candidate.PE = math.inf
            return math.inf

        try:
            pe = sum(abs(p) for p in candidate.Pi) / len(candidate.Pi)
        except (TypeError, ValueError):
            candidate.PE = math.inf
            return math.inf

        candidate.PE = pe
        return pe

    else:
        # Unknown system type
        candidate.PE = math.inf
        return math.inf


def rank_key(candidate: Candidate) -> float:
    """
    Return the ranking key for a candidate (PE only).

    Lower PE is better. None or invalid PE is treated as +inf.

    Args:
        candidate: The Candidate to get the rank key for.

    Returns:
        The PE value, or +inf if PE is None or invalid.
    """
    if candidate.PE is None:
        return math.inf

    try:
        pe = float(candidate.PE)
        if math.isnan(pe):
            return math.inf
        return pe
    except (TypeError, ValueError):
        return math.inf


class TopNContainer:
    """
    Container that maintains the top N candidates with lowest PE scores.
    """

    def __init__(self, N: int) -> None:
        """
        Initialize the container.

        Args:
            N: Maximum number of candidates to keep.
        """
        self._n = N
        self._items: list[Candidate] = []

    def push(self, candidate: Candidate) -> bool:
        """
        Attempt to add a candidate to the container.

        If the container is not full, the candidate is inserted.
        If the container is full, the candidate replaces the entry with
        the maximum PE if candidate.PE < max_PE.

        Args:
            candidate: The Candidate to potentially add.

        Returns:
            True if the candidate was stored, False if discarded.
        """
        if self._n <= 0:
            return False

        candidate_pe = rank_key(candidate)

        if len(self._items) < self._n:
            self._items.append(candidate)
            return True

        # Find the entry with maximum PE
        worst_idx = 0
        worst_pe = rank_key(self._items[0])

        for idx in range(1, len(self._items)):
            pe = rank_key(self._items[idx])
            if pe > worst_pe:
                worst_pe = pe
                worst_idx = idx

        # Replace if candidate is better
        if candidate_pe < worst_pe:
            self._items[worst_idx] = candidate
            return True

        return False

    def sorted(self) -> list[Candidate]:
        """
        Return candidates sorted by PE in ascending order (best first).

        Returns:
            List of candidates sorted by PE.
        """
        return sorted(self._items, key=rank_key)

    def __len__(self) -> int:
        """Return the number of candidates currently stored."""
        return len(self._items)
