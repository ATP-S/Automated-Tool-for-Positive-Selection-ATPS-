"""Color schemes for sequence visualization in phylogenetic trees and alignments."""

from __future__ import annotations

import logging
from collections.abc import Iterable

logger = logging.getLogger(__name__)

# Amino acid color scheme (based on physicochemical properties)
AMINO_ACID_COLORS: dict[str, str] = {
    # Polar positive (basic) - green
    "K": "#109c4b",
    "R": "#109c4b",
    "H": "#ff8e00",  # orange (special case)
    # Polar negative (acidic) - magenta
    "D": "#ea42fc",
    "E": "#109c4b",
    # Polar uncharged - blue
    "S": "#1e67b6",
    "T": "#1e67b6",
    "N": "#ea42fc",
    "Q": "#109c4b",
    # Hydrophobic - red
    "A": "#1e67b6",
    "V": "#d82626",
    "I": "#d82626",
    "L": "#d82626",
    "M": "#d82626",
    # Aromatic - yellow
    "F": "#fed700",
    "W": "#fed700",
    "Y": "#fed700",
    # Special
    "C": "#00a391",  # teal (cysteine)
    "G": "#8d4712",  # brown (glycine)
    "P": "#ffa9e3",  # pink (proline)
    # Gap
    "-": "#000000",
    ".": "#000000",
    "*": "#000000",
}

# Nucleotide color scheme
NUCLEOTIDE_COLORS: dict[str, str] = {
    "A": "#22ca00",  # green (Adenine)
    "T": "#d82626",  # red (Thymine)
    "U": "#d82626",  # red (Uracil - RNA)
    "G": "#fed700",  # yellow (Guanine)
    "C": "#1e67b6",  # blue (Cytosine)
    "N": "#888888",  # gray (unknown)
    "-": "#000000",  # black (gap)
    ".": "#000000",
}

# Default fallback color for unknown characters
DEFAULT_COLOR = "#cccccc"


def get_color_scheme(sequence_type: str = "aa") -> dict[str, str]:
    """Get the color scheme dictionary for a sequence type.

    Args:
        sequence_type: "aa" for amino acids, "nt" for nucleotides.

    Returns:
        Dictionary mapping characters to hex color codes.
    """
    if sequence_type.lower() in ("aa", "amino", "protein"):
        return AMINO_ACID_COLORS.copy()
    elif sequence_type.lower() in ("nt", "dna", "rna", "nucleotide"):
        return NUCLEOTIDE_COLORS.copy()
    else:
        logger.warning("Unknown sequence type '%s'; defaulting to amino acids.", sequence_type)
        return AMINO_ACID_COLORS.copy()


def get_colors(
    sequences: Iterable[str],
    sequence_type: str = "aa",
    color_scheme: dict[str, str] | None = None,
    default_color: str = DEFAULT_COLOR,
) -> list[str]:
    """Generate colors for each character in sequences for visualization.

    Args:
        sequences: Iterable of sequence strings (e.g., list of protein sequences).
        sequence_type: "aa" for amino acids, "nt" for nucleotides.
        color_scheme: Custom color mapping (overrides sequence_type if provided).
        default_color: Color for unknown characters.

    Returns:
        List of hex color codes, one for each character across all sequences.

    Example:
        >>> get_colors(["AKL", "MV-"], sequence_type="aa")
        ['#1e67b6', '#109c4b', '#d82626', '#d82626', '#d82626', '#000000']
    """
    if color_scheme is None:
        color_scheme = get_color_scheme(sequence_type)

    colors = []
    unknown_chars = set()

    for seq in sequences:
        for char in seq:
            upper_char = char.upper()
            if upper_char in color_scheme:
                colors.append(color_scheme[upper_char])
            else:
                colors.append(default_color)
                unknown_chars.add(char)

    if unknown_chars:
        logger.debug("Unknown characters assigned default color: %s", unknown_chars)

    return colors


def get_color(
    char: str,
    sequence_type: str = "aa",
    color_scheme: dict[str, str] | None = None,
    default_color: str = DEFAULT_COLOR,
) -> str:
    """Get the color for a single character.

    Args:
        char: Single character (amino acid or nucleotide).
        sequence_type: "aa" for amino acids, "nt" for nucleotides.
        color_scheme: Custom color mapping.
        default_color: Fallback color for unknown characters.

    Returns:
        Hex color code string.
    """
    if color_scheme is None:
        color_scheme = get_color_scheme(sequence_type)

    return color_scheme.get(char.upper(), default_color)


__all__ = [
    "AMINO_ACID_COLORS",
    "NUCLEOTIDE_COLORS",
    "DEFAULT_COLOR",
    "get_color_scheme",
    "get_colors",
    "get_color",
]
