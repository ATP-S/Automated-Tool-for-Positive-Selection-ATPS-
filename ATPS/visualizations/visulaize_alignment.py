"""Sequence alignment visualization utilities using Bokeh."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# Optional imports
try:
    import numpy as np

    _HAS_NUMPY = True
except ImportError:
    np = None
    _HAS_NUMPY = False

try:
    from bokeh.layouts import gridplot
    from bokeh.models import ColumnDataSource, Range1d, Rect, Text
    from bokeh.plotting import figure, show

    _HAS_BOKEH = True
except ImportError:
    figure = None
    show = None
    ColumnDataSource = None
    Range1d = None
    Rect = None
    Text = None
    gridplot = None
    _HAS_BOKEH = False

try:
    from Bio import AlignIO
    from Bio.Align import MultipleSeqAlignment

    _HAS_BIOPYTHON = True
except ImportError:
    AlignIO = None
    MultipleSeqAlignment = None
    _HAS_BIOPYTHON = False


def _check_dependencies() -> None:
    """Check that all required dependencies are available."""
    missing = []
    if not _HAS_NUMPY:
        missing.append("numpy")
    if not _HAS_BOKEH:
        missing.append("bokeh")
    if missing:
        raise RuntimeError(
            f"Missing required packages: {', '.join(missing)}. "
            f"Install with: pip install {' '.join(missing)}"
        )


def _get_alignment_colors(
    sequences: Sequence[str],
    color_scheme: dict[str, str] | None = None,
) -> list[str]:
    """Get colors for each residue in the alignment.

    Args:
        sequences: List of sequence strings.
        color_scheme: Optional custom color mapping. If None, uses default amino acid colors.

    Returns:
        Flat list of colors for each residue position.
    """
    # Default amino acid color scheme (similar to ClustalX)
    if color_scheme is None:
        color_scheme = {
            # Hydrophobic (blue)
            "A": "#80a0f0",
            "I": "#80a0f0",
            "L": "#80a0f0",
            "M": "#80a0f0",
            "V": "#80a0f0",
            # Aromatic (orange)
            "F": "#f0a000",
            "W": "#f0a000",
            "Y": "#f0a000",
            # Polar (green)
            "N": "#00ff00",
            "Q": "#00ff00",
            "S": "#00ff00",
            "T": "#00ff00",
            # Positive (red)
            "K": "#ff0000",
            "R": "#ff0000",
            "H": "#ff0000",
            # Negative (magenta)
            "D": "#c048c0",
            "E": "#c048c0",
            # Special (yellow/cyan)
            "C": "#f0f000",
            "G": "#f0a000",
            "P": "#00ff00",
            # Gap
            "-": "#ffffff",
            # Unknown
            "X": "#808080",
        }

    colors = []
    for seq in sequences:
        for residue in str(seq).upper():
            colors.append(color_scheme.get(residue, "#808080"))

    return colors


def load_alignment(
    alignment_file: Path | str,
    file_format: str = "fasta",
) -> MultipleSeqAlignment:
    """Load a multiple sequence alignment from file.

    Args:
        alignment_file: Path to alignment file.
        file_format: Alignment format (fasta, clustal, phylip, etc.).

    Returns:
        Biopython MultipleSeqAlignment object.

    Raises:
        RuntimeError: If Biopython is not installed.
        FileNotFoundError: If file doesn't exist.
    """
    if not _HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required. Install with: pip install biopython")

    path = Path(alignment_file)
    if not path.exists():
        raise FileNotFoundError(f"Alignment file not found: {path}")

    alignment = AlignIO.read(str(path), file_format)
    logger.info(
        "Loaded alignment with %d sequences of length %d",
        len(alignment),
        alignment.get_alignment_length(),
    )

    return alignment


def view_alignment(
    alignment: MultipleSeqAlignment | Path | str,
    fontsize: str = "9pt",
    plot_width: int = 800,
    color_scheme: dict[str, str] | None = None,
    title: str | None = None,
    show_plot: bool = True,
    file_format: str = "fasta",
) -> Any:
    """Create an interactive Bokeh visualization of a sequence alignment.

    Creates a two-panel view:
    - Top panel: Overview of entire alignment (zoomable)
    - Bottom panel: Detailed view with residue letters (scrollable)

    Args:
        alignment: Biopython alignment object, or path to alignment file.
        fontsize: Font size for residue text (e.g., "9pt", "12pt").
        plot_width: Width of the plot in pixels.
        color_scheme: Optional custom residue color mapping.
        title: Optional title for the plot.
        show_plot: Whether to display the plot immediately.
        file_format: Format if alignment is a file path.

    Returns:
        Bokeh gridplot object.

    Raises:
        RuntimeError: If required dependencies are not installed.
    """
    _check_dependencies()

    # Load alignment if path provided
    if isinstance(alignment, (str, Path)):
        if not _HAS_BIOPYTHON:
            raise RuntimeError("Biopython required to load alignment files")
        alignment = load_alignment(alignment, file_format)

    # Extract sequences and IDs
    sequences = [str(rec.seq) for rec in alignment]
    seq_ids = [rec.id for rec in alignment]

    # Validate alignment
    if not sequences:
        raise ValueError("Alignment contains no sequences")

    seq_length = len(sequences[0])
    num_seqs = len(sequences)

    logger.info("Visualizing alignment: %d sequences x %d positions", num_seqs, seq_length)

    # Build data for Bokeh
    text = [residue for seq in sequences for residue in seq]
    colors = _get_alignment_colors(sequences, color_scheme)

    # Create coordinate grids
    x = np.arange(1, seq_length + 1)
    y = np.arange(0, num_seqs, 1)
    xx, yy = np.meshgrid(x, y)
    gx = xx.ravel()
    gy = yy.flatten()
    recty = gy + 0.5  # Offset for rectangle centers

    # Create data source
    source = ColumnDataSource(
        {
            "x": gx,
            "y": gy,
            "recty": recty,
            "text": text,
            "colors": colors,
        }
    )

    # Calculate dimensions
    plot_height = num_seqs * 15 + 50
    x_range = Range1d(0, seq_length + 1, bounds="auto")
    view_len = min(100, seq_length)
    view_range = (0, view_len)

    # Panel 1: Overview (entire alignment, zoomable)
    overview = figure(
        title=title,
        width=plot_width,
        height=50,
        x_range=x_range,
        y_range=(0, num_seqs),
        tools="xpan,xwheel_zoom,reset,save",
        min_border=0,
        toolbar_location="below",
    )

    overview_rects = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_color=None,
        fill_alpha=0.6,
    )
    overview.add_glyph(source, overview_rects)
    overview.yaxis.visible = False
    overview.grid.visible = False

    # Panel 2: Detail view (with text, scrollable)
    detail = figure(
        title=None,
        width=plot_width,
        height=plot_height,
        x_range=view_range,
        y_range=seq_ids,
        tools="xpan,reset",
        min_border=0,
        toolbar_location="below",
    )

    # Add residue text
    text_glyph = Text(
        x="x",
        y="y",
        text="text",
        text_align="center",
        text_color="black",
        text_font="monospace",
        text_font_size=fontsize,
    )
    detail.add_glyph(source, text_glyph)

    # Add background rectangles
    detail_rects = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_color=None,
        fill_alpha=0.4,
    )
    detail.add_glyph(source, detail_rects)

    # Style detail panel
    detail.grid.visible = False
    detail.xaxis.major_label_text_font_style = "bold"
    detail.yaxis.minor_tick_line_width = 0
    detail.yaxis.major_tick_line_width = 0

    # Combine panels
    grid = gridplot([[overview], [detail]], toolbar_location="below")

    if show_plot:
        show(grid)

    return grid


def save_alignment_html(
    alignment: MultipleSeqAlignment | Path | str,
    output_file: Path | str,
    **kwargs,
) -> Path:
    """Save alignment visualization to an HTML file.

    Args:
        alignment: Alignment object or path to alignment file.
        output_file: Path for output HTML file.
        **kwargs: Additional arguments passed to view_alignment.

    Returns:
        Path to the saved HTML file.

    Raises:
        RuntimeError: If Bokeh is not installed.
    """
    _check_dependencies()

    from bokeh.io import save as bokeh_save
    from bokeh.resources import CDN

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create plot without showing
    plot = view_alignment(alignment, show_plot=False, **kwargs)

    # Save to file
    bokeh_save(plot, filename=str(output_path), resources=CDN, title="Alignment Viewer")

    logger.info("Saved alignment visualization to %s", output_path)
    return output_path


def get_alignment_statistics(
    alignment: MultipleSeqAlignment | Path | str,
    file_format: str = "fasta",
) -> dict[str, Any]:
    """Calculate basic statistics for an alignment.

    Args:
        alignment: Alignment object or path to alignment file.
        file_format: Format if alignment is a file path.

    Returns:
        Dictionary with alignment statistics.
    """
    if not _HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required. Install with: pip install biopython")

    # Load if path
    if isinstance(alignment, (str, Path)):
        alignment = load_alignment(alignment, file_format)

    sequences = [str(rec.seq) for rec in alignment]
    seq_length = len(sequences[0]) if sequences else 0
    num_seqs = len(sequences)

    # Calculate gap statistics
    total_residues = seq_length * num_seqs
    gap_count = sum(seq.count("-") for seq in sequences)
    gap_fraction = gap_count / total_residues if total_residues > 0 else 0

    # Calculate identity (pairwise)
    identity_sum = 0
    comparisons = 0
    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            matches = sum(1 for a, b in zip(sequences[i], sequences[j]) if a == b and a != "-")
            non_gap = sum(1 for a, b in zip(sequences[i], sequences[j]) if a != "-" and b != "-")
            if non_gap > 0:
                identity_sum += matches / non_gap
                comparisons += 1

    avg_identity = identity_sum / comparisons if comparisons > 0 else 0

    return {
        "num_sequences": num_seqs,
        "alignment_length": seq_length,
        "total_residues": total_residues,
        "gap_count": gap_count,
        "gap_fraction": round(gap_fraction, 4),
        "average_pairwise_identity": round(avg_identity, 4),
        "sequence_ids": [rec.id for rec in alignment],
    }


__all__ = [
    "view_alignment",
    "load_alignment",
    "save_alignment_html",
    "get_alignment_statistics",
]
