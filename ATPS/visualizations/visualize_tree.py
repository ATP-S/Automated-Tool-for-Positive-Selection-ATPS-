"""Phylogenetic tree visualization utilities."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

logger = logging.getLogger(__name__)

# Optional imports
try:
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Clade, Tree
    _HAS_BIOPYTHON = True
except ImportError:
    Phylo = None
    Clade = None
    Tree = None
    _HAS_BIOPYTHON = False

try:
    import matplotlib.pyplot as plt
    _HAS_MATPLOTLIB = True
except ImportError:
    plt = None
    _HAS_MATPLOTLIB = False


def _check_biopython() -> None:
    """Raise error if Biopython is not available."""
    if not _HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required. Install with: pip install biopython")


def _check_matplotlib() -> None:
    """Raise error if matplotlib is not available."""
    if not _HAS_MATPLOTLIB:
        raise RuntimeError("matplotlib is required. Install with: pip install matplotlib")


def load_tree(
    tree_file: Path | str,
    tree_format: str = "newick",
) -> "Tree":
    """Load a phylogenetic tree from file.

    Args:
        tree_file: Path to the tree file.
        tree_format: Format of the tree file (newick, nexus, phyloxml, etc.).

    Returns:
        Biopython Tree object.

    Raises:
        FileNotFoundError: If tree file doesn't exist.
        RuntimeError: If Biopython is not installed.
    """
    _check_biopython()

    tree_path = Path(tree_file)
    if not tree_path.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_path}")

    tree = Phylo.read(str(tree_path), tree_format)
    logger.info("Loaded tree from %s (%s format)", tree_path, tree_format)

    return tree


def draw_ascii_tree(
    tree_file: Path | str,
    tree_format: str = "newick",
    output_file: Optional[Path | str] = None,
) -> str:
    """Draw an ASCII representation of a phylogenetic tree.

    Args:
        tree_file: Path to the tree file.
        tree_format: Format of the tree file.
        output_file: Optional file to write ASCII tree to.

    Returns:
        ASCII tree as a string.
    """
    _check_biopython()

    tree = load_tree(tree_file, tree_format)

    # Capture ASCII output
    import io
    output = io.StringIO()
    Phylo.draw_ascii(tree, file=output)
    ascii_tree = output.getvalue()

    if output_file:
        Path(output_file).write_text(ascii_tree, encoding="utf-8")
        logger.info("Wrote ASCII tree to %s", output_file)

    return ascii_tree


def visualize_tree(
    tree_file: Path | str,
    tree_format: str = "newick",
    interest_species: Optional[str] = None,
    highlight_color: str = "salmon",
    rooted: bool = True,
    output_file: Optional[Path | str] = None,
    figsize: tuple = (12, 8),
    show: bool = True,
    dpi: int = 150,
) -> Optional["Tree"]:
    """Visualize a phylogenetic tree with optional species highlighting.

    Args:
        tree_file: Path to the tree file.
        tree_format: Format of the tree file (newick, nexus, phyloxml).
        interest_species: Species name to highlight in the tree.
        highlight_color: Color for highlighting the species of interest.
        rooted: Whether to treat the tree as rooted.
        output_file: Optional path to save the figure.
        figsize: Figure size as (width, height) in inches.
        show: Whether to display the plot.
        dpi: Resolution for saved figure.

    Returns:
        The loaded Tree object, or None if visualization fails.

    Raises:
        FileNotFoundError: If tree file doesn't exist.
        RuntimeError: If required libraries are not installed.
    """
    _check_biopython()
    _check_matplotlib()

    tree = load_tree(tree_file, tree_format)
    tree.rooted = rooted

    # Highlight species of interest
    if interest_species:
        interest_key = str(interest_species).lower()
        try:
            # Find and highlight the clade containing the species
            mrca = tree.common_ancestor({"name": interest_species})
            mrca.color = highlight_color
            logger.info("Highlighted species '%s' with color '%s'", interest_species, highlight_color)
        except Exception as e:
            logger.warning("Could not highlight species '%s': %s", interest_species, e)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    try:
        Phylo.draw(tree, axes=ax, do_show=False)
        ax.set_title(f"Phylogenetic Tree")

        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
            logger.info("Saved tree visualization to %s", output_path)

        if show:
            plt.show()
        else:
            plt.close(fig)

    except Exception as e:
        logger.error("Failed to draw tree: %s", e)
        plt.close(fig)
        return None

    return tree


def highlight_clades(
    tree: "Tree",
    clade_colors: Dict[str, str],
) -> "Tree":
    """Highlight multiple clades in a tree with different colors.

    Args:
        tree: Biopython Tree object.
        clade_colors: Dictionary mapping species/clade names to colors.

    Returns:
        Modified tree with colored clades.
    """
    _check_biopython()

    for name, color in clade_colors.items():
        try:
            clade = tree.common_ancestor({"name": name})
            clade.color = color
            logger.debug("Colored clade '%s' with '%s'", name, color)
        except Exception as e:
            logger.warning("Could not find clade '%s': %s", name, e)

    return tree


def get_tree_statistics(tree: "Tree") -> Dict[str, Union[int, float, List[str]]]:
    """Get basic statistics about a phylogenetic tree.

    Args:
        tree: Biopython Tree object.

    Returns:
        Dictionary with tree statistics.
    """
    _check_biopython()

    terminals = tree.get_terminals()
    nonterminals = tree.get_nonterminals()

    stats = {
        "num_terminals": len(terminals),
        "num_internal_nodes": len(nonterminals),
        "total_branch_length": tree.total_branch_length(),
        "terminal_names": [t.name for t in terminals if t.name],
        "is_bifurcating": tree.is_bifurcating(),
    }

    # Calculate tree depth if possible
    try:
        depths = tree.depths()
        if depths:
            stats["max_depth"] = max(depths.values())
    except Exception:
        pass

    return stats


def compare_trees(
    tree1_file: Path | str,
    tree2_file: Path | str,
    tree_format: str = "newick",
) -> Dict[str, any]:
    """Compare two phylogenetic trees.

    Args:
        tree1_file: Path to first tree file.
        tree2_file: Path to second tree file.
        tree_format: Format of tree files.

    Returns:
        Dictionary with comparison results.
    """
    _check_biopython()

    tree1 = load_tree(tree1_file, tree_format)
    tree2 = load_tree(tree2_file, tree_format)

    terminals1 = {t.name for t in tree1.get_terminals() if t.name}
    terminals2 = {t.name for t in tree2.get_terminals() if t.name}

    comparison = {
        "tree1_terminals": len(terminals1),
        "tree2_terminals": len(terminals2),
        "shared_terminals": len(terminals1 & terminals2),
        "only_in_tree1": list(terminals1 - terminals2),
        "only_in_tree2": list(terminals2 - terminals1),
        "same_taxa": terminals1 == terminals2,
    }

    return comparison


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def visualization_tree(interest: str) -> None:
    """DEPRECATED: Use `visualize_tree()` instead.

    Args:
        interest: Species of interest to highlight.
    """
    import warnings

    warnings.warn(
        "visualization_tree() is deprecated; use visualize_tree() instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    visualize_tree(
        tree_file="Species_Phylogenetic_tree.txt",
        tree_format="newick",
        interest_species=interest,
        highlight_color="salmon",
        rooted=True,
        show=True,
    )


__all__ = [
    "load_tree",
    "draw_ascii_tree",
    "visualize_tree",
    "highlight_clades",
    "get_tree_statistics",
    "compare_trees",
    # Legacy
    "visualization_tree",
]   