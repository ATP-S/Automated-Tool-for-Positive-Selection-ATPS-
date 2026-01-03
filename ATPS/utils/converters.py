"""File format converters for phylogenetic analysis (FASTA, PHYLIP, Newick)."""
from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Optional Biopython imports
try:
    from Bio import AlignIO, Phylo
    from Bio.Phylo import TreeConstruction
    _HAS_BIOPYTHON = True
except ImportError:
    AlignIO = None
    Phylo = None
    TreeConstruction = None
    _HAS_BIOPYTHON = False

# Default paths (module-relative)
_MODULE_DIR = Path(__file__).resolve().parent
_JMODELTEST_JAR = _MODULE_DIR.parent / "assets" / "tools" / "jmodeltest-2.1.7" / "jModelTest.jar"


def convert_fasta_to_phylip(
    input_file: Path | str,
    output_file: Optional[Path | str] = None,
    jar_path: Optional[Path | str] = None,
    java_cmd: str = "java",
) -> Path:
    """Convert a FASTA file to PHYLIP format using jModelTest.

    Args:
        input_file: Path to the input FASTA file.
        output_file: Path for the output PHYLIP file (default: input + ".phy").
        jar_path: Path to jModelTest.jar (defaults to bundled tool).
        java_cmd: Java executable command.

    Returns:
        Path to the output PHYLIP file.

    Raises:
        FileNotFoundError: If input file or jModelTest.jar is missing.
        RuntimeError: If conversion fails.
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    jar = Path(jar_path) if jar_path else _JMODELTEST_JAR
    if not jar.exists():
        raise FileNotFoundError(f"jModelTest JAR not found: {jar}")

    # jModelTest writes output as input_file.phy
    expected_output = input_path.with_suffix(input_path.suffix + ".phy")
    out_path = Path(output_file) if output_file else expected_output

    cmd = [java_cmd, "-jar", str(jar), "-d", str(input_path), "-getPhylip"]
    logger.info("Converting FASTA to PHYLIP: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error("Conversion failed:\n%s", result.stderr)
        raise RuntimeError(f"FASTA to PHYLIP conversion failed: {result.stderr}")

    # Rename if user specified different output path
    if output_file and expected_output.exists() and out_path != expected_output:
        expected_output.rename(out_path)

    logger.info("Converted to PHYLIP: %s", out_path)
    return out_path


def convert_to_newick(
    alignment_file: Path | str,
    tree_file: Path | str,
    output_newick: Path | str = "Species_Phylogenetic_tree_newick.nwk",
    output_phyloxml: Optional[Path | str] = None,
    alignment_format: str = "phylip-relaxed",
) -> Path:
    """Convert a phylogenetic tree to Newick format with parsimony optimization.

    Uses Biopython to read the alignment and tree, apply NNI parsimony search,
    and output the tree in Newick format.

    Args:
        alignment_file: Path to the alignment file (PHYLIP format).
        tree_file: Path to the input tree file (Newick format).
        output_newick: Path for the output Newick tree file.
        output_phyloxml: Path for intermediate PhyloXML file (optional).
        alignment_format: Format of the alignment file.

    Returns:
        Path to the output Newick tree file.

    Raises:
        RuntimeError: If Biopython is not installed.
        FileNotFoundError: If input files are missing.
    """
    if not _HAS_BIOPYTHON:
        raise RuntimeError(
            "Biopython is required for tree conversion. Install with: pip install biopython"
        )

    align_path = Path(alignment_file)
    tree_path = Path(tree_file)
    newick_path = Path(output_newick)
    phyloxml_path = Path(output_phyloxml) if output_phyloxml else newick_path.with_suffix(".xml")

    if not align_path.exists():
        raise FileNotFoundError(f"Alignment file not found: {align_path}")
    if not tree_path.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_path}")

    logger.info("Reading alignment from %s", align_path)
    aln = AlignIO.read(str(align_path), alignment_format)

    logger.info("Reading tree from %s", tree_path)
    tree = Phylo.read(str(tree_path), "newick")

    # Apply parsimony tree search
    logger.info("Applying NNI parsimony search...")
    scorer = TreeConstruction.ParsimonyScorer()
    searcher = TreeConstruction.NNITreeSearcher(scorer)
    constructor = TreeConstruction.ParsimonyTreeConstructor(searcher, tree)
    pars_tree = constructor.build_tree(aln)

    # Write intermediate PhyloXML
    Phylo.write(pars_tree, str(phyloxml_path), "phyloxml")
    logger.info("Wrote PhyloXML: %s", phyloxml_path)

    # Convert to Newick
    Phylo.convert(str(phyloxml_path), "phyloxml", str(newick_path), "newick")
    logger.info("Converted to Newick: %s", newick_path)

    return newick_path


def convert_alignment_format(
    input_file: Path | str,
    output_file: Path | str,
    input_format: str,
    output_format: str,
) -> Path:
    """Convert between alignment file formats using Biopython.

    Args:
        input_file: Path to the input alignment file.
        output_file: Path for the output alignment file.
        input_format: Input format (e.g., "fasta", "phylip", "clustal").
        output_format: Output format.

    Returns:
        Path to the output file.

    Raises:
        RuntimeError: If Biopython is not installed.
        FileNotFoundError: If input file is missing.
    """
    if not _HAS_BIOPYTHON:
        raise RuntimeError(
            "Biopython is required for alignment conversion. Install with: pip install biopython"
        )

    input_path = Path(input_file)
    output_path = Path(output_file)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    logger.info("Converting %s (%s) -> %s (%s)", input_path, input_format, output_path, output_format)

    alignment = AlignIO.read(str(input_path), input_format)
    AlignIO.write(alignment, str(output_path), output_format)

    logger.info("Conversion complete: %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def convert_fst_phy() -> None:
    """DEPRECATED: Use `convert_fasta_to_phylip(input_file)` instead."""
    import warnings
    warnings.warn("convert_fst_phy() is deprecated; use convert_fasta_to_phylip() instead.", DeprecationWarning, stacklevel=2)

    convert_fasta_to_phylip("Reverse_Translation_Seq.txt-gb1")


def convert_to_newickTree() -> None:
    """DEPRECATED: Use `convert_to_newick(alignment_file, tree_file)` instead."""
    import warnings
    warnings.warn("convert_to_newickTree() is deprecated; use convert_to_newick() instead.", DeprecationWarning, stacklevel=2)

    convert_to_newick(
        alignment_file="Reverse_Translation_Seq.txt-gb1.phy",
        tree_file="Species_Phylogenetic_tree.txt",
    )


__all__ = [
    "convert_fasta_to_phylip",
    "convert_to_newick",
    "convert_alignment_format",
    "convert_fst_phy",
    "convert_to_newickTree",
]
