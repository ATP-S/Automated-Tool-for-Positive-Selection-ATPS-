"""PhyML wrapper for maximum-likelihood phylogenetic tree inference."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

# Default paths (module-relative)
_MODULE_DIR = Path(__file__).resolve().parent
_PHYML_EXE = _MODULE_DIR.parent / "assets" / "tools" / "phyml" / "src" / "phyml"


def run_phyml(
    input_file: Path | str,
    output_dir: Path | str = ".",
    partition: str = "GTR",
    freq: str = "e",
    pinvar: str = "e",
    bootstrap_replicates: int = 100,
    data_type: str = "nt",
    phyml_exe: Path | str | None = None,
    rename_outputs: bool = True,
) -> dict:
    """Run PhyML to build a maximum-likelihood phylogenetic tree.

    Args:
        input_file: Path to the input alignment file (PHYLIP format).
        output_dir: Directory for output files (default current directory).
        partition: Substitution model (e.g., "GTR", "HKY85", "JC69").
        freq: Base frequency estimation ("e" for empirical, or fixed values).
        pinvar: Proportion of invariable sites ("e" for estimated).
        bootstrap_replicates: Number of bootstrap replicates (-b).
        data_type: Data type ("nt" for nucleotides, "aa" for amino acids).
        phyml_exe: Path to PhyML executable (defaults to bundled tool).
        rename_outputs: Rename output files to friendlier names.

    Returns:
        Dict with paths to generated files: tree, boot_trees, stats, boot_stats.

    Raises:
        FileNotFoundError: If input file or PhyML executable is missing.
        RuntimeError: If PhyML returns a non-zero exit code.
    """
    input_path = Path(input_file)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input alignment not found: {input_path}")

    exe = Path(phyml_exe) if phyml_exe else _PHYML_EXE
    if not exe.exists():
        # Try system PATH
        exe = Path("phyml")

    cmd = [
        str(exe),
        "-i",
        str(input_path),
        "-d",
        data_type,
        "-b",
        str(bootstrap_replicates),
        "-m",
        partition,
        "-f",
        freq,
        "-v",
        pinvar,
    ]

    logger.info("Running PhyML: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Write logs for debugging
    (out_dir / "phyml_stdout.txt").write_text(result.stdout)
    (out_dir / "phyml_stderr.txt").write_text(result.stderr)

    if result.returncode != 0:
        logger.error("PhyML failed (rc=%s):\n%s", result.returncode, result.stderr)
        raise RuntimeError(f"PhyML failed (rc={result.returncode}): {result.stderr}")

    logger.info("PhyML completed successfully")

    # PhyML generates files with predictable names based on input
    base_name = input_path.name
    outputs = {
        "tree": input_path.parent / f"{base_name}_phyml_tree.txt",
        "boot_trees": input_path.parent / f"{base_name}_phyml_boot_trees.txt",
        "stats": input_path.parent / f"{base_name}_phyml_stats.txt",
        "boot_stats": input_path.parent / f"{base_name}_phyml_boot_stats.txt",
    }

    # Rename to friendlier names if requested
    renamed = {}
    if rename_outputs:
        rename_map = {
            "tree": out_dir / "Species_Phylogenetic_tree.txt",
            "boot_trees": out_dir / "Species_Phylogenetic_boot_trees.txt",
            "stats": out_dir / "Species_Phylogenetic_stats.txt",
            "boot_stats": out_dir / "Species_Phylogenetic_boot_stats.txt",
        }
        for key, src in outputs.items():
            dest = rename_map[key]
            if src.exists():
                shutil.move(str(src), str(dest))
                renamed[key] = dest
                logger.info("Renamed %s -> %s", src.name, dest.name)
            else:
                logger.warning("Expected output not found: %s", src)
                renamed[key] = None
        return renamed

    return outputs


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def phyml(partition: str, freq: str, pinvar: str, replica: int) -> None:
    """DEPRECATED: Use `run_phyml(input_file, ...)` instead."""
    import warnings

    warnings.warn(
        "phyml() is deprecated; use run_phyml() instead.", DeprecationWarning, stacklevel=2
    )

    run_phyml(
        input_file="Reverse_Translation_Seq.txt-gb1.phy",
        partition=partition,
        freq=freq,
        pinvar=pinvar,
        bootstrap_replicates=replica,
    )


__all__ = ["run_phyml", "phyml"]
