"""PHAST wrapper for phylogenetic conservation analysis and wig-score computation."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

# Optional imports for visualization
try:
    import matplotlib.pyplot as plt
    import pandas as pd

    _HAS_PLOTTING = True
except ImportError:
    pd = None
    plt = None
    _HAS_PLOTTING = False


def run_phylofit(
    tree_file: Path | str,
    alignment_file: Path | str,
    output_file: Path | str = "phyloFit.mod",
    phylofit_cmd: str = "phyloFit",
) -> Path:
    """Run phyloFit to estimate a phylogenetic model from an alignment.

    Args:
        tree_file: Path to the Newick tree file.
        alignment_file: Path to the alignment file.
        output_file: Path for the output model file (default "phyloFit.mod").
        phylofit_cmd: Command/path for phyloFit executable.

    Returns:
        Path to the generated model file.

    Raises:
        FileNotFoundError: If input files are missing.
        RuntimeError: If phyloFit fails.
    """
    tree_path = Path(tree_file)
    align_path = Path(alignment_file)
    out_path = Path(output_file)

    if not tree_path.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_path}")
    if not align_path.exists():
        raise FileNotFoundError(f"Alignment file not found: {align_path}")

    cmd = [phylofit_cmd, "--tree", str(tree_path), str(align_path)]
    logger.info("Running phyloFit: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error("phyloFit failed:\n%s", result.stderr)
        raise RuntimeError(f"phyloFit failed (rc={result.returncode}): {result.stderr}")

    # phyloFit writes to phyloFit.mod by default; rename if needed
    default_output = Path("phyloFit.mod")
    if default_output.exists() and out_path != default_output:
        shutil.move(str(default_output), str(out_path))

    logger.info("phyloFit output: %s", out_path)
    return out_path


def run_phylop(
    model_file: Path | str,
    alignment_file: Path | str,
    output_file: Path | str = "wigscore.txt",
    method: str = "LRT",
    mode: str = "CONACC",
    phylop_cmd: str = "phyloP",
) -> Path:
    """Run phyloP to compute conservation/acceleration scores.

    Args:
        model_file: Path to the phyloFit model file (.mod).
        alignment_file: Path to the alignment file.
        output_file: Path for the wig-scores output (default "wigscore.txt").
        method: Statistical method (default "LRT").
        mode: Analysis mode (default "CONACC").
        phylop_cmd: Command/path for phyloP executable.

    Returns:
        Path to the wig-scores output file.

    Raises:
        FileNotFoundError: If input files are missing.
        RuntimeError: If phyloP fails.
    """
    model_path = Path(model_file)
    align_path = Path(alignment_file)
    out_path = Path(output_file)

    if not model_path.exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    if not align_path.exists():
        raise FileNotFoundError(f"Alignment file not found: {align_path}")

    cmd = [
        phylop_cmd,
        "--wig-scores",
        "--method",
        method,
        "--mode",
        mode,
        str(model_path),
        str(align_path),
    ]
    logger.info("Running phyloP: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error("phyloP failed:\n%s", result.stderr)
        raise RuntimeError(f"phyloP failed (rc={result.returncode}): {result.stderr}")

    out_path.write_text(result.stdout)
    logger.info("phyloP wig-scores written to: %s", out_path)
    return out_path


def plot_wigscores(
    wigscore_file: Path | str,
    output_image: Path | str = "wigscore.png",
    output_excel: Path | str | None = "phast_output.xlsx",
    dpi: int = 100,
) -> Path | None:
    """Parse wig-scores and generate a line plot.

    Args:
        wigscore_file: Path to the wig-scores file from phyloP.
        output_image: Path for the output PNG image.
        output_excel: Path for Excel export (None to skip).
        dpi: Image resolution.

    Returns:
        Path to the saved image, or None if plotting libraries unavailable.

    Raises:
        FileNotFoundError: If wigscore file is missing.
        RuntimeError: If pandas/matplotlib are not installed.
    """
    if not _HAS_PLOTTING:
        logger.warning("pandas/matplotlib not available; skipping plot generation.")
        return None

    wig_path = Path(wigscore_file)
    if not wig_path.exists():
        raise FileNotFoundError(f"Wig-score file not found: {wig_path}")

    img_path = Path(output_image)

    df = pd.read_table(wig_path)
    df.columns.values[0] = "Name"

    # Filter rows containing a period (numeric-like)
    df = df[df["Name"].astype(str).str.contains(r"[.]", regex=True)]
    df = df.apply(pd.to_numeric, errors="coerce").dropna()

    if output_excel:
        excel_path = Path(output_excel)
        df.to_excel(excel_path, sheet_name="WigScores", index=False)
        logger.info("Exported wig-scores to Excel: %s", excel_path)

    # Plot
    fig, ax = plt.subplots()
    df.plot.line(ax=ax, color="#0e6655", legend=False)
    ax.set_facecolor("#d0d3d4")
    ax.set_xlabel("Position")
    ax.set_ylabel("Conservation Score")
    ax.set_title("PhyloP Wig-Scores")

    fig.savefig(img_path, dpi=dpi)
    plt.close(fig)

    logger.info("Saved wig-score plot: %s", img_path)
    return img_path


def run_phast_pipeline(
    gene: str,
    tree_file: Path | str = "Species_Phylogenetic_tree.txt",
    alignment_file: Path | str = "Alignment.ali",
    output_dir: Path | str = ".",
) -> dict:
    """Run the full PHAST pipeline: phyloFit → phyloP → plot.

    Args:
        gene: Gene name (used for naming output files).
        tree_file: Path to the input tree file.
        alignment_file: Path to the alignment file.
        output_dir: Directory for output files.

    Returns:
        Dict with paths to generated files: model, wigscore, image, excel.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Copy tree to .nwk extension if needed
    tree_path = Path(tree_file)
    nwk_path = out_dir / "Species_Phylogenetic_tree.nwk"
    if tree_path.exists():
        shutil.copyfile(tree_path, nwk_path)
    else:
        raise FileNotFoundError(f"Tree file not found: {tree_path}")

    # Run phyloFit
    model_file = out_dir / f"{gene}.mod"
    run_phylofit(nwk_path, alignment_file, output_file=model_file)

    # Rename phyloFit.mod → phyloFit.txt for phyloP compatibility
    phylofit_txt = out_dir / "phyloFit.txt"
    if model_file.exists():
        shutil.copyfile(model_file, phylofit_txt)

    # Run phyloP
    wigscore_file = out_dir / "wigscore.txt"
    run_phylop(phylofit_txt, alignment_file, output_file=wigscore_file)

    # Plot
    image_file = out_dir / "wigscore.png"
    excel_file = out_dir / "phast_output.xlsx"
    plot_wigscores(wigscore_file, output_image=image_file, output_excel=excel_file)

    return {
        "model": model_file,
        "wigscore": wigscore_file,
        "image": image_file,
        "excel": excel_file,
    }


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
# def phast(gene: str) -> None:
#     """DEPRECATED: Use `run_phast_pipeline(gene)` instead."""
#     import warnings
#     warnings.warn("phast() is deprecated; use run_phast_pipeline() instead.", DeprecationWarning, stacklevel=2)

#     original = r"Species_Phylogenetic_tree.txt"
#     target = r"Species_Phylogenetic_tree.nwk"
#     shutil.copyfile(original , target)
#     try:
#         os.system("phyloFit --tree Species_Phylogenetic_tree.nwk Alignment.ali > "+gene+".mod")
#     except:
#         #os.system("sudo apt-get install -y phast")
#         os.system("phyloFit --tree Species_Phylogenetic_tree.nwk Alignment.ali > " + gene + ".mod")
#     try:

#         shutil.move('phyloFit.mod', 'phyloFit.txt')
#     except:
#         #os.system("sudo apt-get install -y phast")
#         os.system("phyloFit --tree Species_Phylogenetic_tree.nwk Alignment.ali > " + gene + ".mod")
#         shutil.move('phyloFit.mod', 'phyloFit.txt')

#     os.system("phyloP --wig-scores --method LRT --mode CONACC phyloFit.txt Alignment.ali > wigscore")
#     file1 = open("wigscore")
#     df = pd.read_table('wigscore')
#     df.to_excel('phast_output.xlsx', 'Sheet1')
#     #df['fixed'] = df['fixed'].str.replace(r'\D', '').astype(float)
#     df.columns.values[0] = "Name"
#     df=df[df.Name.str.contains(r'[.]')]
#     print(df.head())
#     df.head()
#     df=df.astype(float)
#     plt.style.use('seaborn-whitegrid')
#     df=df.astype(float)
#     p=df.plot.line(color="#0e6655")
#     fig1 = plt.gcf()
#     p.set_facecolor('#d0d3d4')
#     #plt.draw()
#     #try:
#     #    plt.show()
#     #except:
#     #    plt.show()
#     #else:
#     fig1.savefig('wigscore.png', dpi=100)


__all__ = [
    "run_phylofit",
    "run_phylop",
    "plot_wigscores",
    "run_phast_pipeline",
]
