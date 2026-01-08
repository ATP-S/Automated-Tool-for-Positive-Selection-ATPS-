"""jModelTest wrapper for nucleotide substitution model selection."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

# Default paths (module-relative)
_MODULE_DIR = Path(__file__).resolve().parent
_JMODELTEST_JAR = _MODULE_DIR.parent / "assets" / "tools" / "jmodeltest-2.1.7" / "jModelTest.jar"


def run_jmodeltest(
    input_file: Path | str,
    output_file: Path | str = "Jmodeltest_output",
    jar_path: Path | str | None = None,
    num_substitution_schemes: int = 11,
    num_rate_categories: int = 4,
    tree_method: str = "BIONJ",
    criterion: str = "AICc",
    include_invariant: bool = True,
    include_freq: bool = True,
    java_cmd: str = "java",
    extra_args: list[str] | None = None,
) -> Path:
    """Run jModelTest to select the best-fit substitution model.

    Args:
        input_file: Path to the input alignment file (PHYLIP format).
        output_file: Path/name for the output file (default "Jmodeltest_output").
        jar_path: Path to jModelTest.jar (defaults to bundled tool).
        num_substitution_schemes: Number of substitution schemes (-s, default 11).
        num_rate_categories: Number of gamma rate categories (-g, default 4).
        tree_method: Tree search method (-t, default "BIONJ").
        criterion: Model selection criterion (e.g., "AIC", "AICc", "BIC").
        include_invariant: Include proportion of invariable sites (-i).
        include_freq: Include models with unequal base frequencies (-f).
        java_cmd: Java executable (default "java").
        extra_args: Additional command-line arguments to pass.

    Returns:
        Path to the output file.

    Raises:
        FileNotFoundError: If input file or jModelTest.jar is missing.
        RuntimeError: If jModelTest returns a non-zero exit code.
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input alignment not found: {input_path}")

    jar = Path(jar_path) if jar_path else _JMODELTEST_JAR
    if not jar.exists():
        raise FileNotFoundError(f"jModelTest JAR not found: {jar}")

    output_path = Path(output_file)

    cmd = [
        java_cmd,
        "-jar",
        str(jar),
        "-d",
        str(input_path),
        "-s",
        str(num_substitution_schemes),
        "-g",
        str(num_rate_categories),
        "-t",
        tree_method,
        f"-{criterion}",
        "-o",
        str(output_path),
    ]

    if include_freq:
        cmd.append("-f")
    if include_invariant:
        cmd.append("-i")
    if extra_args:
        cmd.extend(extra_args)

    logger.info("Running jModelTest: %s", " ".join(cmd))

    # jModelTest needs to run from its own directory to find conf/jmodeltest.conf
    jmodeltest_dir = jar.parent
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(jmodeltest_dir))

    # Write logs for debugging
    log_dir = output_path.parent if output_path.parent.exists() else Path.cwd()
    (log_dir / "jmodeltest_stdout.txt").write_text(result.stdout)
    (log_dir / "jmodeltest_stderr.txt").write_text(result.stderr)

    if result.returncode != 0:
        logger.error("jModelTest failed (rc=%s):\n%s", result.returncode, result.stderr)
        raise RuntimeError(f"jModelTest failed (rc={result.returncode}): {result.stderr}")

    logger.info("jModelTest completed; output: %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def jmodel() -> None:
    """DEPRECATED: Use `run_jmodeltest(input_file)` instead."""
    import warnings

    warnings.warn(
        "jmodel() is deprecated; use run_jmodeltest() instead.", DeprecationWarning, stacklevel=2
    )

    run_jmodeltest("Reverse_Translation_Seq.txt-gb1.phy")


__all__ = ["run_jmodeltest", "jmodel"]
