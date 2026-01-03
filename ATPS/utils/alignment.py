"""Multiple sequence alignment wrappers for MUSCLE, Clustal Omega, and MAFFT."""
from __future__ import annotations

import logging
import subprocess
from enum import Enum
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class Aligner(str, Enum):
    """Supported alignment methods."""
    MUSCLE = "muscle"
    CLUSTALO = "clustalo"
    MAFFT = "mafft"

    @classmethod
    def from_short(cls, code: str) -> "Aligner":
        """Convert short codes to Aligner enum.

        Args:
            code: "mu" for MUSCLE, "cl" for Clustal Omega, "mf" for MAFFT.

        Returns:
            Corresponding Aligner enum value.

        Raises:
            ValueError: If code is not recognized.
        """
        mapping = {"mu": cls.MUSCLE, "cl": cls.CLUSTALO, "mf": cls.MAFFT}
        if code.lower() not in mapping:
            raise ValueError(f"Unknown aligner code '{code}'. Use: {list(mapping.keys())}")
        return mapping[code.lower()]


def run_alignment(
    input_file: Path | str,
    output_file: Path | str = "Alignment.ali",
    aligner: Aligner | str = Aligner.MUSCLE,
    aligner_exe: Optional[str] = None,
    extra_args: Optional[list[str]] = None,
) -> Path:
    """Run multiple sequence alignment using the specified aligner.

    Args:
        input_file: Path to the input sequences file (FASTA format).
        output_file: Path for the output alignment file.
        aligner: Alignment method (Aligner enum or short code: "mu", "cl", "mf").
        aligner_exe: Custom path/command for the aligner executable.
        extra_args: Additional command-line arguments to pass.

    Returns:
        Path to the output alignment file.

    Raises:
        FileNotFoundError: If input file is missing.
        ValueError: If aligner is not recognized.
        RuntimeError: If alignment fails.
    """
    input_path = Path(input_file)
    output_path = Path(output_file)

    if not input_path.exists():
        raise FileNotFoundError(f"Input sequences not found: {input_path}")

    # Convert short code to Aligner enum if needed
    if isinstance(aligner, str):
        if aligner.lower() in ("mu", "cl", "mf"):
            aligner = Aligner.from_short(aligner)
        else:
            try:
                aligner = Aligner(aligner.lower())
            except ValueError:
                raise ValueError(f"Unknown aligner '{aligner}'. Use: {[a.value for a in Aligner]}")

    exe = aligner_exe or aligner.value
    extra = extra_args or []

    # Build command based on aligner
    if aligner == Aligner.MUSCLE:
        # MUSCLE 5 uses -align/-output, MUSCLE 3 uses -in/-out
        # Try MUSCLE 5 syntax first (more common now)
        cmd = [exe, "-align", str(input_path), "-output", str(output_path)] + extra
        use_stdout = False
    elif aligner == Aligner.CLUSTALO:
        cmd = [exe, "-i", str(input_path), "-o", str(output_path)] + extra
        use_stdout = False
    elif aligner == Aligner.MAFFT:
        cmd = [exe] + extra + [str(input_path)]
        use_stdout = True  # MAFFT outputs to stdout
    else:
        raise ValueError(f"Unsupported aligner: {aligner}")

    logger.info("Running %s alignment: %s", aligner.value, " ".join(cmd))

    try:
        if use_stdout:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                output_path.write_text(result.stdout)
        else:
            result = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        raise RuntimeError(f"{aligner.value} executable not found: '{exe}'. Please install it.")

    # Write logs for debugging
    log_dir = output_path.parent if output_path.parent.exists() else Path.cwd()
    (log_dir / f"{aligner.value}_stdout.txt").write_text(result.stdout if not use_stdout else "")
    (log_dir / f"{aligner.value}_stderr.txt").write_text(result.stderr)

    if result.returncode != 0:
        logger.error("%s failed (rc=%s):\n%s", aligner.value, result.returncode, result.stderr)
        raise RuntimeError(f"{aligner.value} failed (rc={result.returncode}): {result.stderr}")

    logger.info("%s completed; output: %s", aligner.value, output_path)
    return output_path


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def diff_aligners(file_path: str, align_type: str) -> None:
    """DEPRECATED: Use `run_alignment(input_file, aligner=...)` instead."""
    import warnings
    warnings.warn("diff_aligners() is deprecated; use run_alignment() instead.", DeprecationWarning, stacklevel=2)

    if align_type.lower() not in ("mu", "cl", "mf"):
        logger.error("Unknown alignment type '%s'. Exiting.", align_type)
        raise ValueError(f"Unknown alignment type: {align_type}")

    run_alignment(file_path, aligner=align_type)


__all__ = ["Aligner", "run_alignment", "diff_aligners"]