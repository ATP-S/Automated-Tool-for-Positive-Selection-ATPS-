"""Gblocks wrapper for trimming poorly aligned regions from sequence alignments."""
from __future__ import annotations

import logging
import os
import platform
import shutil
import subprocess
from pathlib import Path, PurePosixPath

logger = logging.getLogger(__name__)

# Default paths (module-relative)
_MODULE_DIR = Path(__file__).resolve().parent
_GBLOCKS_EXE = _MODULE_DIR.parent / "assets" / "tools" / "Gblocks_0.91b" / "Gblocks"

# Check if running on Windows
_IS_WINDOWS = platform.system() == "Windows"


def _windows_path_to_wsl(win_path: Path) -> str:
    """Convert a Windows path to WSL path format.
    
    E.g., C:\\Users\\name\\file.txt -> /mnt/c/Users/name/file.txt
    """
    # Resolve to absolute path
    abs_path = win_path.resolve()
    path_str = str(abs_path)
    
    # Handle drive letter (C: -> /mnt/c)
    if len(path_str) >= 2 and path_str[1] == ':':
        drive = path_str[0].lower()
        rest = path_str[2:].replace('\\', '/')
        return f"/mnt/{drive}{rest}"
    
    return path_str.replace('\\', '/')


def run_gblocks(
    input_file: Path | str,
    output_suffix: str = "-gb1",
    gblocks_exe: Path | str | None = None,
    skip: bool = False,
    use_wsl: bool | None = None,
) -> Path:
    """Apply Gblocks to trim poorly aligned positions from a codon alignment.

    Args:
        input_file: Path to the input alignment file (e.g., Reverse_Translation_Seq.txt).
        output_suffix: Suffix Gblocks appends to the output (default "-gb1").
        gblocks_exe: Path to Gblocks executable (defaults to bundled tool).
        skip: If True, skip Gblocks and just copy input to output path.
        use_wsl: If True, use WSL to run Gblocks on Windows. 
                 If None (default), auto-detect based on OS.

    Returns:
        Path to the Gblocks output file.

    Raises:
        FileNotFoundError: If input file or Gblocks executable is missing.
        RuntimeError: If Gblocks returns a non-zero exit code.
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input alignment not found: {input_path}")

    output_path = input_path.with_suffix(input_path.suffix + output_suffix)

    if skip:
        logger.info("Skipping Gblocks; copying %s -> %s", input_path, output_path)
        output_path.write_text(input_path.read_text())
        return output_path

    exe = Path(gblocks_exe) if gblocks_exe else _GBLOCKS_EXE
    if not exe.exists():
        raise FileNotFoundError(f"Gblocks executable not found: {exe}")

    # Auto-detect whether to use WSL
    if use_wsl is None:
        use_wsl = _IS_WINDOWS
    
    if use_wsl and _IS_WINDOWS:
        # Use WSL to run the Linux binary on Windows
        wsl_exe = _windows_path_to_wsl(exe)
        wsl_input = _windows_path_to_wsl(input_path)
        
        # Build the command string to run via bash -c
        gblocks_cmd = f'"{wsl_exe}" "{wsl_input}" -t=c -e={output_suffix} -b5=h -d=y -b2=0'
        cmd = ["wsl", "bash", "-c", gblocks_cmd]
        logger.info("Running Gblocks via WSL: %s", gblocks_cmd)
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # WSL may return exit code 1 due to PATH translation warnings, but Gblocks may still succeed
        # Check if output file was created to determine actual success
        if result.returncode != 0:
            if output_path.exists():
                # Gblocks succeeded despite WSL warning
                logger.warning("WSL warning (ignored): %s", result.stderr.strip())
            else:
                logger.error("Gblocks failed:\n%s", result.stderr)
                raise RuntimeError(f"Gblocks failed (rc={result.returncode}): {result.stderr}")
    else:
        cmd = [str(exe), str(input_path), "-t=c", "-e=" + output_suffix, "-b5=h", "-d=y", "-b2=0"]
        logger.info("Running Gblocks: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error("Gblocks failed:\n%s", result.stderr)
            raise RuntimeError(f"Gblocks failed (rc={result.returncode}): {result.stderr}")

    logger.info("Gblocks output: %s", output_path)
    return output_path


def remove_spaces(input_file: Path | str, output_file: Path | str | None = None) -> Path:
    """Remove internal spaces from a Gblocks output file.

    Args:
        input_file: Path to the Gblocks output (e.g., Reverse_Translation_Seq.txt-gb1).
        output_file: Destination path (default: input_file + ".fst").

    Returns:
        Path to the cleaned output file.
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    out_path = Path(output_file) if output_file else input_path.with_suffix(input_path.suffix + ".fst")

    content = input_path.read_text()
    cleaned = content.replace(" ", "")

    out_path.write_text(cleaned)
    logger.info("Removed spaces: %s -> %s", input_path, out_path)
    return out_path


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def Gblocks(state: str) -> None:
    """DEPRECATED: Use `run_gblocks(input_file, skip=state.lower() != 't')` instead."""
    import warnings
    warnings.warn("Gblocks() is deprecated; use run_gblocks() instead.", DeprecationWarning, stacklevel=2)

    input_file = Path("Reverse_Translation_Seq.txt")
    run_gblocks(input_file, skip=(state.lower() != "t"))


def rem_spaces() -> None:
    """DEPRECATED: Use `remove_spaces(input_file)` instead."""
    import warnings
    warnings.warn("rem_spaces() is deprecated; use remove_spaces() instead.", DeprecationWarning, stacklevel=2)

    remove_spaces("Reverse_Translation_Seq.txt-gb1")


__all__ = ["run_gblocks", "remove_spaces", "Gblocks", "rem_spaces"]
    