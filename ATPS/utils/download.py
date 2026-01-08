"""Download and install external tools (jModelTest, Gblocks) for phylogenetic analysis."""

from __future__ import annotations

import logging
import tarfile
import zipfile
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlretrieve

logger = logging.getLogger(__name__)

# Default download URLs
JMODELTEST_URL = "https://github.com/ddarriba/jmodeltest2/releases/download/v2.1.9r20160115/jmodeltest-2.1.7-win32.zip"
GBLOCKS_URL = "https://drive.google.com/uc?export=download&id=1syZCAT748J8_1BOyLOtuRgHXdRK54J2k"

# Default paths
_MODULE_DIR = Path(__file__).resolve().parent
_DEFAULT_TOOLS_DIR = _MODULE_DIR.parent / "assets" / "tools"


def _download_file(url: str, dest: Path, chunk_size: int = 8192) -> Path:
    """Download a file from URL to destination.

    Args:
        url: URL to download from.
        dest: Destination file path.
        chunk_size: Download chunk size in bytes.

    Returns:
        Path to the downloaded file.

    Raises:
        RuntimeError: If download fails.
    """
    logger.info("Downloading %s -> %s", url, dest)
    try:
        urlretrieve(url, str(dest))
    except URLError as e:
        raise RuntimeError(f"Failed to download {url}: {e}") from e

    logger.info("Download complete: %s (%.2f MB)", dest.name, dest.stat().st_size / 1024 / 1024)
    return dest


def _extract_zip(archive: Path, dest_dir: Path) -> Path:
    """Extract a ZIP archive.

    Args:
        archive: Path to the ZIP file.
        dest_dir: Directory to extract into.

    Returns:
        Path to the extraction directory.
    """
    logger.info("Extracting ZIP: %s -> %s", archive, dest_dir)
    with zipfile.ZipFile(archive, "r") as zf:
        zf.extractall(dest_dir)
    return dest_dir


def _extract_tar(archive: Path, dest_dir: Path) -> Path:
    """Extract a TAR archive (including .tar.Z, .tar.gz, .tgz).

    Args:
        archive: Path to the TAR archive.
        dest_dir: Directory to extract into.

    Returns:
        Path to the extraction directory.
    """
    logger.info("Extracting TAR: %s -> %s", archive, dest_dir)
    with tarfile.open(archive, "r:*") as tf:
        tf.extractall(dest_dir)
    return dest_dir


def download_jmodeltest(
    dest_dir: Path | str | None = None,
    url: str = JMODELTEST_URL,
    force: bool = False,
) -> Path:
    """Download and extract jModelTest.

    Args:
        dest_dir: Directory to install into (default: ATPS/assets/tools).
        url: Download URL.
        force: Re-download even if already exists.

    Returns:
        Path to the extracted jModelTest directory.

    Raises:
        RuntimeError: If download or extraction fails.
    """
    tools_dir = Path(dest_dir) if dest_dir else _DEFAULT_TOOLS_DIR
    tools_dir.mkdir(parents=True, exist_ok=True)

    jmodeltest_dir = tools_dir / "jmodeltest-2.1.7"
    archive_path = tools_dir / "jmodeltest-2.1.7.zip"

    if jmodeltest_dir.exists() and not force:
        logger.info("jModelTest already installed: %s", jmodeltest_dir)
        return jmodeltest_dir

    # Download
    _download_file(url, archive_path)

    # Extract
    _extract_zip(archive_path, tools_dir)

    # Cleanup archive
    if archive_path.exists():
        archive_path.unlink()
        logger.debug("Removed archive: %s", archive_path)

    logger.info("jModelTest installed: %s", jmodeltest_dir)
    return jmodeltest_dir


def download_gblocks(
    dest_dir: Path | str | None = None,
    url: str = GBLOCKS_URL,
    force: bool = False,
) -> Path:
    """Download and extract Gblocks.

    Args:
        dest_dir: Directory to install into (default: ATPS/assets/tools).
        url: Download URL.
        force: Re-download even if already exists.

    Returns:
        Path to the extracted Gblocks directory.

    Raises:
        RuntimeError: If download or extraction fails.
    """
    tools_dir = Path(dest_dir) if dest_dir else _DEFAULT_TOOLS_DIR
    tools_dir.mkdir(parents=True, exist_ok=True)

    gblocks_dir = tools_dir / "Gblocks_0.91b"
    archive_path = tools_dir / "Gblocks_Linux64_0.91b.tar.Z"

    if gblocks_dir.exists() and not force:
        logger.info("Gblocks already installed: %s", gblocks_dir)
        return gblocks_dir

    # Download
    _download_file(url, archive_path)

    # Extract
    _extract_tar(archive_path, tools_dir)

    # Cleanup archive
    if archive_path.exists():
        archive_path.unlink()
        logger.debug("Removed archive: %s", archive_path)

    logger.info("Gblocks installed: %s", gblocks_dir)
    return gblocks_dir


def download_all_tools(
    dest_dir: Path | str | None = None,
    force: bool = False,
) -> dict:
    """Download and install all required external tools.

    Args:
        dest_dir: Directory to install tools into.
        force: Re-download even if already exists.

    Returns:
        Dict with paths to installed tools.
    """
    tools_dir = Path(dest_dir) if dest_dir else _DEFAULT_TOOLS_DIR

    results = {
        "jmodeltest": download_jmodeltest(tools_dir, force=force),
        "gblocks": download_gblocks(tools_dir, force=force),
    }

    logger.info("All tools installed to: %s", tools_dir)
    return results


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def download() -> None:
    """DEPRECATED: Use `download_all_tools()` instead."""
    import warnings

    warnings.warn(
        "download() is deprecated; use download_all_tools() instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    download_all_tools()


__all__ = [
    "download_jmodeltest",
    "download_gblocks",
    "download_all_tools",
    "download",
]
