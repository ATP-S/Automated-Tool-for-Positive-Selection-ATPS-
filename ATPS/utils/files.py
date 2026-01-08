"""File and session management utilities for ATPS pipeline."""

from __future__ import annotations

import logging
import shutil
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pipeline Session Manager (recommended approach)
# ---------------------------------------------------------------------------
class PipelineSession:
    """Context manager for pipeline operations with automatic cleanup.

    This class manages a working directory where all intermediate files
    are created. When the session ends, intermediate files are cleaned up
    and only final outputs are preserved.

    Usage:
        with PipelineSession(output_dir="./results", gene_name="ATP6") as session:
            # All tools write to session.work_dir
            run_alignment(input_file, output=session.work_dir / "aligned.fasta")
            run_phyml(session.work_dir / "aligned.phy")

            # Mark files to keep
            session.keep(session.work_dir / "final_tree.nwk")

        # Session ends: intermediate files deleted, kept files moved to output_dir

    Attributes:
        work_dir: Path to the working directory for this session.
        output_dir: Path where final outputs will be saved.
        gene_name: Name of the gene being analyzed.
    """

    # Files that are always considered intermediate (can be deleted)
    INTERMEDIATE_PATTERNS: list[str] = [
        "*.ali",
        "*.htm",
        "*.mod",
        "*_stats.txt",
        "*_boot_*.txt",
        "codeml.ctl",
        "rst",
        "rst1",
        "rub",
        "2NG.dN",
        "2NG.dS",
        "2NG.t",
        "4fold.nuc",
        "lnf",
    ]

    def __init__(
        self,
        output_dir: Path | str = ".",
        gene_name: str | None = None,
        keep_intermediates: bool = False,
        use_temp: bool = False,
    ) -> None:
        """Initialize a pipeline session.

        Args:
            output_dir: Directory for final output files.
            gene_name: Name of the gene (used for subdirectory naming).
            keep_intermediates: If True, don't delete intermediate files.
            use_temp: If True, use system temp directory for work_dir.
        """
        self.output_dir = Path(output_dir).resolve()
        self.gene_name = gene_name
        self.keep_intermediates = keep_intermediates
        self.use_temp = use_temp

        self._work_dir: Path | None = None
        self._temp_dir: tempfile.TemporaryDirectory | None = None
        self._files_to_keep: set[Path] = set()
        self._registered_files: set[Path] = set()

    @property
    def work_dir(self) -> Path:
        """Get the working directory path."""
        if self._work_dir is None:
            raise RuntimeError("Session not started. Use 'with PipelineSession() as session:'")
        return self._work_dir

    def __enter__(self) -> PipelineSession:
        """Start the session and create the working directory."""
        if self.use_temp:
            self._temp_dir = tempfile.TemporaryDirectory(prefix="atps_")
            self._work_dir = Path(self._temp_dir.name)
        else:
            suffix = f"_{self.gene_name}" if self.gene_name else ""
            self._work_dir = self.output_dir / f".atps_session{suffix}"
            self._work_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Started pipeline session in: %s", self._work_dir)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """End the session and perform cleanup."""
        if exc_type is not None:
            logger.warning("Session ended with exception: %s", exc_val)

        if not self.keep_intermediates:
            self._move_kept_files()
            self._cleanup()
        else:
            logger.info("Keeping intermediate files in: %s", self._work_dir)

        if self._temp_dir:
            self._temp_dir.cleanup()

    def register(self, path: Path | str) -> Path:
        """Register a file created during the session.

        Args:
            path: Path to the file (relative paths resolved against work_dir).

        Returns:
            Absolute path to the registered file.
        """
        file_path = Path(path)
        if not file_path.is_absolute():
            file_path = self.work_dir / file_path
        self._registered_files.add(file_path)
        return file_path

    def keep(self, path: Path | str) -> None:
        """Mark a file to be preserved after session cleanup.

        Args:
            path: Path to the file to keep.
        """
        file_path = Path(path)
        if not file_path.is_absolute():
            file_path = self.work_dir / file_path
        self._files_to_keep.add(file_path)
        logger.debug("Marked for keeping: %s", file_path)

    def get_path(self, filename: str) -> Path:
        """Get a path within the working directory.

        Args:
            filename: Name of the file.

        Returns:
            Full path to the file in work_dir.
        """
        return self.work_dir / filename

    def create_subdir(self, name: str) -> Path:
        """Create a subdirectory in the working directory.

        Args:
            name: Name of the subdirectory.

        Returns:
            Path to the created subdirectory.
        """
        subdir = self.work_dir / name
        subdir.mkdir(parents=True, exist_ok=True)
        return subdir

    def _move_kept_files(self) -> None:
        """Move files marked to keep to the output directory."""
        if not self._files_to_keep:
            return

        self.output_dir.mkdir(parents=True, exist_ok=True)

        for src in self._files_to_keep:
            if src.exists():
                dest = self.output_dir / src.name
                shutil.copy2(src, dest)
                logger.info("Saved: %s", dest)

    def _cleanup(self) -> None:
        """Remove the working directory and all contents."""
        if self._work_dir and self._work_dir.exists() and not self.use_temp:
            try:
                shutil.rmtree(self._work_dir)
                logger.info("Cleaned up session directory: %s", self._work_dir)
            except Exception as e:
                logger.warning("Failed to cleanup session directory: %s", e)


@contextmanager
def pipeline_session(
    output_dir: Path | str = ".",
    gene_name: str | None = None,
    keep_intermediates: bool = False,
) -> Iterator[PipelineSession]:
    """Functional interface for PipelineSession.

    Args:
        output_dir: Directory for final output files.
        gene_name: Name of the gene being analyzed.
        keep_intermediates: If True, don't delete intermediate files.

    Yields:
        A PipelineSession instance.
    """
    with PipelineSession(output_dir, gene_name, keep_intermediates) as session:
        yield session


# ---------------------------------------------------------------------------
# Codeml Directory Management
# ---------------------------------------------------------------------------
CODEML_MODELS = ["codeml078", "codeml8a", "codeml2", "codeml2a"]


def create_codeml_dirs(base_dir: Path | str = ".") -> list[Path]:
    """Create directories for codeml model outputs.

    Args:
        base_dir: Base directory where model dirs will be created.

    Returns:
        List of created directory paths.
    """
    base = Path(base_dir)
    created = []

    for model in CODEML_MODELS:
        model_dir = base / model
        model_dir.mkdir(parents=True, exist_ok=True)
        created.append(model_dir)
        logger.debug("Created codeml directory: %s", model_dir)

    return created


def delete_codeml_dirs(base_dir: Path | str = ".") -> None:
    """Remove codeml model directories.

    Args:
        base_dir: Base directory containing model dirs.
    """
    base = Path(base_dir)

    for model in CODEML_MODELS:
        model_dir = base / model
        if model_dir.exists():
            try:
                shutil.rmtree(model_dir)
                logger.debug("Removed codeml directory: %s", model_dir)
            except Exception as e:
                logger.warning("Failed to remove %s: %s", model_dir, e)


# ---------------------------------------------------------------------------
# File Operations
# ---------------------------------------------------------------------------
def delete_files(
    files: list[Path | str],
    base_dir: Path | str | None = None,
    ignore_missing: bool = True,
) -> int:
    """Delete a list of files.

    Args:
        files: List of file paths or names to delete.
        base_dir: Base directory for relative paths (default: current dir).
        ignore_missing: If True, don't raise errors for missing files.

    Returns:
        Number of files successfully deleted.
    """
    base = Path(base_dir) if base_dir else Path.cwd()
    deleted = 0

    for file_path in files:
        path = Path(file_path)
        if not path.is_absolute():
            path = base / path

        if path.exists():
            try:
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
                deleted += 1
                logger.debug("Deleted: %s", path)
            except Exception as e:
                logger.warning("Failed to delete %s: %s", path, e)
        elif not ignore_missing:
            logger.warning("File not found: %s", path)

    return deleted


def delete_by_patterns(
    patterns: list[str],
    base_dir: Path | str | None = None,
) -> int:
    """Delete files matching glob patterns.

    Args:
        patterns: List of glob patterns (e.g., ["*.ali", "*.mod"]).
        base_dir: Directory to search in.

    Returns:
        Number of files deleted.
    """
    base = Path(base_dir) if base_dir else Path.cwd()
    deleted = 0

    for pattern in patterns:
        for path in base.glob(pattern):
            try:
                if path.is_file():
                    path.unlink()
                    deleted += 1
                    logger.debug("Deleted: %s", path)
            except Exception as e:
                logger.warning("Failed to delete %s: %s", path, e)

    return deleted


def cleanup_intermediate_files(base_dir: Path | str | None = None) -> int:
    """Delete common intermediate files from pipeline operations.

    This replaces the old deletion_files() function with pattern-based cleanup.

    Args:
        base_dir: Directory to clean up (default: current directory).

    Returns:
        Number of files deleted.
    """
    patterns = [
        "*.ali",
        "*.htm",
        "*.mod",
        "codeml.ctl",
        "rst",
        "rst1",
        "rub",
        "lnf",
        "2NG.*",
        "4fold.nuc",
    ]

    specific_files = [
        "CodingSequences.fasta",
        "ProteinSequences.fasta",
        "Reverse_Translation_Seq.txt",
        "Reverse_Translation_Seq.txt-gb1",
        "Reverse_Translation_Seq.txt-gb1.fst",
        "Reverse_Translation_Seq.txt-gb1.phy",
        "Reverse_Translation_Seq.txt-gb1.txt",
        "Reverse_Translation_Seq.txt-gb1PS",
        "Species_Phylogenetic_tree.txt",
        "Species_Phylogenetic_tree.nwk",
        "Species_Phylogenetic_tree_newick.nwk",
        "Species_Phylogenetic_tree_newick_nodistances.nwk",
        "Species_Phylogenetic_tree_phyloxl.xml",
        "Species_Phylogenetic_stats.txt",
        "Species_Phylogenetic_boot_stats.txt",
        "Species_Phylogenetic_boot_trees.txt",
        "Jmodeltest_output",
        "interest.txt",
        "BEB.csv",
        "phast_output.xlsx",
        "phyloFit.txt",
        "wigscore",
        "wigscore.png",
        "positive_selection_sites.txt",
    ]

    count = delete_by_patterns(patterns, base_dir)
    count += delete_files(specific_files, base_dir, ignore_missing=True)

    logger.info("Cleaned up %d intermediate files", count)
    return count


def save_gene_results(
    gene_name: str,
    destination: Path | str,
    source_dir: Path | str | None = None,
    exclude_patterns: list[str] | None = None,
) -> Path:
    """Save gene analysis results to a dedicated directory.

    Args:
        gene_name: Name of the gene (used for directory naming).
        destination: Destination directory for saved results.
        source_dir: Source directory containing results (default: cwd).
        exclude_patterns: Glob patterns for files to exclude.

    Returns:
        Path to the saved gene directory.
    """
    source = Path(source_dir) if source_dir else Path.cwd()
    dest = Path(destination)
    gene_dir = dest / f"{gene_name}.gene"

    # Default exclusions (tool directories, scripts)
    default_excludes = {
        "ATPS.py",
        "ATPS_functions.py",
        "Study_Output.csv",
        "Gblocks_0.91b",
        "phyml",
        "jmodeltest-2.1.7",
        ".atps_session*",
    }

    # Collect files to copy
    exclude_set = default_excludes.copy()
    if exclude_patterns:
        for pattern in exclude_patterns:
            for path in source.glob(pattern):
                exclude_set.add(path.name)

    # Create gene directory
    gene_dir.mkdir(parents=True, exist_ok=True)

    # Copy files
    copied = 0
    for item in source.iterdir():
        if item.name in exclude_set:
            continue
        # Also skip if matches any default exclude pattern
        skip = False
        for pattern in default_excludes:
            if item.match(pattern):
                skip = True
                break
        if skip:
            continue

        try:
            dest_path = gene_dir / item.name
            if item.is_dir():
                shutil.copytree(item, dest_path, dirs_exist_ok=True)
            else:
                shutil.copy2(item, dest_path)
            copied += 1
        except Exception as e:
            logger.warning("Failed to copy %s: %s", item, e)

    logger.info("Saved %d items to %s", copied, gene_dir)
    return gene_dir


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def saving_(gene_name: str, save: str) -> None:
    """DEPRECATED: Use `save_gene_results()` instead."""
    import warnings

    warnings.warn(
        "saving_() is deprecated; use save_gene_results() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    save_gene_results(gene_name, save)


def deletion_files() -> None:
    """DEPRECATED: Use `cleanup_intermediate_files()` instead."""
    import warnings

    warnings.warn(
        "deletion_files() is deprecated; use cleanup_intermediate_files() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    cleanup_intermediate_files()


def del_codeml_dir() -> None:
    """DEPRECATED: Use `delete_codeml_dirs()` instead."""
    import warnings

    warnings.warn(
        "del_codeml_dir() is deprecated; use delete_codeml_dirs() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    delete_codeml_dirs()


def creat_codeml_dir() -> None:
    """DEPRECATED: Use `create_codeml_dirs()` instead."""
    import warnings

    warnings.warn(
        "creat_codeml_dir() is deprecated; use create_codeml_dirs() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    create_codeml_dirs()


__all__ = [
    # Session management
    "PipelineSession",
    "pipeline_session",
    # Codeml directories
    "create_codeml_dirs",
    "delete_codeml_dirs",
    "CODEML_MODELS",
    # File operations
    "delete_files",
    "delete_by_patterns",
    "cleanup_intermediate_files",
    "save_gene_results",
    # Legacy
    "saving_",
    "deletion_files",
    "del_codeml_dir",
    "creat_codeml_dir",
]
