"""Parsers for jModelTest, tree files, and codeml BEB output."""
from __future__ import annotations

import logging
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data Classes
# ---------------------------------------------------------------------------
@dataclass
class JModelTestResult:
    """Results from jModelTest parsing.

    Attributes:
        partition: Partition scheme (e.g., "012345").
        frequencies: Nucleotide frequencies as "freqA,freqC,freqG,freqT".
        pinvar: Proportion of invariable sites.
        model_name: Name of the selected model (if available).
    """

    partition: str
    frequencies: str
    pinvar: str
    model_name: Optional[str] = None

    @property
    def freq_list(self) -> List[float]:
        """Parse frequencies string into list of floats."""
        if not self.frequencies:
            return []
        try:
            return [float(f) for f in self.frequencies.split(",") if f]
        except ValueError:
            return []

    @property
    def pinvar_float(self) -> Optional[float]:
        """Parse pinvar as float."""
        try:
            return float(self.pinvar)
        except (ValueError, TypeError):
            return None


@dataclass
class BEBSite:
    """A site under positive selection from BEB analysis.

    Attributes:
        position: Position in the alignment (1-based).
        amino_acid: Amino acid at this position.
        probability: Posterior probability of positive selection.
        omega: Estimated omega (dN/dS) value.
        original_position: Position in original sequence (before Gblocks trimming).
    """

    position: int
    amino_acid: str
    probability: float
    omega: Optional[float] = None
    original_position: Optional[int] = None


# ---------------------------------------------------------------------------
# jModelTest Parsing
# ---------------------------------------------------------------------------
def parse_jmodeltest(
    input_file: Path | str = "Jmodeltest_output",
) -> JModelTestResult:
    """Parse jModelTest output to find the best substitution model.

    Args:
        input_file: Path to jModelTest output file.

    Returns:
        JModelTestResult with partition, frequencies, and pinvar.

    Raises:
        FileNotFoundError: If input file doesn't exist.
        ValueError: If required model parameters not found.
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"jModelTest output not found: {input_path}")

    content = input_path.read_text(encoding="utf-8", errors="replace")
    lines = content.splitlines()

    partition = ""
    pinvar = ""
    freq_a = freq_c = freq_g = freq_t = ""
    model_name = ""
    in_model_section = False

    for line in lines:
        # Detect model selection section
        if "Model selected:" in line:
            in_model_section = True
            # Extract model name
            match = re.search(r"Model selected:\s*(\S+)", line)
            if match:
                model_name = match.group(1)

        if "Tree " in line:
            in_model_section = False

        # Extract partition
        if "partition = " in line:
            match = re.search(r"partition\s*=\s*(\S+)", line)
            if match:
                partition = match.group(1)[:6]

        # Extract pinvar
        if "p-inv = " in line:
            match = re.search(r"p-inv\s*=\s*(\S+)", line)
            if match:
                pinvar = match.group(1)[:6]

        # Extract frequencies (only in model section)
        if in_model_section:
            if "freqA = " in line:
                match = re.search(r"freqA\s*=\s*(\S+)", line)
                if match:
                    freq_a = match.group(1)[:6]
            if "freqC = " in line:
                match = re.search(r"freqC\s*=\s*(\S+)", line)
                if match:
                    freq_c = match.group(1)[:6]
            if "freqG = " in line:
                match = re.search(r"freqG\s*=\s*(\S+)", line)
                if match:
                    freq_g = match.group(1)[:6]
            if "freqT = " in line:
                match = re.search(r"freqT\s*=\s*(\S+)", line)
                if match:
                    freq_t = match.group(1)[:6]

    # Build frequencies string
    frequencies = ",".join(f for f in [freq_a, freq_c, freq_g, freq_t] if f)

    result = JModelTestResult(
        partition=partition,
        frequencies=frequencies,
        pinvar=pinvar,
        model_name=model_name,
    )

    logger.info(
        "Parsed jModelTest: model=%s, partition=%s, pinvar=%s",
        model_name,
        partition,
        pinvar,
    )

    return result


def parse_jmodeltest_fallback(
    input_file: Path | str = "Jmodeltest_output",
) -> JModelTestResult:
    """Parse jModelTest for alternative model if best model lacks data.

    This searches through ranked models to find one with complete frequency data.

    Args:
        input_file: Path to jModelTest output file.

    Returns:
        JModelTestResult from first model with complete data.

    Raises:
        FileNotFoundError: If input file doesn't exist.
        ValueError: If no suitable model found.
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"jModelTest output not found: {input_path}")

    content = input_path.read_text(encoding="utf-8", errors="replace")

    # Find model ranking section
    start_marker = "cumWeight\n-"
    stop_marker = "-\n-lnL"

    start_idx = content.find(start_marker)
    stop_idx = content.find(stop_marker)

    if start_idx == -1 or stop_idx == -1:
        raise ValueError("Could not find model ranking section in jModelTest output")

    # Extract model names from ranking
    ranking_section = content[start_idx + 85:stop_idx - 84]
    model_names = []
    for line in ranking_section.split("\n"):
        if line.strip():
            parts = line.split()
            if parts:
                model_names.append(parts[0])

    logger.debug("Found %d ranked models", len(model_names))

    # Try each model until we find one with complete data
    for model_name in model_names[1:]:  # Skip header
        model_marker = f"   Model = {model_name}"
        model_start = content.find(model_marker)

        if model_start == -1:
            continue

        # Extract model section
        model_end = content.find("\n \n", model_start)
        if model_end == -1:
            model_end = model_start + 1000  # Fallback

        model_section = content[model_start:model_end]

        # Check if this model has frequency data
        if "freqA = " not in model_section:
            continue

        # Parse this model's parameters
        partition = ""
        pinvar = ""
        freq_a = freq_c = freq_g = freq_t = ""

        for line in model_section.split("\n"):
            if "partition =" in line:
                match = re.search(r"partition\s*=\s*(\S+)", line)
                if match:
                    partition = match.group(1)[:6]
            if "p-inv = " in line:
                match = re.search(r"p-inv\s*=\s*(\S+)", line)
                if match:
                    pinvar = match.group(1)[:6]
            if "freqA = " in line:
                match = re.search(r"freqA\s*=\s*(\S+)", line)
                if match:
                    freq_a = match.group(1)[:6]
            if "freqC = " in line:
                match = re.search(r"freqC\s*=\s*(\S+)", line)
                if match:
                    freq_c = match.group(1)[:6]
            if "freqG = " in line:
                match = re.search(r"freqG\s*=\s*(\S+)", line)
                if match:
                    freq_g = match.group(1)[:6]
            if "freqT = " in line:
                match = re.search(r"freqT\s*=\s*(\S+)", line)
                if match:
                    freq_t = match.group(1)[:6]

        frequencies = ",".join(f for f in [freq_a, freq_c, freq_g, freq_t] if f)

        logger.info("Using fallback model: %s", model_name)

        return JModelTestResult(
            partition=partition,
            frequencies=frequencies,
            pinvar=pinvar,
            model_name=model_name,
        )

    raise ValueError("No model with complete frequency data found")


# ---------------------------------------------------------------------------
# Tree File Parsing
# ---------------------------------------------------------------------------
def remove_branch_lengths(
    input_file: Path | str,
    output_file: Optional[Path | str] = None,
) -> Path:
    """Remove branch lengths from a Newick tree file.

    Codeml requires trees without branch lengths for some analyses.

    Args:
        input_file: Path to input Newick tree file.
        output_file: Path for output file. If None, adds "_nodistances" suffix.

    Returns:
        Path to the output file.

    Raises:
        FileNotFoundError: If input file doesn't exist.
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Tree file not found: {input_path}")

    if output_file is None:
        output_path = input_path.with_stem(input_path.stem + "_nodistances")
    else:
        output_path = Path(output_file)

    # Read tree content
    content = input_path.read_text(encoding="utf-8")

    # Remove branch lengths (format: :0.12345 after node names or closing parens)
    # Pattern matches colon followed by digits and optional decimal
    cleaned = re.sub(r":\d+\.?\d*(?:e[+-]?\d+)?", "", content)

    # Also remove any remaining colons that might be orphaned
    cleaned = cleaned.replace(":", "")

    output_path.write_text(cleaned, encoding="utf-8")

    logger.info("Removed branch lengths: %s -> %s", input_path, output_path)
    return output_path


# ---------------------------------------------------------------------------
# Codeml BEB Parsing
# ---------------------------------------------------------------------------
def parse_beb_results(
    codeml_output: Path | str,
) -> List[BEBSite]:
    """Parse Bayes Empirical Bayes (BEB) results from codeml output.

    Args:
        codeml_output: Path to codeml output file (mlc file).

    Returns:
        List of BEBSite objects for sites under positive selection.

    Raises:
        FileNotFoundError: If output file doesn't exist.
        ValueError: If BEB section not found.
    """
    output_path = Path(codeml_output)
    if not output_path.exists():
        raise FileNotFoundError(f"Codeml output not found: {output_path}")

    content = output_path.read_text(encoding="utf-8", errors="replace")

    # Find BEB section
    start_marker = "Bayes Empirical Bayes (BEB)"
    end_marker = "the grid"

    start_idx = content.find(start_marker)
    if start_idx == -1:
        raise ValueError("BEB section not found in codeml output")

    end_idx = content.find(end_marker, start_idx)
    if end_idx == -1:
        end_idx = len(content)

    beb_section = content[start_idx:end_idx]

    # Parse BEB lines
    sites = []
    # Split into blocks and find the data block
    blocks = beb_section.split("\n\n")

    for block in blocks:
        lines = block.strip().split("\n")
        for line in lines:
            line = line.strip()
            if not line:
                continue

            # BEB format: "  123 A  0.987*  1.234 +- 0.123"
            # Match position, amino acid, probability
            match = re.match(
                r"^\s*(\d+)\s+([A-Z*])\s+(\d+\.?\d*)\*?\*?\s*",
                line,
            )
            if match:
                position = int(match.group(1))
                amino_acid = match.group(2)
                probability = float(match.group(3))

                sites.append(
                    BEBSite(
                        position=position,
                        amino_acid=amino_acid,
                        probability=probability,
                    )
                )

    logger.info("Parsed %d BEB sites from %s", len(sites), output_path)
    return sites


def map_beb_to_original_positions(
    beb_sites: List[BEBSite],
    gblocks_html: Path | str,
    output_csv: Optional[Path | str] = None,
) -> List[BEBSite]:
    """Map BEB positions back to original sequence positions.

    Gblocks removes poorly aligned regions, so BEB positions need to be
    mapped back to the original sequence coordinates.

    Args:
        beb_sites: List of BEBSite objects from parse_beb_results.
        gblocks_html: Path to Gblocks HTML output file.
        output_csv: Optional path to write CSV with mapped positions.

    Returns:
        Updated list of BEBSite objects with original_position filled in.

    Raises:
        FileNotFoundError: If Gblocks file doesn't exist.
    """
    gblocks_path = Path(gblocks_html)
    if not gblocks_path.exists():
        raise FileNotFoundError(f"Gblocks file not found: {gblocks_path}")

    content = gblocks_path.read_text(encoding="utf-8", errors="replace")

    # Parse Gblocks flanks
    # Format: "Flanks: [123  456]  [789  1011]"
    flanks_match = re.search(r"Flanks:\s*(.+)", content)
    if not flanks_match:
        logger.warning("Could not find Flanks in Gblocks output")
        return beb_sites

    flanks_str = flanks_match.group(1)
    # Extract numbers
    positions = [int(x) for x in re.findall(r"\d+", flanks_str)]

    # Convert to codon positions (divide by 3) and create ranges
    codon_positions = [p // 3 for p in positions]
    ranges = [
        (codon_positions[i], codon_positions[i + 1])
        for i in range(0, len(codon_positions) - 1, 2)
    ]

    logger.debug("Gblocks ranges: %s", ranges)

    # Map each BEB site to original position
    for site in beb_sites:
        cumulative = 0
        remaining = site.position

        for start, end in ranges:
            range_size = end - start + 1

            if remaining <= cumulative + range_size:
                # Position falls within this range
                offset = remaining - cumulative - 1
                site.original_position = start + offset
                break

            cumulative += range_size
        else:
            logger.warning(
                "Could not map BEB position %d to original sequence",
                site.position,
            )

    # Write CSV if requested
    if output_csv:
        csv_path = Path(output_csv)
        csv_lines = ["position,amino_acid,probability,original_position"]
        for site in beb_sites:
            csv_lines.append(
                f"{site.position},{site.amino_acid},{site.probability},{site.original_position or ''}"
            )
        csv_path.write_text("\n".join(csv_lines), encoding="utf-8")
        logger.info("Wrote BEB mapping to %s", csv_path)

    return beb_sites


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def parsing_jmodeltest() -> Tuple[str, str, str]:
    """DEPRECATED: Use `parse_jmodeltest()` instead."""
    import warnings

    warnings.warn(
        "parsing_jmodeltest() is deprecated; use parse_jmodeltest() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    result = parse_jmodeltest()
    return result.partition, result.frequencies, result.pinvar


def spare_parse() -> Tuple[str, str, str]:
    """DEPRECATED: Use `parse_jmodeltest_fallback()` instead."""
    import warnings

    warnings.warn(
        "spare_parse() is deprecated; use parse_jmodeltest_fallback() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    result = parse_jmodeltest_fallback()
    return result.partition, result.frequencies, result.pinvar


def parsing_treefile() -> None:
    """DEPRECATED: Use `remove_branch_lengths()` instead."""
    import warnings

    warnings.warn(
        "parsing_treefile() is deprecated; use remove_branch_lengths() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    remove_branch_lengths(
        "Species_Phylogenetic_tree_newick.nwk",
        "Species_Phylogenetic_tree_newick_nodistances.nwk",
    )


def parse_BEB(path: str) -> str:
    """DEPRECATED: Use `parse_beb_results()` instead."""
    import warnings

    warnings.warn(
        "parse_BEB() is deprecated; use parse_beb_results() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    sites = parse_beb_results(Path(path) / "codeml078" / "codeml078_mlc.txt")

    # Return in legacy format
    lines = [f"{s.position},{s.amino_acid},{s.probability}" for s in sites]
    return "\n".join(lines)


def Positive_selection_sites(
    BEB_list: str,
    interest: str,
    path: str,
) -> List[int]:
    """DEPRECATED: Use `map_beb_to_original_positions()` instead."""
    import warnings

    warnings.warn(
        "Positive_selection_sites() is deprecated; use map_beb_to_original_positions() instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    # Parse legacy BEB_list format
    sites = []
    for line in BEB_list.strip().split("\n"):
        parts = line.split(",")
        if len(parts) >= 3:
            sites.append(
                BEBSite(
                    position=int(parts[0]),
                    amino_acid=parts[1],
                    probability=float(parts[2]),
                )
            )

    # Copy HTML to TXT (legacy behavior)
    path_obj = Path(path)
    original = path_obj / "Reverse_Translation_Seq.txt-gb1.htm"
    target = path_obj / "Reverse_Translation_Seq.txt-gb1.txt"
    if original.exists():
        shutil.copyfile(original, target)

    mapped = map_beb_to_original_positions(
        sites,
        target,
        path_obj / "BEB.csv",
    )

    return [s.original_position for s in mapped if s.original_position]


__all__ = [
    # Data classes
    "JModelTestResult",
    "BEBSite",
    # jModelTest parsing
    "parse_jmodeltest",
    "parse_jmodeltest_fallback",
    # Tree parsing
    "remove_branch_lengths",
    # BEB parsing
    "parse_beb_results",
    "map_beb_to_original_positions",
    # Legacy
    "parsing_jmodeltest",
    "spare_parse",
    "parsing_treefile",
    "parse_BEB",
    "Positive_selection_sites",
]