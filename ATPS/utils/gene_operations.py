"""Gene sequence operations for reverse translation and alignment processing."""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Optional Biopython import
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    _HAS_BIOPYTHON = True
except ImportError:
    SeqIO = None
    Seq = None
    SeqRecord = None
    _HAS_BIOPYTHON = False


def _check_biopython() -> None:
    """Raise error if Biopython is not available."""
    if not _HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required. Install with: pip install biopython")


def reverse_translate_alignment(
    coding_sequences_file: Path | str,
    protein_alignment_file: Path | str,
    output_file: Path | str,
    interest_species: str | None = None,
) -> Path:
    """Convert protein alignment back to codon-aligned nucleotide sequences.

    This function takes an aligned protein sequence and maps it back to the
    original coding sequences, inserting gaps (---) where the protein alignment
    has gaps (-). The result is a codon-aware nucleotide alignment.

    Args:
        coding_sequences_file: Path to FASTA with original coding sequences.
        protein_alignment_file: Path to FASTA with aligned protein sequences.
        output_file: Path for the output nucleotide alignment.
        interest_species: Species of interest to place first in output.
            If provided, this species will be the first sequence in the output.

    Returns:
        Path to the output file.

    Raises:
        RuntimeError: If Biopython is not installed.
        FileNotFoundError: If input files don't exist.
        ValueError: If sequences don't match between files.
    """
    _check_biopython()

    coding_path = Path(coding_sequences_file)
    protein_path = Path(protein_alignment_file)
    output_path = Path(output_file)

    # Validate input files exist
    if not coding_path.exists():
        raise FileNotFoundError(f"Coding sequences file not found: {coding_path}")
    if not protein_path.exists():
        raise FileNotFoundError(f"Protein alignment file not found: {protein_path}")

    # Load sequences as dictionaries
    coding_dict: dict[str, SeqRecord] = SeqIO.to_dict(SeqIO.parse(str(coding_path), "fasta"))
    protein_dict: dict[str, SeqRecord] = SeqIO.to_dict(SeqIO.parse(str(protein_path), "fasta"))

    # Get species list and optionally reorder with interest species first
    species_list = list(coding_dict.keys())

    if interest_species:
        interest_key = interest_species.lower()
        if interest_key in species_list:
            species_list.remove(interest_key)
            species_list.insert(0, interest_key)
            logger.debug("Placed %s first in output", interest_key)
        else:
            logger.warning("Interest species '%s' not found in sequences", interest_species)

    logger.info("Processing %d species", len(species_list))

    # Process each species
    records = []
    for species in species_list:
        if species not in protein_dict:
            logger.warning("Species %s not found in protein alignment, skipping", species)
            continue

        # Split coding sequence into codons (triplets)
        coding_seq = str(coding_dict[species].seq)
        codons = [coding_seq[i : i + 3] for i in range(0, len(coding_seq), 3)]

        # Get aligned protein sequence
        aligned_protein = str(protein_dict[species].seq)

        # Build reverse-translated alignment
        codon_aligned = _map_codons_to_alignment(codons, aligned_protein)

        # Create output record
        record = SeqRecord(
            Seq(codon_aligned),
            id=species,
            description="",
        )
        records.append(record)

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        SeqIO.write(records, handle, "fasta")

    logger.info("Wrote reverse-translated alignment to %s", output_path)
    return output_path


def _map_codons_to_alignment(codons: list[str], aligned_protein: str) -> str:
    """Map codons to aligned protein sequence, inserting gap codons.

    Args:
        codons: List of codon triplets from original coding sequence.
        aligned_protein: Aligned protein sequence (may contain gaps).

    Returns:
        Codon-aligned nucleotide sequence string.

    Raises:
        ValueError: If codon count doesn't match non-gap amino acid count.
    """
    result = []
    codon_index = 0

    for amino_acid in aligned_protein:
        if amino_acid == "-":
            # Gap in protein alignment -> insert gap codon
            result.append("---")
        else:
            # Amino acid -> use corresponding codon
            if codon_index >= len(codons):
                raise ValueError(
                    f"Codon index {codon_index} exceeds available codons ({len(codons)}). "
                    "Protein alignment may not match coding sequence."
                )
            result.append(codons[codon_index])
            codon_index += 1

    # Verify all codons were used
    if codon_index != len(codons):
        logger.warning(
            "Used %d of %d codons. Some codons may be unused.",
            codon_index,
            len(codons),
        )

    return "".join(result)


def get_codon_positions(
    aligned_sequence: str,
    position: int = 1,
) -> str:
    """Extract specific codon positions from an aligned nucleotide sequence.

    Useful for analyzing specific codon positions (1st, 2nd, or 3rd).

    Args:
        aligned_sequence: Codon-aligned nucleotide sequence.
        position: Codon position to extract (1, 2, or 3).

    Returns:
        String containing only the nucleotides at the specified position.

    Raises:
        ValueError: If position is not 1, 2, or 3.
    """
    if position not in (1, 2, 3):
        raise ValueError(f"Position must be 1, 2, or 3, got {position}")

    # Adjust for 0-based indexing
    pos_index = position - 1

    result = []
    for i in range(pos_index, len(aligned_sequence), 3):
        if i < len(aligned_sequence):
            result.append(aligned_sequence[i])

    return "".join(result)


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of a nucleotide sequence.

    Args:
        sequence: Nucleotide sequence (gaps are ignored).

    Returns:
        GC content as a fraction (0.0 to 1.0).
    """
    # Remove gaps
    seq = sequence.upper().replace("-", "").replace("N", "")

    if not seq:
        return 0.0

    gc_count = seq.count("G") + seq.count("C")
    return gc_count / len(seq)


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def reversedd(interest: str) -> None:
    """DEPRECATED: Use `reverse_translate_alignment()` instead.

    Args:
        interest: The species of interest.
    """
    import warnings

    warnings.warn(
        "reversedd() is deprecated; use reverse_translate_alignment() instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    reverse_translate_alignment(
        coding_sequences_file="CodingSequences.fasta",
        protein_alignment_file="Alignment.ali",
        output_file="Reverse_Translation_Seq.txt",
        interest_species=interest,
    )


__all__ = [
    "reverse_translate_alignment",
    "get_codon_positions",
    "calculate_gc_content",
    # Legacy
    "reversedd",
]
