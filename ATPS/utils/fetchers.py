"""Sequence fetching utilities for retrieving gene/protein data from NCBI or local files."""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Optional imports
try:
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    _HAS_BIOPYTHON = True
except ImportError:
    Entrez = None
    SeqIO = None
    Seq = None
    SeqRecord = None
    _HAS_BIOPYTHON = False

try:
    import pandas as pd

    _HAS_PANDAS = True
except ImportError:
    pd = None
    _HAS_PANDAS = False


def _get_longest_valid_sequence(sequences: list[str]) -> str | None:
    """Get the longest sequence that is a valid coding sequence.

    A valid CDS:
    - Starts with ATG
    - Ends with a stop codon (translated ends with *)
    - Has no internal stop codons

    Args:
        sequences: List of nucleotide sequences.

    Returns:
        The longest valid sequence, or None if none found.
    """
    if not _HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required. Install with: pip install biopython")

    # Sort by length descending
    sorted_seqs = sorted(sequences, key=len, reverse=True)

    for seq_str in sorted_seqs:
        if not seq_str:
            continue
        try:
            translated = str(Seq(seq_str).translate())
            # Check: starts with ATG, ends with stop, no internal stops
            if (
                seq_str.upper().startswith("ATG")
                and translated.endswith("*")
                and "*" not in translated[:-1]
            ):
                return seq_str
        except Exception as e:
            logger.debug("Translation failed for sequence: %s", e)
            continue

    return None


def fetch_sequences_from_ncbi(
    protein: str,
    species_list: list[str],
    email: str,
    interest_species: str | None = None,
) -> tuple[dict[str, str], dict[str, str], str | None]:
    """Fetch gene sequences from NCBI nuccore database.

    Args:
        protein: Protein/gene name to search for.
        species_list: List of species names to search.
        email: Email address for NCBI Entrez (required).
        interest_species: Species of interest for downstream analysis.

    Returns:
        Tuple of (nucleotide_dict, protein_dict, interest_species).
        - nucleotide_dict: {species: nucleotide_sequence}
        - protein_dict: {species: protein_sequence}
        - interest_species: Updated species of interest.

    Raises:
        RuntimeError: If Biopython is not installed.
        ValueError: If email is not provided.
    """
    if not _HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required. Install with: pip install biopython")

    if not email:
        raise ValueError("Email address is required for NCBI Entrez queries.")

    Entrez.email = email
    seq_dict: dict[str, str] = {}
    protein_dict: dict[str, str] = {}
    filter_gene = f"[gene={protein}]"
    interest = interest_species

    for species in species_list:
        logger.info("Fetching sequences for species: %s", species)

        # Build search term
        search_term = (
            f"{species}[Organism] OR {species}[All Fields] AND "
            f"{protein}[Title] AND (biomol_mrna[PROP] AND refseq[filter])"
        )

        try:
            # Search NCBI
            handle = Entrez.esearch(
                db="nuccore",
                term=search_term,
                retmax=4000,
                usehistory="y",
                idtype="acc",
            )
            search_results = Entrez.read(handle)
            handle.close()

            count = int(search_results["Count"])
            if count == 0:
                logger.warning("No sequences found for species: %s", species)
                continue

            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            # Fetch sequences
            fetch_handle = Entrez.efetch(
                db="nuccore",
                rettype="fasta_cds_na",
                retmode="text",
                retmax=count,
                webenv=webenv,
                query_key=query_key,
                idtype="acc",
            )
            data = fetch_handle.read()
            fetch_handle.close()

        except Exception as e:
            logger.error("NCBI fetch failed for %s: %s", species, e)
            continue

        # Parse fetched sequences
        candidate_sequences = []
        for gene in str(data).split("\n\n")[:-1]:
            if filter_gene.lower() in gene.lower():
                # Extract sequence (after the header line)
                seq_start = gene.find("]\n")
                if seq_start != -1:
                    sequence = gene[seq_start + 2 :].replace("\n", "").strip()
                    if sequence:
                        candidate_sequences.append(sequence)

        if not candidate_sequences:
            logger.warning("No matching sequences for %s in species %s", protein, species)
            continue

        # Get the longest valid coding sequence
        valid_seq = _get_longest_valid_sequence(candidate_sequences)
        if valid_seq:
            seq_dict[species] = valid_seq
            logger.info("Found valid sequence for %s (%d bp)", species, len(valid_seq))
        else:
            logger.warning("No valid CDS found for species: %s", species)
            if species == interest:
                # Update interest species if current one not found
                remaining = [s for s in species_list if s in seq_dict]
                if remaining:
                    interest = remaining[0]
                    logger.info("Updated interest species to: %s", interest)

    # Translate to proteins
    for species, nuc_seq in seq_dict.items():
        try:
            protein_dict[species] = str(Seq(nuc_seq).translate())
        except Exception as e:
            logger.error("Translation failed for %s: %s", species, e)

    return seq_dict, protein_dict, interest


def fetch_sequences_from_file(
    gene_path: Path | str,
    protein: str,
) -> dict[str, str]:
    """Load sequences from a local FASTA file.

    Args:
        gene_path: Directory containing FASTA files.
        protein: Protein/gene name (filename without extension).

    Returns:
        Dictionary mapping sequence IDs to sequences.

    Raises:
        FileNotFoundError: If the FASTA file doesn't exist.
        RuntimeError: If Biopython is not installed.
    """
    if not _HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required. Install with: pip install biopython")

    fasta_path = Path(gene_path) / f"{protein}.fasta"
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    seq_dict = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        seq_dict[str(record.id)] = str(record.seq)

    logger.info("Loaded %d sequences from %s", len(seq_dict), fasta_path)
    return seq_dict


def write_sequences_to_fasta(
    sequences: dict[str, str],
    output_file: Path | str,
    description: str = "",
) -> Path:
    """Write sequences to a FASTA file.

    Args:
        sequences: Dictionary mapping IDs to sequences.
        output_file: Output file path.
        description: Optional description for each record.

    Returns:
        Path to the output file.
    """
    if not _HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required. Install with: pip install biopython")

    output_path = Path(output_file)

    records = []
    for seq_id, sequence in sequences.items():
        # Normalize ID: lowercase, underscores for spaces
        normalized_id = seq_id.replace(" ", "_").lower()
        record = SeqRecord(Seq(sequence), id=normalized_id, description=description)
        records.append(record)

    with output_path.open("w") as handle:
        SeqIO.write(records, handle, "fasta")

    logger.info("Wrote %d sequences to %s", len(records), output_path)
    return output_path


def fetch_and_save_sequences(
    protein: str,
    species_list: list[str] | None = None,
    interest_species: str | None = None,
    fetch_from_ncbi: bool = True,
    gene_path: Path | str | None = None,
    email: str | None = None,
    output_dir: Path | str = ".",
) -> tuple[bool, str | None]:
    """Fetch sequences and save to FASTA files.

    This is the main entry point that combines fetching and saving.

    Args:
        protein: Protein/gene name.
        species_list: List of species (required if fetch_from_ncbi=True).
        interest_species: Species of interest for downstream analysis.
        fetch_from_ncbi: If True, fetch from NCBI; if False, load from file.
        gene_path: Path to local FASTA files (required if fetch_from_ncbi=False).
        email: Email for NCBI Entrez (required if fetch_from_ncbi=True).
        output_dir: Directory for output files.

    Returns:
        Tuple of (success, interest_species).
        - success: True if sequences were found and saved.
        - interest_species: The species of interest (may be updated).
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if fetch_from_ncbi:
        if not species_list:
            raise ValueError("species_list is required when fetching from NCBI.")
        if not email:
            raise ValueError("email is required when fetching from NCBI.")

        seq_dict, protein_dict, interest = fetch_sequences_from_ncbi(
            protein=protein,
            species_list=species_list,
            email=email,
            interest_species=interest_species,
        )
    else:
        if not gene_path:
            raise ValueError("gene_path is required when loading from file.")

        seq_dict = fetch_sequences_from_file(gene_path, protein)
        interest = list(seq_dict.keys())[0] if seq_dict else None

        # Translate to proteins
        protein_dict = {}
        for species, nuc_seq in seq_dict.items():
            try:
                protein_dict[species] = str(Seq(nuc_seq).translate())
            except Exception as e:
                logger.error("Translation failed for %s: %s", species, e)

    # Write output files
    if seq_dict:
        write_sequences_to_fasta(seq_dict, out_dir / "CodingSequences.fasta")
    if protein_dict:
        write_sequences_to_fasta(protein_dict, out_dir / "ProteinSequences.fasta")

    success = bool(protein_dict)
    return success, interest


def count_fetched_species(
    input_dir: Path | str = ".",
    output_file: Path | str = "fetched_species.csv",
) -> Path:
    """Count the number of species in each gene's CodingSequences.fasta file.

    Args:
        input_dir: Directory containing .gene subdirectories.
        output_file: Output CSV file path.

    Returns:
        Path to the output CSV file.

    Raises:
        RuntimeError: If pandas is not installed.
    """
    if not _HAS_PANDAS:
        raise RuntimeError("pandas is required. Install with: pip install pandas")

    input_path = Path(input_dir)
    output_path = Path(output_file)

    results = [["Name", "Number of Species"]]

    for gene_dir in input_path.glob("*.gene"):
        fasta_file = gene_dir / "CodingSequences.fasta"
        if fasta_file.exists():
            content = fasta_file.read_text()
            count = content.count(">")
            results.append([gene_dir.name, count])
            logger.debug("Gene %s: %d species", gene_dir.name, count)

    df = pd.DataFrame(results[1:], columns=results[0])
    df.to_csv(output_path, index=False)

    logger.info("Wrote species counts to %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Legacy API (deprecated) — kept for backward compatibility
# ---------------------------------------------------------------------------
def fetchingbyspecies(
    protein: str,
    List_species: list[str] | None = None,
    interest: str | None = None,
    fetch: int | None = None,
    gene_path: str | None = None,
    inp_path: str | None = None,
) -> tuple[bool, str | None]:
    """DEPRECATED: Use `fetch_and_save_sequences()` instead."""
    import warnings

    warnings.warn(
        "fetchingbyspecies() is deprecated; use fetch_and_save_sequences() instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return fetch_and_save_sequences(
        protein=protein,
        species_list=List_species,
        interest_species=interest,
        fetch_from_ncbi=(fetch == 1),
        gene_path=gene_path,
        email="example@gmail.com",  # Legacy default
    )


def _extracted_from_fetchingbyspecies_95(gene_path: str, protein: str) -> dict[str, str]:
    """DEPRECATED: Use `fetch_sequences_from_file()` instead."""
    import warnings

    warnings.warn(
        "_extracted_from_fetchingbyspecies_95() is deprecated; use fetch_sequences_from_file() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return fetch_sequences_from_file(gene_path, protein)


def number_of_fetched_species() -> None:
    """DEPRECATED: Use `count_fetched_species()` instead."""
    import warnings

    warnings.warn(
        "number_of_fetched_species() is deprecated; use count_fetched_species() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    count_fetched_species()


def get_max_str_index(lst: list[str]) -> tuple[int, str]:
    """Get the index and value of the longest string in a list.

    Args:
        lst: List of strings.

    Returns:
        Tuple of (index, string) for the longest string.
    """
    return max(enumerate(lst), key=lambda x: len(x[1]))


__all__ = [
    "fetch_sequences_from_ncbi",
    "fetch_sequences_from_file",
    "write_sequences_to_fasta",
    "fetch_and_save_sequences",
    "count_fetched_species",
    "get_max_str_index",
    # Legacy
    "fetchingbyspecies",
    "number_of_fetched_species",
]
