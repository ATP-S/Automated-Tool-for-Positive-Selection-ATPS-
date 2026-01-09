"""
Test module for fetch_sequences_from_ncbi and related sequence fetching functions.

Tests the function defined in ATPS/utils/fetchers.py.
Use case: Genes TP53 with Species Homo_sapiens, Felis_catus, Rattus_norvegicus
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ATPS.utils.fetchers import (  # noqa: E402
    fetch_and_save_sequences,
    fetch_sequences_from_ncbi,
    write_sequences_to_fasta,
)

# Test configuration
TEST_EMAIL = "marwan.js@live.com"
TEST_GENE = "TP53"
TEST_SPECIES = ["Homo_sapiens", "Felis_catus", "Rattus_norvegicus"]
INTEREST_SPECIES = "Homo_sapiens"

# Path to example test data
EXAMPLE_TEST_DIR = PROJECT_ROOT / "example_tests" / "Fetching_BioPython"
EXPECTED_CODING_SEQS = EXAMPLE_TEST_DIR / "output" / "CodingSequences.fasta"
EXPECTED_PROTEIN_SEQS = EXAMPLE_TEST_DIR / "output" / "ProteinSequences.fasta"


class TestFetchSequencesFromNCBI:
    """Tests for fetch_sequences_from_ncbi function."""

    @pytest.mark.skipif(
        not os.environ.get("ATPS_RUN_NETWORK_TESTS", False),
        reason="Network tests disabled. Set ATPS_RUN_NETWORK_TESTS=1 to enable.",
    )
    def test_fetch_sequences_basic(self):
        """Test fetching sequences from NCBI for a single species."""
        seq_dict, protein_dict, _interest = fetch_sequences_from_ncbi(
            protein=TEST_GENE,
            species_list=[INTEREST_SPECIES],
            email=TEST_EMAIL,
            interest_species=INTEREST_SPECIES,
        )

        assert len(seq_dict) >= 1, "Expected at least one nucleotide sequence"
        assert len(protein_dict) >= 1, "Expected at least one protein sequence"

        for species, seq in seq_dict.items():
            assert all(c in "ATCGatcg" for c in seq), f"Invalid nucleotides in {species}"
            assert seq.upper().startswith("ATG"), f"Sequence for {species} should start with ATG"

    @pytest.mark.skipif(
        not os.environ.get("ATPS_RUN_NETWORK_TESTS", False),
        reason="Network tests disabled. Set ATPS_RUN_NETWORK_TESTS=1 to enable.",
    )
    def test_fetch_sequences_multiple_species(self):
        """Test fetching sequences from NCBI for multiple species."""
        seq_dict, _protein_dict, interest = fetch_sequences_from_ncbi(
            protein=TEST_GENE,
            species_list=TEST_SPECIES,
            email=TEST_EMAIL,
            interest_species=INTEREST_SPECIES,
        )

        assert len(seq_dict) >= 1, "Expected sequences for at least one species"
        assert interest is not None

    def test_fetch_requires_email(self):
        """Test that email is required for NCBI queries."""
        with pytest.raises(ValueError, match="Email.*required"):
            fetch_sequences_from_ncbi(
                protein=TEST_GENE,
                species_list=TEST_SPECIES,
                email="",
            )


class TestWriteSequencesToFasta:
    """Tests for write_sequences_to_fasta function."""

    def test_write_sequences_basic(self, tmp_path):
        """Test writing sequences to a FASTA file."""
        test_sequences = {
            "species_a": "ATGATGATGATG",
            "species_b": "ATGCATGCATGC",
        }

        output_file = tmp_path / "test_output.fasta"
        result = write_sequences_to_fasta(test_sequences, output_file)

        assert result.exists(), "Output file should exist"

        content = output_file.read_text()
        assert ">species_a" in content
        assert ">species_b" in content
        assert "ATGATGATGATG" in content
        assert "ATGCATGCATGC" in content

    def test_write_sequences_normalizes_ids(self, tmp_path):
        """Test that species IDs are normalized (lowercase, underscores)."""
        test_sequences = {
            "Homo Sapiens": "ATGATG",
        }

        output_file = tmp_path / "test_output.fasta"
        write_sequences_to_fasta(test_sequences, output_file)

        content = output_file.read_text()
        assert ">homo_sapiens" in content


class TestFetchAndSaveSequences:
    """Tests for fetch_and_save_sequences function."""

    @pytest.mark.skipif(
        not os.environ.get("ATPS_RUN_NETWORK_TESTS", False),
        reason="Network tests disabled. Set ATPS_RUN_NETWORK_TESTS=1 to enable.",
    )
    def test_fetch_and_save_ncbi(self, tmp_path):
        """Test the complete fetch and save workflow."""
        success, _interest = fetch_and_save_sequences(
            protein=TEST_GENE,
            species_list=[INTEREST_SPECIES],
            interest_species=INTEREST_SPECIES,
            fetch_from_ncbi=True,
            email=TEST_EMAIL,
            output_dir=tmp_path,
        )

        assert success, "Fetching should succeed"

        coding_file = tmp_path / "CodingSequences.fasta"
        protein_file = tmp_path / "ProteinSequences.fasta"

        assert coding_file.exists(), "CodingSequences.fasta should be created"
        assert protein_file.exists(), "ProteinSequences.fasta should be created"

    def test_fetch_from_ncbi_requires_species_list(self):
        """Test that species_list is required when fetching from NCBI."""
        with pytest.raises(ValueError, match="species_list.*required"):
            fetch_and_save_sequences(
                protein=TEST_GENE,
                species_list=None,
                fetch_from_ncbi=True,
                email=TEST_EMAIL,
            )


class TestOutputFormat:
    """Tests to verify output format matches expected examples."""

    def test_expected_output_format(self):
        """Verify expected output file format from example_tests."""
        if not EXPECTED_CODING_SEQS.exists():
            pytest.skip("Expected output file not found")

        content = EXPECTED_CODING_SEQS.read_text()
        lines = content.strip().split("\n")

        headers = [line for line in lines if line.startswith(">")]
        assert len(headers) == 3, "Expected 3 species in example output"

        header_text = " ".join(headers)
        assert "homo_sapiens" in header_text.lower()
        assert "felis_catus" in header_text.lower()
        assert "rattus_norvegicus" in header_text.lower()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
