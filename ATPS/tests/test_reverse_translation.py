"""
Test module for reverse translation function.

Tests the reverse_translate_alignment function defined in ATPS/utils/gene_operations.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ATPS.utils.gene_operations import (  # noqa: E402
    calculate_gc_content,
    get_codon_positions,
    reverse_translate_alignment,
)

# Path to example test data
EXAMPLE_TEST_DIR = PROJECT_ROOT / "example_tests" / "Reverse_Translation_Python"
INPUT_ALIGNMENT = EXAMPLE_TEST_DIR / "input" / "Alignment.ali"
INPUT_CODING_SEQS = EXAMPLE_TEST_DIR / "input" / "CodingSequences.fasta"
EXPECTED_OUTPUT = EXAMPLE_TEST_DIR / "output" / "Reverse_Translation_Seq.txt"


class TestReverseTranslateAlignment:
    """Tests for reverse_translate_alignment function."""

    @pytest.fixture
    def input_files(self) -> tuple:
        """Return paths to input files."""
        if not INPUT_ALIGNMENT.exists():
            pytest.skip(f"Input alignment file not found: {INPUT_ALIGNMENT}")
        if not INPUT_CODING_SEQS.exists():
            pytest.skip(f"Input coding sequences file not found: {INPUT_CODING_SEQS}")
        return INPUT_CODING_SEQS, INPUT_ALIGNMENT

    def test_reverse_translation_basic(self, input_files, tmp_path):
        """Test basic reverse translation functionality."""
        coding_seqs, alignment = input_files
        output_file = tmp_path / "Reverse_Translation_Seq.txt"

        result = reverse_translate_alignment(
            coding_sequences_file=coding_seqs,
            protein_alignment_file=alignment,
            output_file=output_file,
        )

        assert result.exists(), "Output file should be created"
        content = result.read_text()

        assert ">" in content, "Output should be FASTA format"

        lines = [line for line in content.split("\n") if line and not line.startswith(">")]
        for line in lines:
            non_gap_chars = line.replace("-", "").replace(" ", "")
            if non_gap_chars:
                assert all(c in "ATCGatcgN" for c in non_gap_chars), (
                    f"Unexpected character in sequence: {non_gap_chars}"
                )

    def test_reverse_translation_with_interest_species(self, input_files, tmp_path):
        """Test reverse translation with interest species first."""
        coding_seqs, alignment = input_files
        output_file = tmp_path / "Reverse_Translation_Seq.txt"

        result = reverse_translate_alignment(
            coding_sequences_file=coding_seqs,
            protein_alignment_file=alignment,
            output_file=output_file,
            interest_species="homo_sapiens",
        )

        assert result.exists()
        content = result.read_text()

        first_line = content.split("\n")[0]
        assert "homo_sapiens" in first_line.lower(), "Interest species should be first in output"

    def test_missing_coding_sequences_file(self, tmp_path):
        """Test that missing coding sequences file raises error."""
        with pytest.raises(FileNotFoundError):
            reverse_translate_alignment(
                coding_sequences_file=tmp_path / "nonexistent.fasta",
                protein_alignment_file=tmp_path / "alignment.ali",
                output_file=tmp_path / "output.txt",
            )

    def test_missing_alignment_file(self, tmp_path):
        """Test that missing alignment file raises error."""
        coding_file = tmp_path / "coding.fasta"
        coding_file.write_text(">test\nATGATG\n")

        with pytest.raises(FileNotFoundError):
            reverse_translate_alignment(
                coding_sequences_file=coding_file,
                protein_alignment_file=tmp_path / "nonexistent.ali",
                output_file=tmp_path / "output.txt",
            )


class TestGetCodonPositions:
    """Tests for get_codon_positions function."""

    def test_codon_position_1(self):
        """Test extracting first position of codons."""
        sequence = "ATGATGATG"
        result = get_codon_positions(sequence, position=1)
        assert result == "AAA", "Should get first position of each codon"

    def test_codon_position_2(self):
        """Test extracting second position of codons."""
        sequence = "ATGATGATG"
        result = get_codon_positions(sequence, position=2)
        assert result == "TTT", "Should get second position of each codon"

    def test_codon_position_3(self):
        """Test extracting third position of codons."""
        sequence = "ATGATGATG"
        result = get_codon_positions(sequence, position=3)
        assert result == "GGG", "Should get third position of each codon"

    def test_invalid_position(self):
        """Test that invalid position raises ValueError."""
        with pytest.raises(ValueError, match="Position must be 1, 2, or 3"):
            get_codon_positions("ATGATG", position=4)


class TestCalculateGCContent:
    """Tests for calculate_gc_content function."""

    def test_gc_content_50_percent(self):
        """Test GC content for 50% GC sequence."""
        sequence = "ATGC"
        result = calculate_gc_content(sequence)
        assert result == 0.5

    def test_gc_content_0_percent(self):
        """Test GC content for 0% GC sequence."""
        sequence = "ATAT"
        result = calculate_gc_content(sequence)
        assert result == 0.0

    def test_gc_content_100_percent(self):
        """Test GC content for 100% GC sequence."""
        sequence = "GCGC"
        result = calculate_gc_content(sequence)
        assert result == 1.0

    def test_gc_content_ignores_gaps(self):
        """Test that gaps are ignored in GC calculation."""
        sequence = "AT--GC"
        result = calculate_gc_content(sequence)
        assert result == 0.5

    def test_gc_content_empty_sequence(self):
        """Test GC content for empty sequence."""
        result = calculate_gc_content("")
        assert result == 0.0


class TestExpectedOutputFormat:
    """Tests to verify output format matches expected example."""

    def test_expected_output_format(self):
        """Verify expected output file format."""
        if not EXPECTED_OUTPUT.exists():
            pytest.skip("Expected output file not found")

        content = EXPECTED_OUTPUT.read_text()
        lines = content.strip().split("\n")

        headers = [line for line in lines if line.startswith(">")]
        assert len(headers) == 3, "Expected 3 species in output"

        seq_lines = [line for line in lines if line and not line.startswith(">")]
        seq_content = "".join(seq_lines)

        assert "---" in seq_content, "Should have triplet gaps for codon alignment"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
