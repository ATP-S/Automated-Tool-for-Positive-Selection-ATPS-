"""
Test module for jModelTest conversion functions.

Tests the convert_fasta_to_phylip function defined in ATPS/utils/converters.py.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ATPS.utils.converters import (  # noqa: E402
    convert_alignment_format,
    convert_fasta_to_phylip,
)

# Path to example test data
EXAMPLE_TEST_DIR = PROJECT_ROOT / "example_tests" / "jModelTest_Conversion"
INPUT_FILE = EXAMPLE_TEST_DIR / "input" / "Reverse_Translation_Seq.txt-gb1"
EXPECTED_OUTPUT = EXAMPLE_TEST_DIR / "output" / "Reverse_Translation_Seq.txt-gb1.phy"

# jModelTest JAR path
JMODELTEST_JAR = (
    PROJECT_ROOT / "ATPS" / "assets" / "tools" / "jmodeltest-2.1.7" / "jModelTest.jar"
)


def _check_java_available() -> bool:
    """Check if Java is available."""
    return shutil.which("java") is not None


def _check_jmodeltest_available() -> bool:
    """Check if jModelTest JAR is available."""
    return JMODELTEST_JAR.exists()


class TestConvertFastaToPhylip:
    """Tests for convert_fasta_to_phylip function."""

    @pytest.fixture
    def input_file(self) -> Path:
        """Return path to input file."""
        if not INPUT_FILE.exists():
            pytest.skip(f"Input file not found: {INPUT_FILE}")
        return INPUT_FILE

    @pytest.mark.skipif(
        not _check_java_available(),
        reason="Java not installed",
    )
    @pytest.mark.skipif(
        not _check_jmodeltest_available(),
        reason="jModelTest JAR not found",
    )
    def test_convert_fasta_to_phylip(self, input_file, tmp_path):
        """Test converting FASTA to PHYLIP format."""
        test_input = tmp_path / "test_input.txt"
        test_input.write_text(input_file.read_text())

        result = convert_fasta_to_phylip(
            input_file=test_input,
            jar_path=JMODELTEST_JAR,
        )

        assert result.exists(), "Output PHYLIP file should be created"
        content = result.read_text()

        first_line = content.strip().split("\n")[0]
        parts = first_line.strip().split()

        assert len(parts) >= 2, "PHYLIP header should have at least 2 numbers"

    def test_convert_missing_input(self, tmp_path):
        """Test that missing input file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            convert_fasta_to_phylip(
                input_file=tmp_path / "nonexistent.fasta",
            )

    def test_convert_missing_jar(self, tmp_path):
        """Test that missing jModelTest JAR raises FileNotFoundError."""
        input_file = tmp_path / "input.fasta"
        input_file.write_text(">test\nATGATG\n")

        with pytest.raises(FileNotFoundError):
            convert_fasta_to_phylip(
                input_file=input_file,
                jar_path=tmp_path / "nonexistent.jar",
            )


class TestConvertAlignmentFormat:
    """Tests for convert_alignment_format function."""

    def test_convert_fasta_to_clustal(self, tmp_path):
        """Test converting FASTA to Clustal format using Biopython."""
        input_file = tmp_path / "input.fasta"
        input_file.write_text(">seq1\nATGATG---ATG\n>seq2\nATGATGATGATG\n")

        output_file = tmp_path / "output.aln"

        result = convert_alignment_format(
            input_file=input_file,
            output_file=output_file,
            input_format="fasta",
            output_format="clustal",
        )

        assert result.exists()
        content = result.read_text()
        assert "CLUSTAL" in content.upper() or "seq1" in content

    def test_convert_missing_input(self, tmp_path):
        """Test that missing input raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            convert_alignment_format(
                input_file=tmp_path / "nonexistent.fasta",
                output_file=tmp_path / "output.aln",
                input_format="fasta",
                output_format="clustal",
            )


class TestExpectedOutputFormat:
    """Tests to verify output format matches expected example."""

    def test_expected_output_is_phylip(self):
        """Verify expected output file is valid PHYLIP format."""
        if not EXPECTED_OUTPUT.exists():
            pytest.skip("Expected output file not found")

        content = EXPECTED_OUTPUT.read_text()
        lines = [line for line in content.strip().split("\n") if line.strip()]

        first_line = lines[0].strip()
        parts = first_line.split()

        assert len(parts) >= 2, "PHYLIP header should have at least 2 numbers"

        try:
            num_seqs = int(parts[0])
            seq_len = int(parts[1])
            assert num_seqs > 0
            assert seq_len > 0
        except ValueError:
            pytest.fail("PHYLIP header should contain two integers")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
