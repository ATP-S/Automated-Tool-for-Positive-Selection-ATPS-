"""
Test module for Gblocks function.

Tests the run_gblocks function defined in ATPS/models/gblocks.py.
"""

from __future__ import annotations

import platform
import sys
from pathlib import Path

import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ATPS.models.gblocks import remove_spaces, run_gblocks  # noqa: E402

# Path to example test data
EXAMPLE_TEST_DIR = PROJECT_ROOT / "example_tests" / "Gblocks"
INPUT_FILE = EXAMPLE_TEST_DIR / "input" / "Reverse_Translation_Seq.txt"
EXPECTED_OUTPUT = EXAMPLE_TEST_DIR / "output" / "Reverse_Translation_Seq.txt-gb1"

# Gblocks executable path
GBLOCKS_EXE = PROJECT_ROOT / "ATPS" / "assets" / "tools" / "Gblocks_0.91b" / "Gblocks"


def _check_gblocks_available() -> bool:
    """Check if Gblocks executable is available."""
    return GBLOCKS_EXE.exists()


class TestRunGblocks:
    """Tests for run_gblocks function."""

    @pytest.fixture
    def input_file(self) -> Path:
        """Return path to input file."""
        if not INPUT_FILE.exists():
            pytest.skip(f"Input file not found: {INPUT_FILE}")
        return INPUT_FILE

    @pytest.mark.skipif(
        not _check_gblocks_available(),
        reason="Gblocks executable not found",
    )
    @pytest.mark.skipif(
        platform.system() == "Windows",
        reason="Gblocks test not supported on Windows without WSL",
    )
    def test_gblocks_basic(self, input_file, tmp_path):
        """Test basic Gblocks execution."""
        test_input = tmp_path / "Reverse_Translation_Seq.txt"
        test_input.write_text(input_file.read_text())

        result = run_gblocks(
            input_file=test_input,
            gblocks_exe=GBLOCKS_EXE,
        )

        assert result.exists(), "Output file should be created"
        content = result.read_text()
        assert ">" in content, "Output should be FASTA format"

    @pytest.mark.skipif(
        not _check_gblocks_available(),
        reason="Gblocks executable not found",
    )
    def test_gblocks_skip(self, input_file, tmp_path):
        """Test Gblocks with skip=True."""
        test_input = tmp_path / "Reverse_Translation_Seq.txt"
        test_input.write_text(input_file.read_text())

        result = run_gblocks(
            input_file=test_input,
            skip=True,
        )

        assert result.exists()
        assert result.read_text() == test_input.read_text()

    def test_missing_input_file(self, tmp_path):
        """Test that missing input file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            run_gblocks(
                input_file=tmp_path / "nonexistent.txt",
            )

    def test_missing_gblocks_executable(self, tmp_path):
        """Test that missing Gblocks executable raises FileNotFoundError."""
        input_file = tmp_path / "input.txt"
        input_file.write_text(">test\nATGATG\n")

        with pytest.raises(FileNotFoundError):
            run_gblocks(
                input_file=input_file,
                gblocks_exe=tmp_path / "nonexistent_gblocks",
            )


class TestRemoveSpaces:
    """Tests for remove_spaces function."""

    def test_remove_spaces_basic(self, tmp_path):
        """Test removing spaces from Gblocks output."""
        input_file = tmp_path / "test_input.txt"
        content_with_spaces = ">test\nATG ATG ATG\n>test2\nGCA GCA GCA\n"
        input_file.write_text(content_with_spaces)

        result = remove_spaces(input_file)

        assert result.exists()
        content = result.read_text()
        assert " " not in content, "Spaces should be removed"
        assert "ATGATGATG" in content

    def test_remove_spaces_custom_output(self, tmp_path):
        """Test removing spaces with custom output path."""
        input_file = tmp_path / "test_input.txt"
        output_file = tmp_path / "custom_output.txt"
        input_file.write_text(">test\nATG ATG\n")

        result = remove_spaces(input_file, output_file)

        assert result == output_file
        assert output_file.exists()

    def test_remove_spaces_missing_input(self, tmp_path):
        """Test that missing input file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            remove_spaces(tmp_path / "nonexistent.txt")


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
        assert len(seq_lines) > 0, "Should have sequence lines"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
