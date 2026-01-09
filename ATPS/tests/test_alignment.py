"""
Test module for MSA (Multiple Sequence Alignment) functions.

Tests the run_alignment function defined in ATPS/utils/alignment.py.
Tests alignment with MUSCLE, MAFFT, and ClustalO.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ATPS.utils.alignment import Aligner, run_alignment  # noqa: E402

# Path to example test data
EXAMPLE_TEST_DIR = PROJECT_ROOT / "example_tests" / "Alignment_MUSCLE_ClustalO_MAFFT"
INPUT_FILE = EXAMPLE_TEST_DIR / "input" / "ProteinSequences.fasta"
EXPECTED_OUTPUT = EXAMPLE_TEST_DIR / "output" / "Alignment.ali"


def _check_aligner_available(aligner: str) -> bool:
    """Check if an alignment tool is available."""
    if shutil.which(aligner) is not None:
        return True
    common_paths = [
        f"/home/{os.environ.get('USER', 'user')}/miniconda3/bin/{aligner}",
        f"/home/{os.environ.get('USER', 'user')}/anaconda3/bin/{aligner}",
        f"/usr/bin/{aligner}",
        f"/usr/local/bin/{aligner}",
    ]
    return any(Path(p).exists() for p in common_paths)


class TestAligner:
    """Tests for Aligner enum."""

    def test_from_short_muscle(self):
        """Test short code conversion for MUSCLE."""
        assert Aligner.from_short("mu") == Aligner.MUSCLE

    def test_from_short_clustalo(self):
        """Test short code conversion for ClustalO."""
        assert Aligner.from_short("cl") == Aligner.CLUSTALO

    def test_from_short_mafft(self):
        """Test short code conversion for MAFFT."""
        assert Aligner.from_short("mf") == Aligner.MAFFT

    def test_from_short_invalid(self):
        """Test invalid short code raises ValueError."""
        with pytest.raises(ValueError, match="Unknown aligner"):
            Aligner.from_short("invalid")


class TestRunAlignment:
    """Tests for run_alignment function."""

    @pytest.fixture
    def input_fasta(self) -> Path:
        """Return path to input FASTA file."""
        if not INPUT_FILE.exists():
            pytest.skip(f"Input file not found: {INPUT_FILE}")
        return INPUT_FILE

    @pytest.mark.skipif(
        not _check_aligner_available("muscle"),
        reason="MUSCLE not installed",
    )
    def test_muscle_alignment(self, input_fasta, tmp_path):
        """Test alignment using MUSCLE."""
        output_file = tmp_path / "alignment_muscle.ali"

        result = run_alignment(
            input_file=input_fasta,
            output_file=output_file,
            aligner=Aligner.MUSCLE,
        )

        assert result.exists(), "Output file should be created"
        content = result.read_text()

        assert content.startswith(">") or ">" in content
        assert "-" in content, "Alignment should contain gaps"

    @pytest.mark.skipif(
        not _check_aligner_available("mafft"),
        reason="MAFFT not installed",
    )
    def test_mafft_alignment(self, input_fasta, tmp_path):
        """Test alignment using MAFFT."""
        output_file = tmp_path / "alignment_mafft.ali"

        result = run_alignment(
            input_file=input_fasta,
            output_file=output_file,
            aligner=Aligner.MAFFT,
        )

        assert result.exists(), "Output file should be created"
        content = result.read_text()
        assert ">" in content, "Output should be FASTA format"

    @pytest.mark.skipif(
        not _check_aligner_available("clustalo"),
        reason="ClustalO not installed",
    )
    def test_clustalo_alignment(self, input_fasta, tmp_path):
        """Test alignment using ClustalO."""
        output_file = tmp_path / "alignment_clustalo.ali"

        result = run_alignment(
            input_file=input_fasta,
            output_file=output_file,
            aligner=Aligner.CLUSTALO,
        )

        assert result.exists(), "Output file should be created"
        content = result.read_text()
        assert ">" in content, "Output should be FASTA format"

    def test_aligner_from_string(self, input_fasta, tmp_path):
        """Test using string for aligner parameter."""
        output_file = tmp_path / "alignment_string.ali"

        if not _check_aligner_available("muscle") and not _check_aligner_available(
            "mafft"
        ):
            pytest.skip("No aligners available")

        aligner_name = "muscle" if _check_aligner_available("muscle") else "mafft"

        result = run_alignment(
            input_file=input_fasta,
            output_file=output_file,
            aligner=aligner_name,
        )

        assert result.exists()

    def test_missing_input_file(self, tmp_path):
        """Test that missing input file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            run_alignment(
                input_file=tmp_path / "nonexistent.fasta",
                output_file=tmp_path / "output.ali",
            )


class TestExpectedOutputFormat:
    """Tests to verify output format matches expected example."""

    def test_expected_output_exists(self):
        """Verify expected output file exists and has correct format."""
        if not EXPECTED_OUTPUT.exists():
            pytest.skip("Expected output file not found")

        content = EXPECTED_OUTPUT.read_text()
        lines = [line for line in content.strip().split("\n") if line]

        headers = [line for line in lines if line.startswith(">")]
        assert len(headers) == 3, "Expected 3 sequences in alignment"

        sequence_lines = [line for line in lines if not line.startswith(">")]
        full_sequence = "".join(sequence_lines)
        assert "-" in full_sequence, "Aligned sequences should contain gaps"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
