"""
Test module for jModelTest model selection function.

Tests the run_jmodeltest function defined in ATPS/models/jmodeltest.py.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ATPS.models.jmodeltest import run_jmodeltest  # noqa: E402

# Path to example test data
EXAMPLE_TEST_DIR = PROJECT_ROOT / "example_tests" / "jModelTest_Model_Selection"
INPUT_FILE = EXAMPLE_TEST_DIR / "input" / "Reverse_Translation_Seq.txt-gb1.phy"
EXPECTED_OUTPUT = EXAMPLE_TEST_DIR / "output" / "Jmodeltest_output"

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


class TestRunJmodeltest:
    """Tests for run_jmodeltest function."""

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
    def test_jmodeltest_basic(self, input_file, tmp_path):
        """Test basic jModelTest execution."""
        output_file = tmp_path / "Jmodeltest_output"

        result = run_jmodeltest(
            input_file=input_file,
            output_file=output_file,
            jar_path=JMODELTEST_JAR,
        )

        assert result.exists(), "Output file should be created"
        content = result.read_text()
        assert len(content) > 0, "Output should not be empty"

    @pytest.mark.skipif(
        not _check_java_available(),
        reason="Java not installed",
    )
    @pytest.mark.skipif(
        not _check_jmodeltest_available(),
        reason="jModelTest JAR not found",
    )
    def test_jmodeltest_with_options(self, input_file, tmp_path):
        """Test jModelTest with custom options."""
        output_file = tmp_path / "Jmodeltest_output"

        result = run_jmodeltest(
            input_file=input_file,
            output_file=output_file,
            jar_path=JMODELTEST_JAR,
            num_substitution_schemes=7,
            num_rate_categories=4,
            criterion="AIC",
        )

        assert result.exists()

    def test_missing_input_file(self, tmp_path):
        """Test that missing input file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            run_jmodeltest(
                input_file=tmp_path / "nonexistent.phy",
            )

    def test_missing_jar_file(self, tmp_path):
        """Test that missing jModelTest JAR raises FileNotFoundError."""
        input_file = tmp_path / "input.phy"
        input_file.write_text(
            "3 12\nseq1 ATGATGATGAT\nseq2 ATGATGATGAT\nseq3 ATGATGATGAT\n"
        )

        with pytest.raises(FileNotFoundError):
            run_jmodeltest(
                input_file=input_file,
                jar_path=tmp_path / "nonexistent.jar",
            )


class TestExpectedOutputFormat:
    """Tests to verify output format matches expected example."""

    def test_expected_output_format(self):
        """Verify expected output file format."""
        if not EXPECTED_OUTPUT.exists():
            pytest.skip("Expected output file not found")

        content = EXPECTED_OUTPUT.read_text()

        assert len(content) > 0, "Output should not be empty"

        model_keywords = ["model", "lnL", "AIC", "BIC", "AICc"]
        content_lower = content.lower()

        has_model_info = any(kw.lower() in content_lower for kw in model_keywords)
        assert has_model_info, "Output should contain model selection information"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-m", "not slow"])
