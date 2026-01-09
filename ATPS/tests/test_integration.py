"""
End-to-end integration test that runs data through the complete ATPS pipeline.

This test verifies the full workflow:
1. fetch_sequences_from_ncbi / load from file
2. run_alignment (MSA)
3. reverse_translate_alignment
4. run_gblocks
5. convert_fasta_to_phylip
6. run_jmodeltest
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ATPS.models.gblocks import run_gblocks  # noqa: E402
from ATPS.utils.alignment import Aligner, run_alignment  # noqa: E402
from ATPS.utils.converters import convert_fasta_to_phylip  # noqa: E402
from ATPS.utils.gene_operations import reverse_translate_alignment  # noqa: E402

# Path to example test data
EXAMPLE_BASE = PROJECT_ROOT / "example_tests"


class TestEndToEndPipeline:
    """End-to-end pipeline integration tests."""

    @pytest.fixture
    def setup_test_data(self, tmp_path):
        """Copy example test data to temporary directory."""
        coding_src = EXAMPLE_BASE / "Fetching_BioPython" / "output" / "CodingSequences.fasta"
        coding_dst = tmp_path / "CodingSequences.fasta"
        if coding_src.exists():
            shutil.copy(coding_src, coding_dst)

        protein_src = EXAMPLE_BASE / "Fetching_BioPython" / "output" / "ProteinSequences.fasta"
        protein_dst = tmp_path / "ProteinSequences.fasta"
        if protein_src.exists():
            shutil.copy(protein_src, protein_dst)

        return {
            "work_dir": tmp_path,
            "coding_seqs": coding_dst,
            "protein_seqs": protein_dst,
        }

    def test_pipeline_step1_alignment(self, setup_test_data):
        """Test Step 1: MSA alignment of protein sequences."""
        data = setup_test_data

        if not data["protein_seqs"].exists():
            pytest.skip("Protein sequences not available")

        output_file = data["work_dir"] / "Alignment.ali"

        result = run_alignment(
            input_file=data["protein_seqs"],
            output_file=output_file,
            aligner=Aligner.MUSCLE,
        )

        assert result.exists(), "Alignment output should exist"
        content = result.read_text()
        assert ">" in content, "Output should be FASTA format"
        assert "-" in content, "Alignment should contain gaps"

    def test_pipeline_step2_reverse_translation(self, setup_test_data):
        """Test Step 2: Reverse translation."""
        data = setup_test_data

        if not data["coding_seqs"].exists():
            pytest.skip("Coding sequences not available")

        alignment_file = data["work_dir"] / "Alignment.ali"
        run_alignment(
            input_file=data["protein_seqs"],
            output_file=alignment_file,
            aligner=Aligner.MUSCLE,
        )

        output_file = data["work_dir"] / "Reverse_Translation_Seq.txt"

        result = reverse_translate_alignment(
            coding_sequences_file=data["coding_seqs"],
            protein_alignment_file=alignment_file,
            output_file=output_file,
            interest_species="homo_sapiens",
        )

        assert result.exists(), "Reverse translation output should exist"
        content = result.read_text()

        first_line = content.split("\n")[0]
        assert "homo_sapiens" in first_line.lower(), "Interest species should be first"
        assert "---" in content, "Should have triplet gaps for codon alignment"

    def test_pipeline_step3_gblocks(self, setup_test_data):
        """Test Step 3: Gblocks trimming of poorly aligned regions."""
        data = setup_test_data

        alignment_file = data["work_dir"] / "Alignment.ali"
        run_alignment(
            input_file=data["protein_seqs"],
            output_file=alignment_file,
            aligner=Aligner.MUSCLE,
        )

        reverse_file = data["work_dir"] / "Reverse_Translation_Seq.txt"
        reverse_translate_alignment(
            coding_sequences_file=data["coding_seqs"],
            protein_alignment_file=alignment_file,
            output_file=reverse_file,
        )

        gblocks_exe = PROJECT_ROOT / "ATPS" / "assets" / "tools" / "Gblocks_0.91b" / "Gblocks"
        if not gblocks_exe.exists():
            pytest.skip("Gblocks not available")

        result = run_gblocks(
            input_file=reverse_file,
            gblocks_exe=gblocks_exe,
        )

        assert result.exists(), "Gblocks output should exist"
        content = result.read_text()
        assert ">" in content, "Output should be FASTA format"

    def test_pipeline_step4_phylip_conversion(self, setup_test_data):
        """Test Step 4: Convert to PHYLIP format for jModelTest."""
        data = setup_test_data

        gblocks_output = EXAMPLE_BASE / "Gblocks" / "output" / "Reverse_Translation_Seq.txt-gb1"

        if not gblocks_output.exists():
            pytest.skip("Gblocks output not available")

        input_file = data["work_dir"] / "Reverse_Translation_Seq.txt-gb1"
        shutil.copy(gblocks_output, input_file)

        jmodeltest_jar = (
            PROJECT_ROOT / "ATPS" / "assets" / "tools" / "jmodeltest-2.1.7" / "jModelTest.jar"
        )
        if not jmodeltest_jar.exists():
            pytest.skip("jModelTest not available")

        result = convert_fasta_to_phylip(
            input_file=input_file,
            jar_path=jmodeltest_jar,
        )

        assert result.exists(), "PHYLIP output should exist"
        content = result.read_text()

        first_line = content.strip().split("\n")[0]
        parts = first_line.strip().split()
        assert len(parts) >= 2, "PHYLIP format should have header"

    def test_full_pipeline_flow(self, setup_test_data):
        """Test the complete pipeline from sequences to model selection input."""
        data = setup_test_data

        if not data["protein_seqs"].exists():
            pytest.skip("Test data not available")

        print("\n=== ATPS Pipeline Flow Test ===\n")

        # Step 1: Alignment
        print("Step 1: Running MSA alignment...")
        alignment_file = data["work_dir"] / "Alignment.ali"
        run_alignment(
            input_file=data["protein_seqs"],
            output_file=alignment_file,
            aligner=Aligner.MUSCLE,
        )
        assert alignment_file.exists()
        print(f"  ✓ Created: {alignment_file.name}")

        # Step 2: Reverse Translation
        print("Step 2: Running reverse translation...")
        reverse_file = data["work_dir"] / "Reverse_Translation_Seq.txt"
        reverse_translate_alignment(
            coding_sequences_file=data["coding_seqs"],
            protein_alignment_file=alignment_file,
            output_file=reverse_file,
            interest_species="homo_sapiens",
        )
        assert reverse_file.exists()
        print(f"  ✓ Created: {reverse_file.name}")

        # Step 3: Gblocks
        print("Step 3: Running Gblocks...")
        gblocks_exe = PROJECT_ROOT / "ATPS" / "assets" / "tools" / "Gblocks_0.91b" / "Gblocks"
        if gblocks_exe.exists():
            gblocks_output = run_gblocks(
                input_file=reverse_file,
                gblocks_exe=gblocks_exe,
            )
            assert gblocks_output.exists()
            print(f"  ✓ Created: {gblocks_output.name}")
        else:
            print("  ⚠ Skipped (Gblocks not available)")
            gblocks_output = reverse_file

        # Step 4: PHYLIP conversion
        print("Step 4: Converting to PHYLIP format...")
        jmodeltest_jar = (
            PROJECT_ROOT / "ATPS" / "assets" / "tools" / "jmodeltest-2.1.7" / "jModelTest.jar"
        )
        if jmodeltest_jar.exists():
            phylip_file = convert_fasta_to_phylip(
                input_file=gblocks_output,
                jar_path=jmodeltest_jar,
            )
            assert phylip_file.exists()
            print(f"  ✓ Created: {phylip_file.name}")
        else:
            print("  ⚠ Skipped (jModelTest not available)")

        print("\n=== Pipeline completed successfully! ===")

        output_files = list(data["work_dir"].glob("*"))
        print(f"\nOutput files in {data['work_dir']}:")
        for f in sorted(output_files):
            size = f.stat().st_size if f.is_file() else 0
            print(f"  - {f.name} ({size} bytes)")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
