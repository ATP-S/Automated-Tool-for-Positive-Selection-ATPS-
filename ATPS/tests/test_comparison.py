"""
Comparison test that validates function outputs against ground truth in example_tests.

This test compares the actual function outputs with the expected outputs
from the example_tests directory to ensure correctness.
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

# Path to example test data (ground truth)
EXAMPLE_BASE = PROJECT_ROOT / "example_tests"


def parse_fasta(content: str) -> dict[str, str]:
    """Parse FASTA content into a dictionary of {id: sequence}."""
    sequences = {}
    current_id = None
    current_seq: list[str] = []

    for line in content.strip().split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if current_id:
                sequences[current_id] = "".join(current_seq)
            current_id = line[1:].split()[0].lower()
            current_seq = []
        else:
            current_seq.append(line.replace(" ", ""))

    if current_id:
        sequences[current_id] = "".join(current_seq)

    return sequences


def compare_sequences(actual: dict, expected: dict) -> dict:
    """Compare two sequence dictionaries and return comparison results."""
    results: dict = {
        "species_match": False,
        "sequence_similarity": {},
        "length_comparison": {},
        "details": [],
    }

    actual_species = set(actual.keys())
    expected_species = set(expected.keys())

    results["species_match"] = actual_species == expected_species
    results["actual_species"] = list(actual_species)
    results["expected_species"] = list(expected_species)
    results["missing_species"] = list(expected_species - actual_species)
    results["extra_species"] = list(actual_species - expected_species)

    common_species = actual_species & expected_species
    for species in common_species:
        actual_seq = actual[species].upper().replace("-", "").replace(" ", "")
        expected_seq = expected[species].upper().replace("-", "").replace(" ", "")

        results["length_comparison"][species] = {
            "actual": len(actual_seq),
            "expected": len(expected_seq),
            "match": len(actual_seq) == len(expected_seq),
        }

        if len(actual_seq) == len(expected_seq):
            matches = sum(a == b for a, b in zip(actual_seq, expected_seq))
            similarity = matches / len(actual_seq) if actual_seq else 1.0
        else:
            min_len = min(len(actual_seq), len(expected_seq))
            if min_len > 0:
                matches = sum(
                    actual_seq[i] == expected_seq[i] for i in range(min_len)
                )
                similarity = matches / max(len(actual_seq), len(expected_seq))
            else:
                similarity = 1.0 if actual_seq == expected_seq else 0.0

        results["sequence_similarity"][species] = similarity

    return results


class TestCompareWithGroundTruth:
    """Compare function outputs with ground truth from example_tests."""

    def test_compare_alignment_output(self, tmp_path):
        """Compare alignment output with ground truth."""
        print("\n" + "=" * 60)
        print("COMPARISON: Alignment Output vs Ground Truth")
        print("=" * 60)

        input_file = (
            EXAMPLE_BASE / "Alignment_MUSCLE_ClustalO_MAFFT" / "input" / "ProteinSequences.fasta"
        )
        expected_file = (
            EXAMPLE_BASE / "Alignment_MUSCLE_ClustalO_MAFFT" / "output" / "Alignment.ali"
        )

        if not input_file.exists() or not expected_file.exists():
            pytest.skip("Ground truth files not available")

        output_file = tmp_path / "Alignment.ali"
        run_alignment(
            input_file=input_file,
            output_file=output_file,
            aligner=Aligner.MUSCLE,
        )

        actual = parse_fasta(output_file.read_text())
        expected = parse_fasta(expected_file.read_text())

        results = compare_sequences(actual, expected)

        print(f"\n📊 Species Match: {'✓' if results['species_match'] else '✗'}")
        print(f"   Actual species: {results['actual_species']}")
        print(f"   Expected species: {results['expected_species']}")

        print("\n📏 Length Comparison:")
        for species, lengths in results["length_comparison"].items():
            status = "✓" if lengths["match"] else "✗"
            print(f"   {species}: {lengths['actual']} vs {lengths['expected']} {status}")

        print("\n🔬 Sequence Similarity:")
        for species, similarity in results["sequence_similarity"].items():
            status = "✓" if similarity >= 0.95 else "⚠"
            print(f"   {species}: {similarity:.1%} {status}")

        assert results["species_match"], f"Species mismatch: {results}"
        avg_similarity = sum(results["sequence_similarity"].values()) / len(
            results["sequence_similarity"]
        )
        print(f"\n📈 Average Similarity: {avg_similarity:.1%}")

        assert avg_similarity >= 0.90, f"Low sequence similarity: {avg_similarity:.1%}"

    def test_compare_reverse_translation_output(self, tmp_path):
        """Compare reverse translation output with ground truth."""
        print("\n" + "=" * 60)
        print("COMPARISON: Reverse Translation Output vs Ground Truth")
        print("=" * 60)

        coding_file = EXAMPLE_BASE / "Reverse_Translation_Python" / "input" / "CodingSequences.fasta"
        alignment_file = EXAMPLE_BASE / "Reverse_Translation_Python" / "input" / "Alignment.ali"
        expected_file = (
            EXAMPLE_BASE / "Reverse_Translation_Python" / "output" / "Reverse_Translation_Seq.txt"
        )

        if not all(f.exists() for f in [coding_file, alignment_file, expected_file]):
            pytest.skip("Ground truth files not available")

        output_file = tmp_path / "Reverse_Translation_Seq.txt"
        reverse_translate_alignment(
            coding_sequences_file=coding_file,
            protein_alignment_file=alignment_file,
            output_file=output_file,
            interest_species="homo_sapiens",
        )

        actual = parse_fasta(output_file.read_text())
        expected = parse_fasta(expected_file.read_text())

        results = compare_sequences(actual, expected)

        print(f"\n📊 Species Match: {'✓' if results['species_match'] else '✗'}")
        print(f"   Actual species: {results['actual_species']}")
        print(f"   Expected species: {results['expected_species']}")

        print("\n📏 Length Comparison:")
        for species, lengths in results["length_comparison"].items():
            status = "✓" if lengths["match"] else "✗"
            print(f"   {species}: {lengths['actual']} vs {lengths['expected']} {status}")

        print("\n🔬 Sequence Similarity:")
        for species, similarity in results["sequence_similarity"].items():
            status = "✓" if similarity >= 0.99 else "✗"
            print(f"   {species}: {similarity:.1%} {status}")

        avg_similarity = sum(results["sequence_similarity"].values()) / len(
            results["sequence_similarity"]
        )
        print(f"\n📈 Average Similarity: {avg_similarity:.1%}")

        assert results["species_match"], f"Species mismatch: {results}"
        assert avg_similarity >= 0.95, f"Low sequence similarity: {avg_similarity:.1%}"

    def test_compare_gblocks_output(self, tmp_path):
        """Compare Gblocks output with ground truth."""
        print("\n" + "=" * 60)
        print("COMPARISON: Gblocks Output vs Ground Truth")
        print("=" * 60)

        input_file = EXAMPLE_BASE / "Gblocks" / "input" / "Reverse_Translation_Seq.txt"
        expected_file = EXAMPLE_BASE / "Gblocks" / "output" / "Reverse_Translation_Seq.txt-gb1"

        if not input_file.exists() or not expected_file.exists():
            pytest.skip("Ground truth files not available")

        gblocks_exe = PROJECT_ROOT / "ATPS" / "assets" / "tools" / "Gblocks_0.91b" / "Gblocks"
        if not gblocks_exe.exists():
            pytest.skip("Gblocks not available")

        test_input = tmp_path / "Reverse_Translation_Seq.txt"
        shutil.copy(input_file, test_input)

        output_file = run_gblocks(
            input_file=test_input,
            gblocks_exe=gblocks_exe,
        )

        actual = parse_fasta(output_file.read_text())
        expected = parse_fasta(expected_file.read_text())

        results = compare_sequences(actual, expected)

        print(f"\n📊 Species Match: {'✓' if results['species_match'] else '✗'}")
        print(f"   Actual species: {results['actual_species']}")
        print(f"   Expected species: {results['expected_species']}")

        print("\n📏 Length Comparison (after trimming):")
        for species, lengths in results["length_comparison"].items():
            status = "✓" if lengths["match"] else "⚠"
            print(f"   {species}: {lengths['actual']} vs {lengths['expected']} {status}")

        print("\n🔬 Sequence Similarity:")
        for species, similarity in results["sequence_similarity"].items():
            status = "✓" if similarity >= 0.95 else "⚠"
            print(f"   {species}: {similarity:.1%} {status}")

        avg_similarity = sum(results["sequence_similarity"].values()) / len(
            results["sequence_similarity"]
        )
        print(f"\n📈 Average Similarity: {avg_similarity:.1%}")

        assert results["species_match"], f"Species mismatch: {results}"
        assert avg_similarity >= 0.90, f"Low sequence similarity: {avg_similarity:.1%}"

    def test_compare_phylip_conversion_output(self, tmp_path):
        """Compare PHYLIP conversion output with ground truth."""
        print("\n" + "=" * 60)
        print("COMPARISON: PHYLIP Conversion Output vs Ground Truth")
        print("=" * 60)

        input_file = (
            EXAMPLE_BASE / "jModelTest_Conversion" / "input" / "Reverse_Translation_Seq.txt-gb1"
        )
        expected_file = (
            EXAMPLE_BASE / "jModelTest_Conversion" / "output" / "Reverse_Translation_Seq.txt-gb1.phy"
        )

        if not input_file.exists() or not expected_file.exists():
            pytest.skip("Ground truth files not available")

        jmodeltest_jar = (
            PROJECT_ROOT / "ATPS" / "assets" / "tools" / "jmodeltest-2.1.7" / "jModelTest.jar"
        )
        if not jmodeltest_jar.exists():
            pytest.skip("jModelTest not available")

        test_input = tmp_path / "Reverse_Translation_Seq.txt-gb1"
        shutil.copy(input_file, test_input)

        output_file = convert_fasta_to_phylip(
            input_file=test_input,
            jar_path=jmodeltest_jar,
        )

        actual_content = output_file.read_text().strip()
        expected_content = expected_file.read_text().strip()

        actual_lines = actual_content.split("\n")
        expected_lines = expected_content.split("\n")

        actual_header = actual_lines[0].strip().split()
        expected_header = expected_lines[0].strip().split()

        print("\n📊 PHYLIP Header Comparison:")
        print(f"   Actual:   {actual_header[0]} sequences, {actual_header[1]} positions")
        print(f"   Expected: {expected_header[0]} sequences, {expected_header[1]} positions")

        header_match = actual_header == expected_header
        print(f"   Match: {'✓' if header_match else '✗'}")

        actual_data = "\n".join(actual_lines[1:])
        expected_data = "\n".join(expected_lines[1:])

        actual_normalized = " ".join(actual_data.split())
        expected_normalized = " ".join(expected_data.split())

        content_match = actual_normalized == expected_normalized
        print(f"\n🔬 Sequence Data Match: {'✓' if content_match else '✗'}")

        if not content_match:
            min_len = min(len(actual_normalized), len(expected_normalized))
            matches = sum(
                a == b for a, b in zip(actual_normalized[:min_len], expected_normalized[:min_len])
            )
            similarity = matches / max(len(actual_normalized), len(expected_normalized))
            print(f"   Similarity: {similarity:.1%}")
        else:
            similarity = 1.0

        print(f"\n📈 Overall Match: {'✓' if header_match and content_match else '⚠'}")

        assert int(actual_header[0]) == int(expected_header[0]), "Number of sequences mismatch"
        assert similarity >= 0.95, f"Low content similarity: {similarity:.1%}"

    def test_full_pipeline_comparison(self, tmp_path):
        """Run full pipeline and compare all outputs with ground truth."""
        print("\n" + "=" * 70)
        print("FULL PIPELINE COMPARISON WITH GROUND TRUTH")
        print("=" * 70)

        coding_src = EXAMPLE_BASE / "Fetching_BioPython" / "output" / "CodingSequences.fasta"
        protein_src = EXAMPLE_BASE / "Fetching_BioPython" / "output" / "ProteinSequences.fasta"

        if not coding_src.exists() or not protein_src.exists():
            pytest.skip("Ground truth input files not available")

        coding_file = tmp_path / "CodingSequences.fasta"
        protein_file = tmp_path / "ProteinSequences.fasta"
        shutil.copy(coding_src, coding_file)
        shutil.copy(protein_src, protein_file)

        results_summary: dict[str, float] = {}

        # Step 1: Alignment
        print("\n📌 Step 1: MSA Alignment")
        alignment_file = tmp_path / "Alignment.ali"
        run_alignment(protein_file, alignment_file, aligner=Aligner.MUSCLE)

        expected_alignment = (
            EXAMPLE_BASE / "Alignment_MUSCLE_ClustalO_MAFFT" / "output" / "Alignment.ali"
        )
        if expected_alignment.exists():
            actual = parse_fasta(alignment_file.read_text())
            expected = parse_fasta(expected_alignment.read_text())
            comparison = compare_sequences(actual, expected)
            avg_sim = sum(comparison["sequence_similarity"].values()) / len(
                comparison["sequence_similarity"]
            )
            results_summary["alignment"] = avg_sim
            print(f"   Similarity: {avg_sim:.1%} {'✓' if avg_sim >= 0.90 else '⚠'}")

        # Step 2: Reverse Translation
        print("\n📌 Step 2: Reverse Translation")
        reverse_file = tmp_path / "Reverse_Translation_Seq.txt"
        reverse_translate_alignment(
            coding_file, alignment_file, reverse_file, interest_species="homo_sapiens"
        )

        expected_reverse = (
            EXAMPLE_BASE / "Reverse_Translation_Python" / "output" / "Reverse_Translation_Seq.txt"
        )
        if expected_reverse.exists():
            actual = parse_fasta(reverse_file.read_text())
            expected = parse_fasta(expected_reverse.read_text())
            comparison = compare_sequences(actual, expected)
            avg_sim = sum(comparison["sequence_similarity"].values()) / len(
                comparison["sequence_similarity"]
            )
            results_summary["reverse_translation"] = avg_sim
            print(f"   Similarity: {avg_sim:.1%} {'✓' if avg_sim >= 0.95 else '⚠'}")

        # Step 3: Gblocks
        print("\n📌 Step 3: Gblocks")
        gblocks_exe = PROJECT_ROOT / "ATPS" / "assets" / "tools" / "Gblocks_0.91b" / "Gblocks"
        gblocks_output = None
        if gblocks_exe.exists():
            gblocks_output = run_gblocks(reverse_file, gblocks_exe=gblocks_exe)

            expected_gblocks = (
                EXAMPLE_BASE / "Gblocks" / "output" / "Reverse_Translation_Seq.txt-gb1"
            )
            if expected_gblocks.exists():
                actual = parse_fasta(gblocks_output.read_text())
                expected = parse_fasta(expected_gblocks.read_text())
                comparison = compare_sequences(actual, expected)
                avg_sim = sum(comparison["sequence_similarity"].values()) / len(
                    comparison["sequence_similarity"]
                )
                results_summary["gblocks"] = avg_sim
                print(f"   Similarity: {avg_sim:.1%} {'✓' if avg_sim >= 0.90 else '⚠'}")
        else:
            print("   ⚠ Gblocks not available")

        # Step 4: PHYLIP Conversion
        print("\n📌 Step 4: PHYLIP Conversion")
        jmodeltest_jar = (
            PROJECT_ROOT / "ATPS" / "assets" / "tools" / "jmodeltest-2.1.7" / "jModelTest.jar"
        )
        if jmodeltest_jar.exists() and gblocks_output is not None:
            phylip_file = convert_fasta_to_phylip(gblocks_output, jar_path=jmodeltest_jar)

            expected_phylip = (
                EXAMPLE_BASE
                / "jModelTest_Conversion"
                / "output"
                / "Reverse_Translation_Seq.txt-gb1.phy"
            )
            if expected_phylip.exists():
                actual_header = phylip_file.read_text().split("\n")[0].strip().split()
                expected_header = expected_phylip.read_text().split("\n")[0].strip().split()
                header_match = actual_header[0] == expected_header[0]
                results_summary["phylip_conversion"] = 1.0 if header_match else 0.5
                print(f"   Header Match: {'✓' if header_match else '⚠'}")
        else:
            print("   ⚠ jModelTest not available")

        # Summary
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        print("\n📊 Comparison Results:")
        for step, similarity in results_summary.items():
            status = "✓ PASS" if similarity >= 0.90 else "⚠ CHECK"
            print(f"   {step:25s}: {similarity:6.1%} {status}")

        overall = sum(results_summary.values()) / len(results_summary) if results_summary else 0
        print(f"\n📈 Overall Score: {overall:.1%}")

        assert overall >= 0.85, f"Overall similarity too low: {overall:.1%}"
        print("\n✅ All comparisons passed!")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
