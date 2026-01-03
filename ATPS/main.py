#!/usr/bin/env python3
"""ATPS - Automated Tool for Positive Selection Analysis.

This pipeline performs:
1. Sequence fetching from NCBI or local files
2. Multiple sequence alignment
3. Reverse translation to codon alignment
4. Gblocks filtering
5. Model selection with jModelTest
6. Phylogenetic tree construction with PhyML
7. Positive selection analysis with codeml (PAML)
8. Statistical analysis and output generation

Usage:
    python main.py -G gene1,gene2 -S species1,species2 -I interest_species -O output_dir

Arguments:
    -G  : Comma-separated list of gene names (required)
    -S  : Comma-separated list of species names (optional if using -IF)
    -O  : Output directory path (required)
    -I  : Species of interest for branch-site models (optional)
    -IF : Input folder with FASTA files (optional, alternative to -S)
    -A  : Alignment method (muscle, clustalo, mafft) (default: muscle)
    -R  : Bootstrap replicates for PhyML (default: 100)
    -GS : Gblocks stringency (T/F) (default: T)
"""
from __future__ import annotations

import argparse
import csv
import gc
import logging
import sys
from pathlib import Path
from typing import List, Optional, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# Optional imports with graceful fallback
try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    np = None
    _HAS_NUMPY = False

try:
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:
    pd = None
    _HAS_PANDAS = False

try:
    import mne
    _HAS_MNE = True
except ImportError:
    mne = None
    _HAS_MNE = False

# Import ATPS modules
from ATPS.utils.files import (
    PipelineSession,
    cleanup_intermediate_files,
    create_codeml_dirs,
    delete_codeml_dirs,
    save_gene_results,
)
from ATPS.utils.fetchers import fetch_and_save_sequences, count_fetched_species
from ATPS.utils.gene_operations import reverse_translate_alignment
from ATPS.utils.alignment import run_alignment, Aligner
from ATPS.utils.converters import convert_fasta_to_phylip
from ATPS.utils.parsers import (
    parse_jmodeltest,
    parse_jmodeltest_fallback,
    remove_branch_lengths,
    parse_beb_results,
    map_beb_to_original_positions,
)
from ATPS.models.gblocks import run_gblocks, remove_spaces
from ATPS.models.jmodeltest import run_jmodeltest
from ATPS.models.phyml import run_phyml
from ATPS.models.codeml import model078, model8a, model2a, model2, codeml_output, hashing


# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------
def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(
        description="ATPS - Automated Tool for Positive Selection Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "-G", "--genes",
        type=str,
        required=True,
        help="Comma-separated list of gene names",
    )
    parser.add_argument(
        "-S", "--species",
        type=str,
        default="",
        help="Comma-separated list of species names",
    )
    parser.add_argument(
        "-O", "--output",
        type=str,
        required=True,
        help="Output directory path",
    )
    parser.add_argument(
        "-I", "--interest",
        type=str,
        default="",
        help="Species of interest for branch-site models",
    )
    parser.add_argument(
        "-IF", "--input-folder",
        type=str,
        default="",
        help="Input folder containing FASTA files",
    )
    parser.add_argument(
        "-A", "--aligner",
        type=str,
        default="muscle",
        choices=["muscle", "clustalo", "mafft"],
        help="Alignment method (default: muscle)",
    )
    parser.add_argument(
        "-R", "--replicates",
        type=int,
        default=100,
        help="Bootstrap replicates for PhyML (default: 100)",
    )
    parser.add_argument(
        "-GS", "--gblocks-stringency",
        type=str,
        default="T",
        choices=["T", "F"],
        help="Gblocks stringency (T=strict, F=relaxed)",
    )
    parser.add_argument(
        "-E", "--email",
        type=str,
        default="",
        help="Email for NCBI Entrez queries",
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Use faster settings for jModelTest (3 schemes instead of 11)",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep intermediate files for debugging",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Pipeline Steps
# ---------------------------------------------------------------------------
def run_sequence_fetching(
    gene: str,
    species_list: Optional[List[str]],
    interest_species: str,
    input_folder: Optional[Path],
    email: str,
    work_dir: Path,
) -> Tuple[bool, str]:
    """Fetch sequences for a gene.

    Args:
        gene: Gene name.
        species_list: List of species to fetch.
        interest_species: Species of interest.
        input_folder: Folder with local FASTA files (if any).
        email: Email for NCBI.
        work_dir: Working directory.

    Returns:
        Tuple of (success, updated_interest_species).
    """
    logger.info("=" * 60)
    logger.info("Fetching sequences for gene: %s", gene)
    logger.info("=" * 60)

    fetch_from_ncbi = input_folder is None

    success, interest = fetch_and_save_sequences(
        protein=gene,
        species_list=species_list,
        interest_species=interest_species,
        fetch_from_ncbi=fetch_from_ncbi,
        gene_path=input_folder,
        email=email if fetch_from_ncbi else None,
        output_dir=work_dir,
    )

    if not success:
        logger.warning("No sequences found for gene: %s", gene)
        return False, interest_species

    return True, interest.lower() if interest else interest_species


def run_alignment_step(
    work_dir: Path,
    aligner: str,
) -> Path:
    """Run multiple sequence alignment.

    Args:
        work_dir: Working directory.
        aligner: Alignment method name.

    Returns:
        Path to alignment file.
    """
    logger.info("=" * 60)
    logger.info("Running multiple sequence alignment with %s", aligner)
    logger.info("=" * 60)

    input_file = work_dir / "ProteinSequences.fasta"
    output_file = work_dir / "Alignment.ali"

    # Map string to Aligner enum
    aligner_map = {
        "muscle": Aligner.MUSCLE,
        "clustalo": Aligner.CLUSTALO,
        "mafft": Aligner.MAFFT,
    }

    run_alignment(
        input_file=input_file,
        output_file=output_file,
        aligner=aligner_map.get(aligner, Aligner.MUSCLE),
    )

    return output_file


def run_reverse_translation(
    work_dir: Path,
    interest_species: str,
) -> Path:
    """Reverse translate protein alignment to codon alignment.

    Args:
        work_dir: Working directory.
        interest_species: Species of interest.

    Returns:
        Path to reverse-translated alignment.
    """
    logger.info("Reverse translating alignment to codons")

    return reverse_translate_alignment(
        coding_sequences_file=work_dir / "CodingSequences.fasta",
        protein_alignment_file=work_dir / "Alignment.ali",
        output_file=work_dir / "Reverse_Translation_Seq.txt",
        interest_species=interest_species,
    )


def run_filtering_step(
    work_dir: Path,
    stringency: str,
) -> Path:
    """Run Gblocks filtering.

    Args:
        work_dir: Working directory.
        stringency: Gblocks stringency setting ("T" for strict, "F" to skip).

    Returns:
        Path to filtered alignment.
    """
    logger.info("=" * 60)
    logger.info("Filtering alignment with Gblocks")
    logger.info("=" * 60)

    input_file = work_dir / "Reverse_Translation_Seq.txt"

    # Run Gblocks (skip if stringency is "F")
    skip_gblocks = stringency.upper() == "F"
    gblocks_output = run_gblocks(input_file, skip=skip_gblocks)

    # Remove spaces from Gblocks output (FASTA format)
    cleaned_fasta = remove_spaces(gblocks_output)

    # Convert to PHYLIP format
    phylip_output = work_dir / "Reverse_Translation_Seq.txt-gb1.phy"
    convert_fasta_to_phylip(cleaned_fasta, phylip_output)

    return phylip_output


def run_model_selection(
    work_dir: Path,
    phylip_file: Path,
    fast_mode: bool = False,
) -> Tuple[str, str, str]:
    """Run jModelTest for model selection.

    Args:
        work_dir: Working directory.
        phylip_file: Input PHYLIP alignment file.
        fast_mode: Use fewer substitution schemes for faster execution.

    Returns:
        Tuple of (partition, frequencies, pinvar).
    """
    logger.info("=" * 60)
    logger.info("Running jModelTest for model selection%s", " (fast mode)" if fast_mode else "")
    logger.info("=" * 60)

    output_file = work_dir / "Jmodeltest_output"
    
    # Fast mode uses 3 schemes (24 models) instead of 11 (88 models)
    num_schemes = 3 if fast_mode else 11
    run_jmodeltest(phylip_file, output_file, num_substitution_schemes=num_schemes)

    # Parse results
    try:
        result = parse_jmodeltest(output_file)
        partition = result.partition
        freq = result.frequencies
        pinvar = result.pinvar
    except Exception as e:
        logger.warning("Failed to parse primary model, trying fallback: %s", e)
        result = parse_jmodeltest_fallback(output_file)
        partition = result.partition
        freq = result.frequencies
        pinvar = result.pinvar

    # Validate results
    if not partition or not freq:
        logger.warning("Incomplete model parameters, using fallback")
        result = parse_jmodeltest_fallback(output_file)
        partition = result.partition
        freq = result.frequencies
        pinvar = result.pinvar

    if not pinvar:
        pinvar = "e"

    logger.info("Model selection: partition=%s, pinvar=%s", partition, pinvar)
    return partition, freq, pinvar


def run_tree_building(
    work_dir: Path,
    phylip_file: Path,
    partition: str,
    freq: str,
    pinvar: str,
    replicates: int,
) -> Path:
    """Build phylogenetic tree with PhyML.

    Args:
        work_dir: Working directory.
        phylip_file: Input alignment file.
        partition: Substitution model partition.
        freq: Nucleotide frequencies.
        pinvar: Proportion of invariable sites.
        replicates: Bootstrap replicates.

    Returns:
        Path to tree file without branch lengths.
    """
    logger.info("=" * 60)
    logger.info("Building phylogenetic tree with PhyML")
    logger.info("=" * 60)

    outputs = run_phyml(
        input_file=phylip_file,
        output_dir=work_dir,
        partition=partition,
        freq=freq,
        pinvar=pinvar,
        bootstrap_replicates=replicates,
    )

    # PhyML outputs the tree directly - get the tree file
    tree_file = outputs.get("tree")
    if not tree_file or not Path(tree_file).exists():
        raise RuntimeError("PhyML did not produce a tree file")

    # The PhyML tree is already in Newick format, copy to standard name
    newick_file = work_dir / "Species_Phylogenetic_tree_newick.nwk"
    if Path(tree_file) != newick_file:
        import shutil
        shutil.copy2(tree_file, newick_file)

    # Remove branch lengths for codeml
    nodist_file = remove_branch_lengths(newick_file)

    return nodist_file


def run_codeml_analysis(
    work_dir: Path,
    interest_species: str,
    gene: str,
) -> List:
    """Run codeml positive selection analysis.

    Args:
        work_dir: Working directory.
        interest_species: Species of interest.
        gene: Gene name.

    Returns:
        List of model results for CSV output.
    """
    logger.info("=" * 60)
    logger.info("Running codeml positive selection analysis")
    logger.info("=" * 60)

    # Create codeml directories
    create_codeml_dirs(work_dir)

    # Run models
    logger.info("Running model 0, 7, 8 (sites models)")
    model078(base_dir=work_dir)

    # Parse BEB results
    genes_without_beb = []
    try:
        beb_sites = parse_beb_results(work_dir / "codeml078" / "codeml078_mlc.txt")
        if beb_sites:
            gblocks_htm = work_dir / "Reverse_Translation_Seq.txt-gb1.htm"
            if gblocks_htm.exists():
                map_beb_to_original_positions(
                    beb_sites,
                    gblocks_htm,
                    work_dir / "BEB.csv",
                )
    except Exception as e:
        logger.warning("Could not parse BEB results: %s", e)
        genes_without_beb.append(gene)

    # Write genes without BEB
    if genes_without_beb:
        with open(work_dir / "genes_without_BEB.txt", "a") as f:
            f.write("\n".join(genes_without_beb) + "\n")

    logger.info("Running model 8a (null model)")
    model8a(base_dir=work_dir)

    # Branch-site models (require interest species)
    has_branch_site = bool(interest_species)
    if has_branch_site:
        logger.info("Running branch-site models (2a, 2) for species: %s", interest_species)
        hashing(interest_species)

        logger.info("Running model 2a")
        model2a(base_dir=work_dir)

        logger.info("Running model 2")
        model2(base_dir=work_dir)

    # Get model output (state=1 if branch-site models run, 0 otherwise)
    models = codeml_output(state=1 if has_branch_site else 0, protein=gene)

    return models


def apply_multiple_testing_correction(
    output_file: Path,
    save_dir: Path,
) -> Path:
    """Apply Bonferroni correction to p-values.

    Args:
        output_file: Path to Study_Output.csv.
        save_dir: Output directory.

    Returns:
        Path to corrected output file.
    """
    if not _HAS_PANDAS or not _HAS_NUMPY or not _HAS_MNE:
        logger.warning(
            "pandas, numpy, or mne not available. Skipping multiple testing correction."
        )
        return output_file

    logger.info("Applying Bonferroni correction to p-values")

    df = pd.read_csv(output_file)

    # Apply correction to each p-value column
    pvalue_cols = ["p-v7vs8", "p-v8avs8", "p-v2avs2"]
    adj_cols = ["adj78", "adj88a", "adj22a"]

    for pval_col, adj_col in zip(pvalue_cols, adj_cols):
        if pval_col in df.columns:
            try:
                pvals = np.asfarray(df[pval_col])
                _, adjusted = mne.stats.bonferroni_correction(pvals, alpha=0.05)
                df[adj_col] = adjusted
            except Exception as e:
                logger.warning("Could not correct %s: %s", pval_col, e)

    corrected_file = save_dir / "Study_Output_adjPvalue.csv"
    df.to_csv(corrected_file, index=False)

    logger.info("Saved corrected p-values to %s", corrected_file)
    return corrected_file


# ---------------------------------------------------------------------------
# Main Pipeline
# ---------------------------------------------------------------------------
def run_pipeline(args: argparse.Namespace) -> int:
    """Run the complete ATPS pipeline.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Exit code (0 for success, non-zero for failure).
    """
    # Parse arguments
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    species_list = [s.strip() for s in args.species.split(",") if s.strip()] if args.species else None
    output_dir = Path(args.output)
    input_folder = Path(args.input_folder) if args.input_folder else None
    interest_species = args.interest.lower() if args.interest else ""
    aligner = args.aligner.lower()
    replicates = args.replicates
    gblocks_stringency = args.gblocks_stringency
    email = args.email
    fast_mode = args.fast

    # Validate output directory
    if not output_dir.exists():
        logger.error("Output directory does not exist: %s", output_dir)
        return 1

    # If using input folder, get genes from FASTA files
    if input_folder:
        genes = [f.stem for f in input_folder.glob("*.fasta")]
        if not genes:
            logger.error("No FASTA files found in: %s", input_folder)
            return 1
        logger.info("Found %d genes from input folder", len(genes))

    # Validate email for NCBI fetching
    if not input_folder and not email:
        logger.warning("No email provided for NCBI Entrez. Using placeholder.")
        email = "user@example.com"

    # Create output CSV
    study_output = output_dir / "Study_Output.csv"
    with open(study_output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Name", "Model0", "Model7", "Model8", "Model8a",
            "Model2a", "Model2", "LRT(M7_vs_M8)", "LRT(M8a_vs_M8)",
            "LRT(M2a_vs_M2)", "p-v7vs8", "p-v8avs8", "p-v2avs2",
        ])

    # Process each gene
    successful_genes = 0
    for gene in genes:
        logger.info("")
        logger.info("=" * 70)
        logger.info("Processing gene: %s", gene)
        logger.info("=" * 70)

        # Use PipelineSession for automatic cleanup
        with PipelineSession(
            output_dir=output_dir,
            gene_name=gene,
            keep_intermediates=args.keep_intermediates,
        ) as session:
            work_dir = session.work_dir

            try:
                # Step 1: Fetch sequences
                success, interest = run_sequence_fetching(
                    gene=gene,
                    species_list=species_list,
                    interest_species=interest_species,
                    input_folder=input_folder,
                    email=email,
                    work_dir=work_dir,
                )
                if not success:
                    logger.warning("Skipping gene %s: no sequences found", gene)
                    continue

                # Step 2: Multiple sequence alignment
                run_alignment_step(work_dir, aligner)

                # Step 3: Reverse translation
                run_reverse_translation(work_dir, interest)

                # Step 4: Gblocks filtering
                phylip_file = run_filtering_step(work_dir, gblocks_stringency)

                # Step 5: Model selection
                partition, freq, pinvar = run_model_selection(work_dir, phylip_file, fast_mode)

                # Step 6: Tree building
                run_tree_building(
                    work_dir, phylip_file, partition, freq, pinvar, replicates
                )

                # Step 7: Codeml analysis
                models = run_codeml_analysis(work_dir, interest, gene)

                # Write results to CSV
                with open(study_output, "a", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(models)

                # Mark important files to keep
                for pattern in ["*.csv", "*.nwk", "codeml*/"]:
                    for file_path in work_dir.glob(pattern):
                        session.keep(file_path)

                successful_genes += 1
                logger.info("Successfully processed gene: %s", gene)

            except Exception as e:
                logger.error("Failed to process gene %s: %s", gene, e)
                if args.verbose:
                    import traceback
                    traceback.print_exc()
                continue

            finally:
                # Cleanup codeml directories
                try:
                    delete_codeml_dirs(work_dir)
                except Exception:
                    pass

                # Force garbage collection
                gc.collect()

    # Post-processing
    logger.info("")
    logger.info("=" * 70)
    logger.info("Pipeline complete: %d/%d genes processed successfully", successful_genes, len(genes))
    logger.info("=" * 70)

    if successful_genes > 0:
        # Apply multiple testing correction
        apply_multiple_testing_correction(study_output, output_dir)

        # Count fetched species
        try:
            count_fetched_species(output_dir, output_dir / "fetched_species.csv")
        except Exception as e:
            logger.warning("Could not count fetched species: %s", e)

    return 0 if successful_genes > 0 else 1


def main() -> int:
    """Main entry point."""
    args = parse_arguments()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        return run_pipeline(args)
    except KeyboardInterrupt:
        logger.info("Pipeline interrupted by user")
        return 130
    except Exception as e:
        logger.error("Pipeline failed: %s", e)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
