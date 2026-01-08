"""codeml helpers.

This module expects per-model ctl files to already exist in
`ATPS/models/configs/` (e.g. `codeml078.ctl`, `codeml8a.ctl`, `codeml2a.ctl`, `codeml2.ctl`).
Behavior:
- Each `model...()` function copies the appropriate ctl from the `configs`
  directory to `codeml.ctl` in the current working directory and runs
  the `codeml` command (must be available on PATH or provided as full path).
- `codeml_output` parses the model outputs and writes `Gene_Output.csv`.
"""

from __future__ import annotations

import csv
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path

try:
    from scipy.stats import chi2
except Exception:  # pragma: no cover - present at runtime if scipy installed
    chi2 = None

logger = logging.getLogger(__name__)


def _parser_from_fileobj(fileobj) -> list[float]:
    lis = []
    for line in fileobj:
        if "lnL" in line:
            # try to extract the float after '=' if present, otherwise fall back
            m = re.search(r"=\s*([-+]?\d+\.\d+)", line)
            if m:
                lis.append(float(m.group(1)))
            else:
                # original code looked for '):' then whitespace
                start = line.find("):")
                end = line.find("      ")
                try:
                    lis.append(float(line[start + 3 : end].strip()))
                except Exception:
                    continue
    return lis


def codeml_output(state: int, protein: str):
    """Create a CSV row of lnL and LRT/p-values from codeml outputs.

    Args:
        state: 0 or 1 — whether model2/model2a were run and should be parsed.
        protein: name of the protein (written into CSV).

    Returns:
        list containing the same fields written to `Gene_Output.csv`.
    """
    if chi2 is None:
        logger.warning(
            "scipy not available; p-values will be empty. Install scipy to enable p-value computation."
        )
        compute_pvals = False
    else:
        compute_pvals = True

    def parser_path(path: Path) -> list[float]:
        if not path.exists():
            return []
        with path.open("r") as fh:
            return _parser_from_fileobj(fh)

    reader078 = Path("codeml078/codeml078_mlc.txt")
    models078 = parser_path(reader078)

    reader8a = Path("codeml8a/codeml8a_mlc.txt")
    model8a = parser_path(reader8a)

    model2a = [""]
    model2 = [""]
    lrt22a = ""

    if state == 1:
        reader2a = Path("codeml2a/codeml2a_mlc.txt")
        model2a = parser_path(reader2a)

        reader2 = Path("codeml2/codeml2_mlc.txt")
        model2 = parser_path(reader2)

        if model2a and model2:
            lrt22a = 2 * (model2a[0] - model2[0])
        else:
            lrt22a = ""

    # Be defensive about missing indices
    try:
        lrt78 = 2 * (models078[2] - models078[1])
    except Exception:
        lrt78 = ""
    try:
        lrt88a = 2 * (models078[2] - model8a[0])
    except Exception:
        lrt88a = ""

    if compute_pvals:
        chi78 = 1 - chi2.cdf(lrt78, 2) if isinstance(lrt78, (int, float)) else ""
        chi88a = 1 - chi2.cdf(lrt88a, 1) if isinstance(lrt88a, (int, float)) else ""
        chi22a = 1 - chi2.cdf(lrt22a, 1) if state == 1 and isinstance(lrt22a, (int, float)) else ""
    else:
        chi78 = ""
        chi88a = ""
        chi22a = ""

    out_path = Path("Gene_Output.csv")
    with out_path.open("w", newline="") as newfile:
        wr = csv.writer(newfile)
        wr.writerow(
            [
                "Name",
                "Model0",
                "Model7",
                "Model8",
                "Model8a",
                "Model2a",
                "Model2",
                "LRT(M7_vs_M8)",
                "LRT(M8a_vs_M8)",
                "LRT(M2a_vs_M2)",
                "p-v7vs8",
                "p-v8avs8",
                "p-v2avs2",
            ]
        )
        wr.writerow(
            [
                protein,
                models078[0] if len(models078) > 0 else "",
                models078[1] if len(models078) > 1 else "",
                models078[2] if len(models078) > 2 else "",
                model8a[0] if len(model8a) > 0 else "",
                model2a[0] if len(model2a) > 0 else "",
                model2[0] if len(model2) > 0 else "",
                lrt78,
                lrt88a,
                lrt22a,
                chi78,
                chi88a,
                chi22a,
            ]
        )

    return [
        protein,
        models078[0] if len(models078) > 0 else "",
        models078[1] if len(models078) > 1 else "",
        models078[2] if len(models078) > 2 else "",
        model8a[0] if len(model8a) > 0 else "",
        model2a[0] if len(model2a) > 0 else "",
        model2[0] if len(model2) > 0 else "",
        lrt78,
        lrt88a,
        lrt22a,
        chi78,
        chi88a,
        chi22a,
    ]


def hashing(interest: str):
    """Mark the species of interest in the tree and write a new tree file.

    Inserts ' #1' after the first occurrence of the interest name (case-insensitive)
    and writes to `Species_Phylogenetic_tree_newick_Interst#.nwk`.
    """
    src = Path("Species_Phylogenetic_tree_newick_nodistances.nwk")
    dest = Path("Species_Phylogenetic_tree_newick_Interst#.nwk")
    text = src.read_text()
    idx = text.lower().find(interest.lower())
    if idx == -1:
        # nothing found; write original tree to dest and warn
        dest.write_text(text)
        print(f"Interest '{interest}' not found in {src}; wrote original tree to {dest}")
        return
    idx_end = idx + len(interest)
    new_text = text[:idx_end] + " #1" + text[idx_end:]
    dest.write_text(new_text)
    print(f"Wrote tree with interest marker to {dest}")


def _copy_ctl_from_configs(fname: str, dest_dir: Path) -> None:
    """Copy a ctl file from the `configs` directory (module-relative) into `dest_dir` as `codeml.ctl`."""
    module_dir = Path(__file__).resolve().parent
    src = module_dir / "configs" / fname
    if not src.exists():
        raise FileNotFoundError(f"Expected ctl in configs: {src}")
    dest_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dest_dir / "codeml.ctl")


def _run_codeml(codeml_cmd: str = "codeml", cwd: Path | None = None):
    """Run codeml in the provided `cwd` (Path). Writes stdout/stderr into that directory.

    Returns the process return code.
    """
    if cwd is None:
        cwd = Path.cwd()
    try:
        completed = subprocess.run([codeml_cmd], cwd=str(cwd), capture_output=True, text=True)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"codeml executable not found: '{codeml_cmd}'. Provide full path or install PAML."
        ) from e
    # write logs for debugging into the working directory
    (cwd / "codeml_stdout.txt").write_text(completed.stdout)
    (cwd / "codeml_stderr.txt").write_text(completed.stderr)
    return completed.returncode


def model078(codeml_cmd: str = "codeml", base_dir: Path | str | None = None):
    """Run Codeml models 0,7,8.

    The model's ctl is copied into `base_dir/codeml078/codeml.ctl` and
    `codeml` is executed with that folder as the working directory.
    """
    if base_dir is None:
        module_dir = Path(__file__).resolve().parent
        out_dir = module_dir / "codeml078"
    else:
        out_dir = Path(base_dir) / "codeml078"
    _copy_ctl_from_configs("codeml078.ctl", out_dir)
    return _run_codeml(codeml_cmd, cwd=out_dir)


def model8a(codeml_cmd: str = "codeml", base_dir: Path | str | None = None):
    """Run Codeml model 8a in `base_dir/codeml8a/` (defaults to package models dir)."""
    if base_dir is None:
        module_dir = Path(__file__).resolve().parent
        out_dir = module_dir / "codeml8a"
    else:
        out_dir = Path(base_dir) / "codeml8a"
    _copy_ctl_from_configs("codeml8a.ctl", out_dir)
    return _run_codeml(codeml_cmd, cwd=out_dir)


def model2a(codeml_cmd: str = "codeml", base_dir: Path | str | None = None):
    """Run Codeml model 2a in `base_dir/codeml2a/` (defaults to package models dir)."""
    if base_dir is None:
        module_dir = Path(__file__).resolve().parent
        out_dir = module_dir / "codeml2a"
    else:
        out_dir = Path(base_dir) / "codeml2a"
    _copy_ctl_from_configs("codeml2a.ctl", out_dir)
    return _run_codeml(codeml_cmd, cwd=out_dir)


def model2(codeml_cmd: str = "codeml", base_dir: Path | str | None = None):
    """Run Codeml model 2 in `base_dir/codeml2/` (defaults to package models dir)."""
    if base_dir is None:
        module_dir = Path(__file__).resolve().parent
        out_dir = module_dir / "codeml2"
    else:
        out_dir = Path(base_dir) / "codeml2"
    _copy_ctl_from_configs("codeml2.ctl", out_dir)
    return _run_codeml(codeml_cmd, cwd=out_dir)


def path_dir():
    return os.getcwd()


__all__ = [
    "codeml_output",
    "hashing",
    "model078",
    "model8a",
    "model2a",
    "model2",
    "path_dir",
]


# def codeml_output(state, protein):
#     """
#         Function creates a spreadsheet for the p-values

#     Args:
#         state (int 0 or 1): state shows whether model2 and model 2a is required
#         protein (string):  the name of the protein
#     """
#     def parser(file):
#         lis = []
#         lines = file.readlines()
#         for i in lines:
#             if "lnL" in i:
#                 start = i.find("):")
#                 end = i.find("      ")
#                 lis.append(float(i[start+3:end].strip()))
#         return lis

#     reader078 = open("codeml078/codeml078_mlc.txt","r")
#     models078 = parser(reader078)

#     reader8a = open("codeml8a/codeml8a_mlc.txt","r")
#     model8a = parser(reader8a)

#     model2a = [""]
#     model2 = [""]
#     lrt22a = ''

#     if state == 1:
#         reader2a = open("codeml2a/codeml2a_mlc.txt","r")
#         model2a = parser(reader2a)

#         reader2 = open("codeml2/codeml2_mlc.txt","r")
#         model2 = parser(reader2)

#         lrt22a = 2*(model2a[0] - model2[0])

#     lrt78 = 2*(models078[2] - models078[1])
#     lrt88a = 2*(models078[2] - model8a[0])

#     chi78 =  1 - chi2.cdf(lrt78, 2)
#     chi88a =  1 - chi2.cdf(lrt88a, 1)
#     chi22a = 1 - chi2.cdf(lrt22a, 1) if state == 1 else ''
#     with open("Gene_Output.csv", "w") as newfile:
#         wr = csv.writer(newfile)
#         wr.writerow(["Name","Model0", "Model7", "Model8","Model8a", "Model2a", "Model2", "LRT(M7_vs_M8)", "LRT(M8a_vs_M8)", "LRT(M2a_vs_M2)", "p-v7vs8", "p-v8avs8", "p-v2avs2"])
#         wr.writerow([protein, models078[0], models078[1], models078[2], model8a[0], model2a[0], model2[0],lrt78, lrt88a, lrt22a, chi78, chi88a, chi22a])
#         return [protein, models078[0], models078[1], models078[2], model8a[0], model2a[0], model2[0],lrt78, lrt88a, lrt22a, chi78, chi88a, chi22a]


# def hashing(interest):
#     """
#         Function that sets an indicator for the gene of interest for codeml

#     Args:
#         interest (_type_): _description_
#     """
#     newick_open = open("Species_Phylogenetic_tree_newick_nodistances.nwk")
#     newick_create = open("Species_Phylogenetic_tree_newick_Interst#.nwk", 'w')
#     newick_read = newick_open.read()
#     newick_list = newick_read.split(',')
#     result = newick_read.find(str(interest.lower()))
#     print(result)
#     result += len(interest)
#     print(result)
#     newick_read = newick_read[:result] + " #1" + newick_read[result:]
#     newick_open.close()
#     print(newick_read)
#     newick_create.write(newick_read)


# def model078():
#     """
#         Function that Runs Codeml models 0,7,8
#     """
#     f = open("codeml.ctl", "w")
#     m = open("codeml078.ctl" , "r")
#     f.writelines(m)
#     f.close()
#     os.system("codeml")

# def model8a():
#     """
#         Function that Runs Codeml model 8a
#     """
#     f = open("codeml.ctl", "w")
#     m = open("codeml8a.ctl" , "r")
#     f.writelines(m)
#     f.close()
#     os.system("codeml")

# def model2a():
#     """
#         Function that Runs Codeml model 2a
#     """
#     f = open("codeml.ctl", "w")
#     m = open("codeml2a.ctl" , "r")
#     f.writelines(m)
#     f.close()
#     os.system("codeml")

# #sys.exit(0)
# def model2():
#     """
#         Function that Runs Codeml model 2
#     """
#     f = open("codeml.ctl", "w")
#     m = open("codeml2.ctl" , "r")
#     f.writelines(m)
#     f.close()
#     os.system("codeml")

# def path_dir():
#     return os.getcwd()


# def codeml_creating_file():
#     """
#         function that creates configuration files for codeml
#     """


#     c078 = """          seqfile = Reverse_Translation_Seq.txt-gb1.fst * fasta file
#          treefile = Species_Phylogenetic_tree_newick_nodistances.nwk *
#           outfile = codeml078_mlc.txt

#             noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
#           verbose = 1   * 1: detailed output, 0: concise output
#           runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
#                         * 3: StepwiseAddition; (4,5):PerturbationNNI

#           seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
#         CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
#             clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
#             model = 0   * 0 - site models 2 - branch and branch-site models
#                         * models for codons:
#                             * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

#           NSsites = 0 7 8   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
#             icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below

#         fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
#             kappa = 2   * initial or fixed kappa
#         fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
#             omega = 2   * initial or fixed omega, for codons or codon-transltd AAs

#         fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
#             alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
#            Malpha = 0   * different alphas for genes
#             ncatG = 4   * # of categories in the dG or AdG models of rates

#             getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
#      RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
#            method = 0   * 0: simultaneous; 1: one branch at a time
#       fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
#         cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?


#     * Specifications for duplicating results for the small data set in table 1
#     * of Yang (1998 MBE 15:568-573).
#     * see the tree file lysozyme.trees for specification of node (branch) labels"""

#     c8a = """          seqfile = Reverse_Translation_Seq.txt-gb1.fst * fasta file
#          treefile = Species_Phylogenetic_tree_newick_nodistances.nwk *
#           outfile = codeml8a_mlc.txt

#             noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
#           verbose = 1   * 1: detailed output, 0: concise output
#           runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
#                         * 3: StepwiseAddition; (4,5):PerturbationNNI

#           seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
#         CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
#             clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
#             model = 0   * 0 - site models 2 - branch and branch-site models
#                         * models for codons:
#                             * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

#           NSsites = 8  * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
#             icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below

#         fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
#             kappa = 2   * initial or fixed kappa
#         fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate
#             omega = 1   * initial or fixed omega, for codons or codon-transltd AAs

#         fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
#             alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
#            Malpha = 0   * different alphas for genes
#             ncatG = 4   * # of categories in the dG or AdG models of rates

#             getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
#      RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
#            method = 0   * 0: simultaneous; 1: one branch at a time
#       fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
#         cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?


#     * Specifications for duplicating results for the small data set in table 1
#     * of Yang (1998 MBE 15:568-573).
#     * see the tree file lysozyme.trees for specification of node (branch) labels"""

#     c2a = """          seqfile = Reverse_Translation_Seq.txt-gb1.fst * fasta file
#          treefile = Species_Phylogenetic_tree_newick_Interst#.nwk *
#           outfile = codeml2a_mlc.txt

#             noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
#           verbose = 1   * 1: detailed output, 0: concise output
#           runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
#                         * 3: StepwiseAddition; (4,5):PerturbationNNI

#           seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
#         CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
#             clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
#             model = 2   * 0 - site models 2 - branch and branch-site models
#                         * models for codons:
#                             * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

#           NSsites = 2  * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
#             icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below

#         fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
#             kappa = 2   * initial or fixed kappa
#         fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
#             omega = 2   * initial or fixed omega, for codons or codon-transltd AAs

#         fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
#             alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
#            Malpha = 0   * different alphas for genes
#             ncatG = 4   * # of categories in the dG or AdG models of rates

#             getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
#      RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
#            method = 0   * 0: simultaneous; 1: one branch at a time
#       fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
#         cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?


#     * Specifications for duplicating results for the small data set in table 1
#     * of Yang (1998 MBE 15:568-573).
#     * see the tree file lysozyme.trees for specification of node (branch) labels"""

#     c2 = """          seqfile = Reverse_Translation_Seq.txt-gb1.fst * fasta file
#          treefile = Species_Phylogenetic_tree_newick_Interst#.nwk *
#           outfile = codeml2_mlc.txt

#             noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
#           verbose = 1   * 1: detailed output, 0: concise output
#           runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
#                         * 3: StepwiseAddition; (4,5):PerturbationNNI

#           seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
#         CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
#             clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
#             model = 2  * 0 - site models 2 - branch and branch-site models
#                         * models for codons:
#                             * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

#           NSsites = 2   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
#             icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below

#         fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
#             kappa = 2   * initial or fixed kappa
#         fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate
#             omega = 1   * initial or fixed omega, for codons or codon-transltd AAs

#         fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
#             alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
#            Malpha = 0   * different alphas for genes
#             ncatG = 4   * # of categories in the dG or AdG models of rates

#             getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
#      RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
#            method = 0   * 0: simultaneous; 1: one branch at a time
#       fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
#         cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?


#     * Specifications for duplicating results for the small data set in table 1
#     * of Yang (1998 MBE 15:568-573).
#     * see the tree file lysozyme.trees for specification of node (branch) labels"""
#     codeml078 = open("codeml078.ctl" , "w")
#     codeml078.write(c078)
#     codeml8a = open("codeml8a.ctl" , "w")
#     codeml8a.write(c8a)
#     codeml2a = open("codeml2a.ctl","w")
#     codeml2a.write(c2a)
#     codeml2 = open("codeml2.ctl", "w")
#     codeml2.write(c2)


# def codeml_output(state, protein):
#     """
#         Function creates a spreadsheet for the p-values

#     Args:
#         state (int 0 or 1): state shows whether model2 and model 2a is required
#         protein (string):  the name of the protein
#     """
#     def parser(file):
#         lis = []
#         lines = file.readlines()
#         for i in lines:
#             if "lnL" in i:
#                 start = i.find("):")
#                 end = i.find("      ")
#                 lis.append(float(i[start+3:end].strip()))
#         return lis

#     reader078 = open("codeml078/codeml078_mlc.txt","r")
#     models078 = parser(reader078)

#     reader8a = open("codeml8a/codeml8a_mlc.txt","r")
#     model8a = parser(reader8a)

#     model2a = [""]
#     model2 = [""]
#     lrt22a = ''

#     if state == 1:
#         reader2a = open("codeml2a/codeml2a_mlc.txt","r")
#         model2a = parser(reader2a)

#         reader2 = open("codeml2/codeml2_mlc.txt","r")
#         model2 = parser(reader2)

#         lrt22a = 2*(model2a[0] - model2[0])

#     lrt78 = 2*(models078[2] - models078[1])
#     lrt88a = 2*(models078[2] - model8a[0])

#     chi78 =  1 - chi2.cdf(lrt78, 2)
#     chi88a =  1 - chi2.cdf(lrt88a, 1)
#     chi22a = 1 - chi2.cdf(lrt22a, 1) if state == 1 else ''
#     with open("Gene_Output.csv", "w") as newfile:
#         wr = csv.writer(newfile)
#         wr.writerow(["Name","Model0", "Model7", "Model8","Model8a", "Model2a", "Model2", "LRT(M7_vs_M8)", "LRT(M8a_vs_M8)", "LRT(M2a_vs_M2)", "p-v7vs8", "p-v8avs8", "p-v2avs2"])
#         wr.writerow([protein, models078[0], models078[1], models078[2], model8a[0], model2a[0], model2[0],lrt78, lrt88a, lrt22a, chi78, chi88a, chi22a])
#         return [protein, models078[0], models078[1], models078[2], model8a[0], model2a[0], model2[0],lrt78, lrt88a, lrt22a, chi78, chi88a, chi22a]
