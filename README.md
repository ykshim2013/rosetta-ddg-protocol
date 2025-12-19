# Rosetta Cartesian ddG Protocol for Variant Stability Analysis

**A comprehensive guide for computational protein stability predictions using Rosetta cartesian_ddG**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rosetta](https://img.shields.io/badge/Rosetta-2023.49+-blue.svg)](https://www.rosettacommons.org/)

---

## Table of Contents

1. [Overview](#overview)
2. [Software Requirements](#software-requirements)
3. [Structure Preparation](#structure-preparation)
4. [Relaxation Protocol](#relaxation-protocol)
5. [Cartesian ddG Calculation](#cartesian-ddg-calculation)
6. [Output Parsing and Analysis](#output-parsing-and-analysis)
7. [Statistical Analysis](#statistical-analysis)
8. [Interpretation Guidelines](#interpretation-guidelines)
9. [Common Issues and Troubleshooting](#common-issues-and-troubleshooting)
10. [Best Practices for Publication](#best-practices-for-publication)
11. [Complete Workflow Scripts](#complete-workflow-scripts)
12. [References](#references)

---

## Overview

### What is Cartesian ddG?

Rosetta `cartesian_ddG` calculates the change in protein stability (ΔΔG) upon mutation using Cartesian-space energy minimization. Unlike torsional-space methods, Cartesian ddG allows backbone flexibility during optimization, providing more accurate predictions for destabilizing mutations.

### Key Formula

```
ΔΔG = ΔG_mutant - ΔG_wildtype

Where:
  ΔΔG > 0: Destabilizing mutation (positive = less stable)
  ΔΔG < 0: Stabilizing mutation (negative = more stable)
  ΔΔG ≈ 0: Neutral effect
```

### Method Highlights

| Feature | Cartesian ddG | Notes |
|---------|--------------|-------|
| **Backbone sampling** | Yes (limited) | Local minimization around mutation |
| **Sidechain repacking** | Yes | Rotamer optimization |
| **Score function** | REF2015_cart | Cartesian-optimized weights |
| **Typical accuracy** | ~1 kcal/mol RMSE | For soluble proteins |
| **Runtime** | 5-30 min/mutation | Depends on iterations |

### Workflow Overview

```
┌─────────────────────────────────────────────────────────────────┐
│  1. Structure Preparation                                       │
│     - Clean PDB (remove waters, ligands, alternate conformers)  │
│     - Extract target chain(s)                                   │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  2. Rosetta Relax                                               │
│     - Constrained energy minimization                           │
│     - Optimize hydrogen networks                                │
│     - Remove steric clashes                                     │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  3. Cartesian ddG Calculation                                   │
│     - Generate mutation files                                   │
│     - Run multiple iterations (10-20 recommended)               │
│     - Calculate ΔΔG for each variant                            │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  4. Results Analysis                                            │
│     - Parse output files                                        │
│     - Statistical aggregation                                   │
│     - Classification and interpretation                         │
└─────────────────────────────────────────────────────────────────┘
```

---

## Software Requirements

### Rosetta Installation

**Version:** Rosetta 2023.49 or later (academic license required)

**Download:** https://www.rosettacommons.org/software/license-and-download

**Compilation:**
```bash
# Navigate to Rosetta source directory
cd rosetta/source

# Compile with release mode (recommended)
./scons.py -j8 mode=release bin

# Output binaries location:
# Linux:   rosetta/source/bin/relax.linuxgccrelease
#          rosetta/source/bin/cartesian_ddg.linuxgccrelease
# macOS:   rosetta/source/bin/relax.macosclangrelease
#          rosetta/source/bin/cartesian_ddg.macosclangrelease
```

### Required Applications

| Application | Purpose | Binary Name |
|-------------|---------|-------------|
| `relax` | Structure energy minimization | `relax.{platform}release` |
| `cartesian_ddg` | ΔΔG calculation | `cartesian_ddg.{platform}release` |

### Python Dependencies (for analysis)

```bash
pip install numpy pandas scipy matplotlib biopython
```

### System Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU | 4 cores | 8+ cores |
| RAM | 8 GB | 16+ GB |
| Storage | 10 GB | 50+ GB (for large studies) |

---

## Structure Preparation

### Input Structure Requirements

1. **Format:** PDB format (not mmCIF)
2. **Quality:** Resolution ≤ 3.5 Å for experimental structures; pLDDT > 70 for AlphaFold models
3. **Completeness:** No missing residues in regions of interest
4. **Chain ID:** Single chain (chain A recommended) for monomer analysis

### Structure Sources

| Source | Use Case | Quality Check |
|--------|----------|---------------|
| **PDB** | Experimental structures | Resolution, R-free |
| **AlphaFold DB** | No experimental structure | pLDDT scores > 70 |
| **AlphaFold3** | Complexes, custom sequences | Ranking confidence |
| **SWISS-MODEL** | Homology models | QMEAN scores |

### Cleaning PDB Files

```bash
# Remove waters, ligands, and extract ATOM records only
grep "^ATOM" input.pdb > cleaned.pdb

# For multi-chain complexes, extract specific chain
grep "^ATOM.*A " input.pdb > chain_A.pdb

# Verify structure
head -20 cleaned.pdb
```

### Using BioPython for Structure Cleaning

```python
from Bio.PDB import PDBParser, PDBIO, Select

class CleanSelect(Select):
    """Select only protein atoms from specified chain."""
    def __init__(self, chain_id='A'):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        # Exclude water (HOH) and heteroatoms
        return residue.id[0] == ' '

    def accept_atom(self, atom):
        # Exclude alternate conformations (keep only 'A' or ' ')
        return atom.get_altloc() in (' ', 'A')

# Parse and clean structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'input.pdb')

# Save cleaned structure
io = PDBIO()
io.set_structure(structure)
io.save('cleaned.pdb', CleanSelect('A'))

print("Structure cleaned and saved to cleaned.pdb")
```

---

## Relaxation Protocol

### Why Relax?

Rosetta relaxation optimizes the structure within Rosetta's energy function, removing steric clashes and optimizing hydrogen bonding networks. This is **essential** before ddG calculations to ensure:

1. Structure is compatible with Rosetta's scoring function
2. No artificial energy contributions from minor clashes
3. Consistent baseline for wildtype energy calculations

### Rosetta Relax Command

```bash
relax.macosclangrelease \
  -s cleaned.pdb \
  -relax:constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
  -relax:ramp_constraints false \
  -ex1 \
  -ex2 \
  -use_input_sc \
  -nstruct 1 \
  -out:prefix relaxed_ \
  -out:pdb
```

### Parameter Explanation

| Parameter | Description | Purpose |
|-----------|-------------|---------|
| `-s` | Input structure file | Specifies the PDB to relax |
| `-relax:constrain_relax_to_start_coords` | Restrain atoms to starting positions | Prevents large deviations from input |
| `-relax:coord_constrain_sidechains` | Also constrain sidechains | Maintains overall conformation |
| `-relax:ramp_constraints false` | Keep constraints constant | More conservative relaxation |
| `-ex1` | Extra chi1 rotamer sampling | Improved sidechain modeling |
| `-ex2` | Extra chi2 rotamer sampling | Improved sidechain modeling |
| `-use_input_sc` | Include input rotamers in sampling | Preserves native conformations |
| `-nstruct 1` | Number of output structures | Single structure for efficiency |
| `-out:prefix` | Output file prefix | Naming convention |

### Complete Relax Script

```bash
#!/bin/bash
# run_rosetta_relax.sh

# Configuration
ROSETTA_BIN="/path/to/rosetta/source/bin"
RELAX="${ROSETTA_BIN}/relax.macosclangrelease"  # or .linuxgccrelease
INPUT_PDB="cleaned.pdb"
OUTPUT_PREFIX="relaxed_"

# Verify executable exists
if [ ! -f "$RELAX" ]; then
    echo "ERROR: Rosetta relax not found at: $RELAX"
    echo "Update ROSETTA_BIN path to your Rosetta installation"
    exit 1
fi

# Verify input exists
if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: Input PDB not found: $INPUT_PDB"
    exit 1
fi

echo "Starting Rosetta Relax..."
echo "Input: $INPUT_PDB"
echo "Expected time: 30-60 minutes"

# Run Rosetta Relax
$RELAX \
  -s $INPUT_PDB \
  -relax:constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
  -relax:ramp_constraints false \
  -ex1 \
  -ex2 \
  -use_input_sc \
  -nstruct 1 \
  -out:prefix ${OUTPUT_PREFIX} \
  -out:pdb

if [ $? -eq 0 ]; then
    echo "Relax completed successfully!"
    echo "Output: ${OUTPUT_PREFIX}${INPUT_PDB%.pdb}_0001.pdb"
else
    echo "ERROR: Rosetta Relax failed"
    exit 1
fi
```

### Verifying Relaxation Quality

```bash
# Check output file was created
ls -la relaxed_*.pdb

# Compare energies (optional - check score in PDB header)
grep "^pose" relaxed_*.pdb
```

---

## Cartesian ddG Calculation

### Creating Mutation Files

The mutation file specifies which mutations to analyze.

**Format:**
```
total N
1
WT_residue position MUT_residue
```

**Single mutation example (K135R):**
```
total 1
1
K 135 R
```

**Creating mutation files programmatically:**
```python
def create_mutfile(wt_residue, position, mut_residue, output_file="mutfile"):
    """
    Create Rosetta mutation file.

    Args:
        wt_residue: Single-letter wildtype residue code
        position: Residue position number (PDB numbering)
        mut_residue: Single-letter mutant residue code
        output_file: Output filename
    """
    with open(output_file, 'w') as f:
        f.write("total 1\n")
        f.write("1\n")
        f.write(f"{wt_residue} {position} {mut_residue}\n")

    print(f"Created mutation file: {output_file}")
    print(f"Mutation: {wt_residue}{position}{mut_residue}")

# Example: Create K135R mutation file
create_mutfile('K', 135, 'R', 'K135R.mut')
```

**Batch creation for multiple variants:**
```python
import os

# List of variants (format: WT_residue, position, MUT_residue)
variants = [
    ('K', 135, 'R'),
    ('D', 175, 'G'),
    ('R', 138, 'C'),
    ('G', 421, 'R'),
    ('D', 910, 'G'),
]

# Create output directory
os.makedirs('mutation_files', exist_ok=True)

# Generate mutation files
for wt, pos, mut in variants:
    variant_name = f"{wt}{pos}{mut}"
    output_file = f"mutation_files/{variant_name}.mut"

    with open(output_file, 'w') as f:
        f.write("total 1\n")
        f.write("1\n")
        f.write(f"{wt} {pos} {mut}\n")

    print(f"Created: {output_file}")
```

### Cartesian ddG Command

```bash
cartesian_ddg.macosclangrelease \
  -s relaxed_structure.pdb \
  -ddg:mut_file mutfile \
  -ddg:iterations 10 \
  -fa_max_dis 9.0 \
  -ddg:cartesian \
  -score:weights ref2015_cart \
  -ddg:dump_pdbs false \
  -ddg:suppress_checkpointing true
```

### Parameter Reference

| Parameter | Description | Recommended Value |
|-----------|-------------|-------------------|
| `-s` | Input relaxed structure | Required |
| `-ddg:mut_file` | Mutation specification file | Required |
| `-ddg:iterations` | Independent runs per mutation | **10-20** |
| `-fa_max_dis` | Maximum distance for fa terms | 9.0 |
| `-ddg:cartesian` | Enable Cartesian-space optimization | **true** |
| `-score:weights` | Score function | **ref2015_cart** |
| `-ddg:dump_pdbs` | Output mutant structures | false |
| `-ddg:suppress_checkpointing` | Disable checkpoint files | true |

### Advanced Parameters (Optional)

For enhanced accuracy, especially with charged residues:

```bash
# Enhanced electrostatics (recommended for charged mutations)
-score:set_weights fa_elec 1.4

# Beta energy function modifications
-beta_nov16
```

### Rosetta Energy Function (REF2015_cart)

The total energy is computed as:

```
E_total = fa_atr + fa_rep + fa_sol + fa_intra_rep + fa_intra_sol_xover4 +
          lk_ball_wtd + fa_elec + hbond_sr_bb + hbond_lr_bb +
          hbond_bb_sc + hbond_sc + dslf_fa13 + omega + fa_dun +
          p_aa_pp + yhh_planarity + ref + rama_prepro + cart_bonded
```

**Key energy terms:**

| Term | Description |
|------|-------------|
| `fa_atr` | Lennard-Jones attractive interactions |
| `fa_rep` | Lennard-Jones repulsive interactions |
| `fa_sol` | Lazaridis-Karplus solvation energy |
| `fa_elec` | Coulombic electrostatic interactions |
| `hbond_*` | Hydrogen bonding terms (various geometries) |
| `fa_dun` | Dunbrack rotamer energy |
| `rama_prepro` | Ramachandran preferences |
| `cart_bonded` | Cartesian bonded geometry terms |
| `omega` | Omega dihedral penalty |
| `ref` | Reference energies per amino acid |

---

## Output Parsing and Analysis

### Output File Format

Cartesian ddG produces a `mutfile.ddg` (or `ddg_predictions.out`) file:

```
COMPLEX: Round1: WT: 599.176 fa_atr: -6768.495 fa_rep: 1057.577 ...
COMPLEX: Round2: WT: 595.225 fa_atr: -6767.983 fa_rep: 1055.149 ...
...
COMPLEX: Round10: WT: 595.225 fa_atr: -6767.983 fa_rep: 1055.149 ...
COMPLEX: Round1: MUT_175GLY: 598.739 fa_atr: -6760.346 fa_rep: 1053.354 ...
COMPLEX: Round2: MUT_175GLY: 598.692 fa_atr: -6760.365 fa_rep: 1053.396 ...
...
COMPLEX: Round10: MUT_175GLY: 598.692 fa_atr: -6760.365 fa_rep: 1053.396 ...
```

### Python Parsing Script

```python
import re
import numpy as np
from pathlib import Path

def parse_rosetta_ddg(ddg_file):
    """
    Parse Rosetta cartesian_ddg output file.

    Args:
        ddg_file: Path to mutfile.ddg or ddg_predictions.out

    Returns:
        dict: Statistics for WT and mutant energies, plus ΔΔG
    """
    wt_energies = []
    mut_energies = []
    mut_name = None

    with open(ddg_file, 'r') as f:
        for line in f:
            # Skip non-data lines
            if not line.startswith('COMPLEX:'):
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            # Extract total energy (4th field after splitting)
            try:
                energy = float(parts[3])
            except (ValueError, IndexError):
                continue

            # Classify as WT or mutant
            if 'WT:' in line or 'WT_:' in line:
                wt_energies.append(energy)
            elif 'MUT_' in line or 'MUT:' in line:
                mut_energies.append(energy)
                # Extract mutation name
                if mut_name is None:
                    match = re.search(r'MUT_?(\d*[A-Z]{3}|\w+):', line)
                    if match:
                        mut_name = match.group(1)

    if not wt_energies or not mut_energies:
        return None

    # Calculate ΔΔG for each iteration (paired comparison)
    n_pairs = min(len(wt_energies), len(mut_energies))
    ddg_values = np.array(mut_energies[:n_pairs]) - np.array(wt_energies[:n_pairs])

    return {
        'mutation': mut_name,
        'n_iterations': n_pairs,
        'wt_mean': np.mean(wt_energies[:n_pairs]),
        'wt_sd': np.std(wt_energies[:n_pairs], ddof=1),
        'mut_mean': np.mean(mut_energies[:n_pairs]),
        'mut_sd': np.std(mut_energies[:n_pairs], ddof=1),
        'ddg_mean': np.mean(ddg_values),
        'ddg_sd': np.std(ddg_values, ddof=1),
        'ddg_median': np.median(ddg_values),
        'ddg_iqr': np.percentile(ddg_values, 75) - np.percentile(ddg_values, 25),
        'ddg_min': np.min(ddg_values),
        'ddg_max': np.max(ddg_values),
        'ddg_values': ddg_values.tolist()
    }

# Example usage
if __name__ == "__main__":
    results = parse_rosetta_ddg('mutfile.ddg')

    if results:
        print(f"Mutation: {results['mutation']}")
        print(f"Iterations: {results['n_iterations']}")
        print(f"WT Energy:  {results['wt_mean']:.2f} ± {results['wt_sd']:.2f} REU")
        print(f"MUT Energy: {results['mut_mean']:.2f} ± {results['mut_sd']:.2f} REU")
        print(f"ΔΔG Mean:   {results['ddg_mean']:.2f} ± {results['ddg_sd']:.2f} kcal/mol")
        print(f"ΔΔG Median: {results['ddg_median']:.2f} kcal/mol")
    else:
        print("Failed to parse results")
```

### Batch Processing Multiple Variants

```python
import pandas as pd
from pathlib import Path

def process_all_variants(base_dir, pattern="**/mutfile.ddg"):
    """
    Process all cartesian_ddg outputs in a directory tree.

    Args:
        base_dir: Root directory containing variant subdirectories
        pattern: Glob pattern for finding ddg output files

    Returns:
        DataFrame with all results
    """
    results = []

    for ddg_file in Path(base_dir).glob(pattern):
        parsed = parse_rosetta_ddg(ddg_file)

        if parsed:
            # Extract variant name from directory structure
            variant_name = ddg_file.parent.name
            parsed['variant'] = variant_name
            parsed['file_path'] = str(ddg_file)
            results.append(parsed)
        else:
            print(f"Warning: Could not parse {ddg_file}")

    df = pd.DataFrame(results)

    # Reorder columns
    cols = ['variant', 'mutation', 'n_iterations', 'ddg_mean', 'ddg_sd',
            'ddg_median', 'wt_mean', 'wt_sd', 'mut_mean', 'mut_sd']
    df = df[[c for c in cols if c in df.columns] +
            [c for c in df.columns if c not in cols]]

    return df

# Example usage
df = process_all_variants('Raw_Data/')
df.to_csv('rosetta_results_summary.csv', index=False)

print(f"Processed {len(df)} variants")
print(df[['variant', 'ddg_mean', 'ddg_sd']].to_string())
```

---

## Statistical Analysis

### Aggregating Across Multiple Seeds/Replicates

For robust predictions with conformational ensembles:

```python
import pandas as pd
import numpy as np

def aggregate_by_variant(df, groupby_col='variant'):
    """
    Aggregate ΔΔG values across multiple seeds/replicates.

    Args:
        df: DataFrame with per-seed/replicate results
        groupby_col: Column to group by

    Returns:
        DataFrame with aggregated statistics
    """
    aggregated = df.groupby(groupby_col).agg({
        'ddg_mean': ['mean', 'std', 'count'],
        'ddg_sd': 'mean'  # Average within-seed SD
    }).reset_index()

    # Flatten column names
    aggregated.columns = [
        groupby_col,
        'ddg_aggregate_mean',
        'ddg_between_seed_sd',
        'n_seeds',
        'ddg_within_seed_sd'
    ]

    # Calculate combined uncertainty (propagation of errors)
    aggregated['ddg_total_uncertainty'] = np.sqrt(
        aggregated['ddg_between_seed_sd']**2 +
        aggregated['ddg_within_seed_sd']**2
    )

    return aggregated

# Example
agg_results = aggregate_by_variant(df)
print(agg_results)
```

### Statistical Comparison: Pathogenic vs Benign

```python
from scipy import stats
from sklearn.metrics import roc_curve, auc
import numpy as np

def compare_variant_groups(pathogenic_ddg, benign_ddg):
    """
    Compare ΔΔG distributions between pathogenic and benign variants.

    Args:
        pathogenic_ddg: Array of ΔΔG values for pathogenic variants
        benign_ddg: Array of ΔΔG values for benign variants

    Returns:
        dict: Statistical test results
    """
    pathogenic_ddg = np.array(pathogenic_ddg)
    benign_ddg = np.array(benign_ddg)

    # Welch's t-test (unequal variances)
    t_stat, t_pvalue = stats.ttest_ind(pathogenic_ddg, benign_ddg, equal_var=False)

    # Mann-Whitney U test (non-parametric)
    u_stat, u_pvalue = stats.mannwhitneyu(
        pathogenic_ddg, benign_ddg, alternative='two-sided'
    )

    # Effect size (Cohen's d)
    n1, n2 = len(pathogenic_ddg), len(benign_ddg)
    var1, var2 = np.var(pathogenic_ddg, ddof=1), np.var(benign_ddg, ddof=1)
    pooled_sd = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    cohens_d = (np.mean(pathogenic_ddg) - np.mean(benign_ddg)) / pooled_sd

    # ROC-AUC analysis
    y_true = np.concatenate([np.ones(n1), np.zeros(n2)])
    y_score = np.concatenate([pathogenic_ddg, benign_ddg])
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    # Find optimal threshold (Youden's J statistic)
    j_scores = tpr - fpr
    optimal_idx = np.argmax(j_scores)
    optimal_threshold = thresholds[optimal_idx]

    return {
        'pathogenic_n': n1,
        'pathogenic_mean': np.mean(pathogenic_ddg),
        'pathogenic_sd': np.std(pathogenic_ddg, ddof=1),
        'benign_n': n2,
        'benign_mean': np.mean(benign_ddg),
        'benign_sd': np.std(benign_ddg, ddof=1),
        't_statistic': t_stat,
        't_pvalue': t_pvalue,
        'mann_whitney_u': u_stat,
        'mann_whitney_pvalue': u_pvalue,
        'cohens_d': cohens_d,
        'roc_auc': roc_auc,
        'optimal_threshold': optimal_threshold,
        'fpr': fpr,
        'tpr': tpr
    }

# Example usage
results = compare_variant_groups(pathogenic_ddg, benign_ddg)

print("=" * 60)
print("Statistical Comparison: Pathogenic vs Benign")
print("=" * 60)
print(f"\nPathogenic (n={results['pathogenic_n']}):")
print(f"  Mean ΔΔG: {results['pathogenic_mean']:.2f} ± {results['pathogenic_sd']:.2f} kcal/mol")
print(f"\nBenign (n={results['benign_n']}):")
print(f"  Mean ΔΔG: {results['benign_mean']:.2f} ± {results['benign_sd']:.2f} kcal/mol")
print(f"\nStatistical Tests:")
print(f"  Welch's t-test: t = {results['t_statistic']:.3f}, p = {results['t_pvalue']:.2e}")
print(f"  Mann-Whitney U: U = {results['mann_whitney_u']:.0f}, p = {results['mann_whitney_pvalue']:.2e}")
print(f"\nEffect Size:")
print(f"  Cohen's d: {results['cohens_d']:.3f}")
print(f"\nClassification Performance:")
print(f"  ROC-AUC: {results['roc_auc']:.3f}")
print(f"  Optimal threshold: {results['optimal_threshold']:.2f} kcal/mol")
```

---

## Interpretation Guidelines

### ΔΔG Thresholds

| ΔΔG Range (kcal/mol) | Interpretation | Clinical Implication |
|----------------------|----------------|---------------------|
| **> 2.0** | Highly destabilizing | Likely pathogenic |
| **1.0 to 2.0** | Moderately destabilizing | Possibly pathogenic |
| **-1.0 to 1.0** | Neutral | Variant of uncertain significance |
| **-2.0 to -1.0** | Moderately stabilizing | Possibly benign |
| **< -2.0** | Highly stabilizing | Likely benign |

### Important Caveats

1. **Soluble proteins:** Rosetta cartesian_ddG is best validated for globular, soluble proteins

2. **Membrane proteins:** Standard cartesian_ddG does not explicitly model membrane environments
   - Consider using RosettaMP for transmembrane proteins
   - Or acknowledge this limitation in publications

3. **Surface mutations:** Predictions may be less accurate for highly solvent-exposed residues

4. **Protein-specific calibration:** Thresholds should be calibrated when possible against experimentally characterized variants for the specific protein family

### Quality Control Checks

```python
def quality_control(results_df, variance_threshold=5.0):
    """
    Flag potentially unreliable predictions.

    Args:
        results_df: DataFrame with ddg_mean, ddg_sd columns
        variance_threshold: Max acceptable SD (kcal/mol)

    Returns:
        DataFrame with QC flags
    """
    df = results_df.copy()

    # High variance flag
    df['high_variance'] = df['ddg_sd'] > variance_threshold

    # Extreme values (may indicate structural problems)
    df['extreme_ddg'] = df['ddg_mean'].abs() > 20

    # Insufficient iterations
    if 'n_iterations' in df.columns:
        df['low_iterations'] = df['n_iterations'] < 5

    # Overall QC status
    qc_columns = ['high_variance', 'extreme_ddg']
    if 'low_iterations' in df.columns:
        qc_columns.append('low_iterations')

    df['qc_passed'] = ~df[qc_columns].any(axis=1)

    # Summary
    n_passed = df['qc_passed'].sum()
    n_total = len(df)
    print(f"QC Summary: {n_passed}/{n_total} variants passed ({100*n_passed/n_total:.1f}%)")

    if df['high_variance'].any():
        print(f"  - High variance: {df['high_variance'].sum()} variants")
    if df['extreme_ddg'].any():
        print(f"  - Extreme ΔΔG: {df['extreme_ddg'].sum()} variants")

    return df
```

---

## Common Issues and Troubleshooting

### Issue 1: Residue Mismatch Error

**Error message:**
```
ERROR: Wildtype residue at position 135 is ILE, but mutfile specifies LYS
```

**Cause:** The mutation file specifies a wildtype residue that doesn't match the actual residue in the PDB.

**Solution:**
```bash
# Check actual residue in PDB at position 135
grep "^ATOM.*CA.*135 " relaxed_structure.pdb

# Update mutfile to match actual residue
# If position 135 is ILE, not LYS:
echo "total 1
1
I 135 R" > mutfile
```

### Issue 2: Structure Not Relaxed

**Symptom:** Very high or unstable wildtype energies (>1000 REU).

**Cause:** Input structure has clashes or is incompatible with Rosetta scoring.

**Solution:** Always run Rosetta relax before cartesian_ddG:
```bash
relax.macosclangrelease \
  -s input.pdb \
  -relax:constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
  -ex1 -ex2 \
  -use_input_sc \
  -nstruct 1 \
  -out:prefix relaxed_
```

### Issue 3: WT→WT is Undefined

**Situation:** Attempting to calculate ΔΔG for wildtype (no mutation).

**Explanation:** Rosetta cartesian_ddG requires an actual mutation event. WT→WT is conceptually undefined because there is no change to measure.

**Solution:** Skip wildtype systems in Rosetta analysis. This is expected behavior, not an error.

### Issue 4: Missing Atoms

**Error message:**
```
ERROR: Atom CB not found in residue 200
```

**Cause:** Incomplete residues or special residues (Glycine has no CB).

**Solutions:**
1. For Glycine: Mutations at Gly positions are valid; error may indicate other issues
2. For incomplete residues: Use structure repair tools or obtain a more complete structure
3. Check for non-standard residues that may need conversion

### Issue 5: Chain Break Warning

**Warning:**
```
WARNING: Chain break detected between residues 145 and 150
```

**Cause:** Missing residues creating a gap in the chain.

**Solutions:**
1. Model missing loops using Rosetta loop modeling
2. Use a more complete structure
3. Avoid mutations near chain breaks (results may be unreliable)

### Debugging Checklist

| Step | Check | Command |
|------|-------|---------|
| 1 | Structure format | `head -20 input.pdb` |
| 2 | Residue at position | `grep "^ATOM.*CA.*{POS} " input.pdb` |
| 3 | Mutfile format | `cat mutfile` |
| 4 | Rosetta database | `echo $ROSETTA_DATABASE` |
| 5 | Binary exists | `ls -la $(which cartesian_ddg.*)` |
| 6 | Log errors | `tail -50 rosetta.log` |

---

## Best Practices for Publication

### Methods Section Template

> **Protein Stability Prediction Using Rosetta**
>
> Stability effects of missense variants were calculated using Rosetta cartesian_ddG (version 2023.49) with the REF2015_cart score function (Park et al., 2016). Input structures [specify source: PDB ID or AlphaFold prediction] were prepared by constrained energy minimization using Rosetta relax with coordinate constraints to preserve experimental/predicted geometry (-relax:constrain_relax_to_start_coords, -relax:coord_constrain_sidechains, -ex1, -ex2, -use_input_sc). For each variant, cartesian_ddG calculations were performed with N=[10/20] independent iterations using Cartesian-space optimization (-ddg:cartesian, -fa_max_dis 9.0, -score:weights ref2015_cart). ΔΔG values (kcal/mol) were calculated as the difference between mutant and wildtype energies, with positive values indicating destabilization. Results are reported as mean ± standard deviation across iterations. Interpretation thresholds: ΔΔG > 2.0 kcal/mol (destabilizing, likely pathogenic), -2.0 to 2.0 kcal/mol (neutral/uncertain), ΔΔG < -2.0 kcal/mol (stabilizing, likely benign).

### Reporting Checklist

- [ ] Rosetta version specified
- [ ] Score function identified (ref2015_cart)
- [ ] Number of iterations reported
- [ ] Relaxation protocol parameters listed
- [ ] Statistical measures defined (mean ± SD)
- [ ] Interpretation thresholds stated with citations
- [ ] Structure source and quality documented
- [ ] Limitations acknowledged (if applicable)

### Data Availability Statement

> Rosetta cartesian_ddG output files and analysis scripts are available at [repository URL]. Rosetta software is available from RosettaCommons (https://www.rosettacommons.org/) under academic license.

---

## Complete Workflow Scripts

### Master Workflow Script

```bash
#!/bin/bash
# rosetta_ddg_workflow.sh
# Complete workflow for Rosetta cartesian_ddG analysis

set -e  # Exit on error

# ============================================================
# CONFIGURATION - UPDATE THESE PATHS
# ============================================================
ROSETTA_BIN="/path/to/rosetta/source/bin"
PLATFORM="macosclangrelease"  # or "linuxgccrelease"

INPUT_PDB=$1
MUTFILE=$2
OUTPUT_DIR=${3:-"rosetta_output"}
N_ITERATIONS=${N_ITERATIONS:-10}

# Executables
RELAX="${ROSETTA_BIN}/relax.${PLATFORM}"
CARTESIAN_DDG="${ROSETTA_BIN}/cartesian_ddg.${PLATFORM}"

# ============================================================
# VALIDATION
# ============================================================
if [ -z "$INPUT_PDB" ] || [ -z "$MUTFILE" ]; then
    echo "Usage: $0 <input.pdb> <mutfile> [output_dir]"
    echo ""
    echo "Environment variables:"
    echo "  ROSETTA_BIN: Path to Rosetta bin directory"
    echo "  N_ITERATIONS: Number of ddG iterations (default: 10)"
    exit 1
fi

if [ ! -f "$RELAX" ]; then
    echo "ERROR: Rosetta relax not found: $RELAX"
    exit 1
fi

if [ ! -f "$CARTESIAN_DDG" ]; then
    echo "ERROR: Rosetta cartesian_ddg not found: $CARTESIAN_DDG"
    exit 1
fi

# ============================================================
# SETUP
# ============================================================
mkdir -p "$OUTPUT_DIR"
cp "$INPUT_PDB" "$OUTPUT_DIR/input.pdb"
cp "$MUTFILE" "$OUTPUT_DIR/mutfile"
cd "$OUTPUT_DIR"

echo "========================================================"
echo "ROSETTA CARTESIAN DDG WORKFLOW"
echo "========================================================"
echo "Input structure: $INPUT_PDB"
echo "Mutation file: $MUTFILE"
echo "Output directory: $OUTPUT_DIR"
echo "Iterations: $N_ITERATIONS"
echo "========================================================"

# ============================================================
# STEP 1: ROSETTA RELAX
# ============================================================
echo ""
echo "STEP 1: Running Rosetta Relax..."
echo "Expected time: 30-60 minutes"

$RELAX \
  -s input.pdb \
  -relax:constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
  -relax:ramp_constraints false \
  -ex1 \
  -ex2 \
  -use_input_sc \
  -nstruct 1 \
  -out:prefix relaxed_ \
  -out:pdb \
  > relax.log 2>&1

RELAXED_PDB="relaxed_input_0001.pdb"
if [ ! -f "$RELAXED_PDB" ]; then
    echo "ERROR: Relaxation failed. Check relax.log"
    exit 1
fi
echo "Relaxation complete: $RELAXED_PDB"

# ============================================================
# STEP 2: CARTESIAN DDG
# ============================================================
echo ""
echo "STEP 2: Running Cartesian ddG..."
echo "Iterations: $N_ITERATIONS"

$CARTESIAN_DDG \
  -s $RELAXED_PDB \
  -ddg:mut_file mutfile \
  -ddg:iterations $N_ITERATIONS \
  -fa_max_dis 9.0 \
  -ddg:cartesian \
  -score:weights ref2015_cart \
  -ddg:dump_pdbs false \
  -ddg:suppress_checkpointing true \
  > cartesian_ddg.log 2>&1

if [ ! -f "mutfile.ddg" ]; then
    echo "ERROR: Cartesian ddG failed. Check cartesian_ddg.log"
    exit 1
fi
echo "Cartesian ddG complete: mutfile.ddg"

# ============================================================
# STEP 3: PARSE RESULTS
# ============================================================
echo ""
echo "STEP 3: Parsing results..."

python3 << 'PYTHON_SCRIPT'
import re
import numpy as np

wt_energies = []
mut_energies = []
mut_name = None

with open('mutfile.ddg', 'r') as f:
    for line in f:
        if not line.startswith('COMPLEX:'):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            energy = float(parts[3])
        except ValueError:
            continue

        if 'WT:' in line or 'WT_:' in line:
            wt_energies.append(energy)
        elif 'MUT_' in line:
            mut_energies.append(energy)
            if mut_name is None:
                match = re.search(r'MUT_(\d*[A-Z]{3})', line)
                if match:
                    mut_name = match.group(1)

if wt_energies and mut_energies:
    n = min(len(wt_energies), len(mut_energies))
    ddg = np.array(mut_energies[:n]) - np.array(wt_energies[:n])

    # Write summary
    with open('ddg_summary.txt', 'w') as out:
        out.write("=" * 60 + "\n")
        out.write("Rosetta Cartesian_ddG Analysis Summary\n")
        out.write("=" * 60 + "\n\n")
        out.write(f"Mutation: {mut_name}\n")
        out.write(f"Iterations: {n}\n\n")
        out.write(f"Wildtype Energy:  {np.mean(wt_energies[:n]):.2f} ± {np.std(wt_energies[:n], ddof=1):.2f} REU\n")
        out.write(f"Mutant Energy:    {np.mean(mut_energies[:n]):.2f} ± {np.std(mut_energies[:n], ddof=1):.2f} REU\n\n")
        out.write(f"ΔΔG Mean ± SD:    {np.mean(ddg):.2f} ± {np.std(ddg, ddof=1):.2f} kcal/mol\n")
        out.write(f"ΔΔG Median:       {np.median(ddg):.2f} kcal/mol\n")
        out.write(f"ΔΔG Range:        [{np.min(ddg):.2f}, {np.max(ddg):.2f}] kcal/mol\n\n")

        # Interpretation
        mean_ddg = np.mean(ddg)
        if mean_ddg > 2.0:
            interp = "DESTABILIZING (likely pathogenic)"
        elif mean_ddg > 1.0:
            interp = "MILDLY DESTABILIZING"
        elif mean_ddg < -2.0:
            interp = "STABILIZING (likely benign)"
        elif mean_ddg < -1.0:
            interp = "MILDLY STABILIZING"
        else:
            interp = "NEUTRAL (VUS)"

        out.write(f"Interpretation: {interp}\n")

    print(f"Mutation: {mut_name}")
    print(f"ΔΔG = {np.mean(ddg):.2f} ± {np.std(ddg, ddof=1):.2f} kcal/mol")
    print(f"Interpretation: {interp}")
else:
    print("ERROR: Could not parse results from mutfile.ddg")
PYTHON_SCRIPT

echo ""
echo "========================================================"
echo "WORKFLOW COMPLETE"
echo "========================================================"
echo "Results: $OUTPUT_DIR/ddg_summary.txt"
echo "========================================================"
```

### Batch Processing Script

```bash
#!/bin/bash
# batch_rosetta_ddg.sh
# Process multiple variants in parallel

set -e

INPUT_PDB=$1
VARIANT_LIST=$2  # File with variants, one per line (e.g., K135R)
OUTPUT_BASE=${3:-"rosetta_results"}
N_PARALLEL=${N_PARALLEL:-4}
N_ITERATIONS=${N_ITERATIONS:-10}

if [ -z "$INPUT_PDB" ] || [ -z "$VARIANT_LIST" ]; then
    echo "Usage: $0 <input.pdb> <variant_list.txt> [output_base_dir]"
    echo ""
    echo "variant_list.txt format (one variant per line):"
    echo "  K135R"
    echo "  D175G"
    echo "  R138C"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_BASE"

# First, run relax once on the input structure
echo "Running Rosetta Relax on input structure..."
# [Insert relax command here - same as in master workflow]

RELAXED_PDB="${OUTPUT_BASE}/relaxed_structure.pdb"

# Generate mutation files and process
echo "Processing variants from: $VARIANT_LIST"

process_variant() {
    variant=$1
    relaxed=$2
    output_base=$3
    iterations=$4

    # Parse variant (e.g., K135R -> K, 135, R)
    wt_res=${variant:0:1}
    position=$(echo "$variant" | sed 's/^[A-Z]//' | sed 's/[A-Z]$//')
    mut_res=${variant: -1}

    # Create variant directory
    variant_dir="${output_base}/${variant}"
    mkdir -p "$variant_dir"

    # Create mutation file
    cat > "${variant_dir}/mutfile" << EOF
total 1
1
$wt_res $position $mut_res
EOF

    # Run cartesian_ddg
    cartesian_ddg.macosclangrelease \
      -s "$relaxed" \
      -ddg:mut_file "${variant_dir}/mutfile" \
      -ddg:iterations $iterations \
      -fa_max_dis 9.0 \
      -ddg:cartesian \
      -score:weights ref2015_cart \
      -ddg:dump_pdbs false \
      -ddg:suppress_checkpointing true \
      -out:path:all "$variant_dir" \
      > "${variant_dir}/rosetta.log" 2>&1

    echo "Completed: $variant"
}

export -f process_variant

# Run in parallel
cat "$VARIANT_LIST" | xargs -P $N_PARALLEL -I {} \
    bash -c "process_variant {} $RELAXED_PDB $OUTPUT_BASE $N_ITERATIONS"

echo ""
echo "Batch processing complete!"
echo "Results in: $OUTPUT_BASE/"
```

---

## References

### Primary Citations

1. **Rosetta REF2015 Score Function:**

   Park H, Bradley P, Greisen P Jr, et al. Simultaneous Optimization of Biomolecular Energy Functions on Features from Small Molecules and Macromolecules. *J Chem Theory Comput*. 2016;12(12):6201-6212. doi:10.1021/acs.jctc.6b00819

2. **Cartesian ddG Protocol:**

   Frenz B, Lewis SM, King I, et al. Prediction of Protein Mutational Free Energy: Benchmark and Sampling Improvements Increase Classification Accuracy. *Front Bioeng Biotechnol*. 2020;8:558247. doi:10.3389/fbioe.2020.558247

3. **RosettaDDGPrediction:**

   Sora V, Sanchez D, Papaleo E. RosettaDDGPrediction for high-throughput mutational scans: From stability to binding. *Protein Science*. 2023;32(1):e4527. doi:10.1002/pro.4527

4. **Rosetta Relax Protocol:**

   Conway P, Tyka MD, DiMaio F, Konerding DE, Baker D. Relaxation of backbone bond geometry improves protein energy landscape modeling. *Protein Science*. 2014;23(1):47-55. doi:10.1002/pro.2389

### Supporting Literature

5. **Benchmark Validation:**

   Kellogg EH, Leaver-Fay A, Baker D. Role of conformational sampling in computing mutation-induced changes in protein structure and stability. *Proteins*. 2011;79(3):830-838. doi:10.1002/prot.22921

6. **Clinical Variant Interpretation:**

   Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med*. 2015;17(5):405-424. doi:10.1038/gim.2015.30

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2025-12 | Initial release |

---

## License

This protocol documentation is released under the MIT License.

Rosetta software requires a separate academic or commercial license from RosettaCommons (https://www.rosettacommons.org/).

---

## Acknowledgments

This protocol was developed and validated through analysis of FBXO11 and SCN2A variant datasets, incorporating best practices from published literature and the Rosetta community.
