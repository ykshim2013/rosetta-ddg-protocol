# Complete Rosetta Cartesian ddG Protocol

This document provides the complete, detailed protocol for running Rosetta cartesian_ddG calculations.

## Table of Contents

1. [Overview](#overview)
2. [Software Requirements](#software-requirements)
3. [Structure Preparation](#structure-preparation)
4. [Rosetta Relax](#rosetta-relax)
5. [Cartesian ddG Calculation](#cartesian-ddg-calculation)
6. [Output Analysis](#output-analysis)
7. [Statistical Methods](#statistical-methods)

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

### Workflow Diagram

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
cd rosetta/source
./scons.py -j8 mode=release bin

# Output binaries:
# Linux:   rosetta/source/bin/relax.linuxgccrelease
#          rosetta/source/bin/cartesian_ddg.linuxgccrelease
# macOS:   rosetta/source/bin/relax.macosclangrelease
#          rosetta/source/bin/cartesian_ddg.macosclangrelease
```

### Python Dependencies

```bash
pip install numpy pandas scipy matplotlib biopython
```

### System Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU | 4 cores | 8+ cores |
| RAM | 8 GB | 16+ GB |
| Storage | 10 GB | 50+ GB |

---

## Structure Preparation

### Input Requirements

1. **Format:** PDB format (not mmCIF)
2. **Quality:** Resolution ≤ 3.5 Å for experimental; pLDDT > 70 for AlphaFold
3. **Completeness:** No missing residues in regions of interest
4. **Chain ID:** Single chain preferred for monomer analysis

### Cleaning PDB Files

**Command line:**
```bash
# Remove waters, ligands, extract ATOM records
grep "^ATOM" input.pdb > cleaned.pdb

# Extract specific chain
grep "^ATOM.*A " input.pdb > chain_A.pdb
```

**Python (BioPython):**
```python
from Bio.PDB import PDBParser, PDBIO, Select

class CleanSelect(Select):
    def __init__(self, chain_id='A'):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == ' '  # Exclude heteroatoms

    def accept_atom(self, atom):
        return atom.get_altloc() in (' ', 'A')  # Exclude alt confs

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'input.pdb')

io = PDBIO()
io.set_structure(structure)
io.save('cleaned.pdb', CleanSelect('A'))
```

---

## Rosetta Relax

### Purpose

Rosetta relaxation optimizes the structure within Rosetta's energy function:
- Removes steric clashes
- Optimizes hydrogen bonding networks
- Ensures compatibility with Rosetta scoring

### Command

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

### Parameters Explained

| Parameter | Description |
|-----------|-------------|
| `-s` | Input structure |
| `-relax:constrain_relax_to_start_coords` | Restrain atoms to starting positions |
| `-relax:coord_constrain_sidechains` | Also constrain sidechains |
| `-relax:ramp_constraints false` | Keep constraints constant |
| `-ex1 -ex2` | Extra rotamer sampling |
| `-use_input_sc` | Include input rotamers |
| `-nstruct 1` | Single output structure |

### Expected Runtime

- Small proteins (<200 residues): 10-20 minutes
- Medium proteins (200-500 residues): 30-60 minutes
- Large proteins (>500 residues): 1-2 hours

---

## Cartesian ddG Calculation

### Creating Mutation Files

**Format:**
```
total N
1
WT_residue position MUT_residue
```

**Example (K135R):**
```
total 1
1
K 135 R
```

**Python script for batch creation:**
```python
def create_mutfile(wt_res, position, mut_res, output_file="mutfile"):
    with open(output_file, 'w') as f:
        f.write("total 1\n")
        f.write("1\n")
        f.write(f"{wt_res} {position} {mut_res}\n")
```

### Command

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

### Parameters Explained

| Parameter | Description | Recommended |
|-----------|-------------|-------------|
| `-s` | Relaxed input structure | Required |
| `-ddg:mut_file` | Mutation specification | Required |
| `-ddg:iterations` | Independent runs | 10-20 |
| `-fa_max_dis` | Max distance for fa terms | 9.0 |
| `-ddg:cartesian` | Cartesian-space optimization | true |
| `-score:weights` | Score function | ref2015_cart |

### Expected Runtime

- Per variant: 5-30 minutes (depends on protein size and iterations)
- 10 iterations on medium protein: ~15 minutes

---

## Output Analysis

### Output File Format

The `mutfile.ddg` file contains:

```
COMPLEX: Round1: WT: 599.176 fa_atr: -6768.495 fa_rep: 1057.577 ...
COMPLEX: Round2: WT: 595.225 fa_atr: -6767.983 fa_rep: 1055.149 ...
...
COMPLEX: Round1: MUT_175GLY: 598.739 fa_atr: -6760.346 ...
COMPLEX: Round2: MUT_175GLY: 598.692 fa_atr: -6760.365 ...
...
```

### Parsing Script

```python
import re
import numpy as np

def parse_rosetta_ddg(ddg_file):
    wt_energies = []
    mut_energies = []

    with open(ddg_file, 'r') as f:
        for line in f:
            if not line.startswith('COMPLEX:'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            energy = float(parts[3])

            if 'WT:' in line or 'WT_:' in line:
                wt_energies.append(energy)
            elif 'MUT_' in line:
                mut_energies.append(energy)

    n = min(len(wt_energies), len(mut_energies))
    ddg_values = np.array(mut_energies[:n]) - np.array(wt_energies[:n])

    return {
        'ddg_mean': np.mean(ddg_values),
        'ddg_sd': np.std(ddg_values, ddof=1),
        'ddg_median': np.median(ddg_values),
        'n_iterations': n
    }
```

---

## Statistical Methods

### Aggregating Multiple Seeds

```python
def aggregate_by_variant(df):
    return df.groupby('variant').agg({
        'ddg_mean': ['mean', 'std', 'count']
    }).reset_index()
```

### Comparing Groups

```python
from scipy import stats

def compare_groups(pathogenic, benign):
    t_stat, p_value = stats.ttest_ind(pathogenic, benign, equal_var=False)
    return {'t_statistic': t_stat, 'p_value': p_value}
```

### Interpretation Thresholds

| ΔΔG (kcal/mol) | Interpretation |
|----------------|----------------|
| > 2.0 | Destabilizing (likely pathogenic) |
| 1.0 to 2.0 | Mildly destabilizing |
| -1.0 to 1.0 | Neutral (VUS) |
| < -1.0 | Stabilizing (likely benign) |

---

## References

1. Park H, et al. (2016). J Chem Theory Comput. 12(12):6201-6212.
2. Sora V, et al. (2023). Protein Science. 32(1):e4527.
3. Conway P, et al. (2014). Protein Science. 23(1):47-55.
