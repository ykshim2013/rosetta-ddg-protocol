# Rosetta Cartesian ddG Protocol

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rosetta](https://img.shields.io/badge/Rosetta-2023.49+-blue.svg)](https://www.rosettacommons.org/)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.6b00819-green.svg)](https://doi.org/10.1021/acs.jctc.6b00819)

A comprehensive, publication-ready protocol for computational protein stability predictions using Rosetta cartesian_ddG.

## Overview

This repository provides a complete workflow for predicting the stability effects of missense variants using Rosetta's cartesian_ddG application. The protocol has been validated through analysis of FBXO11 and SCN2A variant datasets.

### Workflow

```
┌──────────────────────┐     ┌──────────────────────┐     ┌──────────────────────┐
│  1. Structure Prep   │ ──▶ │   2. Rosetta Relax   │ ──▶ │  3. Cartesian ddG    │
│  - Clean PDB         │     │  - Energy minimize   │     │  - Calculate ΔΔG     │
│  - Extract chain     │     │  - Constrained       │     │  - Multiple iterations│
└──────────────────────┘     └──────────────────────┘     └──────────────────────┘
```

### Key Features

- **Rosetta Relax** - Constrained energy minimization preserving input geometry
- **Cartesian ddG** - Backbone-flexible stability calculations
- **REF2015_cart** - State-of-the-art Rosetta score function
- **Statistical Analysis** - Python scripts for parsing and interpretation
- **Publication-ready** - Methods section templates and citations

## Quick Start

### Prerequisites

- Rosetta 2023.49+ ([Academic License](https://www.rosettacommons.org/software/license-and-download))
- Python 3.8+ with NumPy, Pandas, SciPy

### Installation

```bash
git clone https://github.com/ykshim2013/rosetta-ddg-protocol.git
cd rosetta-ddg-protocol

# Update Rosetta path in scripts
export ROSETTA_BIN=/path/to/rosetta/source/bin
```

### Basic Usage

```bash
# 1. Relax structure
./scripts/01_rosetta_relax.sh input.pdb

# 2. Create mutation file
echo -e "total 1\n1\nK 135 R" > mutfile

# 3. Run cartesian_ddG
./scripts/02_cartesian_ddg.sh relaxed_input.pdb mutfile

# 4. Parse results
python scripts/parse_ddg_results.py mutfile.ddg
```

## Repository Structure

```
rosetta-ddg-protocol/
├── README.md                 # This file
├── LICENSE                   # MIT License
├── CITATION.cff             # Citation information
│
├── docs/
│   ├── PROTOCOL.md          # Complete protocol documentation
│   ├── PARAMETERS.md        # Detailed parameter reference
│   ├── TROUBLESHOOTING.md   # Common issues and solutions
│   └── METHODS_TEMPLATE.md  # Publication methods section
│
├── scripts/
│   ├── 01_rosetta_relax.sh      # Structure relaxation
│   ├── 02_cartesian_ddg.sh      # ΔΔG calculation
│   ├── 03_batch_processing.sh   # Multiple variants
│   ├── parse_ddg_results.py     # Output parsing
│   └── statistical_analysis.py  # Statistical comparison
│
└── examples/
    ├── example_structure.pdb    # Example input
    ├── example_mutfile          # Example mutation file
    └── expected_output/         # Expected results
```

## Protocol Summary

### Step 1: Structure Preparation

Clean your PDB file to remove waters, ligands, and alternate conformations:

```bash
grep "^ATOM" input.pdb > cleaned.pdb
```

### Step 2: Rosetta Relax

Energy minimize with coordinate constraints:

```bash
relax.macosclangrelease \
  -s cleaned.pdb \
  -relax:constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
  -relax:ramp_constraints false \
  -ex1 -ex2 \
  -use_input_sc \
  -nstruct 1 \
  -out:prefix relaxed_
```

### Step 3: Cartesian ddG

Calculate stability change upon mutation:

```bash
cartesian_ddg.macosclangrelease \
  -s relaxed_structure.pdb \
  -ddg:mut_file mutfile \
  -ddg:iterations 10 \
  -fa_max_dis 9.0 \
  -ddg:cartesian \
  -score:weights ref2015_cart
```

### Step 4: Interpretation

| ΔΔG (kcal/mol) | Interpretation | Clinical Implication |
|----------------|----------------|---------------------|
| > 2.0 | Destabilizing | Likely pathogenic |
| 1.0 to 2.0 | Mildly destabilizing | Possibly pathogenic |
| -1.0 to 1.0 | Neutral | VUS |
| < -1.0 | Stabilizing | Likely benign |

## Documentation

| Document | Description |
|----------|-------------|
| [PROTOCOL.md](docs/PROTOCOL.md) | Complete step-by-step protocol |
| [PARAMETERS.md](docs/PARAMETERS.md) | All Rosetta parameters explained |
| [TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) | Common issues and solutions |
| [METHODS_TEMPLATE.md](docs/METHODS_TEMPLATE.md) | Copy-paste methods for papers |

## Citation

If you use this protocol, please cite:

```bibtex
@article{park2016ref2015,
  title={Simultaneous optimization of biomolecular energy functions on features
         from small molecules and macromolecules},
  author={Park, Hahnbeom and Bradley, Philip and Greisen Jr, Per and others},
  journal={Journal of Chemical Theory and Computation},
  volume={12},
  number={12},
  pages={6201--6212},
  year={2016},
  doi={10.1021/acs.jctc.6b00819}
}

@article{sora2023rosettaddg,
  title={RosettaDDGPrediction for high-throughput mutational scans:
         From stability to binding},
  author={Sora, Valentina and Sanchez, David and Papaleo, Elena},
  journal={Protein Science},
  volume={32},
  number={1},
  pages={e4527},
  year={2023},
  doi={10.1002/pro.4527}
}
```

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

Rosetta software requires a separate [academic or commercial license](https://www.rosettacommons.org/software/license-and-download).

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Acknowledgments

This protocol was developed and validated through analysis of FBXO11 and SCN2A variant datasets, incorporating best practices from published literature and the Rosetta community.
