# Rosetta Parameter Reference

Complete reference for all Rosetta parameters used in this protocol.

---

## Rosetta Relax Parameters

### Core Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-s` | string | - | Input PDB structure file |
| `-nstruct` | int | 1 | Number of output structures |
| `-out:prefix` | string | "" | Prefix for output files |
| `-out:pdb` | flag | false | Output in PDB format |

### Constraint Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-relax:constrain_relax_to_start_coords` | flag | false | Add coordinate constraints to starting positions |
| `-relax:coord_constrain_sidechains` | flag | false | Also constrain sidechain atoms |
| `-relax:ramp_constraints` | bool | true | Gradually reduce constraints during relax |

### Rotamer Sampling

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-ex1` | flag | false | Extra chi1 rotamer sampling |
| `-ex2` | flag | false | Extra chi2 rotamer sampling |
| `-ex3` | flag | false | Extra chi3 rotamer sampling |
| `-ex4` | flag | false | Extra chi4 rotamer sampling |
| `-use_input_sc` | flag | false | Include input rotamers in packing |
| `-flip_HNQ` | flag | false | Consider flipping His, Asn, Gln |

### Recommended Relax Command

```bash
relax.macosclangrelease \
  -s input.pdb \
  -relax:constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
  -relax:ramp_constraints false \
  -ex1 -ex2 \
  -use_input_sc \
  -nstruct 1 \
  -out:prefix relaxed_ \
  -out:pdb
```

---

## Cartesian ddG Parameters

### Core Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-s` | string | - | Input relaxed structure |
| `-ddg:mut_file` | string | - | Mutation specification file |
| `-ddg:iterations` | int | 3 | Number of independent runs |

### Cartesian Mode

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-ddg:cartesian` | flag | false | Use Cartesian-space optimization |
| `-fa_max_dis` | real | 9.0 | Maximum distance for fa terms |

### Score Function

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-score:weights` | string | ref2015 | Score function weights file |
| `-score:set_weights` | string | - | Override specific weights |

### Output Control

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-ddg:dump_pdbs` | bool | true | Output mutant structures |
| `-ddg:suppress_checkpointing` | bool | false | Disable checkpoint files |
| `-out:prefix` | string | "" | Prefix for output files |
| `-out:path:all` | string | "." | Output directory |

### Advanced Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-ddg:min` | bool | false | Use minimization protocol |
| `-ddg:mean` | bool | true | Report mean energies |
| `-ddg:legacy` | bool | true | Use legacy protocol |
| `-ddg:bbnbrs` | int | 1 | Backbone neighbors for flexibility |
| `-ddg:frag_nbrs` | int | 2 | Fragment neighbors for sampling |
| `-ddg:score_cutoff` | real | 1.0 | Energy cutoff for acceptance |

### Recommended ddG Command

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

---

## Score Function: REF2015_cart

### Energy Terms

| Term | Weight | Description |
|------|--------|-------------|
| `fa_atr` | 1.0 | Lennard-Jones attractive |
| `fa_rep` | 0.55 | Lennard-Jones repulsive |
| `fa_sol` | 1.0 | Lazaridis-Karplus solvation |
| `fa_intra_rep` | 0.005 | Intra-residue repulsion |
| `fa_intra_sol_xover4` | 1.0 | Intra-residue solvation |
| `lk_ball_wtd` | 1.0 | Weighted LK ball |
| `fa_elec` | 1.0 | Coulombic electrostatics |
| `hbond_sr_bb` | 1.0 | Short-range BB H-bonds |
| `hbond_lr_bb` | 1.0 | Long-range BB H-bonds |
| `hbond_bb_sc` | 1.0 | BB-SC H-bonds |
| `hbond_sc` | 1.0 | SC-SC H-bonds |
| `dslf_fa13` | 1.25 | Disulfide bonds |
| `omega` | 0.4 | Omega dihedral penalty |
| `fa_dun` | 0.7 | Dunbrack rotamer |
| `p_aa_pp` | 0.6 | Backbone geometry |
| `yhh_planarity` | 0.625 | Tyrosine planarity |
| `ref` | 1.0 | Reference energy |
| `rama_prepro` | 0.45 | Ramachandran |
| `cart_bonded` | 0.5 | Cartesian bonded terms |

### Total Energy Equation

```
E_total = fa_atr + fa_rep + fa_sol + fa_intra_rep + fa_intra_sol_xover4 +
          lk_ball_wtd + fa_elec + hbond_sr_bb + hbond_lr_bb +
          hbond_bb_sc + hbond_sc + dslf_fa13 + omega + fa_dun +
          p_aa_pp + yhh_planarity + ref + rama_prepro + cart_bonded
```

---

## Mutation File Format

### Syntax

```
total <N>
<count_for_mutation_1>
<wt_residue> <position> <mut_residue>
[additional mutations...]
```

### Single Mutation Example

```
total 1
1
K 135 R
```

### Multiple Independent Mutations

```
total 3
1
K 135 R
1
D 175 G
1
R 138 C
```

### Residue Codes

| Code | Residue |
|------|---------|
| A | Alanine |
| C | Cysteine |
| D | Aspartate |
| E | Glutamate |
| F | Phenylalanine |
| G | Glycine |
| H | Histidine |
| I | Isoleucine |
| K | Lysine |
| L | Leucine |
| M | Methionine |
| N | Asparagine |
| P | Proline |
| Q | Glutamine |
| R | Arginine |
| S | Serine |
| T | Threonine |
| V | Valine |
| W | Tryptophan |
| Y | Tyrosine |

---

## Environment Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `ROSETTA_BIN` | Path to Rosetta binaries | `/path/to/rosetta/bin` |
| `ROSETTA_DATABASE` | Path to Rosetta database | `/path/to/rosetta/database` |

---

## Platform-Specific Binaries

| Platform | Binary Suffix |
|----------|---------------|
| Linux (GCC) | `.linuxgccrelease` |
| Linux (Clang) | `.linuxclangrelease` |
| macOS (Clang) | `.macosclangrelease` |
| Windows (Cygwin) | `.cygwinrelease` |
