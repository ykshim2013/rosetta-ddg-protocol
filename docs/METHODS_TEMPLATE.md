# Publication Methods Template

Copy-paste ready methods sections for manuscripts.

---

## Short Version (100-150 words)

> Protein stability effects were predicted using Rosetta cartesian_ddG (version 2023.49) with the REF2015_cart score function (Park et al., 2016). Input structures were prepared by constrained energy minimization using Rosetta relax. For each variant, cartesian_ddG calculations were performed with 10 independent iterations using Cartesian-space optimization. ΔΔG values (kcal/mol) were calculated as the difference between mutant and wildtype energies, with positive values indicating destabilization. Results are reported as mean ± SD. Variants with ΔΔG > 2.0 kcal/mol were classified as destabilizing (likely pathogenic), while those with ΔΔG < -1.0 kcal/mol were classified as stabilizing (likely benign).

---

## Standard Version (200-300 words)

> **Protein Stability Prediction**
>
> Stability effects of missense variants were calculated using Rosetta cartesian_ddG (version 2023.49) with the REF2015_cart score function (Park et al., 2016). Input structures [specify source: PDB ID or AlphaFold prediction] were prepared by constrained energy minimization using Rosetta relax with the following parameters: coordinate constraints to preserve experimental/predicted geometry (-relax:constrain_relax_to_start_coords, -relax:coord_constrain_sidechains), constant constraint weights (-relax:ramp_constraints false), and enhanced rotamer sampling (-ex1, -ex2, -use_input_sc).
>
> For each variant, cartesian_ddG calculations were performed with N=10 independent iterations using Cartesian-space optimization (-ddg:cartesian) with the following parameters: maximum distance for full-atom terms of 9.0 Å (-fa_max_dis 9.0) and the REF2015_cart score function (-score:weights ref2015_cart). ΔΔG values (kcal/mol) were calculated as the difference between mutant and wildtype energies (ΔΔG = ΔG_mutant - ΔG_wildtype), with positive values indicating destabilization.
>
> Results are reported as mean ± standard deviation across iterations. Interpretation thresholds were applied as follows: ΔΔG > 2.0 kcal/mol (destabilizing, likely pathogenic), 1.0-2.0 kcal/mol (mildly destabilizing), -1.0 to 1.0 kcal/mol (neutral, variant of uncertain significance), and ΔΔG < -1.0 kcal/mol (stabilizing, likely benign).

---

## Extended Version (300-400 words)

> **Computational Protein Stability Analysis**
>
> Protein stability effects of missense variants were quantified using Rosetta cartesian_ddG (version 2023.49; RosettaCommons) with the REF2015_cart score function, which includes Cartesian-space bonded terms optimized for structure prediction and design applications (Park et al., 2016). This approach allows limited backbone flexibility during energy calculations, improving accuracy for destabilizing mutations compared to fixed-backbone methods.
>
> **Structure Preparation:** Input structures [describe source] were prepared for Rosetta calculations using constrained energy minimization (Rosetta relax application). Structures were minimized with coordinate constraints to the starting geometry (-relax:constrain_relax_to_start_coords) including sidechain atoms (-relax:coord_constrain_sidechains) with constant constraint weights (-relax:ramp_constraints false). Enhanced rotamer sampling was enabled (-ex1, -ex2) with inclusion of input rotamers (-use_input_sc) to preserve experimentally observed conformations where energetically favorable.
>
> **ΔΔG Calculations:** For each variant, the cartesian_ddG application was executed with N=10 independent iterations to ensure statistical robustness. Calculations used Cartesian-space optimization (-ddg:cartesian) with a maximum distance cutoff of 9.0 Å for full-atom energy terms (-fa_max_dis 9.0). The change in folding free energy (ΔΔG) was computed as the difference between mutant and wildtype total energies (ΔΔG = ΔG_mutant - ΔG_wildtype), where positive values indicate destabilization of the native fold.
>
> **Statistical Analysis:** Results are reported as mean ± standard deviation across independent iterations. For variants analyzed with multiple conformational seeds, results were aggregated as the mean across seeds with uncertainty propagated from both within-seed and between-seed variance. Pathogenic versus benign discrimination was assessed using Welch's t-test with effect size (Cohen's d) and receiver operating characteristic (ROC) analysis.
>
> **Interpretation:** Classification thresholds were applied based on established benchmarks: ΔΔG > 2.0 kcal/mol (destabilizing, likely pathogenic), 1.0-2.0 kcal/mol (mildly destabilizing, possibly pathogenic), -1.0 to 1.0 kcal/mol (neutral, variant of uncertain significance), and ΔΔG < -1.0 kcal/mol (stabilizing, likely benign). These thresholds correspond approximately to 10-fold changes in protein stability.

---

## Data Availability Statement

> Rosetta cartesian_ddG output files and analysis scripts are available at [repository URL]. The computational workflow is documented at https://github.com/ykshim2013/rosetta-ddg-protocol. Rosetta software is available from RosettaCommons (https://www.rosettacommons.org/) under academic license.

---

## Required Citations

### Primary (Required)

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
```

### Secondary (Recommended)

```bibtex
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

@article{conway2014relax,
  title={Relaxation of backbone bond geometry improves protein
         energy landscape modeling},
  author={Conway, Patrick and Tyka, Michael D and DiMaio, Frank and
          Konerding, David E and Baker, David},
  journal={Protein Science},
  volume={23},
  number={1},
  pages={47--55},
  year={2014},
  doi={10.1002/pro.2389}
}
```

---

## Supplementary Methods Template

> **Supplementary Methods: Rosetta Cartesian ddG Protocol**
>
> **Software:** Rosetta version 2023.49 (RosettaCommons, https://www.rosettacommons.org/)
>
> **Step 1 - Structure Relaxation:**
> ```
> relax.macosclangrelease \
>   -s input.pdb \
>   -relax:constrain_relax_to_start_coords \
>   -relax:coord_constrain_sidechains \
>   -relax:ramp_constraints false \
>   -ex1 -ex2 \
>   -use_input_sc \
>   -nstruct 1 \
>   -out:prefix relaxed_
> ```
>
> **Step 2 - ΔΔG Calculation:**
> ```
> cartesian_ddg.macosclangrelease \
>   -s relaxed_structure.pdb \
>   -ddg:mut_file mutfile \
>   -ddg:iterations 10 \
>   -fa_max_dis 9.0 \
>   -ddg:cartesian \
>   -score:weights ref2015_cart
> ```
>
> **Mutation file format:**
> ```
> total 1
> 1
> [WT_residue] [position] [MUT_residue]
> ```
>
> **Output parsing:** ΔΔG values were extracted from the mutfile.ddg output file as the difference between MUT and WT total energies for each iteration. Mean and standard deviation were calculated across iterations.
