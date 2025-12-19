# Troubleshooting Guide

Common issues and solutions when running Rosetta cartesian_ddG.

---

## Issue 1: Residue Mismatch Error

### Error Message
```
ERROR: Wildtype residue at position 135 is ILE, but mutfile specifies LYS
```

### Cause
The mutation file specifies a wildtype residue that doesn't match the actual residue in the PDB.

### Solution
```bash
# Check actual residue in PDB at position 135
grep "^ATOM.*CA.*135 " relaxed_structure.pdb

# Update mutfile to match actual residue
# If position 135 is ILE, not LYS:
echo -e "total 1\n1\nI 135 R" > mutfile
```

---

## Issue 2: Structure Not Relaxed

### Symptom
Very high or unstable wildtype energies (>1000 REU).

### Cause
Input structure has steric clashes or is incompatible with Rosetta scoring.

### Solution
Always run Rosetta relax before cartesian_ddG:

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

---

## Issue 3: WT→WT is Undefined

### Situation
Attempting to calculate ΔΔG for wildtype (no mutation).

### Explanation
Rosetta cartesian_ddG requires an actual mutation event. WT→WT is conceptually undefined because there is no change to measure.

### Solution
Skip wildtype systems in Rosetta analysis. This is expected behavior, not an error.

For analysis purposes, WT ΔΔG should be defined as 0.0 by definition.

---

## Issue 4: Missing Atoms

### Error Message
```
ERROR: Atom CB not found in residue 200
```

### Cause
- Incomplete residues in structure
- Glycine residue (has no CB atom)
- Non-standard residues

### Solutions

1. **For Glycine:** This is normal; Gly has no CB. Check if error is about a different atom.

2. **For incomplete residues:**
   ```bash
   # Check for missing atoms
   grep "^ATOM.*200 " input.pdb | wc -l
   # Should be ~8 for most residues, 4 for Gly
   ```

3. **Rebuild missing atoms:** Use structure repair tools or obtain a more complete structure.

---

## Issue 5: Chain Break Warning

### Warning
```
WARNING: Chain break detected between residues 145 and 150
```

### Cause
Missing residues creating a gap in the chain.

### Solutions

1. **Model missing loops:** Use Rosetta loop modeling
2. **Use complete structure:** Obtain structure without gaps
3. **Avoid mutations near breaks:** Results may be unreliable in gap regions

### Check for gaps:
```python
from Bio.PDB import PDBParser

parser = PDBParser()
structure = parser.get_structure('prot', 'input.pdb')

for chain in structure.get_chains():
    residues = list(chain.get_residues())
    for i in range(len(residues)-1):
        curr = residues[i].id[1]
        next_res = residues[i+1].id[1]
        if next_res - curr > 1:
            print(f"Gap between {curr} and {next_res}")
```

---

## Issue 6: Rosetta Binary Not Found

### Error Message
```
-bash: cartesian_ddg.macosclangrelease: command not found
```

### Solution
```bash
# Set Rosetta path
export ROSETTA_BIN=/path/to/rosetta/source/bin

# Use full path
$ROSETTA_BIN/cartesian_ddg.macosclangrelease ...

# Or add to PATH
export PATH=$ROSETTA_BIN:$PATH
```

---

## Issue 7: Database Not Found

### Error Message
```
ERROR: Unable to open database file
```

### Solution
```bash
# Set database path
export ROSETTA_DATABASE=/path/to/rosetta/database

# Or specify in command
cartesian_ddg.macosclangrelease \
  -database /path/to/rosetta/database \
  ...
```

---

## Issue 8: Memory Error

### Error Message
```
std::bad_alloc
```

### Cause
Insufficient RAM for large proteins.

### Solutions

1. **Reduce iterations:**
   ```bash
   -ddg:iterations 5  # Instead of 10 or 20
   ```

2. **Process sequentially:** Don't run multiple jobs simultaneously

3. **Increase swap space:** System administration required

---

## Issue 9: Very Slow Execution

### Symptom
Calculations taking much longer than expected.

### Solutions

1. **Check protein size:** Large proteins (>1000 residues) are slow

2. **Reduce iterations for testing:**
   ```bash
   -ddg:iterations 3
   ```

3. **Disable PDB output:**
   ```bash
   -ddg:dump_pdbs false
   ```

4. **Check disk space:** Full disk slows I/O significantly

---

## Issue 10: Inconsistent Results

### Symptom
Large variation in ΔΔG across iterations.

### Possible Causes

1. **Structure quality issues**
2. **Mutation near flexible region**
3. **Insufficient sampling**

### Solutions

1. **Increase iterations:**
   ```bash
   -ddg:iterations 20
   ```

2. **Use multiple seeds:** Run with different starting structures

3. **Check structure quality:** Verify no clashes or unusual conformations

---

## Debugging Checklist

| Step | Check | Command |
|------|-------|---------|
| 1 | Structure format | `head -20 input.pdb` |
| 2 | Residue at position | `grep "^ATOM.*CA.*{POS} " input.pdb` |
| 3 | Mutfile format | `cat mutfile` |
| 4 | Rosetta database | `echo $ROSETTA_DATABASE` |
| 5 | Binary exists | `ls -la $ROSETTA_BIN/cartesian_ddg.*` |
| 6 | Log errors | `tail -50 rosetta.log` |
| 7 | Disk space | `df -h .` |
| 8 | Memory usage | `free -h` (Linux) or `vm_stat` (macOS) |

---

## Getting Help

1. **Rosetta Forums:** https://www.rosettacommons.org/forum
2. **GitHub Issues:** Open an issue in this repository
3. **Documentation:** https://www.rosettacommons.org/docs/latest/

When reporting issues, include:
- Rosetta version
- Operating system
- Complete error message
- Input files (structure, mutfile)
- Command used
