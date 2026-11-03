# Build Helix from Template Torsions

This directory contains tools to transfer helical torsion patterns from a small oligomer (e.g., tetramer) to larger oligomers (hexamer, octamer, decamer).

## Stand-alone script

`transfer_helix_torsions.py` transfers dihedral angles from a template GRO file to a target GRO file.

### Usage

```bash
conda activate terphenyl
python transfer_helix_torsions.py \
  --template /path/to/helix.gro \
  --target /path/to/target.gro \
  --residue-name MOP \
  --output /path/to/output_helix.gro \
  --template-residue-index 2
```

* `--template` : GRO file with the helical reference (e.g., tetramer helix)
* `--target` : GRO file of the larger oligomer to be remodelled
* `--residue-name` : Residue name to match (e.g., MOP)
* `--output` : Output GRO file with torsions transferred
* `--template-residue-index` : 0-based index of template residue used to define dihedral patterns. Default 0.

The script:
1. Loads both structures with MDAnalysis (topology guessed from GRO)
2. Collects all dihedrals and bond angles fully contained in the reference template residue
3. Uses `terphenyl_simulations.utils.get_torsion_ids` / `get_angle_ids` to map atom names by residue index
4. Copies dihedral angles and bond angles residue-by-residue with periodic repetition
5. Rotates downstream atoms for dihedrals about the B-C bond and for angles about an axis through B perpendicular to the A-B-C plane

## Fixes applied to the library

### `terphenyl_simulations/terphenyl_simulations/edit_conf.py`
- Rewrote `InternalCoordinateEditor.set_torsions_from_template` to correctly:
  - Identify dihedral patterns by filtering universe dihedrals by `resid` instead of `Residue.dihedrals`
  - Map atom names per residue using `get_torsion_ids` without character iteration
  - Select atoms by `resid and name` to avoid cross-residue name collisions
  - Compute one angle per dihedral, not per atom
  - Apply periodic repetition `i % n_template`

### `terphenyl_simulations/terphenyl_simulations/analysis_workflows/md.py`
- Fixed `set_helix_torsions` file naming: now writes `helix_<name>.gro` and centers that file consistently
- Removed unused `em_u` load
- Output paths now match `job.doc["foldamer_gro"]`

These changes allow the existing Signac/Flow workflow to correctly build helical conformations for hexamer/octamer/decamer from the tetramer template.
