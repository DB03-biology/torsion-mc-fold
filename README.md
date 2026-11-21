# Coarse-Grained Torsion Monte Carlo Protein Folding

A lightweight C implementation of a $C_\alpha$ coarse-grained protein folding simulation using Metropolis Monte Carlo, alongside Python scripts for dihedral extraction and energy landscape plotting.

## Folding Trajectory (AsLOV2)

Watch the folding simulation of the *Avena sativa* LOV2 domain (**AsLOV2**) starting from a fully extended linear chain:

[https://github.com/DB03-biology/torsion-mc-fold.git/raw/main/folding_movie.mp4](https://github.com/DB03-biology/torsion-mc-fold/issues/1#issue-6526666341)

*(If the video player does not load inline in your browser, you can view or download [`folding_movie.mp4`](folding_movie.mp4) directly from the repository.)*

## Overview

Starting from an extended polypeptide chain ($3.8\text{ \AA}$ spacing along the $x$-axis), the simulation samples conformational space through random single-residue Cartesian displacements ($T = 1.35$ in reduced units). The potential energy function combines:

- **Bonded terms:** Harmonic bond lengths ($k_b = 100.0$), bond angles ($k_\theta = 20.0$), and periodic pseudo-dihedral potentials ($k_\phi = 15.0$).
- **Native contacts:** 10-12 Lennard-Jones interactions ($\epsilon_{\text{nat}} = 1.5$) for residue pairs separated by $|i - j| \ge 3$ within an $8.0\text{ \AA}$ cutoff in the target structure.
- **Non-native interactions:** Excluded-volume and sequence-dependent hydrophobic Lennard-Jones potentials parameterized using Kyte-Doolittle hydrophobicity scales and residue radii.

### Tested System: AsLOV2
The workflow was tested on **AsLOV2** (sequence stored in `sequence.txt`). Over the trajectory, the chain undergoes a rapid initial hydrophobic collapse followed by slower topological rearrangement, recovering **79%** of native contacts ($Q = 0.79$) relative to the crystallographic reference.

## Repository Contents

- `torsion_fold.c` — Core Monte Carlo folding engine in C.
- `extract.py` — Extracts central backbone pseudo-dihedral angles across saved frames into `torsion_data.csv`.
- `funnel.py` — Generates the 2D folding funnel (`2d_landscape.png`) plotting system energy against the fraction of native contacts ($Q$).
- `sequence.txt` — Input amino acid sequence used for the AsLOV2 run.
- `trajectory.pdb` & `final_str.pdb` — Multi-model folding trajectory and final folded $C_\alpha$ coordinates.
- `analysis.csv` & `torsion_data.csv` — Stepwise energy, $Q$-fraction, and torsion angle logs.
- `2d_landscape.png` & `folding_kinetics.png` — Folding funnel and kinetic progress plots.
- `folding_movie.mp4` — Rendered animation of the AsLOV2 folding trajectory.

## Usage

1. **Compile the C code:**
   ```bash
   gcc -O3 torsion_fold.c -o torsion_fold -lm
   ```

2. **Run a simulation (requires a reference PDB containing Cα atoms matching the sequence length):**
    ```bash
    ./torsion_fold $(cat sequence.txt) target_structure.pdb
    ```

3. **Run post-processing and plotting:**
    ```bash
    pip install -r requirements.txt
    python extract.py
    python funnel.py
    ```
