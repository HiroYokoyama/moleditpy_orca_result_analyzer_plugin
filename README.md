# MoleditPy ORCA Result Analyzer Plugin

A comprehensive plugin for MoleditPy to analyze and visualize results from ORCA quantum chemistry calculations.

![Main UI](img/main.png)

## Features

### 1. SCF Trace
Real-time convergence visualization for SCF energy cycles.
- **Concatenated View**: View a single continuous plot of all SCF cycles throughout the calculation.
- **Interactive Tools**: Full zoom, pan, and save support via an integrated Matplotlib toolbar.

### 2. MO Analysis
- **Levels**: View orbital energies and occupancy.
- **Visualization**: Generate and view 3D Cubes for selected orbitals (Isosurfaces).
- **Limitations**: Currently supports S, P, D, and F shells.

### 3. Optimization / Scan
Analyze **Geometry Optimizations** and **Relaxed Surface Scans**.
- **Interactive Graph**: Plot Energy vs. Step. Click points to update the 3D structure.
- **Animation**: Play/Pause trajectory animations and export as GIFs.

### 4. Forces
Analyze structural forces for the current structure or the entire trajectory.
- **Historical Gradients**: Capture and display force vectors for **every** step of an optimization or scan.
- **Step-by-Step Navigation**: Precise control with `<` and `>` buttons or the trajectory slider.
- **Convergence Tracking**: Multi-line display of **RMS/MAX Gradient** and **RMS/MAX Step**, color-coded (**Green for YES**, **Red for NO**).

### 5. Atomic Charges
- **Populations**: Mulliken, Loewdin, and Hirshfeld populations (if available).
- **3D Coloring**: Color atoms in the 3D viewer based on charge value.

### 6. Dipole Moment
- **Vector Visualization**: Display the total dipole moment vector magnitude and direction in the 3D viewer.

### 7. Frequencies
Visualize vibrational modes and spectra.
- **IR/Raman**: Stick and broadened spectra plots.
- **Visualization**: Animated vibrational modes with vector arrows.

### 8. Thermochemistry
Detailed analysis of thermodynamic properties based on frequency calculations.
- **Broad Summary**: Access Electronic Energy, ZPE, Enthalpy (H), Gibbs Free Energy (G).
- **Detailed Corrections**: Optional breakdown of vibrational, rotational, and translational contributions.

### 9. TDDFT
Analyze electronic excitations and absorption spectra.
- **Spectra**: Absorption and CD spectra with Gaussian broadening and transition analysis.

### 10. NMR
Advanced NMR chemical shielding validation and visualization.
- **Stick Spectrum**: Nucleus-specific stick spectra (1H, 13C, etc.) with experimental reference standards and equivalent atom merging.

## Interface & Usability
- **Logical Workflow**: Buttons are grouped by task: Electronic -> Geometry -> Properties -> Spectroscopy.
- **Concise UI**: Clean, professional labels and modeless windows for side-by-side editing.
- **Keyboard Shortcuts**: `Ctrl+O` (Open), `Ctrl+R` (Reload).

## Installation
Ensure the `orca_result_analyzer` folder is placed in your MoleditPy plugins directory.

## Requirements
- **ORCA Output**: Reads `.out` or `.log` files. Basis set info (`Print[P_Basis] 2`) is required for MO Cube generation.
- **Dependencies**: `matplotlib`, `Pillow` (PIL), and `pyvista` (for 3D vectors). `rdkit` is recommended for high-quality structure generation.

## Required ORCA Keywords
To ensure all features (especially **MO Cube Generation**) work correctly, include the following in your ORCA input:
```text
%output
  Print[P_Basis] 2  # Required for Basis Set parsing (Cube Gen)
  Print[P_Mos] 1    # Ensure MO coefficients are printed
end
```
*Note: Standard output is usually sufficient for Geometry and Energies, but Basis Set information is strictly required for generating cubes.*
