# ORCA Result Analyzer Plugin

A comprehensive plugin for MoleditPy to analyze and visualize results from ORCA quantum chemistry calculations.

## Features

### 1. Overview & Geometry
- **Summary Dashboard**: Quick view of SCF energy, convergence status, and geometry optimization steps.
- **3D Structure**: Visualize the final geometry or any step from a trajectory.

### 2. Trajectory & Scan Analyzer
Analyze **Geometry Optimizations** and **Relaxed Surface Scans**.
- **Interactive Graph**: Plot Energy vs. Step. Click points to update the 3D structure.
- **Animation**: Play/Pause trajectory animations.
- **Exports**:
    - **GIF**: Save animated trajectories with transparency options.
    - **CSV**: Export energy data.
    - **Image**: Save the energy profile graph.

### 3. Frequency Analysis (IR / Raman)
Visualize vibrational modes and spectra.
- **Spectra**:
    - **IR Spectrum**: Transmittance-style plot (inverted Y-axis).
    - **Raman Spectrum**: Normal intensity plot.
- **Visualization**: Animated vibrational modes with vector arrows.
- **GIF Export**: Save animations of specific vibrational modes as GIFs.
- **Scaling**: Adjustable frequency scaling factors.

### 4. Molecular Orbitals (MO)
- **Levels**: View orbital energies and occupancy.
- **Visualization**: Generate and view 3D Cubes for selected orbitals (Isosurfaces).
- **Composition**: Detailed element/atom contributions.

### 5. Electronic & Atomic Properties
- **Atomic Charges**: Mulliken and Loewdin charges.
    - **3D Coloring**: Color atoms in the 3D viewer based on charge (Red=Negative, Blue=Positive).
- **Dipole Moment**: Visualize the total dipole moment vector and magnitude.
- **Forces (Gradients)**: 3D visualization of gradient vectors on atoms.

### 6. Spectroscopic Properties
- **TDDFT/TDA**: UV-Vis absorption spectra with adjustable geometric broadening (Gaussian lineshape).
- **NMR**: Chemical shielding summary validation.

### 7. UI Refinements
- **Color Schemes**: The Charge Viewer now supports custom color schemes (e.g., Red-White-Blue, Blue-White-Red, etc.) with a visual gradient bar legend.

## Installation
Ensure the `orca_result_analyzer` folder is placed in your MoleditPy plugins directory.

## Requirements
- **ORCA Output**: The plugin reads the main output file (`.out` or `.log`) and looks for the `.fchk` file for advanced 3D orbital visualization.
- **Dependencies**: 
    - `matplotlib`: For graphing.
    - `Pillow` (PIL): For GIF generation.

## Required ORCA Keywords
To ensure all features (especially **MO Cube Generation**) work correctly, include the following in your ORCA input:
```text
! NormalPrint
%output
  Print[P_Basis] 2  # Required for Basis Set parsing (Cube Gen)
  Print[P_Mos] 1    # Ensure MO coefficients are printed
end
```
*Note: Standard output is usually sufficient for Geometry and Energies, but Basis Set information is strictly required for generating cubes.*
