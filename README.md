# MoleditPy ORCA Result Analyzer Plugin

A comprehensive plugin for MoleditPy to analyze and visualize results from ORCA quantum chemistry calculations.

![Main UI](img/main.png)

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
- **Performance Note**: For large systems or higher resolution, using ORCA's `orca_plot` utility to generate `.cube` files and the **Cube File Viewer** plugin is recommended for better performance.

### 5. Thermochemistry
Detailed analysis of thermodynamic properties based on frequency calculations.
- **Broad Summary**: Access Electronic Energy, ZPE, Enthalpy (H), Gibbs Free Energy (G), and Entropy (S) terms.
- **Detailed Corrections**: Optional breakdown of vibrational, rotational, and translational contributions to thermal energy and entropy.
- **Export**: Copy table data to clipboard or export to CSV.

### 6. TDDFT Spectrometry
Analyze electronic excitations and absorption spectra.
- **Absorption Spectrum**: View interactive transitions with adjustable Gaussian broadening.
- **Peak List**: Detailed excitation energies, oscillator strengths, and state descriptions.

### 7. Atomic Charges
- **Populations**: Mulliken and Loewdin charge analysis.
- **3D Coloring**: Color atoms in the 3D viewer based on charge (Red=Negative, Blue=Positive).

### 8. Dipole Moment
- **Vector Visualization**: Display the total dipole moment vector magnitude and direction in the 3D viewer.

### 9. Geometry Stability & Forces
Analyze cartesian gradients and structural stability.
- **Forces (Gradients)**: 3D visualization of gradient vectors on atoms.
- **Hessian Analysis**: Manually load ORCA `.hess` files to display force constants and perform stability analysis.

### 10. NMR Analysis
Advanced NMR chemical shielding validation and visualization.
- **Interactive Spectrum**: Stick spectrum (Stem plot) with nucleus-specific filtering (e.g., 1H, 13C).
- **Reference Standards**: Support for experimental reference shielding constants ($\sigma_{ref}$) and chemical shifts ($\delta_{ref}$).
- **Peak Identification**: Click peaks to highlight corresponding atoms in the 3D viewer.
- **Peak Merging**: Group equivalent atoms for simplified spectrum analysis.
- **Persistent Settings**: Manual reference values and merged peaks are saved automatically.

## Interface & Usability
- **Modeless Design**: Analysis windows can stay open alongside the main editor, allowing for side-by-side comparison.
- **Keyboard Shortcuts**: Common actions bound to keys (e.g., `Ctrl+O` to open, `Ctrl+R` to reload).

## Installation
Ensure the `orca_result_analyzer` folder is placed in your MoleditPy plugins directory.

## Requirements
- **ORCA Output**: The plugin reads the main output file (`.out` or `.log`) and retrieves Basis Set information directly for advanced 3D orbital visualization.
- **Dependencies**: 
    - `matplotlib`: For graphing.
    - `Pillow` (PIL): For GIF generation.

## Required ORCA Keywords
To ensure all features (especially **MO Cube Generation**) work correctly, include the following in your ORCA input:
```text
%output
  Print[P_Basis] 2  # Required for Basis Set parsing (Cube Gen)
  Print[P_Mos] 1    # Ensure MO coefficients are printed
end
```
*Note: Standard output is usually sufficient for Geometry and Energies, but Basis Set information is strictly required for generating cubes.*

