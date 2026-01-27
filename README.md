# MoleditPy ORCA Result Analyzer Plugin

A comprehensive plugin for MoleditPy to analyze and visualize results from ORCA quantum chemistry calculations.

![Main UI](img/main.png)

## Features

### 1. Molecular Orbitals (MO)
- **Levels**: View orbital energies and occupancy.
- **Visualization**: Generate and view 3D Cubes for selected orbitals (Isosurfaces).
- **Limitations**: Currently supports S, P, D, and F shells. **G-shells (and higher) are not yet supported** and will be skipped during visualization.
- **Performance Note**: For large systems or higher resolution, using ORCA's `orca_plot` utility to generate `.cube` files and the **Cube File Viewer** plugin is recommended for better performance.

### 2. Frequency Analysis (IR / Raman)
Visualize vibrational modes and spectra.
- **Spectra**:
    - **IR Spectrum**: Transmittance-style plot (inverted Y-axis).
    - **Raman Spectrum**: Normal intensity plot.
- **Visualization**: Animated vibrational modes with vector arrows.
- **GIF Export**: Save animations of specific vibrational modes as GIFs.
- **Scaling**: Adjustable frequency scaling factors.

### 3. Scan & Optimization Results (Trajectory)
Analyze **Geometry Optimizations** and **Relaxed Surface Scans**.
- **Interactive Graph**: Plot Energy vs. Step. Click points to update the 3D structure.
- **Animation**: Play/Pause trajectory animations.
- **Exports**:
    - **GIF**: Save animated trajectories with transparency options.
    - **CSV**: Export energy data.
    - **Image**: Save the energy profile graph.

### 4. Geometry Stability & Forces
Analyze cartesian gradients and structural stability.
- **Forces (Gradients)**: 3D visualization of gradient vectors on atoms.
- **Hessian Analysis**: Manually load ORCA `.hess` files to display force constants and perform stability analysis.

### 5. Thermochemistry
Detailed analysis of thermodynamic properties based on frequency calculations.
- **Broad Summary**: Access Electronic Energy, ZPE, Enthalpy (H), Gibbs Free Energy (G), and Entropy (S) terms.
- **Detailed Corrections**: Optional breakdown of vibrational, rotational, and translational contributions to thermal energy and entropy.
- **Export**: Copy table data to clipboard or export to CSV.

### 6. TDDFT Spectra
Analyze electronic excitations and absorption spectra.
- **Absorption Spectrum**: View interactive transitions with adjustable Gaussian broadening.
- **Spectrum Types**: Support for Absorption (Oscillator Strength) and CD (Rotatory Strength).
- **Peak List**: Detailed excitation energies, oscillator strengths, and state descriptions.

### 7. Dipole Moment
- **Vector Visualization**: Display the total dipole moment vector magnitude and direction in the 3D viewer.
- **Info**: Displays total magnitude and vector components.

### 8. Atomic Charges
- **Populations**: Support for Mulliken, Loewdin, and Hirshfeld charge analysis (if available).
- **3D Coloring**: Color atoms in the 3D viewer based on charge value.
- **Customization**: Customizable color schemes for positive/negative/neutral charges.

### 9. NMR Shielding
Advanced NMR chemical shielding validation and visualization.
- **Interactive Spectrum**: Stick spectrum (Stem plot) with nucleus-specific filtering (e.g., 1H, 13C).
- **Reference Standards**: Support for experimental reference shielding constants ($\sigma_{ref}$) and chemical shifts ($\delta_{ref}$).
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
    - `rdkit` (Optional but recommended): For bond percepton and 3D structure generation.

## Required ORCA Keywords
To ensure all features (especially **MO Cube Generation**) work correctly, include the following in your ORCA input:
```text
%output
  Print[P_Basis] 2  # Required for Basis Set parsing (Cube Gen)
  Print[P_Mos] 1    # Ensure MO coefficients are printed
end
```
*Note: Standard output is usually sufficient for Geometry and Energies, but Basis Set information is strictly required for generating cubes.*
