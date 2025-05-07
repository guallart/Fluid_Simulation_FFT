# Fluid Simulation using FFT in Julia

This repository contains a Julia implementation of a 2D incompressible fluid simulation based on the Fast Fourier Transform (FFT). The simulation models fluid flow on a rectangular grid, applying forces, advection, viscosity, and projection to enforce incompressibility using spectral methods.

---

## Features

- 2D fluid simulation on a rectangular domain using FFT for efficient computation.
- Incorporates viscosity and external forcing.
- Uses Helmholtz decomposition in Fourier space to enforce incompressibility.
- Semi-Lagrangian advection with linear interpolation.
- Visualization of velocity magnitude (kinetic energy) or vorticity using heatmaps.
- Progress bar to track simulation progress.

---

## Requirements

- Julia 1.6 or later
- Packages:
  - `Plots`
  - `FFTW`
  - `LinearAlgebra`
  - `Distributions`
  - `Interpolations`
  - `ProgressMeter`
  - `Printf`

You can install the required packages using Julia's package manager:

```
using Pkg
Pkg.add([
  "Plots",
  "FFTW",
  "LinearAlgebra",
  "Distributions",
  "Interpolations",
  "ProgressMeter",
  "Printf"
])
```


---

## How to Run

1. Clone or download this repository.
2. Open the Julia script containing the simulation code.
3. Run the script in Julia. The simulation will execute for the specified number of time steps (`N_TIME_STEPS`).
4. The velocity magnitude heatmap will be displayed for each timestep.

---

## Code Overview

### Parameters

- `N_X`, `N_Y`: Number of grid points in x and y directions.
- `Lx`, `Ly`: Physical dimensions of the simulation domain.
- `μ`: Fluid viscosity.
- `Δt`: Time step size.
- `N_TIME_STEPS`: Total number of simulation steps.
- `FORCE_STRENGTH`: Magnitude of the external forcing applied to the fluid.

### Key Components

- **Grid and Wavenumbers:** Defines spatial grid and corresponding Fourier wavenumbers.
- **Force Field:** Gaussian-shaped force applied to the fluid, modulated over time.
- **Advection:** Uses semi-Lagrangian backward tracing with interpolation to update velocity fields.
- **Viscosity:** Applied as a spectral filter in Fourier space.
- **Projection:** Enforces incompressibility by projecting velocity onto divergence-free space in Fourier domain.
- **Visualization:** Heatmap of velocity magnitude is plotted each timestep.

---

## Customization

- Modify the grid resolution (`N_X`, `N_Y`) and domain size (`Lx`, `Ly`) to change simulation scale.
- Adjust `μ` to simulate different fluid viscosities.
- Change `FORCE_STRENGTH` or the shape of the forcing functions to experiment with different flow patterns.
- Toggle visualization between velocity magnitude and vorticity by uncommenting/commenting relevant lines in the plotting section.
- Save frames by uncommenting the `savefig` line to generate images for creating animations.

---

## License

This project is provided as-is for educational and research purposes. Feel free to modify and use the code in your projects.

---

## Contact

For questions or suggestions, please open an issue or contact the author.

---

Enjoy exploring fluid dynamics with FFT-based simulation in Julia!
