# 3D Polymer Dynamics Simulation (C++)

I developed this 3D Molecular Dynamics simulation in C++ to model the conformational dynamics and equilibrium statistics of a coarse-grained polymer chain in solution. The system simulates Brownian/Langevin dynamics using the Grønbech-Jensen–Farago (GJF-2GJ) integrator.

---

## Overview

In this project, I model a linear polymer chain as a sequence of spherical beads connected by harmonic springs, subject to excluded volume interactions and immersed in a thermal solvent bath.

To capture the effects of the solvent without explicitly simulating solvent molecules, I employ Langevin dynamics, which introduces both viscous damping (Stokes friction) and stochastic thermal forces satisfying the fluctuation-dissipation theorem.

---

## Physical Interactions & Modeling

### 1. Chain Initialization
I initialize the polymer as a linear chain of `N = 150` beads centered at the origin `(0, 0, 0)` along the x-axis, separated by the equilibrium bond distance `b0 = 2.0`. Initial particle velocities are drawn from a Maxwell-Boltzmann normal distribution with variance proportional to `k_B * T`.

### 2. Harmonic Bonded Potential
Adjacent beads along the chain backbone are connected via a harmonic potential:

```text
U_bond(r) = (1/2) * k_bond * (r - b0)^2
```

- `k_bond = 100.0`: Spring constant
- `b0 = 2.0`: Equilibrium bond length
- `r`: Distance between bead `i` and bead `i + 1`

Forces are applied along the inter-particle displacement vector, with a cutoff guard to prevent numerical divergence if beads overlap.

### 3. Excluded Volume (Weeks-Chandler-Andersen Potential)
To prevent self-intersection and simulate steric hindrance between all non-bonded pairs of beads, I implemented the purely repulsive Weeks-Chandler-Andersen (WCA) potential (a truncated and shifted Lennard-Jones potential):

```text
U_WCA(r) = 4 * epsilon * [(sigma / r)^12 - (sigma / r)^6] + epsilon    for r < r_cut
U_WCA(r) = 0                                                           for r >= r_cut
```

- `sigma = 1.0`: Bead collision diameter
- `epsilon = 1.0`: Potential well depth
- `r_cut = 2^(1/6) * sigma ≈ 1.12246 * sigma`: Minimum of the Lennard-Jones potential

The corresponding pairwise repulsive force magnitude is:

```text
F_WCA(r) = 24 * (epsilon / r) * [2 * (sigma / r)^12 - (sigma / r)^6]
```

### 4. Thermal Bath & Hydrodynamic Friction
Solvent friction is modeled via Stokes drag on spherical monomers:

```text
alpha = 6 * pi * eta * r_bead
```

- Dynamic viscosity: `eta = 1.0`
- Monomer radius: `r_bead = 1.0`
- Friction constant: `alpha ≈ 18.8736`

Thermal fluctuations are injected at each step via independent 3D Gaussian random noise with variance:

```text
Var(beta) = 2 * alpha * k_B * T * dt
```

---

## Grønbech-Jensen–Farago (GJF-2GJ) Integrator

Standard Langevin integrators often introduce systematic errors in kinetic energy or spatial distributions at finite time steps. I chose the Grønbech-Jensen–Farago (GJF) method because it provides robust, statistically sound integration for discrete Langevin systems.

### Integration Coefficients
```text
a = [1 - alpha * dt / (2 * m)] / [1 + alpha * dt / (2 * m)]
b = 1 / [1 + alpha * dt / (2 * m)]
```

### Discrete Step Updates
For step `n` to `n + 1`:

1. **Position Update**:
   ```text
   r(n+1) = r(n) + b * [dt * v(n) + (dt^2 / (2 * m)) * f(n) + (dt / (2 * m)) * beta(n+1)]
   ```

2. **Force Evaluation**:
   ```text
   f(n+1) = F_bond(r(n+1)) + F_WCA(r(n+1))
   ```

3. **Velocity Update**:
   ```text
   v(n+1) = a * v(n) + (dt / (2 * m)) * [a * f(n) + f(n+1)] + (b / m) * beta(n+1)]
   ```

4. **On-site 2GJ Effective Velocity**:
   To ensure the computed kinetic energy matches equipartition without discretization bias:
   ```text
   u_2gj = [r(n+1) - r(n)] / [sqrt(b) * dt]
   ```

---

## Simulation Parameters

| Parameter | Symbol | Value | Description |
| :--- | :--- | :--- | :--- |
| Number of Monomers | `n_part` | `150` | Number of beads in the chain |
| Thermal Energy | `K_Bt` | `1.0` | `k_B * T` (reduced units) |
| Monomer Mass | `m` | `1.0` | Mass of each bead |
| Time Step | `dt` | `0.01` | Integration step size |
| Total Time Steps | `n_steps` | `500,000` | Total simulation duration (`5,000` time units) |
| Sampling Rate | `record_every`| `200` | Sample recorded every 200 steps (2,501 frames) |
| Bond Stiffness | `K_bond` | `100.0` | Harmonic bond constant |
| Equilibrium Bond Length | `b0` | `2.0` | Rest length of springs |
| WCA Diameter | `sigma` | `1.0` | Particle excluded-volume diameter |
| WCA Energy Scale | `epsilon` | `1.0` | Interaction strength |

---

## Measured Observables & Post-Processing

I discard the initial 20% of recorded frames as an equilibration / burn-in period (`burn = n_frames / 5`) before computing equilibrium ensemble averages:

1. **Mean Kinetic Energy per Degree of Freedom**:
   ```text
   <E_k> / (3 * N) = (1 / (3 * N)) * sum_{i} (1/2 * m * |u_2gj,i|^2)
   ```
   I verify this against the equipartition theorem expectation: `(1/2) * k_B * T = 0.5`.

2. **Mean End-to-End Distance**:
   ```text
   R_ee = < |r_{N-1} - r_{0}| >
   ```
   Measures overall chain elongation and spatial extension across snapshots.

3. **Radius of Gyration**:
   ```text
   R_g^2 = (1 / N) * sum_{i=0}^{N-1} |r_i - r_com|^2
   ```
   Characterizes the effective size and compactness of the polymer coil relative to its center of mass `r_com`.

4. **Mean Bond Length**:
   ```text
   <b> = (1 / (N - 1)) * sum_{i=0}^{N-2} |r_{i+1} - r_i|
   ```
   Verifies that the harmonic bond potential successfully maintains the expected bond spacing `b0 = 2.0`.

---

## Building and Running

I implemented the code using standard C++11 / C++17 library components (`<vector>`, `<random>`, `<chrono>`, `<cmath>`). No external dependencies are required.

### Compilation
Compile with an optimizing C++ compiler such as `g++` or `clang++`:

```bash
g++ -O3 -std=c++17 Polymer_Sim_C.cpp -o polymer_sim
```

### Execution
Run the compiled binary:

```bash
./polymer_sim
```

### Example Output
```text
n_steps values is 500000
n_part values is 150
<E_k> per DOF (2GJ velocity) = 0.5002   (theory: 0.5000)
end_end_dist_sq 23.4182
radius of gyration is 7.2145
bond_lengths 1.9998
no of particle is 150
Process completed successfully.
Execution time: 42.1582
```

---

## Future Enhancements
- Spatial cell-linked lists or neighbor lists to reduce WCA computation from O(N^2) to O(N).
- Trajectory exporting to `.xyz` or `.pdb` format for 3D visualization in VMD or Ovito.
- Dynamic structure factor and mean-squared displacement (MSD) calculations to investigate Rouse / Zimm relaxation dynamics.
