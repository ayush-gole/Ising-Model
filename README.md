# 2D Ising Model Simulation

Monte Carlo simulation of the 2D Ising model using the **Metropolis algorithm**, written in C++ with Python for plotting.

Simulates a 2D lattice of magnetic spins and computes thermodynamic observables across temperatures — capturing the ferromagnetic phase transition near the theoretical critical temperature *Tc ≈ 2.269*.

---

## Output

![Thermodynamic Observables](ising_thermo.png)

The plots show energy, magnetization, specific heat and magnetic susceptibility as a function of temperature. The sharp peaks in Cv and χ near T ≈ 2.269 clearly marks the second-order phase transition, which matches theory pretty well for a 30×30 lattice.

---

## Features

- Supports **ferromagnetic** (J = +1) and **antiferromagnetic** (J = −1) coupling
- Both **periodic** and **free** boundary conditions
- Adaptive temperature grid with finer resolution near Tc for cleaner peak detection
- Computes energy per spin, magnetization per spin, specific heat (Cv) and susceptibility (χ)
- Estimates Tc from peaks of Cv and χ
- Results saved to `.txt` file automatically

---

## How to Run

> Both files should be in the **same directory**

**Step 1 — Compile and run the C++ simulation:**
```bash
g++ -O2 -o ising Ising_model_main.cpp
./ising
```
Generates output file like `ising_results_ferro_PB.txt` depending on your settings.

**Step 2 — Plot the results:**
```bash
python Ising_Model_Plotting.py
```
Reads the output file and saves `ising_thermo.png`.

---

## Parameters

All of these are near the top of `Ising_model_main.cpp` and easy to change:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N` | 30 | Lattice size (N × N) |
| `J` | +1.0 | Coupling constant (+1 ferro, −1 anti-ferro) |
| `FREE_BC` | false | Boundary condition (false = periodic) |
| `eq_steps` | 5000 | Equilibration sweeps |
| `measure_steps` | 20000 | Measurement sweeps |

---

## Tech Stack

- **C++** — Metropolis Monte Carlo engine, thermodynamic calculations
- **Python** — Plotting with `numpy` and `matplotlib`

---

## Physics Background

The 2D Ising model is one of the foundational models in statistical mechanics. Each site on the lattice holds a spin *s = ±1* and the Hamiltonian is:

**H = −J Σ sᵢsⱼ**

where the sum runs over nearest-neighbor pairs. Below the critical temperature spins tend to align, producing net magnetization. Above Tc, thermal fluctuations destroy that order. The Metropolis algorithm stochastically samples configurations according to the Boltzmann distribution — which makes it practical for systems where exact enumeration isn't feasible.
