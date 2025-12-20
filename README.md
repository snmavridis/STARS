# Savas' Trajectory And Radar Simulation

## Overview

This repository contains a modified version of Zipfel's CADAC++ software to run simulations of system engagements.
These systems engagements can include air-to-air and air-to-ground scenarios.
This work is a result of Savas' research as published at AIAA SciTech 2026 titled "A Modeling and Simulation Framework for Evaluating Aircraft Survivability".
This repository contains a complete, self-contained workflow for running a C++-based simulation, generating output data, and performing post-processing and analysis using MATLAB. 
It is intended to serve as a single, reproducible “one-stop shop” for simulation, data generation, and analysis.

---

## Repository Structure
```
STARS/
├── src/ # C++ source files for the simulation
├── include/ # C++ header files
├── build/ # Build artifacts (generated; not tracked by Git)
│
├── data/ # Input data and configuration files
├── results/ # Simulation outputs (generated; typically not tracked)
│
├── analysis/ # MATLAB scripts and functions for data analysis
├── scripts/ # Shell scripts for running simulations and automation
├── docs/ # Documentation and supporting material
├── misc/ # Miscellaneous or archival material
│
├── README.md # This file
├── LICENSE
└── .gitignore
```
---

## Directory Descriptions

### `src/`
Contains all C++ implementation files (`.cpp`) for the simulation.
No headers, data files, or scripts should be placed here.

---

### `include/`
Contains all C++ header files (`.h` / `.hpp`) used by the simulation.
The build system includes this directory when compiling.

---

### `build/`
Contains compiled objects and executables generated during the build process.

- This directory is **generated automatically**
- It is **not tracked by Git**
- Typically contains subdirectories such as `obj/` and `bin/`

---

### `data/`
Contains **input files only**, such as:
- ASCII scenario definitions
- Configuration files
- Parameter sets

These files define simulation inputs and are version-controlled so runs can be reproduced.

---

### `results/`
Contains **simulation outputs**, such as:
- Raw output data
- Processed data
- Generated figures

This directory is typically excluded from version control unless example outputs are intentionally retained.

---

### `analysis/`
Contains MATLAB scripts and functions used to analyze simulation outputs.

Typical contents include:
- Data loading scripts
- Signal processing routines
- Plotting and visualization tools

Analysis scripts are written assuming simulation outputs reside in `results/`.

---

### `scripts/`
Contains shell scripts or helper utilities for:
- Running the simulation
- Batch execution
- Cleaning outputs
- Automating end-to-end workflows

These scripts improve reproducibility and ease of use.

---

### `docs/`
Contains documentation such as:
- Build instructions
- File format descriptions
- Theoretical background
- Figures used in reports or papers

---

### `misc/`
Contains miscellaneous, non-critical material such as:
- Archived notes
- Reference material
- Scratch work

Nothing in this directory should be required to build, run, or analyze the simulation.

---

## Typical Workflow

1. **Build the simulation**

- make

2. **Run the simulation**

- ./build/bin/<executable> data/<input-file>

3. **Analyze results**
- Open MATLAB
- Run scripts in `analysis/` to process data in `results/`

---

## Notes

- Generated files (`build/`, `results/`) are excluded from Git to keep the repository clean.
- The directory structure is intentionally modular to support future expansion and collaboration.
