# CD3-TCR Complex Simulation

This repository contains 3 files:

- **`Create_starting_pos.py`**  
  Contains Python script to generate a starting position for the 1 complex configuration.

- **`in.put_cd3_complex.lmp`**  
  Contains LAMMPS script to conduct a coarse-grained simulation of the cytoplasmic domain of the CD3-TCR complex near a membrane.

- **`params.in`**  
  Input file containing the parameters used in `in.put_cd3_complex.lmp`.

---

## Quickstart Guide

### To run 1 complex simulations:

1. Generate the initial configuration using `Create_starting_pos.py`.
2. Define the chosen file name for the initial position in the LAMMPS input script `in.put_cd3_complex.lmp` (line 14).
3. Define the chosen file name for the parameters in `params.in`.
4. Set output file for the trajectory in `in.put_cd3_complex.lmp` (line 17).
5. Set output file for the parameters in `in.put_cd3_complex.lmp` (line 18).
6. Run the LAMMPS script.

### To run 5 complex simulations:

- Same steps as for the 1 complex simulation.  
- The only difference is that `Create_starting_pos.py` must be modified to generate the initial position for 5 complexes.  
- Define the chosen file name for the initial position in the LAMMPS input script `in.put_cd3_complex.lmp` (line 14).

---

The code has been tested and run using **LAMMPS (17 Apr 2024)**, which can be downloaded at:  
👉 [https://www.lammps.org/download.html](https://www.lammps.org/download.html)
