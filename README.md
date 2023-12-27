Please Cite:

Lahkar, S., & Ranganathan, R. (2023). Competing mechanisms govern the thermal rectification behavior in semi-stochastic polycrystalline graphene with graded grain-size distribution. *Carbon*, 118638.

Report any bugs/issues to simantalahkar@hotmail.com / s.lahkar@tue.nl

# Graded-Polygraphene-Thermal-Rectification-Competing-Mechanisms
This repository contains the LAMMPS and python scripts, along with the most important data, that I created to carry out a thorough analysis of the Thermal Rectification (TR) phenomenon in semi-stochastically generated atomistic models of polycrystalline graphene having graded variation of grain sizes - using Molecular Dynamics simulations and phonon calculations based on implementation of the dissipation fluctuation theory. 
![Fig 1. Interplay between TR mechanisms identified in this work](mechanisms3.png)
**Fig 1.** *Possible Interplay between TR mechanisms as identified by us.*

These files have been provided with the intention that the reader is able to reproduce all the main results of our study. Below you will find a set of instructions to use these scripts and data as intended. The BNC.tersoff file contains the potential parameters.

1. The LAMMPS-scripts folder contains all the scripts for various molecular dynamics simulations on the graphene structures:
    - eqlbm_phonon.in.lammps does looped-minimization and NVT equilibration and generates the dynamical matrix during equilibration. *Note that this is done on periodic structures*
    - minimize_eqlbm.in.lammps does looped minimization and NVT equilibration. *Note that this is done on periodic structures*
    - run.min_nemd* scripts carry out minimization, NVT equilibration, followed by an non-equilibrium MD (NEMD) simulation of heat transport (in either of the opposite directions) with the thermal baths near the start and middle of the structure longitudinally. It also calculated the dynamical matrix for the steady state thermal transport simulation. *Note that for this simulation, longitudinally symmetric periodic structures are used*
    -  run.nemd_* scripts carry out NEMD simulation of heat transport (in either of the opposite directions) with the thermal baths at the two ends of the structure longitudinally. *Note that longitudinally non-periodic asymmetric structures are used for these simulations.* Hence, the periodic graded structures post minimization and equilibration generally need to be spliced and re-sequenced, so that the input structures have a *proper unidirectional gradient of grain sizes*. The *\*atomic\** and *\*full\** files read the input structure data file written in the respective formats. *Note that when using the* full *structure, the atomic velocities are taken from the end state of the preceding equilibration and are run not reinitialized before the NEMD simulations.*

2. The Python-scripts contain the different scripts for the pre and post-processing of the different data.
    - renumber_* scipts re-sequence the atomic IDs in the structure file - which is necessary because after the splicing the atomic IDs are no longer continuous.
    - nemd_* scripts analyze the output from the NEMD simulations to calculate the meaningful properties (including the thermal conductivity values) from the raw data and create useful plots. This also creates different plots for monitoring and assessing the system properties throughout the simulation. 
    - spatial_phonon_data.py does various physical phonon calculations (including the phonon energy density maps) using the eigvec.dat file and the structure file (containing the per atom temperature values) for a given simulation. *The eigvec.dat can be generated using the "PHANA" package developed by Prof. Lingti Kong (github repository at lingtikong/phana) from the binary dynamical matrix output file of the corresponding nemd-phonon simulations.* **Copy and replace the contents of the generated spatial_data\*.custom file into either the eqlbm\*.custom or heated\*.custom file (depending on whether it is an NEMD or NVT equilibration run) generated from the LAMMPS simulation and rename the property columns appropriately (the final files for the structures studied are also shared in the *Structures-and-data* folder for reference), which can be easily viewed and analyzed in tools like Ovito.**
    - interprete_phonon_data.py does statistical interpretation of the output data from running the spatial_phonon_data.py script and also plots the participation ratio for the phonon eigen values.
    - mechanism3.py plots a graphical representation of the strengths of the two interacting mechanisms of TR based on one possible inferred relationship of the respective mechanisms to the structural parameters - length and grain size asymmetry.

3. The Structures-and-data folder contains all the main structures and output data obtained in this studied (grouped as per the corresponding Figures in our associated published paper-see Lahkar, S., & Ranganathan, R. (2023), *Carbon*, 118638) that can be used for reproduction of the results, and for comparative reference.