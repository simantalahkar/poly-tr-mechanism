Please Cite:

Lahkar, S., & Ranganathan, R. (2023). Competing mechanisms govern the thermal rectification behavior in semi-stochastic polycrystalline graphene with graded grain-size distribution. *Carbon*, 118638.

Report any bugs/issues to simantalahkar@hotmail.com / s.lahkar@tue.nl

# Graded-Polygraphene-Thermal-Rectification-Competing-Mechanisms
This repository contains the LAMMPS and python scripts, along with the most important data, that I created to carry out a thorough analysis of the Thermal Rectification (TR) phenomenon in semi-stochastically generated atomistic models of polycrystalline graphene having graded variation of grain sizes - using Molecular Dynamics simulations and phonon calculations based on implementation of the dissipation fluctuation theory. 
![Fig 1. Interplay between TR mechanisms identified in this work](mechanisms3.png)
**Fig 1.** *Possible Interplay between TR mechanisms as identified by us.*

These files have been provided with the intention that the reader will be able to reproduce all the main results of our study. Below you will find a set of instructions to use these scripts and data.

1. 

Generate the Dynamical matrix useful for phonon calculations, along with the associated Python analysis scripts that uses the outputs from LAMMPS simulation as well as from the "PHANA" package developed by Prof. Lingti Kong (github repository at lingtikong/phana). 

After running the LAMMPS simulation using the structure and map files created using the latgen tool (input script: in.lammps):

1. Use the phana package to calculate the eigenvectors and eigenvalues from the dynamical matrix obtained from the simulation (eigenvec.dat). The simulation parameters have been optimized for accuracy of phonon calculation and computational efficiency.

2. Run the spatial_phonon_data.py Python script to obtain the necessary information about atomic contributions to propagating phonons and participation ratios in the form of raw data.

3. Run the interprete_phonon_data.py for subsequent postprocessing of the raw data to obtain the participation ratio plots, along with other useful information on the spatial energy density required.

4. Copy and replace the contents of the generated spatial_data*.custom file into the eqlbm*.custom file generated from the LAMMPS simulation and rename the property columns appropriately (a typical final file is also shared for reference), which can be easily viewed and analyzed in tools like Ovito.

5. Use the nemd.py to calculate the thermal conductivity using the temperature profile-data and the energy exchange-data generated from the simulations, which also creates different plots for monitoring and assessing the system properties throughout the simulation. 
