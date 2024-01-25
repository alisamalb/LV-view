# Localized Volume Viewer: creating volume restraints for biased molecular dynamics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7575202.svg)](https://doi.org/10.5281/zenodo.7575202)

## Description

This script will help you creating a volume restraint for gromacs + plumed in the shape of a paraboloid.

Run the code here on Google Colab

This script takes two files as input:

    A .gro file containing the protein and the ligand
    A .ndx file contating information about the atoms of the first file

! You can create the second file using Gromacs make_ndx command, using the first file as input. Do not use the .ndx file of the solvated system.

The script gives two outputs:

    A volume.gro file containing the pseudoatoms distributed in space to represent the chosen paraboloid
    A reference_rotated.gro file containg the protein in .gro format rotated in space so that the protein COM - ligand COM vector lies on the x-axis.

You can use the second file to create a reference for the FIT_TO_TEMPLATE function in plumed2. Use gromacs editconf to convert it to PDB. Use grep to extract the "CA" lines.
