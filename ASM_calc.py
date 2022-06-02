#!/usr/bin/env python

import subprocess
import os
import glob
import shutil
import argparse
import sys
import numpy as np
from datetime import datetime

# List of atoms used to construct fragment A:
Atom_list = [84, 85, 86]

# Atom list of the selected bond   
Atom_bond_list = [61, 84]

# Input settings        
head = ['! RKS BP86 D3BJ def2-TZVP def2/J SlowConv TightSCF Grid6\n\n',
        '%scf\n',
        '  MaxIter  2000\n',
        '  end\n\n',
        '%maxcore  1900\n\n',
        '%pal\n',
        '  nprocs  4\n',
        '  end\n\n',
        '%method\n',
        '  SpecialGridAtoms 53\n',
        '  SpecialGridIntAcc 7\n',
        '  end\n\n',
        '* xyz 0 1\n']

head_solv = ['! RKS BP86 D3BJ def2-TZVP def2/J SlowConv TightSCF Grid6\n\n',
             '%scf\n',
             '  MaxIter  2000\n',
             '  end\n\n',
             '%maxcore  1900\n\n',
             '%pal\n',
             '  nprocs  4\n',
             '  end\n\n',
             '%method\n',
             '  SpecialGridAtoms 53\n',
             '  SpecialGridIntAcc 7\n',
             '  end\n\n',
             '%cpcm\n',
             '  SMD True\n',
             '  SMDsolvent "DIBUTYLETHER"\n',
             '  end\n\n',
             '* xyz 0 1\n']


head_CO2 = ['! RKS BP86 D3BJ def2-TZVP def2/J SlowConv TightSCF Grid6\n\n',
            '%scf\n',
            '  MaxIter  2000\n',
            '  end\n\n',
            '%maxcore  1900\n\n',
            '* xyz 0 1\n']

head_solv_CO2 = ['! RKS BP86 D3BJ def2-TZVP def2/J SlowConv TightSCF Grid6\n\n',
                 '%scf\n',
                 '  MaxIter  2000\n',
                 '  end\n\n',
                 '%maxcore  1900\n\n',
                 '%cpcm\n',
                 '  SMD True\n',
                 '  SMDsolvent "DIBUTYLETHER"\n',
                 '  end\n\n',
                 '* xyz 0 1\n']


# ORCA binaries directory in the computer:
ORCA_dir = "/path/to/ORCA/binaries/orca" # Path to the ORCA binaries


# Functions used in this script
def open_xyz(fname):
    with open(fname, 'r') as f:
        xyz = f.readlines()
        xyz = xyz[2:]
        molecule = []
        for line in xyz:
            molecule.append(line)
    return molecule

def create_fragA(fragA_atoms, full_molecule):
    A = []
    for atom in fragA_atoms:
        A.append(full_molecule[atom])
    return A

def create_fragB(full_molecule, fragment):
    B = []
    for Atom in range(0, len(full_molecule)):
        B.append(full_molecule[Atom])
    for atom in B[:]:
        if atom in fragment:
            B.remove(atom)
    return B

def write_input(inp_settings, structure, name): 
    """
    Takes the XYZ coordinates of the molecule and its name 
    and returns the input of the molecule.
    """
    head = inp_settings
    tail = ['*']
    inp = head + structure + tail
    with open("{}.inp".format(name), 'w+') as f:
        for line in inp:
            f.write(line)

def write_input_solv(inp_settings, structure, name):  
    """
    Takes the XYZ coordinates of the molecule and its name 
    and returns the input of the molecule solvated.
    """
    head = inp_settings
    tail = ['*']
    inp = head + structure + tail
    with open("{}_solv.inp".format(name), 'w+') as f:
        for line in inp:
            f.write(line)

def open_file(file):
    """
    Opens the output file for parsing.
    """    
    with open(file, 'r') as f:
        data = f.readlines()
    return data

def finalE(file):
    """ 
    Takes the last electronic energy value from the output.
    """    
    words = 'FINAL SINGLE POINT ENERGY'
    elec_energy = []
    for line in file:
        if words in line:
            Eel = line.split()
            Eel = float(Eel[-1])
            elec_energy.append(Eel)
            last_electronic_energy = elec_energy[-1]
    return last_electronic_energy

def calculate_distance(molecule, atom_bond_list):
    """
    Opens a XYZ file and calculate the bond distance between two atoms 
    passed as list of atoms in atom_bond_list.
    """
    with open(molecule, 'r') as f:
        xyz = f.readlines()
        xyz = xyz[2:]
        mol = []
        for line in xyz:
            mol.append(line)
            
    x_coord = []
    y_coord = []
    z_coord = []
    for line in range(0, len(mol)):
        data = mol[line].split()
        x_coord.append(data[1])
        y_coord.append(data[2])
        z_coord.append(data[3])
        
    atom1_x_coord, atom1_y_coord, atom1_z_coord = float(x_coord[atom_bond_list[0]]), float(y_coord[atom_bond_list[0]]), float(z_coord[atom_bond_list[0]])
    atom2_x_coord, atom2_y_coord, atom2_z_coord = float(x_coord[atom_bond_list[1]]), float(y_coord[atom_bond_list[1]]), float(z_coord[atom_bond_list[1]])
    x_distance = atom1_x_coord - atom2_x_coord
    y_distance = atom1_y_coord - atom2_y_coord
    z_distance = atom1_z_coord - atom2_z_coord
    bond_length_12 = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
    return bond_length_12


startTime = datetime.now()
parent_dir = os.getcwd()
file_location = os.path.join(parent_dir, "*.xyz")
fnames = glob.glob(file_location)
directories = []
ref_directories = []
files = []
ref_files = []

for fname in range(0, len(fnames)):
    file_name = os.path.basename(fnames[fname])
    if (file_name != "Aref.xyz") and (file_name != "Bref.xyz"):
        files.append(file_name)
        file_name = file_name.split(".")
        directories.append(file_name[0])

for fname in range(0, len(fnames)):
    file_name = os.path.basename(fnames[fname])
    if (file_name == "Aref.xyz") or (file_name == "Bref.xyz"):
        ref_files.append(file_name)
        file_name = file_name.split(".")
        ref_directories.append(file_name[0])

directories = sorted(directories)    
ref_directories = sorted(ref_directories)    
files = sorted(files)   
ref_files = sorted(ref_files)   


# Printing the initial message to the output of this script:
print("""
                  ===================================================
                  |                    ASM_calc                     |           
                  ===================================================
                  |                                                 |
                  | Script written by: Lucas W. de Lima             |
                  | PhD student at the Applied Computational        |
                  | Chemistry Group - Institute of Chemistry - USP  |
                  | Advisor: Prof. Dr. Ataualpa A.C. Braga          |
                  |                                                 |
                  +-------------------------------------------------+
                               Last edition: 03/27/2022

        This script was written to calculate the Activation Strain Model (ASM)
    based on the approach mainly developed in the following articles:

    1 - Chem. Eur. J. 2016, 22, 4431-4439   (ASM for solvated systems using 
        implicity solvation models)
    2 - ACS Catal. 2018, 8, 975-986 (ASM using NEB images and transition state as
        the Minimum Energy Path (MEP) for the calculation)

        The energy expression or ASM in the methodology used is:

             DeltaE(R) = DeltaE_strain(R) + DeltaE_int(R) + DeltaE_solv(R)

    in which DeltaE(R) is the variation of energy in a point R of the MEP in 
    relation to the reference reactants and is decomposed into: 

        i - DeltaE_strain(R): variation of the strain energy in relation to the 
            reactant fragments, i.e.:

                DeltaE_strain(R) = (E(R)_A - Eref_A) + (E(R)_B - Eref_B)

       ii - DeltaE_int(R): the interaction energy between the fragments A and B;
      iii - DeltaE_solv(R): the variation of solvation energy between the solvated
            and vacuum structures, i.e.:

                    DeltaE_solv(R) = solvation_E(R) - solvation_E_ref
    \n
    """)

# Creating directories
print("""
  ------------------------------------
          Creating directories                            
  ------------------------------------   
""")
for i in range(0, len(directories)):
    place = os.path.join(parent_dir, directories[i])
    os.makedirs(place)
    origin = os.path.join(parent_dir, files[i])
    destination = os.path.join(parent_dir, directories[i], files[i])
    shutil.move(origin, destination)

for i in range(0, len(ref_directories)):
    place = os.path.join(parent_dir, ref_directories[i])
    os.makedirs(place)
    origin = os.path.join(parent_dir, ref_files[i])
    destination = os.path.join(parent_dir, ref_directories[i], ref_files[i])
    shutil.move(origin, destination)


print("  Created reference directories:")
for i in range(0, len(ref_directories)):
    print("  {}".format(ref_directories[i]))
print()

print("  Created MEP directories:")
for i in range(0, len(directories)):
    print("  {}".format(directories[i]))
print()
print()


# Moving to directories and performing calculations
print("""
  ------------------------------------
       Starting ASM calculations                            
  ------------------------------------   
""")

print("""
  ************************************
    Reference structure calculations                            
  ************************************   
""")
for i in range(0, len(ref_directories)):
    destination = os.path.join(parent_dir, ref_directories[i])
    os.chdir(destination)   
    print("  Moving to {} directory.\n".format(ref_directories[i])) 
    reference_file = open_xyz("{}".format(ref_files[i]))
    if ref_directories[i] == 'Aref':
        write_input(head_CO2, reference_file, "{}".format(ref_directories[i])) 
        write_input_solv(head_solv_CO2, reference_file, "{}".format(ref_directories[i])) 
        print("  Reference {} inputs created.\n".format(ref_directories[i]))

    elif ref_directories[i] != 'Aref':
        write_input(head, reference_file, "{}".format(ref_directories[i]))  
        write_input_solv(head_solv, reference_file, "{}".format(ref_directories[i])) 
        print("  Reference {} inputs created.\n".format(ref_directories[i]))

    current_dir = os.getcwd()
    ref_input_location = os.path.join(current_dir, "*.inp")
    ref_input_fnames = glob.glob(ref_input_location)
    ref_input_files = []

    for ref_inp_fname in range(0, len(ref_input_fnames)):
        ref_inp_file_name = os.path.basename(ref_input_fnames[ref_inp_fname])
        ref_inp_file_basename = ref_inp_file_name.split('.')
        ref_input_files.append(ref_inp_file_basename[0])

    for ref_inp in range(0, len(ref_input_files)):
        print("  Starting calculation of {}.".format(ref_input_files[ref_inp]))
        subprocess.run("{} {}.inp > {}.out".format(ORCA_dir, ref_input_files[ref_inp], ref_input_files[ref_inp]), shell=True)
        print("  {} calculation done.\n".format(ref_input_files[ref_inp]))
    
    os.chdir(parent_dir)
    
# Reference files energies
Aref_solv_path = os.path.join(os.getcwd(), 'Aref', 'Aref_solv.out')
Aref_solv_outfile = open_file(Aref_solv_path)
Aref_path = os.path.join(os.getcwd(), 'Aref', 'Aref.out')
Aref_outfile = open_file(Aref_path)
Bref_solv_path = os.path.join(os.getcwd(), 'Bref', 'Bref_solv.out')
Bref_solv_outfile = open_file(Bref_solv_path)
Bref_path = os.path.join(os.getcwd(), 'Bref', 'Bref.out')
Bref_outfile = open_file(Bref_path)

Aref_solv_energy = finalE(Aref_solv_outfile)
Aref_energy = finalE(Aref_outfile)
solvation_energy_Aref = Aref_solv_energy - Aref_energy
Bref_solv_energy = finalE(Bref_solv_outfile)
Bref_energy = finalE(Bref_outfile)
solvation_energy_Bref = Bref_solv_energy - Bref_energy

print("""
  ===================================================
             Summary of reference energies                            
  ===================================================   
  Aref energy (solvated): {} a.u.
  Aref energy: {} a.u.
  Aref solvation energy: {:.12f} a.u.

  Bref energy (solvated): {} a.u.
  Bref energy: {} a.u.
  Bref solvation energy: {:.12f} a.u.
 +---------------------------------------------------+

       ***** END OF REFERENCE CALCULATIONS *****

""".format(Aref_solv_energy, Aref_energy, solvation_energy_Aref, Bref_solv_energy, Bref_energy, solvation_energy_Bref))   


print("""
  ************************************
       MEP structures calculations                            
  ************************************   
""")
for i in range(0, len(directories)):
    destination = os.path.join(parent_dir, directories[i])
    os.chdir(destination)   
    print("\n  **** MOVING TO {} DIRECTORY **** \n".format(directories[i])) 
    AB = open_xyz("{}".format(files[i]))
    A = create_fragA(Atom_list, AB) 
    print("  Fragment A of {} created.".format(directories[i]))
    B = create_fragB(AB, A)
    print("  Fragment B of {} created.\n".format(directories[i]))
    write_input(head, AB, "{}_AB".format(directories[i])) 
    write_input_solv(head_solv, AB, "{}_AB".format(directories[i])) 
    write_input(head_CO2, A, "{}_A".format(directories[i])) 
    write_input(head, B, "{}_B".format(directories[i])) 

    current_dir = os.getcwd()
    input_location = os.path.join(current_dir, "*.inp")
    input_fnames = glob.glob(input_location)
    input_files = []

    for inp_fname in range(0, len(input_fnames)):
        inp_file_name = os.path.basename(input_fnames[inp_fname])
        inp_file_basename = inp_file_name.split('.')
        input_files.append(inp_file_basename[0])

    for inp in range(0, len(input_files)):
        print("  Starting calculation of {}.".format(input_files[inp]))
        subprocess.run("{} {}.inp > {}.out".format(ORCA_dir, input_files[inp], input_files[inp]), shell=True)
        print("  {} calculation done.\n".format(input_files[inp]))
    
    # Calculation of bond distance of selected atoms
    xyz_location = os.path.join(current_dir, "*.xyz")
    xyz_filename = glob.glob(xyz_location)
    xyz_fname = os.path.basename(xyz_filename[0])
    bond_distance = calculate_distance(xyz_fname, Atom_bond_list)
    
    # Calculation of energies in relation to the referencies
    AB_solv_output = "{}".format(directories[i]) + "_AB_solv.out"
    AB_solv_output_loc = os.path.join(current_dir, AB_solv_output)
    AB_output = "{}".format(directories[i]) + "_AB.out"
    AB_output_loc = os.path.join(current_dir, AB_output)
    A_output = "{}".format(directories[i]) + "_A.out"
    A_output_loc = os.path.join(current_dir, A_output)
    B_output = "{}".format(directories[i]) + "_B.out"
    B_output_loc = os.path.join(current_dir, B_output)

    AB_solv_outfile = open_file(AB_solv_output_loc) 
    AB_outfile = open_file(AB_output_loc) 
    A_outfile = open_file(A_output_loc)
    B_outfile = open_file(B_output_loc)

    AB_solv_energy = finalE(AB_solv_outfile)
    AB_energy = finalE(AB_outfile)
    A_energy = finalE(A_outfile)
    B_energy = finalE(B_outfile)
    
    # Solvation energy
    solvation_energy_AB = AB_solv_energy - AB_energy
    DeltaE_solv = (solvation_energy_AB - (solvation_energy_Aref + solvation_energy_Bref))*2625.5

    # Strain energy
    A_strain = A_energy - Aref_energy    
    B_strain = B_energy - Bref_energy
    DeltaE_strain = (A_strain + B_strain)*2625.5

    # DeltaE
    DeltaE = (AB_solv_energy - (Aref_solv_energy + Bref_solv_energy))*2625.5

    # Interaction energy
    DeltaE_int = DeltaE - DeltaE_solv - DeltaE_strain

    print("""
      ===================================================
                     Summary of MEP energies           
      ===================================================
      Structure: {}

      DeltaE: {:.2f} kJ·mol**-1
      DeltaE_strain: {:.2f} kJ·mol**-1
      DeltaE_int: {:.2f} kJ·mol**-1
      DeltaE_solv: {:.2f} kJ·mol**-1
      
      Bond distance of selected atoms: {:.3f} angstroms
     +---------------------------------------------------+

                           ********

    """.format(directories[i], DeltaE, DeltaE_strain, DeltaE_int, DeltaE_solv, bond_distance))   

    summary_path = os.path.join(parent_dir, "summary.dat")
    with open(summary_path, "a+") as f:
        f.write("{} {} {} {} {}\n".format(bond_distance, DeltaE, DeltaE_strain, DeltaE_int, DeltaE_solv, DeltaE_solv))
    
    print("  MEP energies written at summary.dat file.\n") 

    os.chdir(parent_dir)


endTime = datetime.now()
totalTime = endTime - startTime
n = totalTime.total_seconds() 

day = n // (24 * 3600) 
n = n % (24 * 3600) 
hour = n // 3600 
 
n %= 3600 
minutes = n // 60 
 
n %= 60 
seconds = n 

print("  Done!")
print("  Total time of execution: {} days, {} hours, {} minutes, {} seconds.".format(day, hour, minutes, seconds))
