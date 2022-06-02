#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt

# Ploting ASM graph
plot_name = os.path.join(os.getcwd(), "ASM_plot.png")
summary_file_location = os.path.join(os.getcwd(), "summary.dat")

data = np.loadtxt(summary_file_location)
data = data[data[:, 0].argsort()]        

bond_distance = data[:, 0]
DeltaE = data[:, 1]
DeltaE_strain = data[:, 2]
DeltaE_int = data[:, 3]
DeltaE_solv = data[:, 4]

plt.plot(bond_distance, DeltaE, label="DeltaE")
plt.plot(bond_distance, DeltaE_strain, label="DeltaE_strain")
plt.plot(bond_distance, DeltaE_int, label="DeltaE_int")
plt.plot(bond_distance, DeltaE_solv, label="DeltaE_solv")
plt.xlim(left=(bond_distance[-1]+0.1), right=(bond_distance[0]-0.1))    
plt.xlabel(r"Bond distance / $\AA$")
plt.ylabel(r"Energy / $\it{kJ·mol^{-1}}$")
plt.legend()
plt.savefig(plot_name)

print("  ASM plot written at ASM_plot.png file.\n") 
print("  Done!")
