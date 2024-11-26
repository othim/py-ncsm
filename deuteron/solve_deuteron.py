'''
    Script to setup the 2N system and diagonalize the Hamiltonian
    in a specific partiel wave.
'''

import numpy as np
import os
from itertools import groupby
from operator import itemgetter

os.sys.path.append("../src/states/")
os.sys.path.append("../src/hamiltonian/")
import setup_Hamiltonian as sH
import basis_states as bs
import load_potential as lp

# From plane wave basis (idaho-n3lo)
# ---------------------
# E =  -2.22458 MeV ***
# _____________________


# *****************************************************************************
# ********************************** ARGUMENTS ********************************
# *****************************************************************************

Omega = 24   # MeV
mN    = 0  # MeV (not needed)
isospin_sym = True
fast        =  False
# Deuteron channel
j  = 1
t  = 0
pi = 1

#interaction_file="../interactions/idaho_n3lo_nmax_40_hw_24.txt"
interaction_file="/Users/toliver/Documents/phd/projects/bayes_dwb/tmp/full_LO.txt"

#Nmax_arr = [0,2,4,6,8,10,12,14,16,20,24,28,30,32,34,36,38,40]
Nmax_arr = [0,2,4,6,8,10,12,14,16,20,24,28]

# *****************************************************************************
# *****************************************************************************

E_arr = []
for Nmax in Nmax_arr:
    
    # Load potential matrix elements
    print(f'Loading potential matrix elements...')
    pot_dict,pot = lp.load_potential_file(interaction_file)
    print(f'Done!')

    # Create NN basis
    alpha_2N = bs.NN_basis_nl(Nmax,False)
    
    # Select only deuteron channel
    grouper = itemgetter("j","t", "pi") 
    JTpi_list, key_list = bs.group_NNN_basis_nl(grouper, alpha_2N, verbose=False)
    print([(j,t,pi) == a for a in key_list])

    # Take the group with correct quantum numbers
    idx = [(j,t,pi) == a for a in key_list].index(True)
    print(f'Selecting block: (j,t,pi)={key_list[idx]}')
    JTpi_block = JTpi_list[idx]
    
    # Print that block
    print(f'JTpi_block: len(JTpi_block)={len(JTpi_block)}')
    for i,s in enumerate(JTpi_block):
        print(s)

    # Setup Hamiltonian and diagonalize
    H_matrix_Gamma_basis = sH.setup_H_alpha_basis(JTpi_block,Omega,mN,pot_dict,\
            isospin_sym,fast)
    eigs,eigv = np.linalg.eigh(H_matrix_Gamma_basis)
    E_arr.append(np.min(eigs))
    print(f'\nnp.min(eigs) = {np.min(eigs):.4f} MeV\n\n')

# Print result
for i,E in enumerate(E_arr):
    print(f'Nmax={Nmax_arr[i]:<4}, E={E:<.4f}') 
