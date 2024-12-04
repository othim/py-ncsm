'''
    Script to setup the 2N system and diagonalize the Hamiltonian
    in a specific (j,t,pi) channel.

    Oliver Thim, 
    Chalmers 2024
'''

import numpy as np
import os
import sys
from itertools import groupby
from operator import itemgetter

# Append the path to py-ncsm
sys.path.append(os.path.abspath('../src/states/'))
sys.path.append(os.path.abspath('../src/hamiltonian/'))
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


# *****************************************************************************
# *****************************************************************************

def diag_2N(Nmax,hw,isospin_sym,fast,j,t,pi,interaction_file,PRINT=True):
    
    # Load potential matrix elements
    if PRINT:
        print(f'Loading potential matrix elements...')
    pot_dict, pot = lp.load_potential_file(interaction_file)

    # Create NN basis
    alpha_2N = bs.NN_basis_nl(Nmax,False)
    
    # Select only deuteron channel
    grouper = itemgetter("j","t", "pi") 
    JTpi_list, key_list = bs.group_NNN_basis_nl(grouper, alpha_2N, verbose=False)
    #if PRINT:
    #    print([(j,t,pi) == a for a in key_list])

    # Take the group with correct quantum numbers
    idx = [(j,t,pi) == a for a in key_list].index(True)
    if PRINT:
        print(f'Selecting block: (j,t,pi)={key_list[idx]}')
    JTpi_block = JTpi_list[idx]
    
    # Print that block
    if PRINT:
        print(f'JTpi_block: len(JTpi_block)={len(JTpi_block)}')
        for i,s in enumerate(JTpi_block):
            print(s)

    # Setup Hamiltonian and diagonalize
    H_matrix_Gamma_basis = sH.setup_H_alpha_basis(JTpi_block,hw,pot_dict,\
            isospin_sym,fast)
    return np.linalg.eigh(H_matrix_Gamma_basis)



if __name__=='__main__':

    interaction_file = "../interactions/idaho_n3lo_nmax_40_hw_24_Np_80_finite.txt"
    Nmax_arr         = [0,2,4,6,8,10,12,14,16,20,24,28,32,36,40]
    hw               = 24   # MeV
    isospin_sym      = True
    fast             = False

    # Deuteron channel
    j  = 1
    t  = 0
    pi = 1

    E_arr = []
    for Nmax in Nmax_arr:
        eigs,eigv = diag_2N(Nmax,hw,isospin_sym,fast,j,t,pi,interaction_file)
        E_arr.append(np.min(eigs))
        print(f'\nnp.min(eigs) = {np.min(eigs):.4f} MeV\n\n')

    # Print result
    for i,E in enumerate(E_arr):
        print(f'Nmax={Nmax_arr[i]:<4}, E={E:<.4f}') 
