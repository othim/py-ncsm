'''
    compute_save_basis.py
    --------------------
    
    This program computes and saves the basis states and transformation matrix.

    Oliver Thim (2024)
'''
import numpy as np
import pywigxjpf as wig
import os
import sys
import json
os.sys.path.append("../basis/")
import basis_states as bs
import antisym_3N as a3n
import gmosh
from itertools import groupby
from operator import itemgetter
from scipy.linalg import block_diag
import time
from datetime import datetime
import argparse

def hms_from_s(t_s):
    ''' 
        Returns hours, minutes and seconds for a 
        given number of seconds t_s (int)
    '''
    T_h = int(t_s/3600)
    T_m = int((t_s-T_h*3600)/60)
    T_s = t_s-T_h*3600-T_m*60
    return T_h,T_m,T_s

# Define a custom argument type for a list of strings
def list_arg(arg):
    return arg.split(',')

def to_int(str_list):
    int_list = []
    for s in str_list:
        int_list.append(int(s))
    return int_list
# *****************************************************************************
# *****************************************************************************
# ********************************** MAIN *************************************
# *****************************************************************************

# This program computes basis states for a given (J,T,pi) block and compute
# the basis change matric going from the alpha to Gamma basis.

# *****************************************************************************
# ******************************** VARIABLES **********************************
verbose = True

J2 = 1
T2 = 1
pi = 1
print_states = True
one_tol = 1e-6

# Specify the directory for the data 
directory_name = f'Nmax={Nmax}_data'
save_data = True

# ***********************************

# Parse and set the NMAX_ARR argument
parser = argparse.ArgumentParser(prog='compute_save_basis')
parser.add_argument('--NMAX',type=list_arg) 
args = parser.parse_args()
NMAX_ARR = to_int(args.NMAX)
print(f'NMAX={NMAX_ARR}')

# *****************************************************************************
# *****************************************************************************

start_prog = time.time()
# Pre-populate arrays in FORTRAN code. This is necessary for calling the function
# gmosh later.
gmosh.angmom_lib.precalculate_binomials() 
gmosh.factorials.precalculate_factorials() 
gmosh.angmom_lib.setup_moshinsky() 

# Setup wigner symbols
wig.wig_table_init(2*100, 9)
wig.wig_temp_init(2*100)

for Nmax in NMAX_ARR:
    # Construct alpha basis (partially antisymmetric) states
    print(f'Constructing states_alpha...')
    states_alpha = bs.NNN_basis_nl(Nmax,verbose)
    print(f'Done! Nmax={Nmax}, len(states_alpha)={len(states_alpha)}')


    # Group states with same J2, T2 and pi (parity)
    print(f'Grouping states with same (J2, T2, pi)...')
    grouper = itemgetter("J2","T2", "pi") 
    JTpi_list, key_list = bs.group_NNN_basis_nl(grouper, states_alpha, verbose=verbose)
    print(f'Done! len={len(JTpi_list)}')

    # Take the group with correct quantum numbers
    idx = [(J2,T2,pi) == a for a in key_list].index(True)
    print(f'Selecting block: (J2,T2,pi)={key_list[idx]}')
    JTpi_block = JTpi_list[idx]

    print(f'JTpi_block: len(JTpi_block)={len(JTpi_block)}')
    if print_states:
        for i,s in enumerate(JTpi_block):
            print(s)


    # Group states according to (N,J2,T2) in this block
    grouper = itemgetter("N","J2","T2") 
    NJT_list, NJT_key_list = bs.group_NNN_basis_nl(grouper, JTpi_block, verbose=False)

    print(f'NJT_list info:')
    print(f'len(NJT_list)={len(NJT_list)}')
    print(f'Matrix sizes:')
    njt_block_size_list = []
    for i,njt_block in enumerate(NJT_list):
        block_size = len(njt_block)
        njt_block_size_list.append(block_size)
        print(f'idx={i}, N={njt_block[0]["N"]}, size={block_size}')
    print("")

        

    # The antisymmetrization matrix is diagonal in these blocks
    eig_Gamma_list = []
    NJTi_states    = []
    alpha_states   = []

    njt_block_time_list = []

    for j,njt_block in enumerate(NJT_list):
        print(f'NJT_block=({j+1}/{len(NJT_list)})')
        #for tmp in njt_block:
        #    print(tmp)
        start = time.time()
        s_A = time.time()
        A = a3n.setup_3N_antisymmetrizer_matrix(njt_block,False)
        e_A = time.time()
        if (a3n.check_A(A) == False):
            print(f'Error in setting up antisymmetrizer, A not projector.')
            sys.exit(1)

        s_e = time.time()
        eigs,eigv = np.linalg.eigh(A)
        e_e = time.time()
        #eigs,eigv = eig(A)
        if (a3n.check_eigs(eigs) == False):
            print(f'Error in setting up antisymmetrizer, eigs.')
            sys.exit(1)
        
        # Check itf eigv is ortogonal
        print(f'eigv ortogonal: {np.max(np.abs(eigv.T-np.linalg.inv(eigv)))}')

        # Keep only eigenvectors with eigenvalue 1
        ones_idx,= np.where(np.abs(eigs-1)<one_tol)
        # eig_Gamma[:,i] is the coefficients for |NJT,i> expressed in the alpha
        # basis
        eig_Gamma = np.copy(eigv[:,ones_idx])
        eig_Gamma_list.append(eig_Gamma)
        # Construct the corresponding states
        for i in range(len(ones_idx)):
            basis_state = {}
            basis_state['N'] = njt_block[i]['N']
            basis_state['J2'] = njt_block[i]['J2']
            basis_state['T2'] = njt_block[i]['T2']
            basis_state['i'] = i
            NJTi_states.append(basis_state)
        
        # Add the alpha states, this is important to get the states in the 
        # correct order
        for i,alpha in enumerate(njt_block):
            alpha_states.append(alpha)

        end = time.time()
        njt_block_time_list.append((end-start)*1000) # ms
        print(f'eig_Gamma shape: {eig_Gamma.shape}')
        print(f'max |A-A^2| = {np.max(np.abs(A-A@A))}')
        print(f'# eigs=1: {np.sum([np.abs(e-1.0)<1e-10 for e in eigs])}')
        print(f'# eigs=0: {np.sum([np.abs(e-0.0)<1e-10 for e in eigs])}')
        print(f'len(eigs)={len(eigs)}\n')
        print(f'Time A   : {(e_A-s_A):<6.3f} s')
        print(f'Time eigs: {(e_e-s_e):<6.3f} s')
        print(f'Time tot : {(end-start):<6.3f} s')
        # Estimate the ammount of time left
        exponent = 2.1
        C=(end-start)*1000/(njt_block_size_list[j])**exponent
        print(f'C={C}\n\n')
        T_left   = int(C*np.sum(np.array(njt_block_size_list[j+1:])**exponent)/1000) # s
        T_left_h,T_left_m,T_left_s = hms_from_s(T_left)
        print(f'Time left={T_left_h:<2.0f}h {T_left_m:<2.0f}m {T_left_s:<2.0f}s\n\n')

    # Construct the full change of-basis matrix Gamma -> alpha
    # The states are Gamma are in the list NJTi_states and the alpha
    # states are in the list JTpi_block
    Gamma_to_alpha = block_diag(*eig_Gamma_list)

    print(f'Gammta_to_alpha.shape={Gamma_to_alpha.shape}')
    print(f'len(NJTi_states)={len(NJTi_states)}')
    print(f'len(JTpi_block)={len(JTpi_block)}')
    if print_states:
        for s in NJTi_states:
            print(s)

        for s in alpha_states:
            print(s)

    if save_data:
        # Create if it fails, instead create a directory where the current date
        # and time is appenden in the name
        try:
            os.mkdir(directory_name)
            print(f"Directory '{directory_name}' created successfully.")
        except FileExistsError:
            print(f"Directory '{directory_name}' already exists.")
            directory_name += str(datetime.now())
            os.mkdir(directory_name)
            print(f"Directory '{directory_name}' created successfully.")

        try:
            with open(directory_name + '/Gamma_basis_Nmax='+str(Nmax)+'.txt','w') as data:
                json.dump(NJTi_states,data)
            with open(directory_name + '/alpha_basis_Nmax='+str(Nmax)+'.txt','w') as data:
                json.dump(alpha_states,data)

            with open(directory_name+'/Gamma_basis_Nmax='+str(Nmax)+'.txt', "r") as fp:
                person_dict = json.load(fp)
            with open(directory_name+'/Gamma_to_alpha_Nmax='+str(Nmax)+'.txt','w') as file:
                np.savetxt(file,Gamma_to_alpha)
                #np.save(file,Gamma_to_alpha)
        except:
            print('Error, data not saved successfully')

end_prog = time.time()

T_left_h,T_left_m,T_left_s = hms_from_s(int(end_prog-start_prog))
print(f'Program ended in: {T_left_h:<2.0f}h {T_left_m:<2.0f}m {T_left_s:<2.0f}s')
