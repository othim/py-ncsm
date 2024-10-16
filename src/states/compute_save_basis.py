import numpy as np
import pywigxjpf as wig
import os
import sys
os.sys.path.append("../basis/")
import basis_states as bs
import antisym_3N as a3n
os.sys.path.append("../external/HO_brackets/")
import gmosh
from itertools import groupby
from operator import itemgetter
from scipy.linalg import block_diag
from scipy.linalg import eig
import time
verbose = False

# Pre-populate arrays in FORTRAN code
gmosh.angmom_lib.precalculate_binomials() 
gmosh.factorials.precalculate_factorials() 
gmosh.angmom_lib.setup_moshinsky() 

# Setup wigner symbols
wig.wig_table_init(2*100, 9)
wig.wig_temp_init(2*100)

# Construct alpha basis (partially antisymmetric) states
Nmax = 40
print(f'Constructing states_alpha...')
states_alpha = bs.NNN_basis_nl(Nmax,verbose)
print(f'Done! Nmax={Nmax}, len(states_alpha)={len(states_alpha)}')


# Group states with same J2, T2 and pi (parity)
print(f'Grouping states with same (J2, T2, pi)...')
grouper = itemgetter("J2","T2", "pi") 
JTpi_list, key_list = bs.group_NNN_basis_nl(grouper, states_alpha, verbose=verbose)
print(f'Done! len={len(JTpi_list)}')

# Take the group with correct quantum numbers
J2 = 1
T2 = 1
pi = 1
print_states = False

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

njt_block_time_list = []

for j,njt_block in enumerate(NJT_list):
    print(f'NJT_block=({j+1}/{len(NJT_list)})')
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
    ones_idx,= np.where(np.abs(eigs-1)<1e-10)
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
    exponent = 2.1
    C=(end-start)*1000/(njt_block_size_list[j])**exponent
    print(f'C={C}\n\n')
    T_left   = int(C*np.sum(np.array(njt_block_size_list[j+1:])**exponent)/1000) # s
    T_left_h = int(T_left/3600)
    T_left_m = int((T_left-T_left_h*3600)/60)
    T_left_s = T_left-T_left_h*3600-T_left_m*60
    print(f'Time left={T_left_h:<2.0f}h {T_left_m:<2.0f}m {T_left_s:<2.0f}s\n\n')

# Construct the full change of-basis matrix Gamma -> alpha
# The states are Gamma are in the list NJTi_states and the alpha
# states are in the list JTpi_block
Gamma_to_alpha = block_diag(*eig_Gamma_list)

print(Gamma_to_alpha.shape)
print(len(NJTi_states))
print(len(JTpi_block))
for s in NJTi_states:
    print(s)

print(njt_block_time_list)
