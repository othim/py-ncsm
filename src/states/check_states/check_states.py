'''
    Script to check if the code generate the corest set of basis states.
'''

import numpy as np
import os
os.sys.path.append("../basis/")
import basis_states as bs

from itertools import groupby
from operator import itemgetter

file = "antisym_out_Andreas.txt"
states_A_txt = np.loadtxt(file,skiprows=18,dtype='str')
print(states_A_txt.shape)



states_A = []
for i in range(len(states_A_txt[:,0])):
    basis_state = {}
    row = states_A_txt[i,:]

    basis_state = {}
    basis_state['n']   = int(row[8])
    basis_state['l']   = int(int(row[9])/2)
    basis_state['s']   = int(int(row[10])/2)
    basis_state['j']   = int(int(row[11])/2)
    basis_state['t']   = int(int(row[12])/2)
    
    basis_state['cN']  = int(row[15])
    basis_state['cL']  = int(int(row[16])/2)
    basis_state['cJ2'] = int(row[18])
    
    basis_state['J2']  = int(row[4])
    basis_state['T2']  = int(row[5])
    basis_state['N']   = int(row[3])
    if row[6] == '+':
        basis_state['pi'] = 1
    else:
        basis_state['pi'] = -1
    
    print(row)
    print(basis_state)
    states_A.append(basis_state)



grouper = itemgetter("N","J2","T2") 
group_A = bs.group_NNN_basis_nl(grouper, states_A, verbose=True)



print(bs.NNN_basis_nl_check_equal(group_A[0],group_A[0], verbose=False))


# Construct states with my method

Nmax = 2 
states = bs.NNN_basis_nl(Nmax,True)
grouper = itemgetter("N","J2","T2") 
NJT_list = bs.group_NNN_basis_nl(grouper, states, verbose=True)


comp1 = group_A[1]
comp2 = NJT_list[9]
print(comp1)
print(comp2)
print(bs.NNN_basis_nl_check_equal(comp1,comp2, verbose=False))


