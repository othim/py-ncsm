'''
    basis_states.py
    ---------------
    
    This file contains some methods to construct and manipulate
    Harmonic Oscillator (HO) states.

    Oliver Thim (2024)
'''
import time
from itertools import groupby


def psi_nlm(n,l,m,r,theta,phi):
    '''
        Computes the single-particle HO wave function as defined in (...)

        Args:


        Returns:
    '''
    return 0


def NN_basis_j(j_min, j_max, mt_min, mt_max):
    '''
        Computes allowed antisymmetrized 2N basis states
        and returns a list of dictionaries containig the basis states.
    Args:
        j_min  : Minimum total angular momentum, (ls)j.
        j_max  : Maximum total angular momentum, (ls)j.
        mt_min : Minimum isospin projection 
        mt_max : Maximum isospin projection
    Returns:
        list(basis_states)
        basis_states : dictionary(l,s,j,t,mt,pi) (pi=parity)
    '''
    basis = []
    for mt in range(mt_min,mt_max+1,1):
        for J in range(j_min,j_max+1,1):
            for S in range(0,2):
                for L in range(abs(J-S),J+S+1,1):
                    for T in range(abs(mt),2,1):
                        if ((L+S+T)%2 != 0):
                            basis_state = {}
                            basis_state['l']  = L
                            basis_state['s']  = S
                            basis_state['j']  = J
                            basis_state['t']  = T
                            basis_state['mt'] = mt
                            basis_state['pi'] = (-1)**L
                            print(basis_state)
                            basis.append(basis_state)

    print('len(basis) = ',len(basis))
    return basis

def triag(L, S, J):
    '''
        if (|L-S| <= J <= |L+S|) returns True, else False.
    '''
    if( J < abs(L - S) or J>L+S):
        return False
    else:
        return True

def test_basis_state(basis_state):
    l = basis_state['l']
    s = basis_state['s']
    t = basis_state['t']
    j = basis_state['j']
    cN = basis_state['cN']
    cL = basis_state['cL']
    cJ2 = basis_state['cJ2']
    J2 = basis_state['J2']
    T2 = basis_state['T2']
    
    # Check 2N Pauli principle
    if (l+s+t)%2 == 0:
        return False

    # Check spin couplings
    if triag(l,s,j) == False:
        return False
    if triag(cL,1/2,cJ2/2) == False:
        return False
    if triag(j,cJ2/2,J2/2) == False:
        return False

    # Check isospin coupling
    if triag(t,1/2,T2/2) == False:
        return False
    
    # If nothing failed, return True
    return True


def NNN_basis_nl(N_max, verbose=False):
    '''
        Computes 3N basis states which are antisymmetrized in particle 1 and 2, 
        i.e. (-1)^{l+s+t} = -1. 

        n,l,s,j,t       : quantum numbers for particle (1,2) subsystem.
        cN,cL,cJ2       : quantum numbers for particle 3 relative to (1,2) subsystem
        J2, T2          : Total J and T for 3N system
        N               : 2*n+l+2*cN+cL 

        NOTE: the 2 in cJ2,J2 and T2 denote that these are two times the 
        respective quantum numbers as this is more convenient to save since 
        these quantum numbers are half integers.
    Args:
        N_max (int)             : Maximum HO energy, N = 2n+l + 2*cN+cL (cN = \mathcal{N})
    Returns:
        basis (list) : list of 'basis_states (dictionary)' that contains keys : 
        dictionary(n,l,s,j,t,mt,cN,cL,cJ2,J2,T2,N)
    '''
    basis = []
    for n in range(0,int(N_max/2)+1,1):  # n = 0,...,int(Nmax/2)
        for l in range(0,N_max-2*n+1,1): # l = 0,...,Nmax-2*n
            for s in range(0,2,1):       # s = 0,1
                for t in range(0,2,1):   # t = 0,1
                    if ((l+s+t)%2 != 0): # continue only for antisym 2N
                        
                        for j in range(abs(l-s),l+s+1,1): # j=|l-s|,...,|l+s|
                            # cN = 0,...,int((Nmax-2n-l)/2)
                            for cN in range(0,int((N_max-2*n-l)/2)+1,1):
                                # cL = 0,...,int(Nmax-2n-l-2*cN)
                                for cL in range(0,N_max-2*n-l-2*cN+1,1):

                                    # The two means that this is two times
                                    # the quantum number. This is more convenient
                                    # since these can be half integers.

                                    # cJ2 = |2*cL-1|,...,|2*cL+1|
                                    for cJ2 in range(abs(2*cL-1),2*cL+2,2):
                                        # J2 = |2*j-cJ2|,...,|2*j+cJ2|
                                        for J2 in range(abs(2*j-cJ2),2*j+cJ2+1,2):
                                            # T2 = |2*t-1|,...,|2*t+1|
                                            for T2 in range(abs(2*t-1),2*t+2,2):
                                                basis_state = {}
                                                basis_state['n']   = int(n)
                                                basis_state['l']   = int(l)
                                                basis_state['s']   = int(s)
                                                basis_state['j']   = int(j)
                                                basis_state['t']   = int(t)
                                                
                                                basis_state['cN']  = int(cN)
                                                basis_state['cL']  = int(cL)
                                                basis_state['cJ2'] = int(cJ2)
                                                
                                                basis_state['J2']  = int(J2)
                                                basis_state['T2']  = int(T2)
                                                basis_state['N']   = \
                                                2*int(n)+int(l)+2*int(cN)+int(cL)
                                                basis_state['pi'] = int((-1)**(l+cL))

                                                basis.append(basis_state)

                                                # Test is basis satisfy the basic
                                                # criteria.
                                                if test_basis_state(basis_state) == False:
                                                    print(f'Error constructing basis states')
                                                    return 0


    if verbose:
        for i,b in enumerate(basis):
            print(b)
        print('len(basis) = ',len(basis))
    return basis

def group_NNN_basis_nl(grouper, basis, verbose=False):
    '''
        Groups the states in the keys defined in items

    Args:
        grouper (Object)  : Grouper for which keys that will be grouped. E.g.
                            itemgetter("N","J2", "T2") 
        basis (list) : list of basis states as computed with 
                       basis_states.NN_basis_nl()

    Returns:
        grouped_basis (list) : List of groups of basis states with the same
                               keys in items
    '''

    grouped_basis = []
    key_list = []
    for key, grp in groupby(sorted(basis, key = grouper), grouper):
        grouped_basis.append(list(grp))
        key_list.append(key)

    if verbose:
        for chn_idx, chn in enumerate(grouped_basis):
            print(f'chn_idx={chn_idx}')
            for state in chn:
                print(state)
            print('')
    return grouped_basis,key_list

def NNN_basis_nl_check_equal(list1,list2, verbose=False):
    ''' 
        Check if two lists of thates contains the same states
    '''

        # Sort both lists of dictionaries for comparison
    sorted_list1 = sorted(list1, key=lambda x: tuple(sorted(x.items())))
    sorted_list2 = sorted(list2, key=lambda x: tuple(sorted(x.items())))

    # Compare the sorted lists
    return sorted_list1 == sorted_list2

##    states_A = copy(
#   for i,sA in enumerate(states_A):
#        for j, sB in enumerate(states_B):
#            if (sA == sB):



# If you run only this file you can check some functions
if __name__ == "__main__":
    start = time.time()
    NNStates = NN_basis_j(0,10,0,0)
    end = time.time()
    print(f'NN_basis_j: {(end-start)*1e3:.4f} ms')
    
    start = time.time()
    Nmax = 30
    NNStates = NNN_basis_nl(Nmax,True)
    end = time.time()
    print(f'NN_basis_nl, Nmax={Nmax}, time={(end-start)*1e3:.4f} ms')
