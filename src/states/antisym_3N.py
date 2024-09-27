'''
    antisym_3N.py
    -------------

    This file contains functions to setup the antisymetrization operator
    in the partially antisymmetric 3N basis to construct a fully antisymmetric
    basis 3N basis.

    Oliver Thim (2024)
'''
import numpy as np
import pywigxjpf as wig
import ctypes
import os
os.sys.path.append("../basis/")
import basis_states as bs


def triag(L, S, J):
    '''
        if (|L-S| <= J <= |L+S|) returns True, else False.
    '''
    if( J < abs(L - S) or J>L+S):
        return False
    else:
        return True

def hat(j):
    return np.sqrt(2*j+1)



def get_hat_factors(bra,ket,L,S2):
    '''
    Args:
        bra (dictionary(n,l,s,j,t,mt,pi,cN,cL,cJ2,J2,T2)) : bra-quantum state
        as defined in basis_states.py.

        ket (dictionary(n,l,s,j,t,mt,pi,cN,cL,cJ2,J2,T2)) : ket-quantum state
        L  (int) : coupled total L 
        S2 (int) : 2*S, where S is the total coupled spin
    Returns:
        All hat-factors in eq. (XX) for the given bra and ket quantum numbers
    '''


    return hat(bra['t'])*hat(ket['t'])*(hat(L)**2)*(hat(S2/2)**2)*hat(bra['j'])*hat(ket['j'])*\
                hat(bra['cJ2']/2)*hat(ket['cJ2']/2)*hat(bra['s'])*hat(ket['s'])

def get_9j_factors(bra,ket,L,S2):
    '''
    Args:
        bra (dictionary(n,l,s,j,t,mt,pi,cN,cL,cJ2,J2,T2)) : bra-quantum state
        as defined in basis_states.py.

        ket (dictionary(n,l,s,j,t,mt,pi,cN,cL,cJ2,J2,T2)) : ket-quantum state
        L (int)  : coupled total L 
        S2 (int) : 2*S, where S is the total coupled spin
    Returns:
        The two 9j symbols in eq. (XX) for the given bra and ket quantum numbers
    '''
    if (bra['J2'] != ket['J2']):
        print("Error, bra.J2 != ket.J2")
        return 0
    
    fac1 = wig.wig9jj(2*bra['l'],  2*bra['s'], 2*bra['j'],\
                      2*bra['cL'], int(2*1/2),   bra['cJ2'],\
                      2*L,      S2,     bra['J2'])
    
    fac2 = wig.wig9jj(2*ket['l'],  2*ket['s'], 2*ket['j'],\
                      2*ket['cL'], int(2*1/2),   ket['cJ2'],\
                      2*L,      S2,     ket['J2'])
    return fac1*fac2

def get_6j_factors(bra,ket,L,S2):
    '''
    Args:
        bra (dictionary(n,l,s,j,t,mt,pi,cN,cL,cJ2,J2,T2)) : bra-quantum state
        as defined in basis_states.py.

        ket (dictionary(n,l,s,j,t,mt,pi,cN,cL,cJ2,J2,T2)) : ket-quantum state
        L  (int) : coupled total L 
        S2 (int) : 2*S, where S is the total coupled spin
    Returns:
        The two 6j symbols in eq. (XX) for the given bra and ket quantum numbers
    '''
    if (bra['T2'] != ket['T2']):
        print("bra.T2 != ket.T2")
        return 0
    facT = wig.wig6jj(int(2*1/2),  int(2*1/2), 2*bra['t'],\
                      int(2*1/2), bra['T2'], 2*ket['t'])
    
    facS = wig.wig6jj(int(2*1/2),  int(2*1/2), 2*bra['s'],\
                      int(2*1/2),  S2,   2*ket['s'])
    return facT*facS

def get_HO_braket(bra,ket,L,S2):
    '''
    Args:
        bra (dictionary(n,l,s,j,t,mt,pi,cN,cL,cJ2,J2,T2)) : bra-quantum state
        as defined in basis_states.py.

        ket (dictionary(n,l,s,j,t,mt,pi,cN,cL,cJ2,J2,T2)) : ket-quantum state
        L  (int) : coupled total L 
        S2 (int) : 2*S, where S is the total coupled spin
    Returns:
        The general HO braket
    '''
    

def setup_3N_antisymmetrizer_matrix(basis,verbose=False):
    '''
        Computes the antisymmetrizer matrix as defined in eq. (XX) in the notes.

    Args:

        basis (list) : list of basis states as computed with 
                       basis_states.NN_basis_nl()
        j2max (int)  : 2*j_max for the wigner symbols

    Returns:
        antisym_matrix (np.array(len(basis),len(basis)))
    '''
       
    A = np.zeros((len(basis),len(basis)),dtype=complex)
    
    for i,bra in enumerate(basis):
        for j,ket in enumerate(basis):
            matrix_el=0
            
            braN  = 2*bra['n']+bra['l'] + 2*bra['cN']+bra['cL']
            ketN = 2*ket['n']+ket['l'] + 2*ket['cN']+ket['cL']

            # Implement delta_{braN,ketN} and check if bra.(J,T) == ket.(J,T)
            if (braN == ketN and bra['J2'] == ket['J2'] and bra['T2'] == ket['T2']): 
                # Sum over L and S
                for L in range(abs(bra['l']-bra['cL']),bra['l']+bra['cL']+1):
                    # Implement other restrictions for L
                    if (triag(ket['l'],ket['cL'],L)):
                        
                        # Note that S2 = 2*S where S is the total spin
                        for S2 in range(abs(2*bra['s']-1),2*bra['s']+2,2):
                            # Implement other restrictions for S
                            if (triag(bra['s'],1/2,S2/2) and triag(bra['J2']/2,L,S2/2)): 

                                # All the hat-factors
                                fac1 = get_hat_factors(bra,ket,L,S2)
                                
                                # Exponential of L
                                fac2 = (-1)**(L)
                                #fac2 = (-1)**(1+int(S2/2)+int(bra['T2']/2))
                                
                                # 9j factors
                                fac3 = get_9j_factors(bra,ket,L,S2)

                                # 6j factors
                                fac4 = get_6j_factors(bra,ket,L,S2)

                                # General HO bracket for two particles
                                fac5 = 1 #get_HO_braket(bra,ket,L,S2)
                        
                                matrix_el += fac1*fac2*fac3*fac4*fac5

            else:
                matrix_el = 0

            # Set calculated matrix element
            #A[i,j] = matrix_el

            A[i,j] = (-2/3)*matrix_el
            if (i==j):
                A[i,j] += 1/3
    
    if verbose:
        print(f'A.shape={A.shape}')
        print(f'A.nonzero={np.count_nonzero(A)}')
    return A


if __name__=="__main__":
    
    # Setup wigner symbols
    wig.wig_table_init(2*100, 9)
    wig.wig_temp_init(2*100)
    # Construct a basis
    Nmax = 0
    basis = bs.NNN_basis_nl(Nmax,-1,1,False)

    # Compute antisummetrizer matrix
    A1 = setup_3N_antisymmetrizer_matrix(basis,True)
    print(A1)
    
    
