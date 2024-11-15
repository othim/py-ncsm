'''
    antisym_3N.py
    -------------

    This file contains functions to setup the antisymetrization operator
    in the partially antisymmetric 3N basis to construct a fully antisymmetric
    basis 3N basis.

    This file requires the module f2py-created module 'gmosh'. 

    Oliver Thim (2024)
'''
import numpy as np
import pywigxjpf as wig
import os
os.sys.path.append("../basis/")
import basis_states as bs
os.sys.path.append("../external/HO_brackets/")
import gmosh

from itertools import groupby
from operator import itemgetter

#import test_mymod as tm

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

    return hat(bra['t'])*hat(ket['t'])*(hat(L)**2)*(hat(S2/2.0)**2)*hat(bra['j'])*hat(ket['j'])*\
                hat(bra['cJ2']/2.0)*hat(ket['cJ2']/2.0)*hat(bra['s'])*hat(ket['s'])

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
    
    fac1 = wig.wig9jj(2*bra['l'],  2*bra['s'],  2*bra['j'],\
                      2*bra['cL'], int(2*1/2),  bra['cJ2'],\
                      2*L,         S2,          bra['J2'])
    
    fac2 = wig.wig9jj(2*ket['l'],  2*ket['s'],  2*ket['j'],\
                      2*ket['cL'], int(2*1/2),  ket['cJ2'],\
                      2*L,         S2,          ket['J2'])
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

def setup_3N_antisymmetrizer_element(bra,ket,verbose=False):
    '''
        Computes the antisymmetrizer matrix elements as defined in 
        Eq. (31) in Phys. Rev. C 61 044001.

    Args:

        bra (basis_state) : basis state as computed with 
                            basis_states.NN_basis_nl()
        
        ket (basis_state) : basis state as computed with 
                            basis_states.NN_basis_nl()
    Returns:
        matrix_el (float) 
    '''
    matrix_el = 0

    for S2 in [1,3]:
        for LL in range(abs(int((bra['J2']-S2))),int((bra['J2']+S2))+1,2):
            L = int(LL/2)
            # All the hat-factors
            fac1 = get_hat_factors(bra,ket,L,S2)
            
            # Exponential of L
            fac2 = (-1)**(L) 
            
            # 9j factors
            fac3 = get_9j_factors(bra,ket,L,S2)

            # 6j factors
            fac4 = get_6j_factors(bra,ket,L,S2)
            
            # Mass ratio
            d=1/3.0

            # ***********************************************
            # Arguments:
            # (n,l,nc,lc,n1,l1,n2,l2,lr,d) 
            # All are (int) except d which is a float.
            # For documentation, see the file my_braket.f90
            # This bracket is defined as in Trjilfaj (1972)
            # but gmosh takes '1/d' as argument, hence the
            # argument passed here is 1/3 which differs
            # compared to Eq. (31) in Phys. Rev. C 61 044001.
            # More details of the implementation is provided
            # in the report.
            # ***********************************************
            fac5 = gmosh.angmom_lib.gmosh(bra['n'],bra['l'],\
                    bra['cN'],bra['cL'],ket['cN'],ket['cL'],\
                    ket['n'],ket['l'],L,d)
            
            matrix_el += fac1*fac2*fac3*fac4*fac5

    
    return matrix_el

def check_A(A):
    '''
        Funtion that returns true if A is a projector.
    '''

    if (np.max(np.abs(A-A@A))>1e-6):
        return False
    else:
        return True

def check_eigs(eigs):
    '''
        Function that checks if all eigenvalues are either zero or one.
    '''

    num1 = np.sum([np.abs(e-1.0)<1e-6 for e in eigs])
    num2 = np.sum([np.abs(e-0.0)<1e-6 for e in eigs])
    if (num1+num2 != len(eigs)):
        print(eigs)
        return False
    else:
        return True

def setup_3N_antisymmetrizer_matrix(basis,verbose=False):
    '''
        Computes the antisymmetrizer matrix <bra|A|ket>.

    Args:

        basis (list) : list of basis states as computed with 
                       basis_states.NN_basis_nl()
    Returns:

        A (np.array(len(basis),len(basis)))
    '''
       
    A = np.zeros((len(basis),len(basis)))
    
    for i,bra in enumerate(basis):
        for j,ket in enumerate(basis):

            if verbose:
                print(f'\ni={i}, j={j}')   
            
            braN = 2*bra['n']+bra['l'] + 2*bra['cN']+bra['cL']
            ketN = 2*ket['n']+ket['l'] + 2*ket['cN']+ket['cL']

            matrix_el = 0
            # Check that (N,J,T) are the same for the states
            if (braN == ketN and bra['J2'] == ket['J2'] and bra['T2'] == ket['T2']): 
                matrix_el = setup_3N_antisymmetrizer_element(bra,ket,verbose)
                #matrix_el_fortran = tm.antisym_me_mymod_braket(bra,ket)
            else:
                matrix_el = 0

            # Set calculated matrix element
            A[i,j] = (-2.0/3.0)*matrix_el
            if (i==j):
                A[i,j] += 1.0/3.0
            
    if verbose:
        print(f'A.shape={A.shape}')
    return A

# *****************************************************************************
# OLD function
# *****************************************************************************
def setup_3N_antisymmetrizer_element_2(bra,ket,verbose=False):
    '''
        Computes the antisymmetrizer matrix elements as defined in eq. (XX) 
        in the notes.

    Args:

        bra (basis_state) : basis state as computed with 
                            basis_states.NN_basis_nl()
        
        ket (basis_state) : basis state as computed with 
                            basis_states.NN_basis_nl()
    Returns:
        matrix_el (float) 
    '''
    CASE = 0
    matrix_el = 0

    print(f'\nGetting A element')
    print(f'bra={bra}')
    print(f'ket={ket}')

    L_range  = range(abs(bra['l']-bra['cL']),bra['l']+bra['cL']+1,1)
    S2_range = range(abs(2*bra['s']-1),2*bra['s']+2,2)

    print(f'L-range: {np.array(L_range)}')
    print(f'S-range: {np.array(S2_range)/2}')
    # Sum over L and S
    for L in L_range:
        # Implement the ket  restriction for L
        if (bs.triag(ket['l'],ket['cL'],L)):
            
            # Note that S2 = 2*S where S is the total spin
            for S2 in S2_range:
                # Implement other restrictions for S
                if (bs.triag(ket['s'],1/2,S2/2) and bs.triag(bra['J2']/2,L,S2/2)): 
                    print(f'L = {L}, S2 = {S2}')
                    # All the hat-factors
                    fac1 = get_hat_factors(bra,ket,L,S2)
                    
                    # 9j factors
                    fac3 = get_9j_factors(bra,ket,L,S2)

                    # 6j factors
                    fac4 = get_6j_factors(bra,ket,L,S2)
                    
                    # As Navratil expression, works for Nmax = 0,1
                    if CASE == 0:
                        # Exponential of L
                        fac2 = (-1)**(L) 
                        
                        # Mass ratio
                        d=1/3.0

                        # *********************************************
                        # Arguments:
                        # (n,l,nc,lc,n1,l1,n2,l2,lr,d) 
                        # All are (int) except d which is a float.
                        # For documentation, see the file my_braket.f90
                        # *********************************************
                        fac5 = gmosh.angmom_lib.gmosh(bra['n'],bra['l'],\
                                bra['cN'],bra['cL'],ket['cN'],ket['cL'],\
                                ket['n'],ket['l'],L,d)
                    
                    elif CASE == 1:

                        #fac2 = (-1)**(L+bra['s']+ket['s']+bra['t']+\
                        #       ket['t']-bra['l']-bra['cL'])
                        
                        fac2 = (-1)**(L+ket['l']+bra['cL']) # equal!
                        
                        f = (-1)**L
                        if (fac2-f)>1e-8:
                            print('fac2 dissagree')
                            print(bra)
                            print(ket)
                            print(f'L={L}, S = {S2/2}')
                        d = 1/3.0
                        #fac5 = gmosh.angmom_lib.gmosh(bra['n'],bra['l'],\
                        #        bra['cN'],bra['cL'],ket['n'],ket['l'],\
                        #         ket['cN'],ket['cL'],L,d)
                        #fac5 = gmosh.angmom_lib.gmosh(bra['n'],bra['l'],\
                        #        bra['cN'],bra['cL'],ket['cN'],ket['cL'],\
                        #        ket['n'],ket['l'],L,d)
                        fac5 = gmosh.angmom_lib.gmosh(bra['cN'],bra['cL'],\
                                bra['n'],bra['l'],ket['n'],ket['l'],\
                                ket['cN'],ket['cL'],L,d)
                        
                        #print(f'{fac5}, {fac5_*(-1)**(bra["l"]-L)}')
                        
                    if verbose:
                        #print(f'fac2={fac2}')
                        print(f'L={L}, S={S2/2:<.1f}')
                        print(f'gmosh={fac5}')           
                    
                    # Check some relations that must hold
                    sumN = bra['n']+bra['cN']+ket['n']+ket['cN']
                    sumL = (bra['l']+bra['cL']-ket['l']-ket['cL'])/2
                    #fac6 = (-1)**(sumN)
                    fac6 = 1
                    #if (((-1)**sumL-(-1)**sumN)>1e-8):
                    #    print("ERROR sumL != sumN")
                    #    print(bra)
                    #    print(ket)
                    #if ((-1)**(bra['l']+bra['cL']) - (-1)**(ket['l']+ket['cL']))>1e-8:
                    #    print("ERROR")
                    
                    #sumN2 = bra['cN']+ket['n']+ket['cN']
                    #fac6 = (-1)**sumN2
                    print(f'fac1={fac1}')
                    print(f'fac2={fac2}')
                    print(f'fac3={fac3}')
                    print(f'fac4={fac4}')
                    print(f'fac5={fac5}')
                    print(f'fac6={fac6}')
                    matrix_el += fac1*fac2*fac3*fac4*fac5*fac6
    
    print(f'matrix_el={matrix_el}')

    return matrix_el
