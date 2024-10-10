'''
    setup_Hamiltonian.py
    --------------------
    
    This file contains some methods to construct Hamiltonian matrix in a
    HO basis.

    Oliver Thim (2024)
'''

import numpy as np





def get_T_ME_HO_alpha_basis(bra,ket,Omega,mN):
    '''
        Computes the matrix element of the kinetic energy
        <bra|T_{12}|ket> = <n cN \alpha | T_{12} | n' cN' \alpha'>. 
        Eq. (XX) in XX.

        Args:
            bra (basis_state) : basis state as computed with 
                                basis_states.NN_basis_nl()
            ket (basis_state) : basis state as computed with 
                                basis_states.NN_basis_nl()
            Omega (float)     : HO frquency (MeV).
            mN (float)        : nucleon mass (MeV).
        Returns:
            <bra|T_{12}|ket> (float) : Kinetic energy ME in the HO basis (MeV).
    '''
    matrix_el = 0
    n  = bra['n']
    np = ket['n']
    lp = ket['l']
    # Implement delta-functions
    if (    bra['cN']  != ket['cN']  or\
            bra['cL']  != ket['cL']  or\
            bra['cJ2'] != ket['cJ2'] or\
            bra['l']   != ket['l']   or\
            bra['s']   != ket['s']   or\
            bra['j']   != ket['j']   or\
            bra['t']   != ket['t']):
        return 0
    elif (bra['n'] == ket['n']+1):
        matrix_el = ((-1)**(n+np))*np.sqrt((np+1)*(np+lp+3/2))
    elif (bra['n'] == ket['n']):
        matrix_el = ((-1)**(n+np))*(2*np+lp+3/2)
    elif (bra['n'] == ket['n']-1):
        matrix_el = ((-1)**(n+np))*np.sqrt(np*(np+lp+1/2))
    else:
        matrix_el = 0
    
    return matrix_el*(Omega/2.0)*(1/(2*mN))

def get_T_M_HO_Gamma_basis(Gbasis,Omega,mN):
    '''
        Computes 3N kinetic energy operator in Gamma basis compused of matrix
        elements <\Gamma|T_{c.m.}|\Gamma'>
    '''
    
    # Compute ME in the alpha basis

    # Change basis to the Gamma basis

    return 0

def get_V_ME_HO_alpha_basis(bra,ket,Omega,mN,V_pp_p_func):



def setup_H_ME_HO_Gamma_basis():
