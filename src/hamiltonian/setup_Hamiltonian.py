'''
    setup_Hamiltonian.py
    --------------------
    
    This file contains some methods to construct Hamiltonian matrix.

    Oliver Thim (2024)
'''

import numpy as np
import json
import sys
import time
import load_potential as lp
import copy

def get_T_me_alpha_basis(bra,ket,Omega,mN):
    '''
        Computes the matrix element of the kinetic energy
        <bra|T_{12}|ket> = <n cN \alpha | T_{12} | n' cN' \alpha'>. 
        Eq. (XX) in XX.

        NOTE: right now there is a ugly hack implemented to get the kinetic
        energy right for the 2N and 3N systems that just checks if the bra
        and ket contains certain quantum numbers. This should be changed.

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
    n   = bra['n']
    n_p = ket['n']
    l_p = ket['l']
    # Implement these delta functions if we have a 3N state
    if ('N' in bra):
        if (    bra['pi']  != ket['pi']   or\
                bra['J2'] != ket['J2']   or\
                bra['T2'] != ket['T2']):
            print(f'Error, alpha_states should only contain a fixed (J,T,pi)')
            sys.exit(1)
            return 0
        if (    bra['cN']  != ket['cN']  or\
                bra['cL']  != ket['cL']  or\
                bra['cJ2'] != ket['cJ2']):
            return 0
    # Implement delta function over 2N subsystem.
    if(     bra['l']   != ket['l']   or\
            bra['s']   != ket['s']   or\
            bra['j']   != ket['j']   or\
            bra['t']   != ket['t']):
        return 0
    # Note that there is a plus sign on the off-diagonal components.
    # this sign depends if the (-1)^n factor in the HO basis functions 
    # sit in position- or momentum-space. In our case we assume that the 
    # (-1)^n factor is in momentum space, and thus we get +.
    elif (bra['n'] == ket['n']+1):
        matrix_el = +np.sqrt((n_p+1)*(n_p+l_p+3.0/2.0)) # NOTE +
    elif (bra['n'] == ket['n']):
        matrix_el = (2*n_p+l_p+3.0/2.0)
    elif (bra['n'] == ket['n']-1):
        matrix_el = +np.sqrt(n_p*(n_p+l_p+1.0/2.0))     # NOTE +
    else:
        matrix_el = 0

    # If 3N system
    if 'N' in bra:
        matrix_el = matrix_el*Omega/3
        return matrix_el
    # If 2N system
    else:
        matrix_el = matrix_el*Omega/2
        return matrix_el



def get_V_HO_alpha_basis(bra,ket,Omega,mN):
    '''
        Computes matrix elements of the Harmonic Oscillator. This function
        can be used as a test to see that the HO eigenvalues for the 2N system 
        is correct.
    '''
    matrix_el = 0
    n   = bra['n']
    n_p = ket['n']
    l_p = ket['l']
    # Implement delta-functions
    if ('N' in bra):
        if (    bra['cN']  != ket['cN']  or\
                bra['cL']  != ket['cL']  or\
                bra['cJ2'] != ket['cJ2']):
            return 0
    
    if(     bra['l']   != ket['l']   or\
            bra['s']   != ket['s']   or\
            bra['j']   != ket['j']   or\
            bra['t']   != ket['t']):
        return 0
    elif (n == n_p+1):
        matrix_el = -np.sqrt((n_p+1)*(n_p+l_p+3/2))
    elif (n == n_p):
        matrix_el = 2*n_p+l_p+3/2
    elif (n == n_p-1):
        matrix_el = -np.sqrt(n_p*(n_p+l_p+1/2))
    else:
        matrix_el=0
    return matrix_el*Omega/2
    

def get_V_HO_matrix_alpha_basis(alpha_basis_list,Omega,mN):
    
    T_matrix = np.zeros((len(alpha_basis_list),len(alpha_basis_list)))
    for i,bra in enumerate(alpha_basis_list):
        for j,ket in enumerate(alpha_basis_list):
            el = get_V_HO_alpha_basis(bra,ket,Omega,mN)
            T_matrix[i,j] = el
    print(T_matrix.shape)
    return T_matrix




def get_two_body_HO_potential_el(n,l,n_p,l_p,s,j,t,pot_dict,isospin_sym):
    el = 0
    if not isospin_sym:
        # NOTE: This implementation is specifically for ^3H!!
        if (t==0):
            m_t = 0
            # np element
            key = (n,l,n_p,l_p,s,j,m_t)
            el = pot_dict[key]
        else:
            # Isospin average of nn and np
            key1 = (n,l,n_p,l_p,s,j,0)
            el += pot_dict[key1]*(1.0/3.0)
            key2 = (n,l,n_p,l_p,s,j,1)
            el += pot_dict[key2]*(2.0/3.0)
    else:
        m_t = 0
        # np element
        key = (n,l,n_p,l_p,s,j,m_t)
        el = pot_dict[key]
    return el

def get_V_me_alpha_basis(bra,ket,Omega,mN,pot_dict,isospin_sym):
    '''
        Computes the potential matrix element in the alpha basis.
    '''
    
    # If in 3N case.
    if 'N' in bra:
        # Check delta function over 3N variables, otherwise stop program
        if (    bra['pi']  != ket['pi']   or\
                bra['J2'] != ket['J2']   or\
                bra['T2'] != ket['T2']):
            print(f'Error, alpha_states should only contain a fixed (J,T,pi)')
            sys.exit(1)
            return 0

        # Implement delta-functions over third particle
        if (    bra['cN']  != ket['cN']   or\
                bra['cL']  != ket['cL']   or\
                bra['cJ2']   != ket['cJ2']):
            return 0
    
    # Implement delta-functions over two-body system
    if (    bra['s']  != ket['s']   or\
            bra['j']  != ket['j']   or\
            bra['t']  != ket['t']):
        return 0
    
    # Else, compute the potential matrix element

    # NOTE that conventions might eneter here for several reasons.
    # 
    # 1. If you potential matrix elements in the interaction file do not follow
    #    the Machleidt convention you might need to add a factor i^(l-l') here
    #    to compensate. See the report for more discussion.
    # 2. This code assumes that the HO ME are computed with the (-1)^n factor
    #    in the momentum space wave functions. If you have the opposite 
    #    convention, you need to:  add a factor (-1)^(n+n') here.
    #
    # I think that these corrections are OK, but I cannot guarantee that the 
    # code must not be changed at other paces too.
    #fac = np.real((1j)**(+bra['l']-ket['l']))
    el = get_two_body_HO_potential_el(bra['n'],bra['l'],ket['n'],ket['l']
            ,bra['s'],bra['j'],bra['t'],pot_dict,isospin_sym)
    return el



def get_V_matrix_alpha_basis(alpha_basis_list,Omega,mN,pot_dict,isospin_sym,fast):
    ''' 
        Construct the potential matrix from the potential elements.
    '''
    
    # Setup empty matrix
    V_matrix = np.zeros((len(alpha_basis_list),len(alpha_basis_list)))
    for i,bra in enumerate(alpha_basis_list):
        for j,ket in enumerate(alpha_basis_list):
            el = get_V_me_alpha_basis(bra,ket,Omega,mN,pot_dict,isospin_sym)
            V_matrix[i,j] = el
    '''
    if fast:
        for i,bra in enumerate(alpha_basis_list):
            for j,ket in enumerate(alpha_basis_list):
                if (j<i) : continue
                el = get_V_me_alpha_basis(bra,ket,Omega,mN,pot_dict,isospin_sym)
                V_matrix[i,j] = el
                V_matrix[j,i] = np.conj(el)
    else:
        for i,bra in enumerate(alpha_basis_list):
            for j,ket in enumerate(alpha_basis_list):
                el = get_V_me_alpha_basis(bra,ket,Omega,mN,pot_dict,isospin_sym)
                V_matrix[i,j] = el
    '''
    return V_matrix

def get_T_matrix_alpha_basis(alpha_basis_list,Omega,mN,fast):
    '''
        Construct the kinetic energy matrix from the kinetic energy elements.
    '''
    
    T_matrix = np.zeros((len(alpha_basis_list),len(alpha_basis_list)))
    for i,bra in enumerate(alpha_basis_list):
        for j,ket in enumerate(alpha_basis_list):
            el = get_T_me_alpha_basis(bra,ket,Omega,mN)
            T_matrix[i,j] = el
    '''
    # If fast=True, use that V^\dagger = V and only compute half of the elements
    if fast:
        for i,bra in enumerate(alpha_basis_list):
            for j,ket in enumerate(alpha_basis_list):
                if (j<i) : continue
                el = get_T_me_alpha_basis(bra,ket,Omega,mN)
                T_matrix[i,j] = el
                T_matrix[j,i] = np.conj(el)
    else:
        for i,bra in enumerate(alpha_basis_list):
            for j,ket in enumerate(alpha_basis_list):
                el = get_T_me_alpha_basis(bra,ket,Omega,mN)
                T_matrix[i,j] = el
    '''
    return T_matrix

def setup_H_alpha_basis(alpha_basis_list,Omega,mN,pot_dict,isospin_sym,fast):
    '''
       Construct the Hamiltonian matrix H = T+V in the alpha basis.
    '''
    
    # Compute kinetic energy operator
    T_alpha_basis = get_T_matrix_alpha_basis(alpha_basis_list,Omega,mN,fast)
    
    # Compute potential operator
    V_alpha_basis = get_V_matrix_alpha_basis(alpha_basis_list,Omega,mN,pot_dict,\
            isospin_sym,fast)
    return T_alpha_basis + V_alpha_basis

def alpha_to_Gamma_basis(Gamma_to_alpha_matrix,M_alpha_basis):
    '''
        Converts a matrix in the alpha (partially antisymmertric) basis to the 
        Gamma basis (fully antisymmetric).
    '''
    return Gamma_to_alpha_matrix.T@M_alpha_basis@Gamma_to_alpha_matrix
    

def setup_H_Gamma_basis(Gamma_to_alpha_matrix,alpha_basis_list,Gamma_basis_list,\
        Omega,mN,pot_dict,isospin_sym,fast):
    
    # Compute kinetic energy operator.
    T_alpha_basis = 3*get_T_matrix_alpha_basis(alpha_basis_list,Omega,mN,fast)
    T_Gamma_basis = alpha_to_Gamma_basis(Gamma_to_alpha_matrix,T_alpha_basis)
    
    # Compute potential operator
    V_alpha_basis = get_V_matrix_alpha_basis(alpha_basis_list,Omega,mN,pot_dict,\
            isospin_sym,fast)
    V_Gamma_basis = 3*alpha_to_Gamma_basis(Gamma_to_alpha_matrix,V_alpha_basis)
    
    return T_Gamma_basis + V_Gamma_basis


# Some test code to test functionality 
if __name__=="__main__":
    # Set varibles
    Omega = 24   # MeV
    mN    = 938  # MeV
    isospin_sym = True
    #interaction_file="../../interactions/vn3lo500_nmax30_jrelmax10_hw24.dat"
    #interaction_file="../../interactions/nmax8_halfmn.txt"
    #interaction_file="../../interactions/nmax16_np100.txt"
    #interaction_file="../../interactions/nmax16_cdbonn.txt"
    #interaction_file="../../interactions/test.txt"
    #interaction_file="../../interactions/nmax8_mnhalf.txt"
    #interaction_file="../../interactions/nmax8_nmhalf_rev.txt"
    #interaction_file="../../interactions/nmax8_mnhalf_non.txt"
    #interaction_file="../../interactions/nmax30_sqb.txt"
    interaction_file="../../interactions/nmax36_sqb_Np100.txt"
    #interaction_file="../../interactions/nmax30_sqb_cdbonn.txt"
    #interaction_file="../../interactions/nmax30_inc_T.txt"
    
    Nmax_arr = [0,2,4,6,8,10,12,14,16,20,24,28,30]
    #Nmax_arr = [0,2,4,6,8,10,12,14,16]
    #Nmax_arr = [0]
    #Nmax_arr = [0,2]
    #Nmax_arr = [0,2,4,6,8]
    E_arr = []
    len_Gamma_arr = []
    for Nmax in Nmax_arr:
        basis_state_dir=f'../states/Nmax={Nmax}_data'
        # Read in the basis states and change-of-basis matrix
        print(f'Reading states from \'{basis_state_dir}\'')
        with open(basis_state_dir+'/Gamma_basis_Nmax='+str(Nmax)+'.txt', "r") as file:
            Gamma_basis_list = json.load(file)
        with open(basis_state_dir + '/alpha_basis_Nmax='+str(Nmax)+'.txt','r') as file:
            alpha_basis_list = json.load(file)
        with open(basis_state_dir+'/Gamma_to_alpha_Nmax='+str(Nmax)+'.txt','r') as file:
            Gamma_to_alpha = np.loadtxt(file)

        # Check Nmax=0 special case and make sure that it is a matrix
        if (Gamma_to_alpha.shape == (2,)):
            Gamma_to_alpha = Gamma_to_alpha.reshape(2,1)

        print(f'Done!')


        print(f'Gammta_to_alpha.shape={Gamma_to_alpha.shape}')
        print(f'len(Gamma_states)={len(Gamma_basis_list)}')
        print(f'len(alpha_states)={len(alpha_basis_list)}\n')
        len_Gamma_arr.append(len(Gamma_basis_list))
        # Load potential matrix elements
        print(f'Loading potential matrix elements...')
        pot_dict,pot = lp.load_potential_file(interaction_file)
        print(f'Done!')
        

        # Setup the hamiltonian
        print(f'Setting up Hamiltonian...')
        start = time.time()
        H_matrix_Gamma_basis = setup_H_Gamma_basis(Gamma_to_alpha,\
                alpha_basis_list,Gamma_basis_list,Omega,mN,pot_dict,isospin_sym)
        print(H_matrix_Gamma_basis.shape)
        end = time.time()
        print(f'Done! Time={(end-start)*1000:.0f} ms')

        # Diagonalize the Hamiltonian
        print(f'Diagonalize Hamiltonian...')
        start = time.time()
        eigs,eigv = np.linalg.eigh(H_matrix_Gamma_basis)
        end = time.time()
        print(f'Done! Time={(end-start)*1000:.3f} ms')
        #print(eigs)
        print(f'\nnp.min(eigs) = {np.min(eigs):.3f} MeV\n\n')
        E_arr.append(np.min(eigs))

    print(f'Nmax \t E \t  dim') 
    for i,E in enumerate(E_arr):
        print(f'{Nmax_arr[i]:<8} {E:<8.3f} {len_Gamma_arr[i]:<8} ') 

