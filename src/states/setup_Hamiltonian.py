'''
    setup_Hamiltonian.py
    --------------------
    
    This file contains some methods to construct Hamiltonian matrix in a
    HO basis.

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
    # Implement delta-functions
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
    
    if(     bra['l']   != ket['l']   or\
            bra['s']   != ket['s']   or\
            bra['j']   != ket['j']   or\
            bra['t']   != ket['t']):
        return 0
    elif (bra['n'] == ket['n']+1):
        matrix_el = +np.sqrt((n_p+1)*(n_p+l_p+3.0/2.0)) # NOTE +
    elif (bra['n'] == ket['n']):
        matrix_el = (2*n_p+l_p+3.0/2.0)
    elif (bra['n'] == ket['n']-1):
        matrix_el = +np.sqrt(n_p*(n_p+l_p+1.0/2.0))     # NOTE +
    else:
        matrix_el = 0
    if 'N' in bra:
        matrix_el = matrix_el*Omega/3
        #print(f'matrix_el={matrix_el}')
        return matrix_el
    else:
        matrix_el = matrix_el*Omega/2
        #print(f'matrix_el={matrix_el}')
        return matrix_el



def get_V_HO_alpha_basis(bra,ket,Omega,mN):
    '''
        Computes matrix elements of the Harmonic Oscillator
        potential in the alpha basis
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

def get_T_matrix_alpha_basis(alpha_basis_list,Omega,mN):
    
    T_matrix = np.zeros((len(alpha_basis_list),len(alpha_basis_list)))
    for i,bra in enumerate(alpha_basis_list):
        for j,ket in enumerate(alpha_basis_list):
            el = get_T_me_alpha_basis(bra,ket,Omega,mN)
            T_matrix[i,j] = el
    return T_matrix



def get_two_body_HO_potential_el(n,l,n_p,l_p,s,j,t,pot_dict):
    el = 0
    if (t==0):
        m_t = 0
        key = (n,l,n_p,l_p,s,j,m_t)
        el = pot_dict[key]
    else:
        key1 = (n,l,n_p,l_p,s,j,0)
        el += pot_dict[key1]*(1.0/3.0)
        key2 = (n,l,n_p,l_p,s,j,1)
        el += pot_dict[key2]*(2.0/3.0)
    return el

def get_V_me_alpha_basis(bra,ket,Omega,mN,pot_dict):
    
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
    #fac = np.real((1j)**(+bra['l']-ket['l']))
    #print(np.imag((1j)**(+bra['l']-ket['l'])))
    #fac = (-1)**(+bra['n']+ket['n'])
    fac = 1
    el = fac*get_two_body_HO_potential_el(bra['n'],bra['l'],ket['n'],ket['l']
            ,bra['s'],bra['j'],bra['t'],pot_dict)
    return el



def get_V_matrix_alpha_basis(alpha_basis_list,Omega,mN,pot_dict):
    
    # Setup empty matrix
    V_matrix = np.zeros((len(alpha_basis_list),len(alpha_basis_list)))
    for i,bra in enumerate(alpha_basis_list):
        for j,ket in enumerate(alpha_basis_list):
            el = get_V_me_alpha_basis(bra,ket,Omega,mN,pot_dict)
            V_matrix[i,j] = el

    return V_matrix


def alpha_to_Gamma_basis(Gamma_to_alpha_matrix,M_alpha_basis):
    return Gamma_to_alpha_matrix.T@M_alpha_basis@Gamma_to_alpha_matrix
    

def setup_H_Gamma_basis(Gamma_to_alpha_matrix,alpha_basis_list,Gamma_basis_list,\
        Omega,mN,pot_dict):
    
    # Compute kinetic energy operator
    T_alpha_basis = 3*get_T_matrix_alpha_basis(alpha_basis_list,Omega,mN)
    T_Gamma_basis = alpha_to_Gamma_basis(Gamma_to_alpha_matrix,T_alpha_basis)
    
    # Compute potential operator
    V_alpha_basis = get_V_matrix_alpha_basis(alpha_basis_list,Omega,mN,pot_dict)
    #V_alpha_basis = get_V_HO_matrix_alpha_basis(alpha_basis_list,Omega,mN)
    V_Gamma_basis = 3*alpha_to_Gamma_basis(Gamma_to_alpha_matrix,V_alpha_basis)
    #print(f'V_alpha={V_alpha_basis}')
    #print(f'T_alpha={T_alpha_basis}')
    #print(f'T_gamma={T_Gamma_basis}')
    #print(f'Gtoa={Gamma_to_alpha_matrix}')
    return T_Gamma_basis + V_Gamma_basis

def setup_H_alpha_basis(alpha_basis_list,Omega,mN,pot_dict):
    
    # Compute kinetic energy operator
    T_alpha_basis = get_T_matrix_alpha_basis(alpha_basis_list,Omega,mN)
    
    # Compute potential operator
    V_alpha_basis = get_V_matrix_alpha_basis(alpha_basis_list,Omega,mN,pot_dict)
    return T_alpha_basis + V_alpha_basis

if __name__=="__main__":
    # Set varibles
    Omega = 24   # MeV
    mN    = 938  # MeV
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
        basis_state_dir=f'Nmax={Nmax}_data'
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
                alpha_basis_list,Gamma_basis_list,Omega,mN,pot_dict)
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

