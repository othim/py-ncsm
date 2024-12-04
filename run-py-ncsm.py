'''
    run-py-ncsm.py
    --------------
    - Main script to run the code, see README.txt

    Oliver Thim,
    Chalmers (2024)
'''
# Import packages
import os
import sys
import numpy as np
import scipy.sparse.linalg as sl
import configparser
import json
import time
import git # Can be removed, just to print version for extra control

# Append the path to py-ncsm
dir_path   = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path + '/src/states/')
sys.path.append(dir_path + '/src/hamiltonian/')
import load_potential as lp
import setup_Hamiltonian as sh


# *****************************************************************************
# ************************** DEFINE DEFAULT ARGUMENTS *************************
# *****************************************************************************
def print_git_version():
    # Find the hash for the version of this repo
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
    print("*******************************************************************************")
    print(f'***** Running py-ncsm. Version: {sha} ******')
    print("*******************************************************************************\n")


def get_default_args():
    '''
        This function returns a set of default arguments
        in case no infile is provided.
    '''
    args = {}
    #args['nmax_arr']         = [0,2,4,6,8]
    args['nmax_arr']         = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]
    args['hbar_omega']       = 24
    args['isospin_sym']      = False
    args['fast_comp']        = False
    args['dim_lanczos']      = 2000
    args['interaction_file'] = "interactions/idaho_n3lo_nmax_40_hw_24_Np_80_finite.txt"
    args['output_file']      = "none"
    return args

def get_args_from_infile(file):
    config.read(file)
    args = dict(config.items('settings'))
    args['nmax_arr'] = json.loads(args['nmax_arr']) # Convert to dictionary
    args['hbar_omega']    = int(args['hbar_omega'])
    args['isospin_sym']   = (args['isospin_sym'] == 'True')
    args['fast_comp']     = (args['fast_comp'] == 'True')
    args['dim_lanczos']   = int(args['dim_lanczos']) 
    return args

def load_states_and_transform_matrix(basis_state_dir,Nmax):
    
    print(f'Reading states from \'{basis_state_dir}\'... ',end='')
    with open(basis_state_dir+'/Gamma_basis_Nmax='+str(Nmax)+'.txt', "r") as file:
        Gamma_basis_list = json.load(file)
    with open(basis_state_dir + '/alpha_basis_Nmax='+str(Nmax)+'.txt','r') as file:
        alpha_basis_list = json.load(file)
    with open(basis_state_dir+'/Gamma_to_alpha_Nmax='+str(Nmax)+'.txt','r') as file:
        Gamma_to_alpha = np.loadtxt(file)

    # Check Nmax=0 special case and make sure that it is a matrix
    if (Gamma_to_alpha.shape == (2,)):
        Gamma_to_alpha = Gamma_to_alpha.reshape(2,1)
    
    print('Done!')
    return alpha_basis_list,Gamma_basis_list,Gamma_to_alpha

def diag_3N(Nmax,hbar_omega,isospin_sym,fast_comp,dim_lanczos,interaction_file):
    time_d = {}
    
    basis_state_dir=dir_path+f'/src/states/Nmax={Nmax}_data'
    # Read in the basis states and change-of-basis matrix
    alpha_basis_list,Gamma_basis_list,Gamma_to_alpha = \
            load_states_and_transform_matrix(basis_state_dir,Nmax)


    print(f'Gammta_to_alpha.shape={Gamma_to_alpha.shape}')
    print(f'len(Gamma_states)={len(Gamma_basis_list)}')
    print(f'len(alpha_states)={len(alpha_basis_list)}\n')

    # Load potential matrix elements
    print('Loading potential matrix elements... ',end='')
    pot_dict,pot = lp.load_potential_file(interaction_file)
    print('Done!')
    

    # Setup the hamiltonian
    print('Setting up Hamiltonian... ',end='')
    start = time.time()
    H_matrix_Gamma_basis = sh.setup_H_Gamma_basis(Gamma_to_alpha,\
            alpha_basis_list,Gamma_basis_list,hbar_omega,pot_dict,\
            isospin_sym,fast_comp)
    print(f'shape={H_matrix_Gamma_basis.shape}. ',end='')
    end = time.time()
    time_d['set_H'] = (end-start)*1000.0
    print(f'Done! Time={(end-start)*1000:.0f} ms')

    # Diagonalize the Hamiltonian
    print('Diagonalizing Hamiltonian... ',end='')
    start = time.time()
    # If the dimension is large enough Lanczos diagonalization is used
    method = ""
    if H_matrix_Gamma_basis.shape[0]<dim_lanczos:
        print('exact diag... ',end='')
        eigs, eigv = np.linalg.eigh(H_matrix_Gamma_basis)
        method = "exact"
    else:
        print('Lanczos diag... ',end='')
        eigs, eigv  = sl.eigsh(H_matrix_Gamma_basis,k=1,tol=0,which='SA')
        method = "Lanczos"
    end = time.time()
    time_d['diag'] = (end-start)*1000.0
    print(f'Done! Time={(end-start)*1000:.3f} ms')
    return eigs,eigv,time_d,len(Gamma_basis_list),len(alpha_basis_list),method


# *****************************************************************************
# *****************************************************************************
# ********************************* MAIN CODE *********************************
# *****************************************************************************
# *****************************************************************************
if __name__=='__main__':

    print_git_version()

    # *****************************************************************************
    # ******************************* READ ARGUMENTS ******************************
    # *****************************************************************************
    config = configparser.ConfigParser(inline_comment_prefixes="#")

    DEFAULT = False
    if len(sys.argv)<2:
        args = get_default_args()
        DEFAULT = True
    else:
        args = get_args_from_infile(sys.argv[1])

    # If an output file is specified, write to that file
    if args['output_file']!='none':
        print(f'Writing output to file: {args["output_file"]}')
        original = sys.stdout
        sys.stdout = open(args['output_file'], 'w')
        print_git_version()

    print("-------------------------------------------------------------------------------")
    print("                              *** PARAMETERS ***")
    print("-------------------------------------------------------------------------------")

    if DEFAULT:
        print(f'input_file=none, reading default arguments\n')
    else:
        print(f'input_file = {sys.argv[1]}\n')

    print('---------------------')
    print(f'        args')
    print('---------------------')
    for key,val in args.items():
        print(f'{key:<16} = {val}')
    print("\n\n")


    print("-------------------------------------------------------------------------------")
    print("                             *** INFO FROM RUN ***")
    print("-------------------------------------------------------------------------------")
    # *****************************************************************************
    # ******************************* LOOP OVER NMAX ******************************
    # *****************************************************************************
    E_arr         = []
    len_Gamma_arr = []
    len_alpha_arr = []
    method_arr    = []
    time_arr      = []
    for Nmax in args['nmax_arr']:

        eigs, eigv, time_d,len_Gamma,len_alpha, method = diag_3N(Nmax,args['hbar_omega'],args['isospin_sym']\
                ,args['fast_comp'],args['dim_lanczos'],args['interaction_file'])
        
        time_arr.append(time_d)
        print(f'E={np.min(eigs):.5f} MeV')
        print('----------------------------------\n\n')
        E_arr.append(np.min(eigs))
        len_Gamma_arr.append(len_Gamma)
        len_alpha_arr.append(len_alpha)
        method_arr.append(method)


    # *****************************************************************************
    # ******************************** PRINT RESULT *******************************
    # *****************************************************************************
    print("")
    print("-------------------------------------------------------------------------------")
    print("                                *** OUTPUT ***")
    print("-------------------------------------------------------------------------------")
    print(f'Nmax  E                    dim (Gamma)  dim (alpha)  Comp. H (ms)  Diag. H (ms)  method') 
    for i,E in enumerate(E_arr):
        print(f'{args["nmax_arr"][i]:<5} {E:<20.16f} {len_Gamma_arr[i]:<12} {len_alpha_arr[i]:<12} {time_arr[i]["set_H"]:<13.0f} {time_arr[i]["diag"]:<13.0f} {method_arr[i]:<8}') 

    if args['output_file']!='none':
        sys.stdout = original
    print(f'Done printing to file, program ended.')
