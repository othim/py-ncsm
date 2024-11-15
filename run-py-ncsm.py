'''
    run-py-ncsm.py
    --------------

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

# Import py-ncsm-file
os.sys.path.append("src/states/")
os.sys.path.append("src/hamiltonian/")
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
    args['isospin_sym']      = True
    args['fast_comp']        = False
    args['dim_lanczos']      = 2000
    args['interaction_file'] = "interactions/nmax36_sqb_Np100.txt"
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

# *****************************************************************************
# *****************************************************************************
# ********************************* MAIN CODE *********************************
# *****************************************************************************
# *****************************************************************************

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
    sys.stdout = open(args['output_file'], 'a')
    print_git_version()

print("-------------------------------------------------------------------------------")
print("                              *** PARAMETERS ***")
print("-------------------------------------------------------------------------------")

if DEFAULT:
    print(f'input_file=none, reading default arguments\n')
else:
    print(f'input_file = {sys.argv[1]}\n')

mN = 0
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
    time_d = {}

    basis_state_dir=f'src/states/Nmax={Nmax}_data'
    # Read in the basis states and change-of-basis matrix
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


    print(f'Gammta_to_alpha.shape={Gamma_to_alpha.shape}')
    print(f'len(Gamma_states)={len(Gamma_basis_list)}')
    print(f'len(alpha_states)={len(alpha_basis_list)}\n')
    len_Gamma_arr.append(len(Gamma_basis_list))
    len_alpha_arr.append(len(alpha_basis_list))

    # Load potential matrix elements
    print('Loading potential matrix elements... ',end='')
    pot_dict,pot = lp.load_potential_file(args['interaction_file'])
    print('Done!')
    

    # Setup the hamiltonian
    print('Setting up Hamiltonian... ',end='')
    start = time.time()
    H_matrix_Gamma_basis = sh.setup_H_Gamma_basis(Gamma_to_alpha,\
            alpha_basis_list,Gamma_basis_list,args['hbar_omega'],mN,pot_dict,\
            args['isospin_sym'],args['fast_comp'])
    print(f'shape={H_matrix_Gamma_basis.shape}. ',end='')
    end = time.time()
    time_d['set_H'] = (end-start)*1000.0
    print(f'Done! Time={(end-start)*1000:.0f} ms')

    # Diagonalize the Hamiltonian
    print('Diagonalizing Hamiltonian... ',end='')
    start = time.time()
    # If the dimension is large enough Lanczos diagonalization is used
    if H_matrix_Gamma_basis.shape[0]<args['dim_lanczos']:
        print('exact diag... ',end='')
        eigs, eigv = np.linalg.eigh(H_matrix_Gamma_basis)
        method_arr.append('exact')
    else:
        print('Lanczos diag... ',end='')
        eigs, eigv  = sl.eigsh(H_matrix_Gamma_basis,k=1,tol=0,which='SA')
        method_arr.append('Lanczos')
    end = time.time()
    time_d['diag'] = (end-start)*1000.0
    time_arr.append(time_d)
    print(f'Done! Time={(end-start)*1000:.3f} ms')
    
    print(f'E={np.min(eigs):.5f} MeV')
    print('----------------------------------\n\n')
    E_arr.append(np.min(eigs))


# *****************************************************************************
# ******************************** PRINT RESULT *******************************
# *****************************************************************************
print("")
print("-------------------------------------------------------------------------------")
print("                                *** OUTPUT ***")
print("-------------------------------------------------------------------------------")
print(f'Nmax  E          dim (Gamma)  dim (alpha)  Comp. H (ms)  Diag. H (ms)  method') 
for i,E in enumerate(E_arr):
    print(f'{args["nmax_arr"][i]:<5} {E:<10.5f} {len_Gamma_arr[i]:<12} {len_alpha_arr[i]:<12} {time_arr[i]["set_H"]:<13.0f} {time_arr[i]["diag"]:<13.0f} {method_arr[i]:<8}') 

if args['output_file']!='none':
    sys.stdout = original
print(f'Done printing to file, program ended.')
