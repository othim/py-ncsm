'''
    run-py-ncsm.py
    --------------

    TODO: 
    - Fix output file generation
    - Fix ISOSPIN_SYM flag. -> fix it in the setup H file also.
    - Add Lancos diag
'''
# Import packages
import os
import sys
import numpy as np
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
def get_default_args():
    '''
        This function returns a set of default arguments
        in case no infile is provided.
    '''
    print(f'input_file=none, reading default arguments\n')
    args = {}
    args['nmax_arr']         = [0,2,4,6,8]
    args['hbar_omega']       = 24
    args['isospin_sym']      = True
    args['interaction_file'] = "interactions/nmax36_sqb_Np100.txt"
    args['output_file']      = "none"
    return args

def get_args_from_infile(file):
    print(f'input_file = {file}\n')
    config.read(file)
    args = dict(config.items('settings'))
    args['nmax_arr'] = json.loads(args['nmax_arr']) # Convert to dictionary
    args['hbar_omega']    = int(args['hbar_omega'])
    args['isospin_sym']   = (args['isospin_sym'] == 'True')
    return args

# *****************************************************************************
# *****************************************************************************
# ********************************* MAIN CODE *********************************
# *****************************************************************************
# *****************************************************************************

# Find the hash for the version of this repo
repo = git.Repo(search_parent_directories=True)
sha = repo.head.object.hexsha
print("*******************************************************************************")
print(f'***** Running py-ncsm. Version: {sha} ******')
print("*******************************************************************************\n")
print("-------------------------------------------------------------------------------")
print("                              *** PARAMETERS ***")
print("-------------------------------------------------------------------------------")

# *****************************************************************************
# ******************************* READ ARGUMENTS ******************************
# *****************************************************************************
config = configparser.ConfigParser(inline_comment_prefixes="#")

if len(sys.argv)<2:
    args = get_default_args()
else:
    args = get_args_from_infile(sys.argv[1])

# If an output file is specified, write to that file
if args['output_file']!='none':
    print(f'Writing output to file: {args["output_file"]}')
    original = sys.stdout
    sys.stdout = open(args['output_file'], 'a')

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
E_arr = []
len_Gamma_arr = []
len_alpha_arr = []
for Nmax in args['nmax_arr']:
    
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
            alpha_basis_list,Gamma_basis_list,args['hbar_omega'],mN,pot_dict,args['isospin_sym'])
    print(f'shape={H_matrix_Gamma_basis.shape}. ',end='')
    end = time.time()
    print(f'Done! Time={(end-start)*1000:.0f} ms')

    # Diagonalize the Hamiltonian
    print('Diagonalizing Hamiltonian... ',end='')
    start = time.time()
    eigs,eigv = np.linalg.eigh(H_matrix_Gamma_basis)
    end = time.time()
    print(f'Done! Time={(end-start)*1000:.3f} ms')
    
    print(f'E={np.min(eigs):.3f} MeV')
    print('----------------------------------\n\n')
    E_arr.append(np.min(eigs))


# *****************************************************************************
# ******************************** PRINT RESULT *******************************
# *****************************************************************************
print("")
print("-------------------------------------------------------------------------------")
print("                                *** OUTPUT ***")
print("-------------------------------------------------------------------------------")
print(f'Nmax \t E \t  dim (Gamma) \t dim (alpha)') 
for i,E in enumerate(E_arr):
    print(f'{args["nmax_arr"][i]:<8} {E:<8.3f} {len_Gamma_arr[i]:<8} \t {len_alpha_arr[i]:<8}') 

if args['output_file']!='none':
    sys.stdout = original
print(f'Done printing to file, program ended.')
