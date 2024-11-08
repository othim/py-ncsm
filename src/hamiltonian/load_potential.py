'''
'''


import numpy as np
import time



def load_potential_file(file_name):
    '''
        Loads a potential file containing matrix-elements.
    '''

    # Load potential into array
    pot = np.loadtxt(file_name)

    # Convert to dictionary with tuple of quantum numbers as
    # key to enable O(1) lookup time for matrix elements.
    data_dict = {tuple(row[:7].astype(int)): row[7] for row in pot}

    return data_dict, pot




if __name__=="__main__":

    # Load potential
    fname = '../../interactions/vn3lo500_nmax30_jrelmax10_hw24.dat'
    data_dict,pot = load_potential_file(fname)

    for i,key in enumerate(data_dict):
        print(f'key={str(key):<30}, data={data_dict[key]}')
    
    start = time.time()
    for i,key in enumerate(data_dict):
        a = data_dict[key]
    end = time.time()
    print(f'Total lookup time         = {(end-start)*1e3:.3f} ms')
    print(f'Total lookup per element  = {(end-start)*1e6/pot.shape[0]:.3f} mu s')
    print(f'pot.shape={pot.shape}')
