'''
    load_potential.py
    -----------------

    This file contains function to load potential matrix elements from files.
'''
import numpy as np
import time

def get_data(file_name):
    with open(file_name, 'r') as f:                                                  
        for index, line in enumerate(f):                                        
            # search string                                                     
            if 'DATA:' in line:  
                # Load potential into array
                pot = np.loadtxt(file_name,skiprows=index+1)
                return pot


def load_potential_file(file_name):
    '''
        Loads a potential file containing matrix elements (ME).

        The file is assumed to have the form:
        
        n   l   np  lp  s   j   mt  ME   
        0   0   0   0   0   0   -1  -7.966724483064686                           
        0   0   1   0   0   0   -1  -2.652044568838252                           
        0   0   2   0   0   0   -1  2.89434840507413                             
        0   0   3   0   0   0   -1  7.911560485250315 

        where the first line conaining the desciption of the quantum numbers
        is not included in the file.

        Retiurns:
            Dictionary with the quantum numbers as a tuple as key for
            each potential ME.
    '''
    #pot = get_data(file_name)
    pot = np.loadtxt(file_name)
    #for row in pot:
    #    print(row)
    # Convert to dictionary with tuple of quantum numbers as
    # key to enable O(1) lookup time for matrix elements.
    data_dict = {tuple(row[:7].astype(int)): row[7] for row in pot}

    return data_dict, pot


# Test the loading and lookup time
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
