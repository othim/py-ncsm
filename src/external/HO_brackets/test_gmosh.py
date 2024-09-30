'''
    Test the module that is built with f2py.
'''

import numpy as np
import gmosh



# *********************************************
# Arguments:
# (n,l,nc,lc,n1,l1,n2,l2,lr,d) 
# All are (int) except d which is a float.
# For documentation, see the file my_braket.f90
# *********************************************
print(gmosh.angmom_lib.gmosh(1,2,1,3,1,2,1,3,3,2.0))
