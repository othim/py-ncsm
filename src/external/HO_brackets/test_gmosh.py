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
n = 2
l = 3
cN = 1
cL = 2
n1 = 2
l1 = 2
n2 = 1
l2 = 3
L = 5
d = 1.0

print(gmosh.angmom_lib.gmosh(n,l,cN,cL,n1,l1,n2,l2,L,d))


# TMB(EE,LL,ER,LR,E1,L1,E2,L2,LM,D)


print(gmosh.tmb_brackets.tmb(int(2*n+l),int(l),int(2*cN+cL),int(cL),int(2*n1+l1)\
        ,int(l1),int(2*n2+l2),int(l2),int(L),d))
