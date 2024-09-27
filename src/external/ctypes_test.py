import ctypes

#REAL(DP) FUNCTION gmosh &                                                     
#318        (n,l,nc,lc,n1,l1,n2,l2,lr,d)                                             
#319     IMPLICIT NONE                                                               
#320     INTEGER, INTENT(IN) :: n,l,nc,lc,n1,l1,n2,l2,lr                             
#321     REAL(DP), INTENT(IN) :: d     

#print(lib.__angmom_lib_MOD_gmosh)
#lib.__angmom_lib_MOD_gmosh.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_int,\
#         ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,\
#         ctypes.c_int,ctypes.c_double]
#lib.__angmom_lib_MOD_gmosh.restype = ctypes.c_double
#print(lib.__angmom_lib_MOD_gmosh(int(1),int(0),int(1),int(0),int(0),int(0),\
#int(0),int(0),int(0),3.0))


#lib = ctypes.CDLL('./lib.so')
#print(lib.__angmom_lib_MOD_hat)
#lib.__angmom_lib_MOD_hat.argtype = [ctypes.c_long]
#lib.__angmom_lib_MOD_hat.restype = ctypes.c_double
#print(lib.__angmom_lib_MOD_hat(0))



lib = ctypes.CDLL('./test.so')
print(lib)
lib.__angmom_lib_MOD_gmosh.argtypes = [ctypes.c_double,ctypes.c_double]
lib.__angmom_lib_MOD_gmosh.restype = ctypes.c_double
print(lib.__angmom_lib_MOD_gmosh(2.0,3.0))
