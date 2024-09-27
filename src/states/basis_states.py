'''
    basis_states.py
    ---------------
    
    This file contains some methods to construct and manipulate
    Harmonic Oscillator (HO) states.

    Oliver Thim (2024)
'''
import time


def psi_nlm(n,l,m,r,theta,phi):
    '''
        Computes the single-particle HO wave function as defined in (...)

        Args:


        Returns:
    '''
    return 0


def NN_basis_j(j_min, j_max, mt_min, mt_max):
    '''
        Computes allowed antisymmetrized 2N basis states
        and returns a list of dictionaries containig the basis states.
    Args:
        j_min  : Minimum total angular momentum, (ls)j.
        j_max  : Maximum total angular momentum, (ls)j.
        mt_min : Minimum isospin projection 
        mt_max : Maximum isospin projection
    Returns:
        list(basis_states)
        basis_states : dictionary(l,s,j,t,mt,pi) (pi=parity)
    '''
    basis = []
    for mt in range(mt_min,mt_max+1,1):
        for J in range(j_min,j_max+1,1):
            for S in range(0,2):
                for L in range(abs(J-S),J+S+1,1):
                    for T in range(abs(mt),2,1):
                        if ((L+S+T)%2 != 0):
                            basis_state = {}
                            basis_state['l']  = L
                            basis_state['s']  = S
                            basis_state['j']  = J
                            basis_state['t']  = T
                            basis_state['mt'] = mt
                            basis_state['pi'] = (-1)**L
                            print(basis_state)
                            basis.append(basis_state)

    print('len(basis) = ',len(basis))
    return basis

def NNN_basis_nl(N_max, mt_min, mt_max,verbose=False):
    '''
        Computes 3N basis states which are antisummetrized in particle 1 and 2, 
        i.e. (-1)^{l+s+t} = -1. 

        n,l,s,t,mt,pi : quantum numbers for particle (1,2) subsystem.
        cN,cL,cJ2     : quantum numbers for particle 3 relative to (1,2) subsystem
        J2, T2        : Total J and T for 3N system

        NOTE: the 2 in cJ2,J2 and T2 denote that these are two times the 
        respective quantum numbers as this is more convenient to save since 
        these quantum numbers are half integers.
    Args:
        N_max (int)             : Maximum HO energy, N = 2n+l + 2*cN+cL (cN = \mathcal{N})
        mt_min -1,0,1           : Minimum isospin projection 
        mt_max -1,0,1, > mt_min : Maximum isospin projection

    Returns:
        basis (list) : list of 'basis_states (dictionary)' that contains keys : 
        dictionary(n,l,s,j,t,mt,pi,cN,cL,cJ2,J2,T2) (pi=parity)
    '''
    basis = []
    for mt in range(mt_min,mt_max+1,1):
        for n in range(0,int(N_max/2)+1,1):
            for l in range(0,N_max-2*n+1,1):
                for s in range(0,2):
                    for j in range(abs(l-s),l+s+1,1):
                        for t in range(abs(mt),2,1):
                            if ((l+s+t)%2 != 0): # continue only for antisym 2N
                            
                                for cN in range(0,int((N_max-2*n-l)/2)+1,1):
                                    for cL in range(0,N_max-2*n-l-2*cN+1,1):

                                        # The two means that this is two times
                                        # the quantum number. This is more convenient
                                        # since these can be half integers.
                                        for cJ2 in range(abs(2*cL-1),2*cL+2,2):
                                            for J2 in range(abs(2*j-cJ2),2*j+cJ2+1,2):
                                                for T2 in range(abs(2*t-1),2*t+2,2):
                                                    basis_state = {}
                                                    basis_state['n']  = int(n)
                                                    basis_state['l']  = int(l)
                                                    basis_state['s']  = int(s)
                                                    basis_state['j']  = int(j)
                                                    basis_state['t']  = int(t)
                                                    basis_state['mt'] = int(mt)
                                                    basis_state['pi'] = int((-1)**l)
                                                    
                                                    basis_state['cN'] = int(cN)
                                                    basis_state['cL'] = int(cL)
                                                    basis_state['cJ2'] = int(cJ2)
                                                    
                                                    basis_state['J2'] = int(J2)
                                                    basis_state['T2'] = int(T2)

                                                    basis.append(basis_state)

    if verbose:
        for b in enumerate(basis):
            print(b)
        print('len(basis) = ',len(basis))
    return basis


# If you run only this file you can check some functions
if __name__ == "__main__":
    start = time.time()
    NNStates = NN_basis_j(0,10,0,0)
    end = time.time()
    print(f'NN_basis_j: {(end-start)*1e3:.4f} ms')
    
    start = time.time()
    Nmax = 12
    NNStates = NNN_basis_nl(Nmax,-1,1,True)
    end = time.time()
    print(f'NN_basis_nl, Nmax={Nmax}, time={(end-start)*1e3:.4f} ms')
