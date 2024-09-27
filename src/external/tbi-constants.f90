! file with a single module
! named constants. This module contains the definition 
! of _all_ mathematical/ physical/ numerical/ flow ...
! parameters that are used throughout the program

MODULE constants
  
  IMPLICIT NONE
  
  ! PROGRAM VERSION ('do not change')
  CHARACTER(LEN=120), PARAMETER :: version = '0.99 (120820)'
  
  ! directories for storing precalculated quantities
  CHARACTER(LEN=80), PARAMETER, PUBLIC :: PRECALCDIR='./precalc'
  
  ! kinds
  !
  INTEGER, PARAMETER, PUBLIC :: dp  = KIND(1.0D0)
  
  ! misc
  !
  INTEGER, PARAMETER, PUBLIC  :: irrel       = -99
  REAL(DP), PARAMETER, PUBLIC :: zero        = 1.d-6 
  REAL(DP), PARAMETER, PUBLIC :: ho_norm_tol = 1.d-6
  REAL(DP), PARAMETER, PUBLIC :: global_eps  = 1.d-6

  ! Mathematical constants
  !
  REAL(DP), PUBLIC, PARAMETER :: pi = 3.141592653589793_dp
    
  ! Physical constants
  !
  REAL(DP), PUBLIC, PARAMETER :: speed_of_light     = 2.99792458_dp    !xE+8  ! m/s
  REAL(DP), PUBLIC, PARAMETER :: hbar               = 1.054571726_dp   !xE-34 ! J*s
  REAL(DP), PUBLIC, PARAMETER :: electron_charge    = 1.602176565_dp   !xE-19 ! Coulomb = 1 eV
  REAL(DP), PUBLIC, PARAMETER :: electric_constant  = 8.854187817620_dp!xE−12! F·m−1
  REAL(DP), PUBLIC, PARAMETER :: electron_mass      = 0.510998928D0    ! MeV
  ! Derived constants
  !
  REAL(DP), PUBLIC, PARAMETER :: hbarc = hbar * speed_of_light/electron_charge * 100.0_dp
  REAL(DP), PUBLIC, PARAMETER :: alpha_fsc = electron_charge**2/(4.0_dp*pi*electric_constant*hbar*speed_of_light)
  ! Numerical constants
  !
  REAL(DP), PUBLIC, PARAMETER :: pi_squared     = pi * pi 
  REAL(DP), PUBLIC, PARAMETER :: sqrt2          = sqrt(2.0_dp)
  REAL(DP), PUBLIC, PARAMETER :: sqrt6          = sqrt(6.0_dp)
  REAL(DP), PUBLIC, PARAMETER :: one_over_sqrt2 = 1.0_dp/sqrt(2.0_dp)
  REAL(DP), PUBLIC, PARAMETER :: one_third      = 1.0_dp/3.0_dp
  REAL(DP), PUBLIC, PARAMETER :: two_thirds     = 2.0_dp/3.0_dp
  REAL(DP), PUBLIC, PARAMETER :: sqrt2_thirds   = sqrt2/3.0_dp
  REAL(DP), PUBLIC, PARAMETER :: one_over_sqrt2_pow3  = (1.0_dp/sqrt(2.0_dp))**3.0_dp
  REAL(DP), PUBLIC, PARAMETER :: hbarc_over_sqrt2_pow3 = (hbarc/sqrt(2.0_dp))**3.0_dp
  REAL(DP), PUBLIC, PARAMETER :: rad2deg        = 180.0_dp/pi
  REAL(DP), PUBLIC, PARAMETER :: deg2rad        = pi/180.0_dp
  COMPLEX*16, PUBLIC, PARAMETER :: cmplxi = (0.0D0,1.0D0) 
  COMPLEX*16, PUBLIC, PARAMETER :: two_i = (0.0D0,2.0D0) 
  COMPLEX*16, PUBLIC, PARAMETER :: one_over_two_i = (0.0D0,-0.5D0) 

  ! Physical parameter, set in SUBROUTINE initialize
  ! masses and coupling constants differ sometimes
  ! depending on the year and definition of the 
  ! interaction.
  REAL(DP), PUBLIC :: gA 
  REAL(DP), PUBLIC :: fpi
  !
  REAL(DP), PUBLIC :: proton_mass               ! MeV/c2 (depends on twobody model parameters)
  REAL(DP), PUBLIC :: neutron_mass              ! MeV/c2 (depends on twobody model parameters)
  REAL(DP), PUBLIC :: nucleon_mass_z(-3:3)      ! 2*Tz = -3, -2, -1, 0, 1, 2, 3
  REAL(DP), PUBLIC :: oscillator_length_z(-3:3) ! hbarc/SQRT((nucleon_mass)*(hbar_omega))
  REAL(DP), PUBLIC :: pion_pm_mass, pion_0_mass, pion_mass
    
  ! input file parameters
  !
  ! explanatory title of calculation
  CHARACTER(LEN=240) :: title
  
  ! flag for ouput to unit=SCR (see below)
  LOGICAL, PUBLIC :: verbose
  ! determines wether the program should run
  ! in twobody or threebody mode. 
  CHARACTER(LEN=80), PUBLIC :: nbody_basis
  ! flag for just carrying out a memcheck
  LOGICAL, PUBLIC :: memcheck

  ! normal order truncation control parameter
  CHARACTER(LEN=80), PUBLIC :: normal_order
  LOGICAL, PUBLIC :: NORMALORDER

  ! desired reference frame for the result
  CHARACTER(LEN=80), PUBLIC :: reference_frame
  
  ! angular momentum coupled result
  CHARACTER(LEN=80), PUBLIC :: jcoupling
  ! which threebody scheme to compute in
  CHARACTER(LEN=80), PUBLIC :: scheme_mode
  ! angular momentum projection limits doubled value
  INTEGER, PUBLIC :: Mproj_min, Mproj_max

  ! flag to diagonalize the four/three/two-body hamiltonian
  ! for all N<=Nmax
  CHARACTER(LEN=80), PUBLIC :: compute_Nmax_series
  CHARACTER(LEN=80), PUBLIC :: compute_OBDM
  CHARACTER(LEN=80), PUBLIC :: compute_isomix
  REAL(DP), PUBLIC :: lawson_beta

  ! oscillator parameter in MeV
  REAL(DP), PUBLIC :: hbar_omega
  
  ! dimension of Hilbert space of basis states
  INTEGER, PUBLIC :: nmax_lab, lmax_lab, NPmax, lmax_3lab
  INTEGER, PUBLIC :: twobody_P_space
  LOGICAL, PUBLIC :: SQUARE, TRIANGULAR, USER_NMAX, USER_LMAX, USER_L3MAX

  ! coulomb inclusion flag
  CHARACTER(LEN=80), PUBLIC :: include_coulomb
  
  CHARACTER(LEN=80), PUBLIC :: CIB
 
  REAL(DP)         , PUBLIC :: T_fac
  REAL(DP)         , PUBLIC :: vNN_fac
  REAL(DP)         , PUBLIC :: v3N_c1_fac
  REAL(DP)         , PUBLIC :: v3N_c3_fac
  REAL(DP)         , PUBLIC :: v3N_c4_fac
  REAL(DP)         , PUBLIC :: v3N_cD_fac
  REAL(DP)         , PUBLIC :: v3N_cE_fac
  REAL(DP)         , PUBLIC :: v3N_2pi_fac
  REAL(DP)         , PUBLIC :: v3N_2pi_1pi_fac
  REAL(DP)         , PUBLIC :: v3N_2pi_cont_fac
  REAL(DP)         , PUBLIC :: v3N_rings_fac
  REAL(DP)         , PUBLIC :: v3N_relCS_fac
  REAL(DP)         , PUBLIC :: v3N_relCT_fac
  REAL(DP)         , PUBLIC :: v3N_rel2pi_fac
   
  REAL(DP)         , PUBLIC :: v3N_c1
  REAL(DP)         , PUBLIC :: v3N_c3
  REAL(DP)         , PUBLIC :: v3N_c4
  REAL(DP)         , PUBLIC :: v3N_cD
  REAL(DP)         , PUBLIC :: v3N_cE
  REAL(DP)         , PUBLIC :: v3N_2pi
  REAL(DP)         , PUBLIC :: v3N_2pi_1pi
  REAL(DP)         , PUBLIC :: v3N_2pi_cont
  REAL(DP)         , PUBLIC :: v3N_rings
  REAL(DP)         , PUBLIC :: v3N_relCS
  REAL(DP)         , PUBLIC :: v3N_relCT
  REAL(DP)         , PUBLIC :: v3N_rel2pi
  
  ! twobody interaction model
  CHARACTER(LEN=900), PUBLIC :: ARG_WOUT_FILENAME = 'none'
  CHARACTER(LEN=900), PUBLIC :: v2
  CHARACTER(LEN=900), PUBLIC :: v2_file
  CHARACTER(LEN=900), PUBLIC :: v3_file
  CHARACTER(LEN=900), PUBLIC :: v3_file_c1
  CHARACTER(LEN=900), PUBLIC :: v3_file_c3
  CHARACTER(LEN=900), PUBLIC :: v3_file_c4
  CHARACTER(LEN=900), PUBLIC :: v3_file_cD
  CHARACTER(LEN=900), PUBLIC :: v3_file_cE
  CHARACTER(LEN=900), PUBLIC :: v3_file_2pi
  CHARACTER(LEN=900), PUBLIC :: v3_file_2pi_1pi
  CHARACTER(LEN=900), PUBLIC :: v3_file_2pi_cont
  CHARACTER(LEN=900), PUBLIC :: v3_file_rings
  CHARACTER(LEN=900), PUBLIC :: v3_file_relCS
  CHARACTER(LEN=900), PUBLIC :: v3_file_relCT
  CHARACTER(LEN=900), PUBLIC :: v3_file_rel2pi

  CHARACTER(LEN=900), PUBLIC :: v3_c1_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_c1_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_c1_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_c1_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_c1_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_c1_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_c1_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_c1_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_c1_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_c1_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_c3_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_c3_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_c3_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_c3_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_c3_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_c3_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_c3_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_c3_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_c3_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_c3_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_c4_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_c4_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_c4_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_c4_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_c4_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_c4_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_c4_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_c4_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_c4_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_c4_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_cD_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_cD_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_cD_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_cD_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_cD_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_cD_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_cD_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_cD_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_cD_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_cD_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_cE_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_cE_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_cE_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_cE_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_cE_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_cE_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_cE_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_cE_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_cE_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_cE_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_2pi_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_2pi_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_2pi_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_2pi_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_2pi_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_2pi_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_2pi_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_2pi_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_2pi_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_2pi_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_2pi1pi_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_2pi1pi_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_2pi1pi_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_2pi1pi_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_2pi1pi_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_2pi1pi_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_2pi1pi_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_2pi1pi_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_2pi1pi_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_2pi1pi_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_rings_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_rings_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_rings_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_rings_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_rings_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_rings_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_rings_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_rings_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_rings_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_rings_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_2picont_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_2picont_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_2picont_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_2picont_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_2picont_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_2picont_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_2picont_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_2picont_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_2picont_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_2picont_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_relCS_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_relCS_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_relCS_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_relCS_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_relCS_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_relCS_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_relCS_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_relCS_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_relCS_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_relCS_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_relCT_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_relCT_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_relCT_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_relCT_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_relCT_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_relCT_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_relCT_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_relCT_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_relCT_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_relCT_chn10
  CHARACTER(LEN=900), PUBLIC :: v3_rel2pi_chn1
  CHARACTER(LEN=900), PUBLIC :: v3_rel2pi_chn2
  CHARACTER(LEN=900), PUBLIC :: v3_rel2pi_chn3
  CHARACTER(LEN=900), PUBLIC :: v3_rel2pi_chn4
  CHARACTER(LEN=900), PUBLIC :: v3_rel2pi_chn5
  CHARACTER(LEN=900), PUBLIC :: v3_rel2pi_chn6
  CHARACTER(LEN=900), PUBLIC :: v3_rel2pi_chn7
  CHARACTER(LEN=900), PUBLIC :: v3_rel2pi_chn8
  CHARACTER(LEN=900), PUBLIC :: v3_rel2pi_chn9
  CHARACTER(LEN=900), PUBLIC :: v3_rel2pi_chn10

  INTEGER          , PUBLIC :: v2_chiral_order
  REAL(DP)         , PUBLIC :: v2_cutoff
  LOGICAL          , PUBLIC :: TWOBODY_MOMENTUM_SPACE
  
  ! relative j min and max for twobody interaction
  INTEGER, PUBLIC :: v2_jmin, v2_jmax

  ! large space used for solving the twobody problem
  INTEGER, PUBLIC :: N2max

  ! threebody interaction model
  CHARACTER(LEN=80), PUBLIC :: v3  
  
  ! relative j min and max for threebody interaction
  INTEGER, PUBLIC :: v3_jmin, v3_jmax 

  ! total isospin and isospin projection min and max for threebody configs
  INTEGER, PUBLIC :: v3_tmin, v3_tmax 
  INTEGER, PUBLIC :: v3_tzmin, v3_tzmax 
  
  ! total parity
  CHARACTER(LEN=2), PUBLIC :: v3_total_parity

  ! large space used for solving the threebody problem
  INTEGER, PUBLIC :: N3max
  INTEGER, PUBLIC :: hoeft_N3max

  ! FOURBODY MODEL SPACE PARAMETERS
  INTEGER, PUBLIC :: N4max
  INTEGER, PUBLIC :: J4
  INTEGER, PUBLIC :: T4
  ! total parity
  CHARACTER(LEN=2), PUBLIC :: fourbody_parity
  
  ! mass of current nucleus
  INTEGER, PUBLIC :: mass_nucleus
  
  ! fermi vaccum definition for normal ordering
  INTEGER, PUBLIC :: number_of_protons
  INTEGER, PUBLIC :: number_of_neutrons

  ! model space geometry
  CHARACTER(LEN=80), PUBLIC :: model_space
  
  ! isospin projection limits
  ! doubled
  
  ! SCATTERING PARAMETERS SWITCH
  CHARACTER(LEN=80), PUBLIC :: phase_shifts
  CHARACTER(LEN=80), PUBLIC :: scattering_observables
  CHARACTER(LEN=80), PUBLIC :: scatter_chi
  
  ! file containing the computed phase shifts
  CHARACTER(LEN=80), PUBLIC :: pp_phasefile
  CHARACTER(LEN=80), PUBLIC :: pn_phasefile
  CHARACTER(LEN=80), PUBLIC :: nn_phasefile
  
  ! file containing the computed scattering observables
  CHARACTER(LEN=80), PUBLIC :: pp_obsfile
  CHARACTER(LEN=80), PUBLIC :: pn_obsfile

  ! file containing the experimental scattering observables
  CHARACTER(LEN=80), PUBLIC :: scattfile
  ! file containing the computed chi squares
  CHARACTER(LEN=80), PUBLIC :: chifile

  ! file with single particle basis
  CHARACTER(LEN=80), PUBLIC :: spfile

  ! file with 2/3(/4) body confs in rel/lab
  CHARACTER(LEN=80), PUBLIC :: confs_file
  
  ! file with nbody matrix elements
  CHARACTER(LEN=80), PUBLIC :: output_format
  LOGICAL          , PUBLIC :: CCFORM, SMFORM, SMFORMTE, STFORM, SMFORMAE
  CHARACTER(LEN=80), PUBLIC :: incl_T0
  CHARACTER(LEN=80), PUBLIC :: incl_U0
  CHARACTER(LEN=80), PUBLIC :: incl_VR
  CHARACTER(LEN=80), PUBLIC :: subtract_HCM
  CHARACTER(LEN=80), PUBLIC :: nbmefile
  
  ! precalculation files
  CHARACTER(LEN=80), PUBLIC :: hocfps3_IN
  CHARACTER(LEN=80), PUBLIC :: hocfps3_OUT
  CHARACTER(LEN=80), PUBLIC :: hocfps4_IN
  CHARACTER(LEN=80), PUBLIC :: hocfps4_OUT
  CHARACTER(LEN=900), PUBLIC :: hocfps3_filename
  CHARACTER(LEN=900), PUBLIC :: hocfps4_filename
  CHARACTER(LEN=300), PUBLIC :: v2_IN
  CHARACTER(LEN=300), PUBLIC :: v2_OUT
  CHARACTER(LEN=300), PUBLIC :: v3_IN
  CHARACTER(LEN=300), PUBLIC :: v3_OUT
  CHARACTER(LEN=300), PUBLIC :: v4_IN
  CHARACTER(LEN=300), PUBLIC :: v4_OUT

  CHARACTER(LEN=300), PUBLIC :: v3_IN_rk
  CHARACTER(LEN=300), PUBLIC :: v3_OUT_rk
  
  CHARACTER(LEN=180), PUBLIC :: V3ci
  CHARACTER(LEN=180), PUBLIC :: V3cD
  CHARACTER(LEN=180), PUBLIC :: V3cE
  CHARACTER(LEN=180), PUBLIC :: V4ci
  CHARACTER(LEN=180), PUBLIC :: V4cD
  CHARACTER(LEN=180), PUBLIC :: V4cE
  REAL(DP), PUBLIC           :: cD_IN
  REAL(DP), PUBLIC           :: cE_IN

  LOGICAL, PUBLIC :: READ_THREEBODY_HOCFPS
  LOGICAL, PUBLIC :: WRITE_THREEBODY_HOCFPS
  LOGICAL, PUBLIC :: READ_FOURBODY_HOCFPS
  LOGICAL, PUBLIC :: WRITE_FOURBODY_HOCFPS
  LOGICAL, PUBLIC :: READ_TWOBODY_INTERACTION
  LOGICAL, PUBLIC :: WRITE_TWOBODY_INTERACTION
  LOGICAL, PUBLIC :: READ_THREEBODY_INTERACTION
  LOGICAL, PUBLIC :: WRITE_THREEBODY_INTERACTION
  LOGICAL, PUBLIC :: READ_FOURBODY_INTERACTION
  LOGICAL, PUBLIC :: WRITE_FOURBODY_INTERACTION
  LOGICAL, PUBLIC :: WRITE_INTERACTION
  LOGICAL, PUBLIC :: READ_INTERACTION
  
  LOGICAL, PUBLIC :: READ_THREEBODY_INTERACTION_rk
  LOGICAL, PUBLIC :: WRITE_THREEBODY_INTERACTION_rk

  LOGICAL, PUBLIC :: READ_THREEBODY_ci_INTERACTION
  LOGICAL, PUBLIC :: READ_THREEBODY_cD_INTERACTION
  LOGICAL, PUBLIC :: READ_THREEBODY_cE_INTERACTION
  LOGICAL, PUBLIC :: READ_FOURBODY_ci_INTERACTION 
  LOGICAL, PUBLIC :: READ_FOURBODY_cD_INTERACTION 
  LOGICAL, PUBLIC :: READ_FOURBODY_cE_INTERACTION 

  ! runfile
  CHARACTER(LEN=80), PUBLIC :: runfile
  
  ! chinn5 IO file
  CHARACTER(LEN=20), PARAMETER :: chinn5file = 'chinn5IO.dat'

  ! ho mesh sizes, using the same value for 
  ! all Tz within twobody and threebody
  INTEGER, PUBLIC :: p3_ho_mesh_size ! threebody k-mesh
  INTEGER, PUBLIC :: r3_ho_mesh_size ! threebody r-mesh
  INTEGER, PUBLIC :: p2_ho_mesh_size ! twobody   k-mesh
  INTEGER, PUBLIC :: r2_ho_mesh_size ! twobody   r-mesh
  ! number of meshpoints used in threenucleon 
  ! momentum integrals
  INTEGER, PUBLIC :: q_mesh_size ! bessel integrals
  INTEGER, PUBLIC :: u_mesh_size ! legendre integrals
  ! switch for setting meshes and cutoffs manually
  LOGICAL, PUBLIC :: automesh
  ! radial cutoffs
  REAL(DP), PUBLIC :: r3_cutoff, p3_cutoff, p2_cutoff, r2_cutoff
  ! ...with isospin dependence (set in code). Mainly for future compatibility.
  REAL(DP), PUBLIC :: r3_cutoff_z(-3:3), p3_cutoff_z(-3:3), p2_cutoff_z(-2:2), r2_cutoff_z(-2:2)
  
  ! fortran units
  !
  ! screen
  INTEGER, PARAMETER :: SCR        = 6
  ! error unit (= 6 for screen)
  INTEGER, PARAMETER :: ERR        = 6
  INTEGER, PARAMETER :: FRUN       = 7
  INTEGER, PARAMETER :: FSP        = 10
  INTEGER, PARAMETER :: FNO1B      = 30
  INTEGER, PARAMETER :: FPP        = 31
  INTEGER, PARAMETER :: FPN        = 32
  INTEGER, PARAMETER :: FNN        = 33
  INTEGER, PARAMETER :: FPPP       = 34
  INTEGER, PARAMETER :: FPPN       = 35
  INTEGER, PARAMETER :: FPNN       = 36
  INTEGER, PARAMETER :: FNNN       = 37
  INTEGER, PARAMETER :: FIO        = 38
  INTEGER, PARAMETER :: FCHINN5OUT = 39
  INTEGER, PARAMETER :: FPHASEPP   = 51
  INTEGER, PARAMETER :: FPHASEPN   = 52
  INTEGER, PARAMETER :: FPHASENN   = 53
  INTEGER, PARAMETER :: FSCATTPP   = 54
  INTEGER, PARAMETER :: FSCATTPN   = 55
  INTEGER, PARAMETER :: FCHIDATA   = 56
  INTEGER, PARAMETER :: FCHIFILE   = 57
  INTEGER, PARAMETER :: FCONFS     = 60
  INTEGER, PARAMETER :: FCONFS_CC  = 61
  INTEGER, PARAMETER :: F3HOCFPS   = 62
  INTEGER, PARAMETER :: F4HOCFPS   = 63
  INTEGER, PARAMETER :: F2V        = 64
  INTEGER, PARAMETER :: F3V        = 65
  INTEGER, PARAMETER :: F3Vrk      = 66
  INTEGER, PARAMETER :: F4V        = 67
  INTEGER, PARAMETER :: FKAI       = 68
  INTEGER, PARAMETER :: FKAIv2     = 69
  INTEGER, PARAMETER :: FKAI_c1    = 70
  INTEGER, PARAMETER :: FKAI_c3    = 71
  INTEGER, PARAMETER :: FKAI_c4    = 72
  INTEGER, PARAMETER :: FKAI_cD    = 73
  INTEGER, PARAMETER :: FKAI_cE    = 74
  INTEGER, PARAMETER :: FKAI_2pi   = 75
  INTEGER, PARAMETER :: FKAI_2pi_1pi = 76
  INTEGER, PARAMETER :: FKAI_2pi_cont = 77
  INTEGER, PARAMETER :: FKAI_rings = 78
  INTEGER, PARAMETER :: FKAI_relCS = 79
  INTEGER, PARAMETER :: FKAI_relCT = 80
  INTEGER, PARAMETER :: FKAI_rel2pi = 81
  
  ! type sizes
  !
  INTEGER, PARAMETER, PUBLIC :: sizeof_int     = 4
  INTEGER, PARAMETER, PUBLIC :: sizeof_double  = 8
  INTEGER, PARAMETER, PUBLIC :: sizeof_float   = 4
  INTEGER, PARAMETER, PUBLIC :: sizeof_logical = 4
  
END MODULE constants
