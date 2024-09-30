!
!
! auxiliary functions for number represenations
! and logical checks etc
!
MODULE aux
  USE constants

  IMPLICIT NONE
  !
  ! defines equidistant points in given limits
  ! This construction is very handy since various
  ! loops over quantum numbers l,j,s,... occur 
  ! all the time
  TYPE, PUBLIC :: limits_type
     SEQUENCE
     INTEGER :: min = 0, max = 0, step = 0
  END TYPE limits_type
  !
  ! pair of values
  TYPE, PUBLIC :: pair
     SEQUENCE
     INTEGER :: v1 = 0, v2 = 0
  END TYPE pair

  INTERFACE check_max_set
     MODULE PROCEDURE check_max_set_int, check_max_set_double
  END INTERFACE

  TYPE(limits_type), PUBLIC, PARAMETER :: zero_limits = limits_type(0,0,0)
  INTEGER, PARAMETER, PUBLIC :: sizeof_limits_type = 3*sizeof_int
  !
CONTAINS
  !
  ! transform integer number into character representation
  FUNCTION i2a (i)
    INTEGER, INTENT (IN) :: i
    CHARACTER(LEN=20)    :: i2a
    WRITE (i2a, FMT = '(I)') i
    i2a = ADJUSTL(i2a)
  END FUNCTION i2a
  
  ! transform real number into character representation
  FUNCTION r2a (r)
    REAL(DP), INTENT (IN) :: r
    CHARACTER(LEN=30)     :: r2a
    WRITE (r2a, FMT = '(F5.2)') r
    r2a = ADJUSTL(r2a)
  END FUNCTION r2a

  FUNCTION l2a(l)
    LOGICAL, INTENT (IN) :: l
    CHARACTER(LEN=1)     :: l2a
    IF (l)       WRITE (l2a, FMT = '(A1)') 'Y'
    IF (.NOT. l) WRITE (l2a, FMT = '(A1)') 'N'
  END FUNCTION l2a
    
  FUNCTION write_limits (a) RESULT (res)
    CHARACTER(LEN=100)       :: res
    TYPE(limits_type), INTENT(IN) :: a
    res = '['//TRIM(ADJUSTL(i2a(a%min)))//', '//TRIM(ADJUSTL(i2a(a%max)))//'] step '//TRIM(ADJUSTL(i2a(a%step)))
  END FUNCTION write_limits
  ! swaps integers
  SUBROUTINE int_swap(a,b)
    INTEGER, INTENT(INOUT) :: a,b
    INTEGER        :: temp
    temp=a
    a=b
    b=temp
  END SUBROUTINE int_swap
  
  ! phase factor of type (-1)**arg
  FUNCTION minus_power (arg) RESULT (res)
    INTEGER              :: res
    INTEGER, INTENT (IN) :: arg
    INTEGER :: exponent
    
    exponent = ABS(arg) ! (-1)**N = (-1)**(-N)
    
    SELECT CASE (MOD(exponent,2))
    CASE (0)
       res =  1
    CASE (1)
       res = -1
    END SELECT

  END FUNCTION minus_power

  ! IF    a>b set b=a
  ! ELSE   do nothing
  SUBROUTINE check_max_set_int(a,b)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)    :: a
    INTEGER, INTENT(INOUT) :: b

    IF (a>b) b=a

  END SUBROUTINE check_max_set_int

  SUBROUTINE check_max_set_double(a,b)
    
    USE constants

    IMPLICIT NONE
    
    REAL(DP), INTENT(IN)    :: a
    REAL(DP), INTENT(INOUT) :: b

    IF (a>b) b=a

  END SUBROUTINE check_max_set_double

  subroutine progress(string,ndone,ntotal)
    implicit none
    character*(*) string
    character*255 prog,oldprog
    !double precision oldtime,hires_time,tl, now
    real oldtime, tl, now
    integer ndone,ntotal,i
    save oldprog,oldtime
    
    !if (ndone.eq.0) oldtime=hires_time()
    if (ndone.eq.0) call cpu_time(oldtime)
    call cpu_time(now)
    tl=now-oldtime
    if (tl.lt.0) tl=0
    if (ndone.gt.0) tl=(1.0*ntotal/ndone)*tl-tl
    
    ! When finished, print the total time taken rather
    ! than just 00:00!
    if (ndone.eq.ntotal) tl=now-oldtime
    
    write(prog,'(a25,1x,''['')') string
    do i=1,40
       prog(27+i:27+i)=' '
    enddo
    if (ndone<ntotal) write(prog(43:51),'(f7.1,''%'')') 100.0*ndone/ntotal
    do i=1,40
       if ((1.0*ndone/ntotal).gt.(1.0*i/40)) then
          if (prog(27+i:27+i).eq.' ') prog(27+i:27+i)='='
       endif
    enddo
           
    prog(67:67)=']'
    
    write(prog(68:72),'(i4.4,'':'')')int(tl/3600)
    write(prog(73:75),'(i2.2,'':'')')int((tl-int(tl/3600)*3600)/60)
    write(prog(76:77),'(i2.2)')int((tl-int(tl/60)*60))
    if (prog.ne.oldprog) write(0,'(a,a,$)') prog(1:77),char(13)
    oldprog=prog
    if (ndone.eq.ntotal) write(0,*)
    return
  end subroutine progress
    
END MODULE aux


MODULE factorials
  USE constants
  IMPLICIT NONE

  INTEGER, PARAMETER             :: maxvalue = 140
  INTEGER, PARAMETER             :: dbl_maxvalue = 300
  REAL(DP), DIMENSION(0:maxvalue), PUBLIC :: factorial
   REAL(DP), DIMENSION(0:dbl_maxvalue), PUBLIC :: double_factorial
    
CONTAINS
  
  SUBROUTINE precalculate_factorials
    INTEGER :: i
    
    factorial(0) = 1.0_dp
    DO i = 1, maxvalue
       factorial(i) = i * factorial(i-1)
    ENDDO
    
    double_factorial(0) = 1.0_dp
    double_factorial(1) = 1.0_dp
    
    DO i = 2, dbl_maxvalue
       double_factorial(i) = i * double_factorial(i-2)
    END DO

END SUBROUTINE precalculate_factorials
  
  DOUBLE PRECISION FUNCTION  dfac(m)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: m
    INTEGER :: i
    
    IF (MOD(m,2).NE.1) STOP 'wrong argument to dfac'
    dfac = 0.0D0
    IF (m == 1)RETURN
    DO i=3,m,2
       dfac=dfac+LOG(FLOAT(i))
    ENDDO
  END FUNCTION  dfac

  DOUBLE PRECISION FUNCTION  fac(m)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: m
    INTEGER :: i
    
    fac = 0.0D0
    IF(m == 0) RETURN
    DO i=1,m
       fac=fac+LOG(FLOAT(i))
    ENDDO
    
  END FUNCTION fac

END MODULE factorials
!           
!    Angular momentum library
!
MODULE angmom_lib
  USE constants
  USE factorials
  IMPLICIT NONE
  INTEGER, PRIVATE, PARAMETER :: maxjj = 140
  REAL(DP), PRIVATE :: f_mb(maxjj),g_mb(maxjj),w_mb(maxjj)
  REAL(DP), PUBLIC :: binom(0:maxjj,0:maxjj)

CONTAINS

  SUBROUTINE precalculate_binomials
  
    INTEGER :: i,j
    
    DO i = 0, maxjj
       DO j = 0, i
          binom(i,j) = fbinom(i,j)
       END DO
    END DO
    
  END SUBROUTINE precalculate_binomials

  ! Checks triangular relation
  ! |a-b| <= ab <= a+b TRUE
  ! ELSE: FALSE
  LOGICAL FUNCTION triag(a, b, ab)
    INTEGER :: a, b, ab
    triag = .FALSE.
    IF ( ab < ABS (a - b) ) RETURN
    IF ( ab > a + b ) RETURN
    triag = .TRUE.
  END FUNCTION triag

  REAL(DP) FUNCTION pladder(jj,mm) 
    
    INTEGER, INTENT(IN) :: jj,mm
    REAL(DP) :: j,m
    
    pladder = 0.0_dp
    
    j=REAL(jj,kind=dp)*0.5_dp
    m=REAL(mm,kind=dp)*0.5_dp
    
    ! stepping out of space
    IF(m>=j) RETURN

    pladder = SQRT(j*(j+1.0_dp) - m*(m+1.0_dp))
    
  END FUNCTION pladder

  REAL(DP) FUNCTION mladder(jj,mm) 
    
    INTEGER, INTENT(IN) :: jj,mm
    REAL(DP) :: j,m
    
    mladder = 0.0_dp
    
    j=REAL(jj,kind=dp)*0.5_dp
    m=REAL(mm,kind=dp)*0.5_dp
    
    ! stepping out of space
    IF(m+j<=0) RETURN

    mladder = SQRT(j*(j+1.0_dp) - m*(m-1.0_dp))
    
  END FUNCTION mladder
    
  REAL(DP) FUNCTION hat(x)
    
    INTEGER, INTENT(IN) :: x
    REAL(DP) :: xx

    xx = REAL(x,kind=dp)+1.0_dp
    
    hat = SQRT(xx)
        
  END FUNCTION hat

   SUBROUTINE setup_moshinsky
    IMPLICIT NONE
    INTEGER :: i
    REAL(DP) :: a
  
    ! Moshinsky brackets
    f_mb(1)=0.
    g_mb(1)=LOG(0.5D0)
    w_mb(1)=0.
    DO i=2,maxjj
       a=i-1
       f_mb(i)=f_mb(i-1)+LOG(a)
       g_mb(i)=g_mb(i-1)+LOG(a+0.5D0)
       w_mb(i)=LOG(a+a+1.)
    ENDDO

  END SUBROUTINE setup_moshinsky

  ! Library of angular momentum coupling coefficient routines in fortran 90
  ! Paul Stevenson, Oxford University/Oak Ridge National Laboratory.
  ! spaul@mail.phy.ornl.gov
  !  integer, parameter :: rk = selected_real_kind(p=15)
  !contains
  
  ! calculates a clebsch-gordan coefficient < j1/2 m1/2 j2/2 m2/2 | j/2 m/2 >
  ! arguments are integer and twice the true value.
  FUNCTION cleb (j1, m1, j2, m2, j, m)
    INTEGER, INTENT(IN) :: j1, m1, j2, m2, j, m
    INTEGER             :: par, z, zmin, zmax
    REAL(DP)            :: cleb, factor, sum
        
    ! some checks for validity (let's just return zero for bogus arguments)  
    IF (2*(j1/2)-INT(2*(j1/2.0)) /= 2*(ABS(m1)/2)-int(2*(ABS(m1)/2.0)) .OR. &
         2*(j2/2)-INT(2*(j2/2.0)) /= 2*(ABS(m2)/2)-int(2*(ABS(m2)/2.0)) .OR. &
         2*(j/2)-INT(2*(j/2.0)) /= 2*(ABS(m)/2)-int(2*(ABS(m)/2.0)) .OR. &
         j1<0 .OR. j2<0 .or. j<0 .OR. ABS(m1)>j1 .OR. ABS(m2)>j2 .OR.&
         ABS(m)>j .or. j1+j2<j .OR. ABS(j1-j2)>j .OR. m1+m2/=m) THEN
       cleb= 0.0_dp
       
    ELSE
       factor = 0.0_dp
       factor = binom(j1,(j1+j2-j)/2) / binom((j1+j2+j+2)/2,(j1+j2-j)/2)
       factor = factor * binom(j2,(j1+j2-j)/2) / binom(j1,(j1-m1)/2)
       factor = factor / binom(j2,(j2-m2)/2) / binom(j,(j-m)/2)
       factor = SQRT(factor)
       
       zmin = MAX(0,j2+(j1-m1)/2-(j1+j2+j)/2,j1+(j2+m2)/2-(j1+j2+j)/2)
       zmax = MIN((j1+j2-j)/2,(j1-m1)/2,(j2+m2)/2)
       
       sum=0.0_dp
       DO z = zmin,zmax
          par=1
          IF (2*(z/2)-INT(2*(z/2.0)) /= 0) par=-1
          sum=sum+par*binom((j1+j2-j)/2,z)*binom((j1-j2+j)/2,(j1-m1)/2-z)*&
               binom((-j1+j2+j)/2,(j2+m2)/2-z)
       ENDDO
       
       cleb = factor*sum
    ENDIF
  END FUNCTION cleb
  
  RECURSIVE FUNCTION fbinom(n,r) RESULT(res)
    INTEGER, INTENT(IN)  :: n,r
    REAL(DP)             :: res
    IF (n==r .OR. r==0) THEN
       res = 1.0_dp
    ELSEIF (r==1) THEN
       res = REAL(n,dp)
    ELSE
       res = REAL(n,dp) / REAL(n-r,dp) * fbinom(n-1,r)
    ENDIF
  END FUNCTION fbinom
  
  FUNCTION threej(ja, ma, jb, mb, jc, mc) 
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: ja,ma,jb,mb,jc,mc
    REAL(DP) :: threej, cg
    
    cg = cleb(ja,ma,jb,mb,jc,-mc)
    threej = cg*((-1)**((ja-jb-mc)/2))/SQRT(jc+1.0_dp)

  END FUNCTION threej

  function sixj(a,b,c,d,e,f)
    implicit none
    integer, intent(in) :: a,b,c,d,e,f
    real(dp) :: sixj
    integer :: nlo, nhi, n
    real(dp) :: outfactors, sum, sumterm
    
    ! calculates a Wigner 6-j symbol. Argument a-f are integer and are
    ! twice the true value of the 6-j's arguments, in the form
    ! { a b c }
    ! { d e f }
    ! Calculated using binomial coefficients to allow for (reasonably) high
    ! arguments.

    ! First check for consistency of arguments:
    sixj=0.0_dp
    if(mod(a+b,2)/=mod(c,2)) return
    if(mod(c+d,2)/=mod(e,2)) return
    if(mod(a+e,2)/=mod(f,2)) return
    if(mod(b+d,2)/=mod(f,2)) return
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(c-d)>e .or. c+d<e) return
    if(abs(a-e)>f .or. a+e<f) return
    if(abs(b-d)>f .or. b+d<f) return

    outfactors = angdelta(a,e,f)/angdelta(a,b,c)
    outfactors = outfactors * angdelta(b,d,f)*angdelta(c,d,e)

    nlo = max( (a+b+c)/2, (c+d+e)/2, (b+d+f)/2, (a+e+f)/2 )
    nhi = min( (a+b+d+e)/2, (b+c+e+f)/2, (a+c+d+f)/2)

    sum=0.0_dp
    do n=nlo,nhi
       sumterm = (-1)**n
       sumterm = sumterm * binom(n+1,n-(a+b+c)/2)
       sumterm = sumterm * binom((a+b-c)/2,n-(c+d+e)/2)
       sumterm = sumterm * binom((a-b+c)/2,n-(b+d+f)/2)
       sumterm = sumterm * binom((b-a+c)/2,n-(a+e+f)/2)
       sum=sum+sumterm
    end do

    sixj = sum * outfactors

  end function sixj
  
  function angdelta(a,b,c)
    implicit none
    integer :: a,b,c
    real(dp)    :: angdelta, scr1
    ! calculate the function delta as defined in varshalovich et al. for
    ! use in 6-j symbol:
    scr1= factorial((a+b-c)/2)
    scr1=scr1/factorial((a+b+c)/2+1)
    scr1=scr1*factorial((a-b+c)/2)
    scr1=scr1*factorial((-a+b+c)/2)
    angdelta=sqrt(scr1)
  end function angdelta
  
  function ninej(a,b,c,d,e,f,g,h,i)
    implicit none
    integer  :: a,b,c,d,e,f,g,h,i
    real(dp) :: ninej, sum
    integer  :: xlo, xhi
    integer  :: x
    ! calculate a 9-j symbol. The arguments are given as integers twice the
    ! value of the true arguments in the form
    ! { a b c }
    ! { d e f }
    ! { g h i }

    ninej=0.0_dp
    ! first check for bogus arguments (and return zero if so)
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(d-e)>f .or. d+e<f) return
    if(abs(g-h)>i .or. g+h<i) return
    if(abs(a-d)>g .or. a+d<g) return
    if(abs(b-e)>h .or. b+e<h) return
    if(abs(c-f)>i .or. c+f<i) return
    
    xlo = max(abs(b-f),abs(a-i),abs(h-d))
    xhi = min(b+f,a+i,h+d)
    
    sum=0.0_dp
    do x=xlo,xhi,2
       sum=sum+(-1)**x*(x+1)*sixj(a,b,c,f,i,x)*sixj(d,e,f,b,x,h)*&
            sixj(g,h,i,x,a,d)
    end do
    ninej=sum

  end function ninej
  
  !     This routine calculates the moshinsky vector bracket      
  !     Note that D=mass1/mass2                                   
  !     Ref  m.sotona and m.gmitro  comp.phys.comm 3(1972)53      
  !
  REAL(DP) FUNCTION gmosh &
       (n,l,nc,lc,n1,l1,n2,l2,lr,d) 
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n,l,nc,lc,n1,l1,n2,l2,lr
    REAL(DP), INTENT(IN) :: d
    INTEGER :: ip,ixf,ix, iyi, iyf, j1f,j2,k1i,k1f,m1f,iy,m2f,k2, &
         m2,m2i,m1,j1,k2f,k2i,k1
    REAL(DP) :: dl, d1l, bb, ba, anorm, y, p, bc, cfac, bm , &
         sm, s, sxy, bxy

    gmosh=0.
    
    IF(n+n+nc+nc+l+lc-n1-n1-n2-n2-l1-l2 /= 0 ) RETURN
    IF(l+lc-lr < 0 ) RETURN
    IF(l1+l2-lr < 0 ) RETURN
    IF(ABS(l-lc)-lr > 0 ) RETURN
    IF(ABS(l1-l2)-lr > 0 ) RETURN

    DL=LOG(D)
    D1L=LOG(D+1.)
    bb=f_mb(n1+1)+f_mb(n2+1)+f_mb(n+1)-f_mb(nc+1)+ &
         g_mb(n1+l1+1)+g_mb(n2+l2+1) &
         -g_mb(n+l+1)-g_mb(nc+lc+1)
    ba=w_mb(l1+1)+w_mb(l2+1)+w_mb(lc+1)+w_mb(l+1)+ &
         f_mb(l1+l2-lr+1)+f_mb(l+lc+lr+2) &
         +f_mb(l+lc-lr+1)+f_mb(lc+lr-l+1)+ &
         f_mb(lr+l-lc+1)-f_mb(l1+l2+lr+2) &
         -f_mb(l1+lr-l2+1)-f_mb(l2+lr-l1+1)-DBLE(l)*d1l
    ip=lr+n+n1+n2
    p=1+2*(ip/2*2-ip)
    anorm=p*EXP(0.5D0*(bb+ba))
    y=0.
    j1f=l+1
    DO j1=1,j1f
       j2=l+2-j1
       k1i=ABS(l1-j1+1)+1
       k1f=l1+j1
       DO k1=k1i,k1f,2
          m1f=n1-(j1+k1-l1)/2+2
          IF(m1f-1 < 0 )  CYCLE
          k2i=MAX(ABS(l2-j2+1),ABS(lc-k1+1))+1
          k2f=MIN(l2+j2,lc+k1)
          IF(k2i-k2f > 0 ) CYCLE
          DO k2=k2i,k2f,2
             m2f=n2-(j2+k2-l2)/2+2
             IF(m2f-1 < 0 )  CYCLE
             ip=j2-1+(l1+k1+j1+l2+k2+j2)/2
             p=1+2*(ip/2*2-ip)
             bc=0.5D0*(DBLE(k1+j2-2)*dl-DBLE(k1+k2-2)*d1l) &
                  +f_mb(k1+l1-j1+1)+f_mb(k1+k2-lc-1)+ &
                  f_mb(k2+l2-j2+1)-f_mb(k1+l1+j1)-f_mb(k1+k2+lc)- &
                  f_mb(k2+l2+j2)+w_mb(k1)+w_mb(k2)+f_mb((k1+l1+j1)/2)+ &
                  f_mb((k1+k2+lc)/2)+f_mb((k2+l2+j2)/2)- &
                  f_mb((k1+l1-j1)/2+1)-f_mb((l1+j1-k1)/2+1)- &
                  f_mb((j1+k1-l1)/2)-f_mb((k1+k2-lc)/2)- &
                  f_mb((k2+lc-k1)/2+1)-f_mb((lc+k1-k2)/2+1) &
                  -f_mb((k2+l2-j2)/2+1)-f_mb((l2+j2-k2)/2+1)- &
                  f_mb((j2+k2-l2)/2)
             cfac=p*EXP(bc)
             sxy=0.
             ixf=MIN(k1+k1,k1+k2-lc)-1
             DO ix=1,ixf
                iyi=MAX(1,ix+j1+l2-k1-lr)
                iyf=MIN(l2+l2+1,l1+l2-lr+1,l2+lc+ix-k1-j2+2)
                IF(iyi-iyf > 0 ) CYCLE
                DO iy=iyi,iyf
                   ip=ix+iy
                   p=1+2*(ip/2*2-ip)
                   bxy=f_mb(k1+k1-ix)+f_mb(l2+l2-iy+2)+ &
                        f_mb(k2+lc-k1+ix)+f_mb(l1+lr-l2+iy) &
                        -f_mb(ix)-f_mb(iy)-f_mb(k1+k2-lc-ix)- &
                        f_mb(l1+l2-lr-iy+2)-f_mb(k1-l2+lr-j1+iy-ix+1)- &
                        f_mb(l2-k1+lc-j2+ix-iy+3)
                   sxy=sxy+p*EXP(bxy)
                ENDDO
             ENDDO
             s=cfac*sxy
             sm=0.
             DO m1=1,m1f
                m2i=MAX(1,nc-m1-(k1+k2-lc)/2+3)
                IF(m2i-m2f > 0 ) CYCLE
                DO m2=m2i,m2f
                   ip=m1+m2
                   p=1+2*(ip/2*2-ip)
                   bm=DBLE(m1-1)*DL-DBLE(m1+m2-2)*d1l+g_mb(1) &
                        +g_mb(m1+m2+(k1+k2+lc)/2-2)-g_mb(k1+m1-1)- &
                        g_mb(k2+m2-1)+f_mb(m1+m2+(k1+k2-lc)/2-2)- &
                        f_mb(m1)-f_mb(m2)-f_mb(n1-m1-(j1+k1-l1)/2+3)- &
                        f_mb(n2-m2-(j2+k2-l2)/2+3) &
                        -f_mb(m1+m2-nc+(k1+k2-lc)/2-2)
                   sm=sm+p*EXP(bm)
                ENDDO
             ENDDO
             y=y+s*sm
          ENDDO
       ENDDO
    ENDDO
    gmosh=anorm*y

  END FUNCTION gmosh

  !recursive function factorial(n) result(res)
  !  implicit none
  !  integer  :: n
  !  real(dp) :: res
  !  
  !  if (n==0 .or. n==1) then
  !     res=1.0_dp
  !  else
  !     res=n*factorial(n-1)
  !  end if
  !end function factorial
 
END MODULE angmom_lib

! HARMONIC OSCILLATOR BRACKETS CALCULATION PROGRAM 
! Version 1.0: May 2001 
! E-mail: Gintautas_Kamuntavicius@fc.vdu.lt
! Reference: nucl-th/0105009
! WWW: http://www.nuclear.physics.vdu.lt  
!
! F95 VERSION: andreas ekstrom
!
! CALL TMB(EE,LL,ER,LR,E1,L1,E2,L2,LM,D)
!
! EE-centre of mass energy, LL-centre of mass angular moment
! ER-relative energy, LR-relative angular moment
! E1-energy of the first particle, L1-angular moment of the first particle
! E2-energy of the second particle, L2-angular moment of the second particle
! LM-total angular moment   

!Output: TMB-HARMONIC OSCILLATOR BRACKET!
MODULE tmb_brackets

  USE angmom_lib
  USE constants

  IMPLICIT NONE
  
  REAL(DP) :: BIN(0:99,0:99)
  REAL(DP) :: TIN(0:99,0:99,0:99)
  REAL(DP) :: T1, skirtm
  
CONTAINS
  
  SUBROUTINE INITIALIZE_TMB_BRACKETS
    CALL BINOM_TMB
    CALL TRINOM
  END SUBROUTINE INITIALIZE_TMB_BRACKETS
  
  FUNCTION TRI(I,J,K) RESULT(RES)
    !     TRIADIC CONDITION FOR MOMENTS I/2,J/2,K/2:
    !     I+J>=K, I+K>=J, J+K>=I,
    !     I/2+J/2+K/2 = INTEGER.
    !     TRI=1, WHEN TRIADIC CONDITION IS FULFILLED, TRI=0 OTHERWISE	 
    IMPLICIT NONE
    INTEGER :: I,J,K,L, RES
    RES=0
    L=I+J+K
    IF(L/2*2.NE.L) RETURN
    L=L/2
    IF((L-I)*(L-J)*(L-K).LT.0) RETURN
    RES=1
    
  END FUNCTION TRI

  FUNCTION GG(E1,L1,EA,LA,EB,LB) RESULT(RES)
    IMPLICIT NONE
    REAL(DP) RES
    INTEGER E1,L1,EA,LA,EB,LB
    
    RES=KL0(LA,LB,L1)*DSQRT(DFLOAT((2*LA+1)*(2*LB+1))*&
         TIN(E1-L1,EA-LA,EB-LB)*TIN(E1+L1+1,EA+LA+1,EB+LB+1))
  END FUNCTION GG
  
  FUNCTION C6J(I,J,K,L,M,N) RESULT(RES)
    !     6J - COEFFICIENT
    !     ( I/2  J/2  K/2 )
    !     ( L/2  M/2  N/2 ) 
    !     [JB 65] (22.1.4)
    IMPLICIT NONE
    INTEGER :: I,J,K,L,M,N,I1,I2,I3,I4,I5,IZ,JZ,KZ
    REAL(DP)  :: T,DZ,RES
    
    RES = 0.D0
    IF(TRI(I,J,K)*TRI(I,M,N)*TRI(J,L,N)*TRI(K,L,M).EQ.0) RETURN
    I1=(I+J+K)/2
    I2=(I+M+N)/2
    I3=(J+L+N)/2
    I4=(K+L+M)/2
    I5=(I+K+L+N)/2
    T=DSQRT(DFLOAT((I1+1)*(I4+1))/DFLOAT((I2+1)*(I3+1))* & 
         BIN(I2,I)*BIN(I,I2-N)*BIN(I4,L)/ & 
         BIN(I3,N)*BIN(L,I4-K)*BIN(I1,K)*BIN(K,I1-J)/&
         BIN(N,I3-L))/DFLOAT(I2-N+1)
    JZ=MAX0(0,(I+L-J-M)/2)
    DZ=1.D0
    IF((JZ+I5)/2*2.NE.(JZ+I5)) DZ=-1.D0
    KZ=MIN0(I2-M,I3-J,I5-K)
    DO IZ=JZ,KZ
       RES=RES+DZ*T*BIN(I2-M,IZ)/BIN(I2,I5-K-IZ)*&
            BIN(I2-I,I3-J-IZ)/BIN(I4-I3+J+IZ+1,I2-N+1)
       DZ=-DZ
    END DO
  
  END FUNCTION C6J
  
  FUNCTION C9J(J1,J2,J3,L1,L2,L3,K1,K2,K3) RESULT(RES)
    !     9J COEFICIENT
    !     (J1/2 J2/2 J3/2)
    !     (L1/2 L2/2 L3/2)
    !	(K1/2 K2/2 K3/2)  	 
    !     [JB 65] (24.33)
    USE angmom_lib
    
    IMPLICIT NONE
    REAL(DP)  :: RES
    INTEGER :: J1,J2,J3,L1,L2,L3,K1,K2,K3,I,J,K,L
    RES=0.D0
    L=TRI(J1,J2,J3)*TRI(L1,L2,L3)*TRI(K1,K2,K3)*&
         TRI(J1,L1,K1)*TRI(J2,L2,K2)*TRI(J3,L3,K3)
    IF(L.EQ.0) RETURN
    J=MAX0(IABS(J1-K3),IABS(J2-L3),IABS(L1-K2))
    K=MIN0(J1+K3,J2+L3,L1+K2)
    DO I=J,K,2
       RES=RES+DFLOAT(I+1)*C6J(J1,J2,J3,L3,K3,I)*&
            C6J(L1,L2,L3,J2,I,K2)*C6J(K1,K2,K3,I,J1,L1)
    END DO
    IF(J/2*2.NE.J) RES=-RES
    RETURN
  END FUNCTION C9J
  
  FUNCTION KL0(I,J,K) RESULT(RES)
    !	KLEBS-GORDAN COEFFICIENT WITH ZERO PROJECTIONS OF MOMENTA
    !	(I, J, K)
    !	(0, 0, 0)  
    !	I,J,K - MOMENTA = INTEGER NUMBERS
    !	[JB,65] (15.10)
    IMPLICIT NONE
    REAL(DP)  :: T, RES
    INTEGER :: I,J,K,L,M
    
    RES=0.D0
    IF(TRI(I,J,K).EQ.0) RETURN
    L=(I+J+K)/2
    M=L-K
    T=1.D0
    IF(M/2*2.NE.M) T=-1.D0
    RES=T*BIN(K,L-J)*BIN(L,K)/&
         DSQRT(BIN(2*K,2*(L-J))*BIN(2*L+1,2*K+1))
  
  END FUNCTION KL0
  
  ! Modified by Andreas:
  ! included (-)*(NN+NR+N1+N2) to comply with our
  ! definition of the radial wave function.
  ! DONT use the C9J function
  ! ninej in angmom_lib has superior
  ! precision, at the cost of a factor 
  ! of two in speed.... 

  FUNCTION TMB(EE,LL,ER,LR,E1,L1,E2,L2,LM,D) RESULT(RES)
    
    !     	   TALMI-MOSHINSKY BRACKET
    !	    (EE,LL;ER,LR:LM/E1,L1;E2,L2:LM)D
    IMPLICIT NONE
    REAL(DP) :: S,D,T,RES,PHASE
    INTEGER :: EE,LL,ER,LR,E1,L1,E2,L2,LM
    INTEGER :: M,ED,LD,EB,LB,EC,LC,EA,LA
    INTEGER :: NN, NR, N1, N2

    RES=0.D0
    IF(EE+ER.NE.E1+E2) RETURN
    IF(TRI(2*LL,2*LR,2*LM)*TRI(2*L1,2*L2,2*LM).EQ.0) RETURN
    T=DSQRT((D**(E1-ER))/((1.D0+D)**(E1+E2)))
    M=MIN0(ER,E2)
    S=1.D0
    DO ED=0,M
       EB=ER-ED
       EC=E2-ED
       EA=E1-ER+ED
       DO LD=ED,0,-2 
          DO LB=EB,0,-2
             IF(TRI(LD,LB,LR).EQ.0) CYCLE
             DO LC=EC,0,-2
                IF(TRI(LD,LC,L2).EQ.0) CYCLE
                DO LA=EA,0,-2
                   IF((TRI(LA,LB,L1).EQ.0).OR.(TRI(LA,LL,LC).EQ.0)) CYCLE
                   RES=RES+S*T* &
                        ninej(2*LA,2*LB,2*L1,2*LC,2*LD,2*L2,2*LL,2*LR,2*LM)* &
                        GG(E1,L1,EA,LA,EB,LB)*GG(E2,L2,EC,LC,ED,LD)* &
                        GG(EE,LL,EA,LA,EC,LC)*GG(ER,LR,EB,LB,ED,LD)
                END DO
             END DO
          END DO
       END DO
       S=S*(-D)
    END DO

    NN = (EE-LL)/2
    NR = (ER-LR)/2
    N1 = (E1-L1)/2
    N2 = (E2-L2)/2
    PHASE = (-1.0_dp)**(NN+NR+N1+N2)
    RES = RES*PHASE

    RETURN
  END FUNCTION TMB

  SUBROUTINE BINOM_TMB
    !	THE ARRAY OF BINOMIAL COEFFICIENTS
    !     BIN(I,J)= = I!/J!/(I-J)! 
    IMPLICIT NONE
    
    INTEGER :: I,K
    
    DO I=0,99
       BIN(I,0)=1.D0
       BIN(I,I)=1.D0
       DO K=1,I/2
          BIN(I,K)=DNINT(BIN(I,K-1)/DFLOAT(K)*DFLOAT(I-K+1))
          BIN(I,I-K)=BIN(I,K)
       END DO
    END DO
    RETURN
  END SUBROUTINE BINOM_TMB
  
  SUBROUTINE TRINOM
    !	THE ARRAY OF TRINOMIAL COEFFICIENTS
    !	TIN(I,J,K)=I!!/J!!/K!!
    IMPLICIT NONE
    
    INTEGER :: I,J,K,M,N
    TIN(0,0,0)=1.D0
    TIN(1,1,1)=1.D0
    DO I=2,99
       M=I-I/2*2
       TIN(I,I,M)=1.D0
       TIN(I,M,I)=1.D0
       N=M+2
       DO J=I,N,-2
          DO K=N,J,2
             TIN(I,J,K)=TIN(I,J,K-2)/DFLOAT(K)
             TIN(I,K,J)=TIN(I,J,K)
          END DO
          TIN(I,J-2,M)=TIN(I,J,M)*DFLOAT(J)
          TIN(I,M,J-2)=TIN(I,J-2,M)
       END DO
    END DO
    RETURN
  END SUBROUTINE TRINOM
  
END MODULE tmb_brackets
!
!     Special function's library
!     CONTAINS:
!       * radial H.O. wave functions, and its 1st+2nd derivatives
!         obtained using automatic differentiation
!
MODULE special_functions
  USE constants
  USE factorials
    
CONTAINS
  !
  !     H.O. functions using generalized Laguerres function
  !  
  REAL(DP) FUNCTION rnl_laguerre(n, l, z, oscl)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l,n
    REAL(DP),INTENT(IN) :: z, oscl
    REAL(DP)  :: cx(0:200), factor, zz,xp, ph
    
    ! MOMENTUM SPACE WF PHASE-CONVENTION
    !ph=(-1.D0)**n
    ! COORDINATE SPACE WF PHASE-CONVENTION
    ph = 1.0_dp
    factor = 0.5D0*((n+l+2)*LOG(2.D0)+fac(n)-dfac(2*n+2*l+1)-0.5D0*LOG(pi))
    factor = EXP(factor)
    
    zz= z*z/(oscl*oscl)
    
    CALL laguerre_general( n, l+0.5D0, zz, cx )
    
    xp = EXP(-zz*0.5D0)*((z/oscl)**l)*cx(n)
        
    rnl_laguerre = xp*ph*factor/(oscl**1.5_dp)
    
  END FUNCTION rnl_laguerre
  !
  !     H.O. functions using Kummers function   
  !     return value dimensionless, needs multiplication
  !     with oscl parameter**(-1.5)
  !
  ! input: n-radial number
  ! input: l-ang mom
  ! input: z-coordinate/oscl_parameter
  REAL(DP) FUNCTION rnl_kummer(n, l, z)
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: l, n
    REAL(DP), INTENT(IN) :: z
    INTEGER              :: lll, nn
    REAL(DP)             :: y, dl, gamfaa, dfll, gamfab, dfnn

    rnl_kummer=0.0_dp ; y=0.5_dp*z*z
    IF(y > 60.0_dp) RETURN
    dl = l
    IF((ABS(z) < 1.0d-6) .AND. (l == 0)) rnl_kummer = 1.0_dp
    IF( ABS(z) > 1.0d-6) rnl_kummer = (z**l) * EXP(-y) * hypkum(n,dl+1.5_dp,z*z)
    gamfaa = 0.5_dp * SQRT(pi)
    IF(l /= 0) THEN
       DO lll = 1, l
          dfll = lll - 1
          gamfaa = gamfaa * (dfll + 1.5_dp)
       ENDDO
    ENDIF
    gamfab = gamfaa
    IF(n /= 0) THEN
       dfll = dl + 0.5_dp
       DO nn = 1, n
          dfnn = nn
          gamfab = gamfab * ((dfnn + dfll) / dfnn)
       ENDDO
    ENDIF
   
    rnl_kummer = rnl_kummer * (SQRT(2.0_dp * gamfab) / gamfaa)
   
  END FUNCTION rnl_kummer
  !
  !     Kummers function, Abramowitz & Stegun   
  !     exp. 13.1.2. a(there) equals (-n)       
  !  
  REAL(DP) FUNCTION hypkum(n, b, z)
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: b, z
    INTEGER              :: nmax, nf
    REAL(DP)             :: af, bf, zf, term, dfnf, xadd, sum


    IF(n < 0) WRITE (*,*)' error exit in hypkum ',  n,b,z
    hypkum = 1.0_dp
    IF(n == 0) RETURN
    nmax = n ; af = - n ; bf = b ; zf = z ; sum = 1.0 ; term = 1.0_dp
    DO nf = 1, nmax
       dfnf = nf
       xadd = dfnf - 1.0_dp
       term = term * ((af + xadd) / (bf + xadd)) * (zf / dfnf)
       IF(ABS(term) <  1.0d-12) EXIT
       sum = sum + term
    ENDDO
    hypkum = sum

  END FUNCTION hypkum
  !
  !  This function sets up the recursive relation
  !  for the associated Legendre polynomials
  !
  REAL(DP) FUNCTION legendre_polynomials(l, m, x)
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: l, m
    REAL(DP), INTENT(IN) :: x
    REAL(DP)             :: fact,pll,pmm,pmmp1,somx2
    INTEGER              :: i,ll
    

    !  check whether m, l and x are ok
    IF((M < 0).OR.(M > L).OR.(ABS(X) > 1.0_dp)) THEN
       WRITE(*,*) 'legendre_polynomials: bad arguments', m, l, x; RETURN
    ENDIF

    !  calculate now pmm as starting point for iterations
    pmm=1.0
    IF (m > 0) THEN
       somx2=SQRT((1.0-x)*(1.0+x))
       fact=1.0_dp;
       DO i=1, m
          pmm = -fact*somx2*pmm
          fact = fact+2.0_dp
       ENDDO
    ENDIF

    !  if l == m we do not need to use recursion relation
    IF (l == m) THEN
       legendre_polynomials=pmm

       !  recursive relation for associated Legendre polynomials
    ELSE
       pmmp1=x*(2*m+1)*pmm

       !  analytical formula for the case l == m+1
       IF (l == (m+1)) THEN
          legendre_polynomials=pmmp1
       ELSE
          DO ll=m+2, l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          legendre_polynomials= pll
       ENDIF
    ENDIF
    
  END FUNCTION legendre_polynomials
  
  SUBROUTINE laguerre_general( n, alpha, x, cx )
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL (dp ) ::  alpha
    REAL ( dp ) :: cx(0:n)
    INTEGER :: i
    REAL ( dp ), INTENT(IN) ::  x

    IF ( alpha <= -1.0D+00 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'LAGUERRE_GENERAL - Fatal error!'
       WRITE ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
       WRITE ( *, '(a)' ) '  but ALPHA must be greater than -1.'
       STOP
    END IF
    IF ( n < 0 ) THEN
       RETURN
    END IF
    cx(0) = 1.0D+00
    IF ( n == 0 ) THEN
       RETURN
    END IF
    cx(1) = 1.0D+00 + alpha - x
    DO i = 2, n
       cx(i) = ( ( REAL ( 2 * i - 1, kind = 8 ) + alpha - x ) * cx(i-1)   &
            + ( REAL (   - i + 1, kind = 8 ) - alpha     ) * cx(i-2) ) &
            / REAL (     i,     kind = 8 )
    END DO

  END SUBROUTINE laguerre_general
  
  SUBROUTINE legendre_second_kind_polynomials (n, x, qn, qd)
    !
    !       ====================================================
    !       Purpose: Compute Legendre functions Qn(x) their derivatives,  Qn'(x)
    !       Input :  x  --- Argument of Qn(x)
    !                n  --- Degree of Qn(x)  ( n = 0,1,2,...)
    !       Output:  QN(n) --- Qn(x)
    !                QD(n) --- Qn'(x)
    !       ====================================================
    !
    !       Examples:     x1 = 0.50,    x2 = 2.50
    !
    !       n      Qn(x1)        Qn'(x1)       Qn(x2)          Qn'(x2)
    !     ----------------------------------------------------------------
    !       0     .54930614    1.33333333   .42364893D+00  -.19047619D+00
    !       1    -.72534693    1.21597281   .59122325D-01  -.52541546D-01
    !       2    -.81866327    -.84270745   .98842555D-02  -.13109214D-01
    !       3    -.19865477   -2.87734353   .17695141D-02  -.31202687D-02
    !       4     .44017453   -2.23329085   .32843271D-03  -.72261513D-03
    !       5     .55508089    1.08422720   .62335892D-04  -.16437427D-03
    !       ===============================================================
    !
    
    IMPLICIT NONE
    
    REAL(DP), INTENT(IN) :: X
    INTEGER,  INTENT(IN) :: N
    REAL(DP), INTENT(INOUT), DIMENSION(0:N) :: QN, QD
    
    INTEGER  :: J, K, L, NL
    REAL(DP) :: EPS, X2, Q0, Q1, QC2, QC1
    REAL(DP) :: QR, QF, QF0, QF1, QF2
    
    ! this method doesn't handle case N=0 correctly
    IF (N == 0) STOP 'N = 0 in legendre_second_kind_polynomials.'
    
    EPS=1.0D-14
    IF (DABS(X).EQ.1.0D0) THEN
       DO K=0,N
          QN(K)=1.0D+300
          QD(K)=1.0D+300
       ENDDO
    ENDIF

    IF (X .LE. 1.021_dp) THEN
       X2 = DABS ((1.0_dp + X)/(1.0_dp-X))
       Q0 = .5_dp * DLOG(X2)
       Q1 = X*Q0 - 1.0_dp
       QN(0) = Q0
       QN(1) = Q1
       QD(0) = 1.0D0/(1.0D0-X*X)
       QD(1) = QN(0)+X*QD(0)

       DO K=2,N
          QF=((2.0D0*K-1.0D0)*X*Q1-(K-1.0D0)*Q0)/K
          QN(K)=QF
          QD(K)=(QN(K-1)-X*QF)*K/(1.0D0-X*X)
          Q0=Q1
          Q1=QF
       ENDDO

    ELSE
       
       QC2=1.0D0/X
       DO J=1,N
          QC2=QC2*J/((2.0*J+1.0D0)*X)
          IF (J.EQ.N-1) QC1=QC2
       ENDDO
       
       DO L=0,1
          NL=N+L
          QF=1.0D0
          QR=1.0D0
          DO K=1,500
             QR=QR*(0.5D0*NL+K-1.0D0)*(0.5D0*(NL-1)+K)/((NL+K-0.5D0)*K*X*X)
             QF=QF+QR
             IF (DABS(QR/QF).LT.EPS) EXIT
          ENDDO
          
          IF (L.EQ.0) THEN
             QN(N-1)=QF*QC1
          ELSE
             QN(N)=QF*QC2
          ENDIF

       ENDDO

       QF2=QN(N)
       QF1=QN(N-1)

       DO K=N,2,-1
          QF0=((2*K-1.0D0)*X*QF1-K*QF2)/(K-1.0D0)
          QN(K-2)=QF0
          QF2=QF1
          QF1=QF0
       ENDDO

       QD(0)=1.0D0/(1.0D0-X*X)
       DO K=1,N
          QD(K)=K*(QN(K-1)-X*QN(K))/(1.0D0-X*X)
       ENDDO

    ENDIF
    RETURN

  END SUBROUTINE legendre_second_kind_polynomials
  !
  !
  !      This routine calculates gauss-legendre mesh points and weights      
  !      input:                                                              
  !      x1   : lower limit of the integration interval                      
  !      x2   : upper limit ---------- "" -------------                      
  !      n    : the desired number of mesh points                            
  !      output :                                                            
  !      x     : gauss-legendre mesh points on the interval (x1,x2)          
  !      w     : the corresponding weights                                   
  !      From  : Numerical recipes
  !      F90 version : M. Hjorth-Jensen
  !
  SUBROUTINE gauss_legendre(x1, x2, x, w, n)
    IMPLICIT NONE
    INTEGER,                  INTENT(IN)    :: n
    REAL(DP),                 INTENT(IN)    :: x1, x2
    REAL(DP), DIMENSION(1:n), INTENT(INOUT) :: x(n), w(n)
    INTEGER                                 :: i, j, m
    REAL(DP)                                :: p1,p2,p3,pp,xl,xm,z,z1
    REAL(DP), PARAMETER                     :: eps = 3.D-14
    
    m=(n+1)/2
    xm=0.5_dp*(x2+x1)
    xl=0.5_dp*(x2-x1)
    DO i=1,m
       z1=0.
       z=COS(pi*(i-.25_dp)/(n+.5_dp))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.0_dp
          p2=0.0_dp
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.0_dp*j-1.)*z*p2-(j-1.0_dp)*p3)/j
          ENDDO
          pp=n*(z*p1-p2)/(z*z-1.)
          z1=z
          z=z-p1/pp
       ENDDO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO

  END SUBROUTINE gauss_legendre

  recursive function lacz_gamma(a) result(g)
    
    IMPLICIT NONE
    
    REAL(DP), INTENT(IN) :: a
    REAL(DP) :: g
    
    REAL(DP), parameter :: pi = 3.14159265358979324
    INTEGER, parameter :: cg = 7
    
    ! these precomputed values are taken by the sample code in Wikipedia,
    ! and the sample itself takes them from the GNU Scientific Library
    REAL(DP), dimension(0:8), parameter :: p = &
         (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
         771.32342877765313, -176.61502916214059, 12.507343278686905, &
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)
    
    REAL(DP) :: t, w, x
    INTEGER :: i
    
    x = a
    
    if ( x < 0.5_dp ) then
       g = pi / ( sin(pi*x) * lacz_gamma(1.0_dp-x) )
    else
       x = x - 1.0_dp
       t = p(0)
       do i=1, cg+1
          t = t + p(i)/(x+real(i,kind=8))
       end do
       w = x + real(cg) + 0.5_dp
       g = sqrt(2.0_dp*pi) * w**(x+0.5_dp) * exp(-w) * t
    end if
  end function lacz_gamma
  
  FUNCTION gammln(XX) RESULT(res)
    
   
    IMPLICIT NONE
    
    REAL*8 COF(6),STP,HALF,ONE,FPF,X,XX,TMP,SER
    INTEGER J
    REAL*8 res

    DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0, &
         -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
    DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/

    X=XX-ONE
    TMP=X+FPF
    TMP=(X+HALF)*LOG(TMP)-TMP
    SER=ONE
    DO J=1,6
       X=X+ONE
       SER=SER+COF(J)/X
    END DO
    
    res=TMP+LOG(STP*SER)
    
  END FUNCTION gammln
  
  ! faster but less precise spherical bessel calculator
  REAL*8 FUNCTION BESL(RU,L)
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION BES(60,2)   
    MM=L+1
    XT=RU**MM-0.01D0
    IF(XT.LT.0.0D0.OR.XT.EQ.0.0D0)THEN
       IF(L.EQ.0)THEN
          BESL=1.0D0
          RETURN
       ELSE
          BESL=RU**L
          KFC=2*L+1
          DNO=1
          DO J=1,KFC,2
             FJ=J
             DNO=DNO*FJ
          ENDDO
          BESL=BESL/DNO
          RETURN
       ENDIF
    ELSEIF(XT.GT.0.0)THEN
       I=2
       BES(1,I)=DSIN(RU)
       BES(2,I)=DSIN(RU)/RU-DCOS(RU)
       M=L+1
       KKK= MAX0(M,3)
       IF(KKK.LT.3)THEN
          BJ=BES(L+1,2)
          BESL=BJ/RU
          RETURN
       ELSE
          DO LL=3,KKK
             A=LL-1
             L1=LL-1
             L2=LL-2
             BES(LL,I)=(2.*A-1.D0)/RU*BES(L1,I)-BES(L2,I)
          ENDDO
       ENDIF
       BJ=BES(L+1,2)
       BESL=BJ/RU
    ENDIF
    RETURN
  END FUNCTION BESL
  
  !---------------------------------------------------------------------
  SUBROUTINE SBESJY  (X,LMAX, J,Y,JP,YP, IFAIL )
    !---------------------------------------------------------------------
    !   REAL SPHERICAL BESSEL FUNCTIONS AND X DERIVATIVES
    !            j , y , j', y'                    FROM L=0 TO L=LMAX
    !        FOR REAL X > SQRT(ACCUR) (E.G. 1D-7)    AND INTEGER LMAX
    ! 
    !  J (L)  =      j/L/(X) STORES   REGULAR SPHERICAL BESSEL FUNCTION:
    !  JP(L)  = D/DX j/L/(X)            j(0) =  SIN(X)/X
    !  Y (L)  =      y/L/(X) STORES IRREGULAR SPHERICAL BESSEL FUNCTION:
    !  YP(L)  = D/DX y/L/(X)            y(0) = -COS(X)/X
    !                                                
    !    IFAIL = -1 FOR ARGUMENTS OUT OF RANGE
    !          =  0 FOR ALL RESULTS SATISFACTORY
    ! 
    !   USING LENTZ-THOMPSON EVALUATION OF CONTINUED FRACTION CF1,
    !   AND TRIGONOMETRIC FORMS FOR L = 0 SOLUTIONS.
    !   LMAX IS LARGEST L NEEDED AND MUST BE <= MAXL, THE ARRAY INDEX.
    !   MAXL CAN BE DELETED AND ALL THE ARRAYS DIMENSIONED (0:*)
    !   SMALL IS MACHINE DEPENDENT, ABOUT SQRT(MINIMUM REAL NUMBER),
    !         SO 1D-150 FOR DOUBLE PRECISION ON VAX, PCS ETC.
    !   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
    !   IN OSCILLATING REGION X .GE.  [ SQRT{LMAX*(LMAX+1)} ]
    !   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8
    !   IS THE SMALLEST NUMBER WITH 1+ACC8.NE.1 FOR OUR WORKING PRECISION
    !   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
    !   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
    !   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADILY WORSE:
    !   THE VARIABLE ERR IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
    !
    !   NOTE: FOR X=1 AND L=100  J = 7.4 E-190     Y = -6.7+E186    1.4.94
    !---------------------------------------------------------------------
    !   AUTHOR :   A.R.BARNETT       MANCHESTER    12 MARCH 1990/95
    !                                AUCKLAND      12 MARCH 1991
    !---------------------------------------------------------------------  
    IMPLICIT    NONE
    INTEGER     LIMIT,         MAXL,       LMAX, IFAIL, NFP, L
    PARAMETER ( LIMIT = 20000, MAXL = 250 )
    DOUBLE PRECISION  J(0:MAXL), Y(0:MAXL), JP(0:MAXL), YP(0:MAXL)
    DOUBLE PRECISION  ZERO,ONE,TWO,THREE,SMALL, ACCUR, TK,SL, ERR
    DOUBLE PRECISION  X,XINV, CF1,DCF1, DEN, C,D, OMEGA, TWOXI
    PARAMETER ( ZERO  = 0.0D0  , ONE   = 1.0D0 , TWO = 2.0D0 )
    PARAMETER ( SMALL = 1.D-150, THREE = 3.0D0 )
    COMMON /STEDE/    ERR,NFP       ! not required in code        
    !-------

    J(:)  = 0.0D0
    Y(:)  = 0.0D0
    JP(:) = 0.0D0
    YP(:) = 0.0D0

    ACCUR = 1.D-14                  ! suitable for double precision
    IFAIL = -1                      ! user to check on exit
    IF (X .LT. DSQRT(ACCUR) )       GOTO 50
    !-------TAKES CARE OF NEGATIVE X ... USE REFLECTION FORMULA
    !-------BEGIN CALCULATION OF CF1 UNLESS LMAX = 0, WHEN SOLUTIONS BELOW
    XINV  = ONE / X
    IF (LMAX .GT. 0) THEN
       TWOXI =     XINV + XINV
       SL  =  REAL(LMAX)* XINV     ! used also in do loop 3
       TK  =  TWO * SL  + XINV * THREE     
       CF1 =  SL                   ! initial value of CF1
       DEN =  ONE                  ! unnormalised j(Lmax,x)
       IF ( ABS(CF1) .LT. SMALL ) CF1 = SMALL
       C   = CF1                   ! inverse ratio of A convergents
       D   = ZERO                  ! direct  ratio of B convergents   
       DO L = 1,LIMIT
          C   = TK - ONE / C
          D   = TK - D
          IF ( ABS(C) .LT. SMALL ) C = SMALL
          IF ( ABS(D) .LT. SMALL ) D = SMALL
          D   = ONE / D
          DCF1= D   * C
          CF1 = CF1 * DCF1
          IF ( D .LT. ZERO ) DEN = - DEN
          IF ( ABS(DCF1 - ONE) .LE. ACCUR ) GOTO 20
          TK   = TK + TWOXI
          NFP  = L                 ! ie number in loop
       END DO
       GOTO 50             ! error exit, no convergence
20     CONTINUE
       ERR = ACCUR * DSQRT(DBLE(NFP))    ! error estimate
       J (LMAX) = DEN      ! lower-case j's  really
       JP(LMAX) = CF1 * DEN                         
       !------ DOWNWARD RECURSION TO L=0  AS SPHERICAL BESSEL FUNCTIONS
       DO L =  LMAX , 1, -1
          J (L-1)  = (SL + XINV) * J(L)   + JP(L)
          SL  =  SL - XINV
          JP(L-1)  =  SL * J(L-1)          - J(L)
       END DO
       DEN = J(0)
    ENDIF                           ! end loop for Lmax GT 0
    !------ CALCULATE THE L=0 SPHERICAL BESSEL FUNCTIONS DIRECTLY
    J (0)   =  XINV * DSIN(X)
    Y (0)   = -XINV * DCOS(X)
    JP(0)   = -Y(0) - XINV * J(0)
    YP(0)   =  J(0) - XINV * Y(0)
    IF (LMAX .GT. 0) THEN
       OMEGA  =  J(0) / DEN
       SL  = ZERO
       DO L = 1 , LMAX
          J (L) = OMEGA * J (L)
          JP(L) = OMEGA * JP(L)
          Y (L) = SL * Y(L-1)   -   YP(L-1)
          SL  = SL + XINV
          YP(L) = Y(L-1)  -  (SL + XINV) * Y(L)
       END DO
    ENDIF
    IFAIL = 0                       ! calculations successful
    RETURN
    !---------------------------------------------------------------------
    !       ERROR TRAPS
    !---------------------------------------------------------------------
50  IF (X .LT. ZERO) THEN
       WRITE(6,1000) X
    ELSEIF (X .EQ. ZERO) THEN
       IFAIL = 0
       J(0) = ONE
       DO L = 1, LMAX
          J(L) = ZERO     ! remaining arrays untouched
       END DO
    ELSE                          ! x .le. sqrt(accur), e.g. 1D-7
       WRITE(6,1001) X
    ENDIF
1000 FORMAT(' X NEGATIVE !',1PE15.5,'    USE REFLECTION FORMULA'/)
1001 FORMAT(' WITH X = ',1PE15.5,'    TRY SMALL-X SOLUTIONS',&
          X  '    j/L/(X)  ->   X**L / (2L+1)!!          AND', &
          X  '    y/L/(X)  ->  -(2L-1)!! / X**(L+1)'/)
    RETURN
  END SUBROUTINE SBESJY
  !---------------------------------------------------------------------
  !       END OF SUBROUTINE SBESJY 
  !---------------------------------------------------------------------
  !
  !  function interface to SBESJY
  !  SPHERICAL BESSEL FUNCTION
  !
  FUNCTION sbessj(N,X) RESULT(res)

    USE constants
    USE factorials

    IMPLICIT NONE
    
    INTEGER , INTENT(IN) :: N
    REAL(DP), INTENT(IN) :: X
    REAL(DP)             :: res
    INTEGER              :: IFAIL
    REAL(DP), DIMENSION(0:250) :: bess_j, bess_y, bess_jp, bess_yp
    
    IF (X<=1.0E-006) THEN
       IF (2*N+1 > dbl_maxvalue) THEN
          res = 0.0_dp
       ELSE
          res = X**N/double_factorial(2*N+1)
       END IF
       IF (isnan(res)) THEN
          WRITE(SCR,*) X, N, double_factorial(2*N+1), X**N
          WRITE(ERR,*) 'error(sbessj): detected NaN'
          STOP
       END IF
       RETURN
    END IF

    IF (X<=1.0_dp .AND. N>20) THEN
       IF (2*N+1 > dbl_maxvalue) THEN
          res = 0.0_dp
       ELSE
          res = X**N/double_factorial(2*N+1)
       END IF
       IF (isnan(res)) THEN
          WRITE(SCR,*) X, N, double_factorial(2*N+1), X**N
          WRITE(ERR,*) 'error(sbessj): detected NaN'
          STOP
       END IF
       RETURN
    END IF    

    IF (X<=2.0_dp .AND. N>100) THEN
       IF (2*N+1 > dbl_maxvalue) THEN
          res = 0.0_dp
       ELSE
          res = X**N/double_factorial(2*N+1)
       END IF
       IF (isnan(res)) THEN
          WRITE(SCR,*) X, N, double_factorial(2*N+1), X**N
          WRITE(ERR,*) 'error(sbessj): detected NaN'
          STOP
       END IF
       RETURN
    END IF    
    
    CALL SBESJY(X, N, bess_j, bess_y, bess_jp, bess_yp, IFAIL)
    IF (IFAIL == -1) THEN
       WRITE(ERR,*) 'error(sbessj): arguments out of range ',  IFAIL
       STOP
    ELSE
       res = bess_j(N)
       IF (isnan(res)) THEN
          WRITE(ERR,*) 'error(sbessj): spherical bessel Not a Number', res
          WRITE(ERR,*) 'asymptotic form limit might be set too low!'
          WRITE(SCR,*) X, N
          STOP
       END IF
    END IF
    
  END FUNCTION sbessj
  !
  ! FIRST KIND BESSEL FUNCTION
  !
  FUNCTION bessj (N,X) RESULT(res)
    
    USE constants
    
    IMPLICIT NONE
    
    INTEGER , INTENT(IN) :: N
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: res
    REAL(DP) TOX,BJM,BJ,BJP,SUM
    
    INTEGER, PARAMETER :: IACC = 40
    REAL(DP),PARAMETER :: BIGNO = 1.D10
    REAL(DP),PARAMETER :: BIGNI = 1.D-10
    
    INTEGER :: J, M, JSUM 

    IF (N.EQ.0) THEN
       res = bessj0(X)
       RETURN
    ENDIF
    
    IF (N.EQ.1) THEN
       res = bessj1(X)
       RETURN
    ENDIF

    IF (X.EQ.0.) THEN
       res = 0.
       RETURN
    ENDIF
    TOX = 2./X
    IF (X.GT.FLOAT(N)) THEN
       BJM = BESSJ0(X)
       BJ  = BESSJ1(X)
       DO J = 1,N-1
          BJP = J*TOX*BJ-BJM
          BJM = BJ
          BJ  = BJP
       END DO
       
       res = BJ
    ELSE
       M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
       res = 0.
       JSUM = 0
       SUM = 0.
       BJP = 0.
       BJ  = 1.
       DO J = M,1,-1
          BJM = J*TOX*BJ-BJP
          BJP = BJ
          BJ  = BJM
          IF (ABS(BJ).GT.BIGNO) THEN
             BJ  = BJ*BIGNI
             BJP = BJP*BIGNI
             res = res*BIGNI
             SUM = SUM*BIGNI
          ENDIF
          IF (JSUM.NE.0) SUM = SUM+BJ
          JSUM = 1-JSUM
          IF (J.EQ.N) res = BJP
       END DO
       SUM = 2.*SUM-BJ
      res = res/SUM
   ENDIF
   RETURN
   
 END FUNCTION bessj
 !
 !     This subroutine calculates the First Kind Bessel Function of
 !     order 1, for any real number X. The polynomial approximation by
 !     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
 !     REFERENCES:
 !     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
 !     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
 !     VOL.5, 1962. 
 FUNCTION BESSJ1 (X) RESULT(res)
   
   USE constants 
   
   IMPLICIT NONE
   REAL(DP), INTENT(IN) :: x
   REAL(DP) :: res
   REAL(DP) AX,FR,FS,Z,FP,FQ,XX
   
   REAL(DP) Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
        ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
   DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
        .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
   DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
        .8449199096D-5,-.88228987D-6,.105787412D-6 /
   DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, & 
        242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
   DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
        18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /
   
   AX = ABS(X)
   IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      res = X*(FR/FS)
   ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-2.35619491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      res = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
   ENDIF
   RETURN
 END FUNCTION BESSJ1
 !     This subroutine calculates the First Kind Bessel Function of
 !     order 0, for any real number X. The polynomial approximation by
 !     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
 !     REFERENCES:
 !     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
 !     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
 !     VOL.5, 1962.
 FUNCTION bessj0(X) RESULT(RES)
   
   USE constants
   
   IMPLICIT NONE
   
   REAL(DP), INTENT(IN) :: X
   REAL(DP) :: RES
   REAL(DP) :: AX,FR,FS,Z,FP,FQ,XX
   
   REAL(DP) Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
        ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
   DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
        -.2073370639D-5,.2093887211D-6 /
   DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
        -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
   DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
        651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
   DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
        9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
   
   RES = 1.0D0
   
   IF(X.EQ.0.D0) RETURN
   AX = ABS (X)
   
   IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      RES = FR/FS
   ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-.785398164
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      RES = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
   ENDIF
   
   RETURN
   
 END FUNCTION bessj0

 !----------------------------------------------------------------------
      SUBROUTINE COUL90(X, ETA, XLMIN,LRANGE, FC,GC,FCP,GCP, KFN,IFAIL)
!----------------------------------------------------------------------
!
!  COULOMB & BESSEL FUNCTION PROGRAM-- COUL90 -- USING STEED'S METHOD
!
!  COUL90 RETURNS ARRAYS FC = F, GC = G, FCP = (D/DX) F, GCP = (D/DX) G
!   FOR REAL X .GT. 0. ,REAL ETA (INCLUDING 0.), AND REAL XLMIN .GT.-1.
!   FOR (LRANGE+1) INTEGER-SPACED LAMBDA VALUES.
!   IT HENCE GIVES POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER
!   EQUATION, TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF
!   THE DIRAC EQUATION.    BY SETTING ETA = 0.0 AND RENORMALISING
!   SPHERICAL & CYLINDRICAL BESSEL FUNCTIONS ARE COMPUTED INSTEAD.
!----------------------------------------------------------------------
!   CALLING VARIABLES; ALL REALS ARE DOUBLE PRECISION (REAL*8)
!
!   X       - REAL ARGUMENT FOR COULOMB FUNCTIONS > 0.0 
!             [ X > SQRT(ACCUR) : ACCUR IS TARGET ACCURACY 1.0D-14 ]
!   ETA     - REAL SOMMERFELD PARAMETER, UNRESTRICTED > = < 0.0
!   XLMIN   - REAL MINIMUM LAMBDA-VALUE (L-VALUE OR ORDER),
!             GENERALLY IN RANGE 0.0 - 1.0 AND MOST USUALLY 0.0
!   LRANGE  - INTEGER NUMBER OF ADDITIONAL L-VALUES : RESULTS RETURNED
!             FOR L-VALUES XLMIN TO XLMIN + LRANGE INCLUSIVE
!   FC ,GC  - REAL VECTORS F,G OF REGULAR, IRREGULAR COULOMB FUNCTIONS
!   FCP,GCP - REAL VECTORS FOR THE X-DERIVATIVES OF  F,G
!             THESE VECTORS TO BE OF LENGTH AT LEAST MINL + LRANGE
!             STARTING ELEMENT MINL = MAX0( IDINT(XLMIN+ACCUR),0 )
!   KFN     - INTEGER CHOICE OF FUNCTIONS TO BE COMPUTED :
!           = 0         REAL COULOMB FUNCTIONS AND DERIVATIVES F & G
!           = 1    SPHERICAL BESSEL      "      "     "        j & y
!           = 2  CYLINDRICAL BESSEL      "      "     "        J & Y
!
!   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
!   IN OSCILLATING REGION X .GE. [ETA + SQRT{ETA**2 + XLM*(XLM+1)}]
!   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8 IS
!   THE SMALLEST NUMBER WITH 1.+ACC8.NE.1. FOR OUR WORKING PRECISION.
!   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
!   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
!   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADILY WORSE:
!   THE VARIABLE PACCQ IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
!----------------------------------------------------------------------
!   ERROR RETURNS                THE USER SHOULD TEST IFAIL ON EXIT
!
!   IFAIL ON INPUT IS SET TO 0                        LIMIT = 20000
!   IFAIL IN OUTPUT =  0 : CALCULATIONS SATISFACTORY
!                   =  1 : CF1 DID NOT CONVERGE AFTER LIMIT ITERATIONS
!                   =  2 : CF2 DID NOT CONVERGE AFTER LIMIT ITERATIONS
!                   = -1 : X < 1D-7 = SQRT(ACCUR)
!                   = -2 : INCONSISTENCY IN ORDER VALUES (L-VALUES) 
!----------------------------------------------------------------------
!  MACHINE-DEPENDENT PARAMETERS:    ACCUR - SEE ABOVE
!           SMALL - OFFSET FOR RECURSION = APPROX SQRT(MIN REAL NO.)
!           IE 1D-30 FOR IBM REAL*8,    1D-150 FOR DOUBLE PRECISION
!----------------------------------------------------------------------
!  PROGRAMMING HISTORY AND BACKGROUND: CPC IS COMPUTER PHYSICS COMMUN.
!  ORIGINAL PROGRAM  RCWFN       IN    CPC  8 (1974) 377-395
!                 +  RCWFF       IN    CPC 11 (1976) 141-142
!  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314
!  REVISED STANDARD  COULFG      IN    CPC 27 (1982) 147-166
!  BACKGROUND MATERIAL IN J. COMP. PHYSICS 46 (1982) 171-188         
!  CURRENT PROGRAM   COUL90  (FORTRAN77) SUPERCEDES THESE EARLIER ONES
!  (WHICH WERE WRITTEN IN FORTRAN 4) AND ALSO BY INCORPORATING THE NEW
!  LENTZ-THOMPSON ALGORITHM FOR EVALUATING THE FIRST CONTINUED FRACTION
!  ..SEE ALSO APPENDIX TO J. COMP. PHYSICS 64 (1986) 490-509     1.4.94
!----------------------------------------------------------------------
!  AUTHOR: A. R. BARNETT           MANCHESTER  MARCH   1981/95
!                                  AUCKLAND    MARCH   1991
!----------------------------------------------------------------------
      IMPLICIT         NONE
      INTEGER          LRANGE, KFN, IFAIL
      DOUBLE PRECISION X, ETA, XLMIN
      DOUBLE PRECISION FC (0:*), GC (0:*), FCP(0:*), GCP(0:*)
!----- ARRAYS INDEXED FROM 0 INSIDE SUBROUTINE: STORAGE FROM MINL
      DOUBLE PRECISION ACCUR,ACCH,SMALL, ONE,ZERO,HALF,TWO,TEN2, RT2DPI
      DOUBLE PRECISION XINV,PK,CF1,C,D,PK1,ETAK,RK2,TK,DCF1,DEN,XLM,XLL
      DOUBLE PRECISION EL,XL,RL,SL, F,FCMAXL,FCMINL,GCMINL,OMEGA,WRONSK
      DOUBLE PRECISION WI, A,B, AR,AI,BR,BI,DR,DI,DP,DQ, ALPHA,BETA
      DOUBLE PRECISION E2MM1, FJWKB,GJWKB, P,Q,PACCQ, GAMMA,GAMMAI
      INTEGER          IEXP, NFP, NPQ, L, MINL,MAXL, LIMIT
      LOGICAL          ETANE0, XLTURN
      PARAMETER      ( LIMIT = 20000, SMALL = 1.0D-150 )
      COMMON  /STEED/  PACCQ,NFP,NPQ,IEXP,MINL    !not required in code
      COMMON  /DESET/  CF1,P,Q,F,GAMMA,WRONSK     !information only
!----------------------------------------------------------------------
!     COUL90 HAS CALLS TO: DSQRT,DABS,MAX0,IDINT,DSIGN,DFLOAT,DMIN1
!----------------------------------------------------------------------
      DATA ZERO,ONE,TWO,TEN2,HALF /0.0D0, 1.0D0, 2.0D0, 1.0D2, 0.5D0/
      DATA RT2DPI /0.797884560802865D0/ 
!Q    DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 Q0/
!-----THIS CONSTANT IS  DSQRT(TWO / PI):
!-----USE Q0 FOR IBM REAL*16: D0 FOR REAL*8 AND DOUBLE PRECISION
!----------------CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
                        ACCUR = 1.0D-14
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO                                
      ACCH  = DSQRT(ACCUR)
!-----   TEST RANGE OF X, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
                IF( X .LE. ACCH )                GO TO 100
      IF( KFN.EQ.2 )   THEN
         XLM = XLMIN - HALF                                  
        ELSE
         XLM = XLMIN                                                     
        ENDIF
      IF( XLM.LE.-ONE .OR. LRANGE.LT.0 )         GO TO 105 
      E2MM1  = XLM * XLM + XLM
      XLTURN = X * (X -  TWO * ETA) .LT. E2MM1
      E2MM1  = E2MM1  +  ETA * ETA
      XLL    = XLM + DFLOAT(LRANGE)
!-----  LRANGE IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
!-----  XLL  IS MAX LAMBDA VALUE [ OR 0.5 SMALLER FOR J,Y BESSELS ]
!-----  DETERMINE STARTING ARRAY ELEMENT (MINL) FROM XLMIN
      MINL  = MAX0( IDINT(XLMIN + ACCUR),0 )     ! index from 0
      MAXL  = MINL + LRANGE
!-----   EVALUATE CF1  =  F   =  DF(L,ETA,X)/DX   /   F(L,ETA,X)
      XINV = ONE / X
      DEN  = ONE                       ! unnormalised F(MAXL,ETA,X)
      PK   = XLL + ONE
      CF1  = ETA / PK  +  PK * XINV                                             
           IF( DABS(CF1).LT.SMALL )    CF1 = SMALL
      RK2  = ONE
         D = ZERO
         C = CF1
!----- BEGIN CF1 LOOP ON PK = K STARTING AT LAMBDA + 1: LENTZ-THOMPSON
      DO 10 L =  1 , LIMIT             ! abort if reach LIMIT (20000)    
          PK1 = PK + ONE
          IF( ETANE0 ) THEN
                ETAK = ETA / PK
                RK2  = ONE + ETAK * ETAK
                 TK  = (PK + PK1) * (XINV + ETAK / PK1)
             ELSE
                 TK  = (PK + PK1) * XINV
             ENDIF
          D   =  TK - RK2 * D          ! direct  ratio of B convergents    
          C   =  TK - RK2 / C          ! inverse ratio of A convergents
            IF( DABS(C).LT.SMALL ) C = SMALL
            IF( DABS(D).LT.SMALL ) D = SMALL
          D   = ONE / D
          DCF1=   D * C
          CF1 = CF1 * DCF1
              IF( D.LT.ZERO )    DEN = -DEN
          PK  = PK1
          IF( DABS(DCF1-ONE).LT.ACCUR )     GO TO  20 ! proper exit
   10 CONTINUE
                                            GO TO 110 ! error exit 
   20       NFP = PK - XLL - 1                        ! number of steps
              F = CF1                                 ! need DEN later
!----DOWNWARD RECURRENCE TO LAMBDA = XLM; ARRAYS GC, GCP STORE RL, SL
      IF( LRANGE.GT.0 )       THEN
          FCMAXL    = SMALL  * DEN 
          FCP(MAXL) = FCMAXL * CF1
          FC (MAXL) = FCMAXL
                    XL = XLL                   
                    RL = ONE
          DO 30 L =  MAXL, MINL+1, -1
             IF( ETANE0 )  THEN
                    EL = ETA / XL                
                    RL = DSQRT( ONE + EL * EL )
                    SL = XL * XINV  + EL
                    GC (L) = RL                  ! storage
                    GCP(L) = SL
                ELSE
                    SL = XL * XINV
                ENDIF
             FC (L-1)  = ( FC(L)   * SL  +  FCP(L) ) / RL
             FCP(L-1)  =   FC(L-1) * SL  -  FC (L) * RL
             XL    =  XL - ONE                   ! end value is XLM
   30     CONTINUE
         IF( DABS(FC(MINL)).LT.ACCUR*SMALL )  FC(MINL) = ACCUR * SMALL
          F   = FCP(MINL) / FC(MINL)             ! F'/F at min L-value
          DEN = FC (MINL)                        ! normalisation
      ENDIF
!---------------------------------------------------------------------
!-----   NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
!-----   EVALUATE CF2 = P + I.Q  USING STEED'S ALGORITHM (NO ZEROS)
!---------------------------------------------------------------------
      IF( XLTURN ) CALL JWKB( X,ETA,DMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP )
      IF( IEXP.GT.1 .OR. GJWKB.GT.(ONE / (ACCH*TEN2)) ) THEN
          OMEGA = FJWKB
          GAMMA = GJWKB * OMEGA
          P     = F
          Q     = ONE
        ELSE                                     ! find cf2                               
          XLTURN = .FALSE.
          PK =  ZERO
          WI =  ETA + ETA
          P  =  ZERO
          Q  =  ONE - ETA * XINV
          AR = -E2MM1
          AI =  ETA
          BR =  TWO * (X - ETA)
          BI =  TWO
          DR =  BR / (BR * BR + BI * BI)
          DI = -BI / (BR * BR + BI * BI)
          DP = -XINV * (AR * DI + AI * DR)
          DQ =  XINV * (AR * DR - AI * DI)
          DO 40 L = 1, LIMIT
             P  = P  + DP
             Q  = Q  + DQ
             PK = PK + TWO
             AR = AR + PK
             AI = AI + WI                                                   
             BI = BI + TWO                                                  
             D  = AR * DR - AI * DI + BR                                        
             DI = AI * DR + AR * DI + BI                                        
             C  = ONE / (D * D + DI * DI)                                         
             DR =  C * D                                                      
             DI = -C * DI                                                     
             A  = BR * DR - BI * DI - ONE                                       
             B  = BI * DR + BR * DI                                             
             C  = DP * A  - DQ * B
             DQ = DP * B  + DQ * A                                              
             DP = C
      IF( DABS(DP)+DABS(DQ).LT.(DABS(P)+DABS(Q)) * ACCUR ) GO TO 50
   40     CONTINUE
                                              GO TO 120 ! error exit
   50     NPQ   = PK / TWO                              ! proper exit
          PACCQ = HALF * ACCUR / DMIN1( DABS(Q),ONE )
          IF( DABS(P).GT.DABS(Q) ) PACCQ = PACCQ * DABS(P)
!---------------------------------------------------------------------
!    SOLVE FOR FCMINL = F AT LAMBDA = XLM AND NORMALISING FACTOR OMEGA
!---------------------------------------------------------------------
          GAMMA   = (F - P) / Q
          GAMMAI  = ONE / GAMMA
          IF( DABS(GAMMA) .LE. ONE )  THEN 
                 OMEGA  = DSQRT( ONE  +  GAMMA * GAMMA )
            ELSE
                 OMEGA  = DSQRT( ONE  +  GAMMAI* GAMMAI) * DABS(GAMMA)
            ENDIF 
          OMEGA  = ONE / ( OMEGA * DSQRT(Q) )
          WRONSK = OMEGA
        ENDIF   
!--------------------------------------------------------------------- 
!    RENORMALISE IF SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
!---------------------------------------------------------------------
      IF( KFN.EQ.1 )       THEN         !   spherical Bessel functions
                 ALPHA = XINV
                 BETA  = XINV
        ELSEIF( KFN.EQ.2 ) THEN         ! cylindrical Bessel functions
                 ALPHA = HALF * XINV
                 BETA  = DSQRT( XINV ) * RT2DPI
        ELSE                            ! kfn = 0,   Coulomb functions
                 ALPHA = ZERO     
                 BETA  = ONE
        ENDIF
      FCMINL = DSIGN( OMEGA,DEN ) * BETA
      IF( XLTURN )   THEN
                        GCMINL =   GJWKB * BETA
        ELSE
                        GCMINL =  FCMINL * GAMMA
        ENDIF
      IF( KFN.NE.0 )    GCMINL = -GCMINL         ! Bessel sign differs
      FC (MINL) = FCMINL
      GC (MINL) = GCMINL
      GCP(MINL) = GCMINL * (P - Q * GAMMAI - ALPHA) 
      FCP(MINL) = FCMINL * (F - ALPHA)
      IF( LRANGE.EQ.0 )                          RETURN
!---------------------------------------------------------------------
!    UPWARD RECURRENCE FROM GC(MINL),GCP(MINL) STORED VALUES ARE RL,SL
!    RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
!      XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
!---------------------------------------------------------------------
      OMEGA = BETA * OMEGA / DABS(DEN)
                 XL = XLM
                 RL = ONE 
      DO 60  L = MINL+1 , MAXL                   ! indexed from 0
                 XL = XL + ONE
          IF( ETANE0 ) THEN
                 RL = GC (L)
                 SL = GCP(L)
            ELSE 
                 SL =  XL * XINV
            ENDIF
          GC (L)  = ( (SL - ALPHA) * GC(L-1) - GCP(L-1) ) / RL
          GCP(L)  =    RL *  GC(L-1)  -  (SL + ALPHA) * GC(L)
          FCP(L)  = OMEGA * ( FCP(L)  -  ALPHA * FC(L) )
          FC (L)  = OMEGA *   FC (L)
   60 CONTINUE
      RETURN
!------------------   ERROR MESSAGES
  100 IFAIL = -1
      WRITE(6,1000) X,ACCH
 1000 FORMAT(' FOR X = ',1PD12.3,'     TRY SMALL-X  SOLUTIONS,', &
     ' OR X IS NEGATIVE'/ ,' SQUARE ROOT (ACCURACY) =  ',D12.3/)
                     RETURN
  105 IFAIL = -2                                                        
      WRITE (6,1005) LRANGE,XLMIN,XLM                                    
 1005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES: LRANGE, XLMIN, XLM = ',    &
     I10,1P2D15.6/)                                                        
                     RETURN                                   
  110 IFAIL =  1                                                        
      WRITE (6,1010) LIMIT, CF1,DCF1, PK,ACCUR                              
 1010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',I10,' ITERATIONS',&
     ' CF1,DCF1,PK,ACCUR =  ',1P4D12.3/)                               
                     RETURN                                       
  120 IFAIL =  2                                                        
      WRITE (6,1020) LIMIT,P,Q,DP,DQ,ACCUR
 1020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',I7,' ITERATIONS',& 
     ' P,Q,DP,DQ,ACCUR =  ',1P4D17.7,D12.3/)
                     RETURN                                              
      END SUBROUTINE COUL90                                                              
!---------------------------------------------------------------------                                                                       
      SUBROUTINE  JWKB   (X,ETA,XL, FJWKB,GJWKB, IEXP)            
      DOUBLE PRECISION    X,ETA,XL, FJWKB,GJWKB, DZERO                      
!----------------------------------------------------------------------
!-----COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS  FOR XL .GE. 0.
!-----AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
!-----CALCULATED IN SINGLE, RETURNED IN DOUBLE PRECISION VARIABLES
!-----CALLS DMAX1, SQRT, ALOG, EXP, ATAN2, FLOAT, INT     
!     AUTHOR:    A.R.BARNETT   FEB 1981    LAST UPDATE MARCH 1991
!----------------------------------------------------------------------
      REAL    ZERO,HALF,ONE,SIX,TEN,RL35,ALOGE
      REAL    GH2,XLL1,HLL,HL,SL,RL2,GH,PHI,PHI10
      INTEGER IEXP, MAXEXP
      PARAMETER  ( MAXEXP = 300 )
      DATA  ZERO,HALF,ONE,SIX,TEN  /0.0E0, 0.5E0, 1.0E0, 6.0E0, 1.0E1/
      DATA DZERO,RL35,ALOGE /0.0D0, 35.0E0, 0.4342945E0 /  
!----------------------------------------------------------------------
!CHOOSE MAXEXP NEAR MAX EXPONENT RANGE E.G. 1.D300 FOR DOUBLE PRECISION
!----------------------------------------------------------------------
      GH2   =  X * (ETA + ETA - X)                                         
      XLL1  = DMAX1( XL * XL + XL, DZERO )                                   
      IF( GH2 + XLL1 .LE. ZERO )                 RETURN
      HLL  = XLL1 + SIX / RL35                                           
      HL   = SQRT(HLL)                                                 
      SL   = ETA / HL + HL / X                                             
      RL2  = ONE + ETA * ETA / HLL                                         
      GH   = SQRT(GH2 + HLL) / X                                         
      PHI  = X*GH - HALF*( HL*ALOG((GH + SL)**2 / RL2) - ALOG(GH) )      
      IF ( ETA.NE.ZERO ) PHI = PHI - ETA * ATAN2(X*GH,X - ETA)         
      PHI10 = -PHI * ALOGE                                                
      IEXP  =  INT(PHI10)                                               
      IF ( IEXP.GT.MAXEXP ) THEN
           GJWKB = TEN**(PHI10 - FLOAT(IEXP))               
      ELSE
           GJWKB = EXP(-PHI)                                
           IEXP  = 0                                        
      ENDIF
      FJWKB = HALF / (GH * GJWKB)                                           
      RETURN                                                            
      END SUBROUTINE JWKB                                                           
!---------------------------------------------------------------------
!     END OF CONTINUED-FRACTION COULOMB & BESSEL PROGRAM  COUL90
!---------------------------------------------------------------------
      
      SUBROUTINE CPSI(X,Y,PSR,PSI)
!
!       =============================================
!       Purpose: Compute the psi function for a
!                complex argument
!       Input :  x   --- Real part of z
!                y   --- Imaginary part of z
!       Output:  PSR --- Real part of psi(z)
!                PSI --- Imaginary part of psi(z)
!       =============================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(8)
        DATA A/-.8333333333333D-01,.83333333333333333D-02, &
            -.39682539682539683D-02,.41666666666666667D-02, &
            -.75757575757575758D-02,.21092796092796093D-01, &
            -.83333333333333333D-01,.4432598039215686D0/
        CPSI_PI=3.141592653589793D0
        IF (Y.EQ.0.0D0.AND.X.EQ.INT(X).AND.X.LE.0.0D0) THEN
           PSR=1.0D+300
           PSI=0.0D0
        ELSE
           IF (X.LT.0.0D0) THEN
              X1=X
              Y1=Y
              X=-X
              Y=-Y
           ENDIF
           X0=X
           IF (X.LT.8.0D0) THEN
              N=8-INT(X)
              X0=X+N
           ENDIF
           IF (X0.EQ.0.0D0.AND.Y.NE.0.0D0) TH=0.5D0*CPSI_PI
           IF (X0.NE.0.0D0) TH=DATAN(Y/X0)
           Z2=X0*X0+Y*Y
           Z0=DSQRT(Z2)
           PSR=DLOG(Z0)-0.5D0*X0/Z2
           PSI=TH+0.5D0*Y/Z2
           DO 10 K=1,8
              PSR=PSR+A(K)*Z2**(-K)*DCOS(2.0D0*K*TH)
10            PSI=PSI-A(K)*Z2**(-K)*DSIN(2.0D0*K*TH)
           IF (X.LT.8.0D0) THEN
              RR=0.0D0
              RI=0.0D0
              DO 20 K=1,N
                 RR=RR+(X0-K)/((X0-K)**2.0D0+Y*Y)
20               RI=RI+Y/((X0-K)**2.0D0+Y*Y)
              PSR=PSR-RR
              PSI=PSI+RI
           ENDIF
           IF (X1.LT.0.0D0) THEN
              TN=DTAN(CPSI_PI*X)
              TM=DTANH(CPSI_PI*Y)
              CT2=TN*TN+TM*TM
              PSR=PSR+X/(X*X+Y*Y)+CPSI_PI*(TN-TN*TM*TM)/CT2
              PSI=PSI-Y/(X*X+Y*Y)-CPSI_PI*TM*(1.0D0+TN*TN)/CT2
              X=X1
              Y=Y1
           ENDIF
        ENDIF
        RETURN
      END SUBROUTINE CPSI
 
END MODULE special_functions
!
!   Gauss-Legendre mesh types and initialization routines
!   Now also extended to include setup_gauss_laguerre_mesh
!
MODULE gauss_legendre_mesh
  
  USE aux
  USE constants
  
  IMPLICIT NONE
  
  ! basic mesh info
  TYPE, PUBLIC :: mesh_info
     INTEGER  :: amount
     REAL(DP) :: xmin, xmax
  END TYPE mesh_info
  
  ! GAUSS-LEGENDRE MESH POINTS
  ! Gauss-Legendre mesh point x, corresponding integration weight w and corresponding x*x*w-value
  TYPE, PUBLIC :: gauleg_mesh_point
     SEQUENCE
     REAL(DP) :: x, w, xxw
  END TYPE gauleg_mesh_point
  
  ! mesh points and weights in momentum space
  TYPE, PUBLIC :: gauleg_mesh
     TYPE(mesh_info)                                    :: info
     TYPE(gauleg_mesh_point), DIMENSION(:), ALLOCATABLE :: pnt
  END TYPE gauleg_mesh

  ! PLAIN MESH POINTS
  ! equidistant mesh points, constructed according to mesh_info values.
  ! auxilary to interpolation procedures
  TYPE, PUBLIC :: evaluation_mesh
     TYPE(mesh_info)                     :: info
     REAL(DP), DIMENSION(:), ALLOCATABLE :: x
  END TYPE evaluation_mesh
  
CONTAINS

  ! constructs mesh of mesh%info%amount equidistant points, from xmin to xmax 
  SUBROUTINE setup_evaluation_mesh (mesh)
    TYPE(evaluation_mesh), INTENT(INOUT) :: mesh
    REAL(DP)                             :: dx
    INTEGER                              :: i
    
    CALL destroy_evaluation_mesh (mesh)
    ALLOCATE(mesh%x(1:mesh%info%amount))
    
    IF (mesh%info%amount == 1) THEN 
       mesh%x(1) = (mesh%info%xmax + mesh%info%xmin)/2.0_dp
       RETURN
    END IF
    
    dx = (mesh%info%xmax - mesh%info%xmin) / (mesh%info%amount - 1)
    DO i = 1, mesh%info%amount 
       mesh%x(i) = mesh%info%xmin + (i-1) * dx
    ENDDO
    
  END SUBROUTINE setup_evaluation_mesh

  ! constructs mesh of mesh%info%amount equidistant points, from xmin to xmax 
  SUBROUTINE allocate_evaluation_mesh (mesh)
    TYPE(evaluation_mesh), INTENT(INOUT) :: mesh

    CALL destroy_evaluation_mesh (mesh)
    ALLOCATE(mesh%x(1:mesh%info%amount))
    
  END SUBROUTINE allocate_evaluation_mesh

  ! constructs mesh of mesh%info%amount equidistant points, from xmin to xmax 
  ! In this modified version, there are more points in the beginning...
  SUBROUTINE setup_evaluation_mesh_mod (mesh)
    TYPE(evaluation_mesh), INTENT(INOUT) :: mesh
    REAL(DP)                             :: dxF, dx, F, xF 
    INTEGER                              :: i, NF, N
    
    CALL destroy_evaluation_mesh (mesh)
    ALLOCATE(mesh%x(1:mesh%info%amount))
    
    F  = 0.20_dp
    xF = 4.0_dp
    
    IF (xF>=1.0_dp/F) xF = xF-1
    
    NF = INT(F*mesh%info%amount-1)
    NF = xF*NF
    N = mesh%info%amount-1-NF

    dx  = (1.0_dp-F)*(mesh%info%xmax - mesh%info%xmin) / N
    dxF = (F)*(mesh%info%xmax - mesh%info%xmin) / NF
    
    DO i = 1, mesh%info%amount 
       IF (i<=NF) THEN
          mesh%x(i) = mesh%info%xmin + (i-1) * dxF
       ELSEIF (i>NF .AND. i<mesh%info%amount) THEN
          mesh%x(i) = mesh%x(NF) + (i-NF) * dx
       ELSEIF(i==mesh%info%amount) THEN
          mesh%x(i) = mesh%info%xmax
       END IF
    ENDDO
    
  END SUBROUTINE setup_evaluation_mesh_mod
   
  SUBROUTINE destroy_evaluation_mesh (mesh)
    TYPE(evaluation_mesh), INTENT(INOUT) :: mesh
    IF(ALLOCATED(mesh%x)) DEALLOCATE(mesh%x)
  END SUBROUTINE destroy_evaluation_mesh
  
  !
  !
  !      This routine calculates gauss-legendre mesh points and weights      
  !      INPUT:                                                              
  !      mesh%info%xmin     : lower limit of the integration interval                      
  !      mesh%info%xmax     : upper limit ---------- "" -------------                      
  !      mesh%info%amount   : the desired number of mesh points                            
  !      OUTPUT:                                                            
  !      mesh%pnt(:)%x      : gauss-legendre mesh points on the interval (x1,x2)          
  !      mesh%pnt(:)%w      : the corresponding weights                                   
  !      FROM               : Numerical recipes
  !      F90 version        : M. Hjorth-Jensen
  !      Object interface   : M. Kartamyshev
  !
  SUBROUTINE setup_gauleg_mesh (mesh)
    USE constants
    CHARACTER(*), PARAMETER                 :: my_name = 'setup_gauleg_mesh'
    TYPE(gauleg_mesh), INTENT(INOUT)        :: mesh
    INTEGER                                 :: i, j, m, n
    REAL(DP)                                :: x1, x2
    REAL(DP), DIMENSION(1:mesh%info%amount) :: x, w
    REAL(DP)                                :: p1,p2,p3,pp,xl,xm,z,z1
    REAL(DP), PARAMETER                     :: EPS = 3.D-14
    
    ! allocate points and weights storages
    CALL destroy_gauleg_mesh (mesh)
    
    IF (mesh%info%amount <= 0 .OR. mesh%info%xmin >= mesh%info%xmax) THEN
       WRITE(*,*) my_name//': incorrect mesh info', mesh%info ; STOP
    ENDIF
    
    ALLOCATE( mesh%pnt( 1:mesh%info%amount ) )
    mesh%pnt(:) = gauleg_mesh_point(0.0_dp, 0.0_dp, 0.0_dp)
    
    ! set values of local variables
    x1 = mesh%info%xmin ; x2 = mesh%info%xmax; n = mesh%info%amount
    
    m=(n+1)/2
    xm=0.5_dp*(x2+x1)
    xl=0.5_dp*(x2-x1)
    DO i=1,m
       z1=0.0_dp
       z=COS(pi*(i - 0.25_dp)/(n + 0.5_dp))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.0_dp
          p2=0.0_dp
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
          END DO
          pp=n*(z*p1-p2)/(z*z-1.0_dp)
          z1=z
          z=z-p1/pp
       END DO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO
    
    ! set return values
    mesh%pnt(:)%x = x(:) ; mesh%pnt(:)%w = w(:) ; mesh%pnt(:)%xxw = x(:) * x(:) * w(:)

  END SUBROUTINE setup_gauleg_mesh
  !
  !
  !      This routine calculates gauss-LAGUERRE mesh points and weights      
  !      NOTE : the weights are multiplied with exp(x)
  !      INPUT:                                                              
  !      mesh%info%xmin     : 0
  !      mesh%info%xmax     : infinity
  !      mesh%info%amount   : the desired number of mesh points                            
  !      alpha              : exponent alpha in laguerre-integrand § x^[alpha] exp[-x] f(x) dx
  !      OUTPUT:                                                            
  !      mesh%pnt(:)%x      : gauss-laguerre mesh points on the interval (0,infinity)          
  !      mesh%pnt(:)%w      : the corresponding weights                                   
  !      FROM               : Numerical recipes
  !      F90 version        : M. Hjorth-Jensen
  !      Object interface   : M. Kartamyshev
  !      laguerre           : A. Ekstrom
  !
  SUBROUTINE setup_gauss_laguerre_mesh (mesh, alpha)
    
    USE constants
    USE special_functions, ONLY: gammln
    
    IMPLICIT NONE
    
    CHARACTER(*), PARAMETER                 :: my_name = 'setup_gauss_laguerre_mesh'
    TYPE(gauleg_mesh), INTENT(INOUT)        :: mesh
    REAL(DP), INTENT(IN)                    :: alpha
    INTEGER                                 :: i, its, j, n, maxit
    REAL(DP), DIMENSION(1:mesh%info%amount) :: x, w
    REAL(DP)                                :: p1,p2,p3,pp,z,z1,ai
    REAL(DP), PARAMETER                     :: EPS = 3.D-14
    
    ! allocate points and weights storages
    CALL destroy_gauleg_mesh (mesh)
    
    IF (mesh%info%amount <= 0 .OR. mesh%info%xmin /= 0.0_dp .OR. mesh%info%xmax /= 0.0_dp) THEN
       WRITE(*,*) my_name//': incorrect mesh info', mesh%info ; STOP
    ENDIF
    
    ALLOCATE( mesh%pnt( 1:mesh%info%amount ) )
    mesh%pnt(:) = gauleg_mesh_point(0.0_dp, 0.0_dp, 0.0_dp)

    maxit = 20
    z1 = 0.0_dp
    n = mesh%info%amount

    DO i=1,n
       IF(i==1) THEN
          z=(1.+alpha)*(3.+.92*alpha)/(1.+2.4*n+1.8*alpha)
       ELSE IF(i==2) THEN
          z=z+(15.+6.25*alpha)/(1.+.9*alpha+2.5*n)
       ELSE
          ai=i-2
          z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alpha/(1.+3.5*ai))* & 
               (z-x(i-2))/(1.+.3*alpha)
       END IF
       its = 0
       DO WHILE ( ABS(z-z1) >  EPS)
          p1=1.d0
          p2=0.d0
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2*j-1+alpha-z)*p2-(j-1+alpha)*p3)/j
          END DO
          pp=(n*p1-(n+alpha)*p2)/z
          z1=z
          z=z1-p1/pp
          its = its + 1
          IF( its > MAXIT) THEN
             WRITE(ERR,*) 'warning(setup_gauss_laguerre_mesh): to many iteration, exiting cycle'
             RETURN
          END IF
       END DO
       x(i)=z
       w(i)=-exp(gammln(alpha+n)-gammln(REAL(n,kind=dp)))/(pp*n*p2) * exp(x(i))
    END DO
    
    ! set return values
    mesh%pnt(:)%x = x(:) ; mesh%pnt(:)%w = w(:) ; mesh%pnt(:)%xxw = x(:) * x(:) * w(:)

  END SUBROUTINE setup_gauss_laguerre_mesh
  
  SUBROUTINE destroy_gauleg_mesh (mesh)
    TYPE(gauleg_mesh), INTENT(INOUT) :: mesh
    IF (ALLOCATED(mesh%pnt) ) DEALLOCATE (mesh%pnt)
  END SUBROUTINE destroy_gauleg_mesh
     
END MODULE gauss_legendre_mesh

MODULE matrix_storage
  
  USE constants

  ! Type to hold the twobody interaction potential
  ! matrix for a specific channel.
  TYPE, PUBLIC :: PS_vmatrix_type
     REAL(DP), DIMENSION(:),ALLOCATABLE :: mtx
  END TYPE PS_vmatrix_type

  ! Type to hold the twobody interaction potential
  ! matrix for a specific channel.
  TYPE, PUBLIC :: vmatrix_type
     REAL(DP), DIMENSION(:,:),ALLOCATABLE :: mtx
  END TYPE vmatrix_type
  
  TYPE, PUBLIC :: interaction_type
     ! each channels holds a <nl|v|n'l'> 
     !matrix-block of interaction elements
     !
     TYPE(vmatrix_type),DIMENSION(:), ALLOCATABLE :: channel ! J, Tz, Pi
     INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: sjtzpi2channel! for given S,J, Tz, Pi : obtain channel number
  END TYPE interaction_type

  TYPE, PUBLIC :: PS_interaction_type
     ! each channels holds a <nl|v|n'l'> 
     !matrix-block of interaction elements
     !
     TYPE(PS_vmatrix_type),DIMENSION(:), ALLOCATABLE :: channel ! J, Tz, Pi
     INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: sjtzpi2channel! for given S,J, Tz, Pi : obtain channel number
  END TYPE PS_interaction_type

  TYPE, PUBLIC :: uncoupled_interaction_type
     ! each channels holds a <nl|v|n'l'> 
     !matrix-block of interaction elements
     !
     TYPE(vmatrix_type),DIMENSION(:), ALLOCATABLE :: channel ! M, Tz, Pi
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: tzpi2channel! for given Tz, Pi : obtain channel number
  END TYPE uncoupled_interaction_type
  
  TYPE, PUBLIC :: coulomb_type
     REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: mtx
  END TYPE coulomb_type

END MODULE matrix_storage

MODULE matrix_stuff

  USE constants

CONTAINS
  
  FUNCTION trace(a) RESULT(res)
    REAL(DP), INTENT(IN) :: a(:,:)
    REAL(DP) :: res
    INTEGER :: nrows, ncols,i
    
    ncols = SIZE(a(:,1))
    nrows = SIZE(a(1,:))
    
    res = 0.0_dp
    
    IF (ncols /= nrows) STOP 'error(trace): matrix not square'
    
    DO i = 1, nrows
       res = res + a(i,i)
    END DO
    
  END FUNCTION trace
  
  FUNCTION PStrace(a) RESULT(res)
    REAL(DP), INTENT(IN) :: a(:)
    REAL(DP) :: res
    INTEGER :: nrows, ncols,i, dim
    
    dim = SIZE(a(:))
    nrows = -0.5_dp + SQRT(0.25_dp + 2.0_dp*dim)
    ncols = nrows
    
    res = 0.0_dp
    
    DO i = 1, ncols
       res = res + a(i*(i-1)/2 + i)
    END DO
    
  END FUNCTION PStrace
  
  SUBROUTINE truncate_matrix (a, new_dim)
    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: a(:,:)
    INTEGER , INTENT(IN) :: new_dim
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: tmp_a
    INTEGER              :: i, k
    
    ALLOCATE(tmp_a(1:new_dim,1:new_dim))
    
    tmp_a(:,:) = 0.0_dp
    
    DO i = 1, new_dim
       DO k = 1, new_dim
          
          tmp_a(i,k) = a(i,k)
          
       ENDDO ! k
    ENDDO ! i
    
    DEALLOCATE(a)
    ALLOCATE(a(1:new_dim,1:new_dim))
    a(:,:) = tmp_a(:,:)
    DEALLOCATE(tmp_a)
    
  END SUBROUTINE truncate_matrix

  SUBROUTINE ps_truncate_matrix (a, new_dim)
    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: a(:)
    INTEGER , INTENT(IN) :: new_dim
    REAL(DP), DIMENSION(:), ALLOCATABLE :: tmp_a
    INTEGER              :: i, k
    
    ALLOCATE(tmp_a(1:new_dim*(new_dim+1)/2))
    
    tmp_a(:) = 0.0_dp
    
    DO i = 1, new_dim
       DO k = 1, new_dim
          
          tmp_a(i*(i-1)/2 + k) = a(i*(i-1)/2 + k)
          
       ENDDO ! k
    ENDDO ! i
    
    DEALLOCATE(a)
    ALLOCATE(a(1:new_dim*(new_dim+1)/2))
    a(:) = tmp_a(:)
    DEALLOCATE(tmp_a)
    
  END SUBROUTINE ps_truncate_matrix

  SUBROUTINE write_matrix (unt,a)
    REAL(DP), INTENT(IN) :: a(:,:)
    INTEGER              :: i, k
    INTEGER, INTENT(IN)  :: unt
    
    DO i = 1, SIZE(a(:,1))
       DO k = 1, SIZE(a(1,:))
          
          WRITE(unt,fmt='(F20.16)',ADVANCE='NO') a(i,k)
          
       ENDDO ! k
       WRITE(unt,*)
    ENDDO ! i
    
  END SUBROUTINE write_matrix

  SUBROUTINE write_matrix_packed_storage (unt,n,a)
    
    USE aux

    REAL(DP), INTENT(IN) :: a(n*(n+1)/2)
    INTEGER, INTENT(IN)  :: n
    INTEGER              :: i, k, idx, row, col
    INTEGER, INTENT(IN)  :: unt
    
    DO row = 1, n
       DO col = 1, n
       
          i = row
          k = col
          
          IF (i<=k)  idx = (k*(k-1)/2 + i)
          IF (i>k )  idx = (i*(i-1)/2 + k)
          
          WRITE(unt,fmt='(F10.4)',ADVANCE='NO') a(idx)
          
       ENDDO ! k
       WRITE(unt,*)
    ENDDO ! i
    
  END SUBROUTINE write_matrix_packed_storage
  
  SUBROUTINE write_matrix_packed (unt,a)
    REAL(DP), INTENT(IN) :: a(:,:)
    INTEGER              :: i, k
    INTEGER, INTENT(IN)  :: unt
    
    DO i = 1, SIZE(a(:,1))
       DO k = 1, SIZE(a(1,:))
          
          WRITE(unt,fmt='(F12.4)',ADVANCE='NO') a(i,k)
          
       ENDDO ! k
       WRITE(unt,*)
    ENDDO ! i
    
  END SUBROUTINE write_matrix_packed

  SUBROUTINE compute_eigenvalues_to_file(unt, dim, H)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unt
    INTEGER, INTENT(INOUT) :: dim
    REAL*8, INTENT(INOUT) :: H(dim,dim)
    REAL*8 ::  E(dim)
    INTEGER :: i
    LOGICAL :: UNT_IS_OPEN
    
    INQUIRE(unt, opened=UNT_IS_OPEN)
    IF (.NOT. UNT_IS_OPEN) RETURN

    CALL dsyev_interface(dim, H, E)
    
    WRITE(unt,"(A)") 'EIGENVALUES:'
    DO i=1, dim
       WRITE(unt,"(I5,F20.9)") i, E(i)
    END DO
    
  END SUBROUTINE compute_eigenvalues_to_file
  
  SUBROUTINE ddot_interface(N, DX, DY, SCALAR)

    IMPLICIT NONE
    INTEGER , INTENT(INOUT) :: N
    DOUBLE PRECISION  , INTENT(INOUT) :: DX(N), DY(N)
    DOUBLE PRECISION  , INTENT(INOUT) :: SCALAR
    INTEGER :: INCX, INCY
    DOUBLE PRECISION  :: DDOT

    INCX = 1
    INCY = 1

    SCALAR = DDOT(N, DX, INCX, DY, INCY)

  END SUBROUTINE ddot_interface

  SUBROUTINE dspmv_interface(N, alpha, A, X, beta, Y)

    IMPLICIT NONE
    
    INTEGER , INTENT(INOUT) :: N
    DOUBLE PRECISION  , INTENT(INOUT) :: alpha
    DOUBLE PRECISION  , INTENT(INOUT) :: beta
    DOUBLE PRECISION  , INTENT(INOUT) :: A(N)
    DOUBLE PRECISION  , INTENT(INOUT) :: X(N)
    DOUBLE PRECISION  , INTENT(INOUT) :: Y(N)
    INTEGER :: INCX, INCY

    INCX  = 1
    INCY  = 1

    CALL DSPMV('U', N, alpha, A, X, INCX, beta, Y, INCY)

  END SUBROUTINE dspmv_interface


  ! get order of square matrix stored
  ! in packed format
  FUNCTION get_sqpsmat_dim(mat) RESULT(res)
    
    USE matrix_storage
    
    IMPLICIT NONE

    TYPE(PS_vmatrix_type), INTENT(IN) :: mat
    REAL(DP) :: dim
    INTEGER :: res
    
    dim = -0.5 +SQRT(0.25_dp+2.0_dp*SIZE(mat%mtx))
    
    res = INT(dim)
    
  END FUNCTION get_sqpsmat_dim

  !         co      co+ca
  !    _________________
  !   |     co    co+ca |
  !   | ro   ________   |
  !   |     |        |  |
  !   |     |        |  |
  !   |     |________|  |
  !   | ro+ra           |
  !   |                 |
  !   |_________________|
  !
  !   Upper triangle of packed storage matrix (psmat)
  !   is filled with values from full format matrix mat.
  !
  !   values below diagonal are ignored
  !
  SUBROUTINE block_fill_ps_matrix(row_offset, row_amount, col_offset, col_amount, psmat, mat)
    
    USE matrix_storage

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: row_offset, row_amount
    INTEGER, INTENT(IN) :: col_offset, col_amount
    TYPE(PS_vmatrix_type), INTENT(INOUT) :: psmat
    TYPE(vmatrix_type), INTENT(IN) :: mat
    INTEGER :: r,c,cidx,ridx

    DO c=1, col_amount
       cidx = c+col_offset-1
       DO r=1, row_amount
          ridx = r+row_offset-1
          IF (ridx>cidx) CYCLE
          psmat%mtx(cidx*(cidx-1)/2 + ridx) = mat%mtx(r,c)
       END DO
    END DO
    
  END SUBROUTINE block_fill_ps_matrix

  ! full format 2 packed format
  SUBROUTINE dtrttp_interface(N, MTR, MTP)
    
    IMPLICIT NONE
    INTEGER , INTENT(INOUT) :: N
    REAL*8  , INTENT(INOUT) :: MTP(N*(N+1)/2)
    REAL*8  , INTENT(INOUT) :: MTR(1:N,1:N)
    INTEGER :: INFO
    
    CALL DTRTTP('U',N,MTR,N,MTP,INFO)
    
  END SUBROUTINE dtrttp_interface
  
  ! packed format 2 full format
  SUBROUTINE dtpttr_interface(N, MTP, MTR)
    
    IMPLICIT NONE
    INTEGER , INTENT(INOUT) :: N
    REAL*8  , INTENT(INOUT) :: MTP(N*(N+1)/2)
    REAL*8  , INTENT(INOUT) :: MTR(1:N,1:N)
    INTEGER :: INFO, i
    
    CALL DTPTTR('U',N,MTP,MTR,N,INFO)
    
    ! uncomment to get the full matrix
    ! not just the Upper/Lower sector.
    MTR = MTR + TRANSPOSE(MTR)
    
    DO i=1,N
       MTR(i,i) = MTR(i,i)-MTP(i*(i-1)/2+i)
    END DO
    
  END SUBROUTINE dtpttr_interface 

  SUBROUTINE dspev_V_interface(N, AP, Z, W)
    
    IMPLICIT NONE
    INTEGER , INTENT(INOUT) :: N
    REAL*8  , INTENT(INOUT) :: AP(N*(N+1)/2)
    REAL*8  , INTENT(INOUT) :: W(N)   ! eigenvalues
    REAL*8  , INTENT(INOUT) :: Z(N,N) ! eigenvectors
    INTEGER :: LDZ, INFO
    REAL*8, DIMENSION(:)  , ALLOCATABLE ::  WORK

    LDZ   = N
    ALLOCATE(WORK(3*N))
    CALL DSPEV('V','Upper',N,AP,W,Z,LDZ,WORK,INFO)
    DEALLOCATE(WORK)
    
  END SUBROUTINE dspev_V_interface

  SUBROUTINE dspev_N_interface(N, AP, W)
    

    IMPLICIT NONE
    INTEGER , INTENT(INOUT) :: N
    REAL*8  , INTENT(INOUT) :: AP(N*(N+1)/2)
    REAL*8  , INTENT(INOUT) :: W(N)   ! eigenvalues
    INTEGER :: LDZ, INFO, Z
    REAL*8, DIMENSION(:)  , ALLOCATABLE ::  WORK
    
    LDZ   = N
    ALLOCATE(WORK(3*N))
    CALL DSPEV('N','Upper',N,AP,W,Z,LDZ,WORK,INFO)
    
  END SUBROUTINE dspev_N_interface

  SUBROUTINE dsyev_interface(N, A, W)
    
    IMPLICIT NONE
    INTEGER , INTENT(INOUT) :: N
    REAL*8  , INTENT(INOUT) :: A(N,N)
    REAL*8  , INTENT(INOUT) :: W(N)
    
    INTEGER :: LDA, LWORK
    INTEGER :: INFO
    REAL*8, DIMENSION(:)  , ALLOCATABLE ::  WORK
    
    LDA   = N
    LWORK = -1
    ALLOCATE(WORK(1))
    CALL DSYEV('Vectors','Upper',N,A,LDA,W,WORK,LWORK,INFO)
    LWORK = INT(WORK(1))
    DEALLOCATE(WORK)
    ALLOCATE(WORK(LWORK))
    CALL DSYEV('Vectors','Upper',N,A,LDA,W,WORK,LWORK,INFO)
  END SUBROUTINE dsyev_interface

  SUBROUTINE dgetrx_interface(N, A)
    
    IMPLICIT NONE
    INTEGER , INTENT(INOUT) :: N
    REAL*8  , INTENT(INOUT) :: A(N,N)
    INTEGER :: LDA, LWORK
    INTEGER :: INFO
    REAL*8, DIMENSION(:)  , ALLOCATABLE ::  WORK
    REAL*8, DIMENSION(:)  , ALLOCATABLE ::  IPIV
    
    LDA   = N
    
    ALLOCATE(IPIV(N))
    CALL DGETRF(N,N,A,LDA,IPIV,INFO)
  
    LWORK = -1
    ALLOCATE(WORK(1))
    
    CALL DGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
    LWORK = INT(WORK(1))
    DEALLOCATE(WORK)
    ALLOCATE(WORK(LWORK))
    
    CALL DGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
    
  END SUBROUTINE dgetrx_interface
  
  FUNCTION is_unit_matrix(M,tol) RESULT(res)
    
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: M(:,:)
    REAL(DP), INTENT(IN) :: tol ! tolerance
    LOGICAL :: res
    INTEGER :: row, col

    res = .FALSE.
    
    IF (SIZE(M,1) /= SIZE(M,2)) THEN
       WRITE(ERR,"(A)") 'error(is_unit_matrix): non-square matrix'
       STOP
    END IF
    
    DO row=1, SIZE(M,1)
       DO col=1, SIZE(M,2)
          IF (row/=col .AND. ABS(M(row,col)) > tol) RETURN
       END DO
    END DO
    
    res = .TRUE.

  END FUNCTION is_unit_matrix

  FUNCTION is_symmetric_matrix(M,tol) RESULT(res)
    
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: M(:,:)
    REAL(DP), INTENT(IN) :: tol ! tolerance
    LOGICAL :: res
    INTEGER :: row, col

    res = .FALSE.
    
    IF (SIZE(M,1) /= SIZE(M,2)) THEN
       WRITE(ERR,"(A)") 'error(is_unit_matrix): non-square matrix'
       STOP
    END IF
    
    DO row=1, SIZE(M,1)
       DO col=1, SIZE(M,2)
          IF (ABS(M(row,col)-M(col,row)) > tol) RETURN
       END DO
    END DO
    
    res = .TRUE.

  END FUNCTION is_symmetric_matrix

  FUNCTION is_zero_matrix(M,tol) RESULT(res)
    
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: M(:,:)
    REAL(DP), INTENT(IN) :: tol ! tolerance
    LOGICAL :: res
    INTEGER :: row, col

    res = .FALSE.
    
    DO row=1, SIZE(M,1)
       DO col=1, SIZE(M,2)
          IF (ABS(M(row,col)) > tol) RETURN
       END DO
    END DO
    
    res = .TRUE.

  END FUNCTION is_zero_matrix

  FUNCTION compare_matrices(A, B, tol) RESULT(res)
    IMPLICIT NONE
    REAL(DP), INTENT(IN)    :: A(:,:), B(:,:)
    REAL(DP), INTENT(IN)    :: tol ! tolerance
    INTEGER :: res
    INTEGER :: row, col
    
    res = 0

    IF (SIZE(A,1) /= SIZE(B,1)) THEN
       res = -1
       RETURN
    END IF

    IF (SIZE(A,2) /= SIZE(B,2)) THEN
       res = -2
       RETURN
    END IF

    DO row=1, SIZE(A,1)
       DO col=1, SIZE(A,2)
          IF (ABS(A(row,col)-B(row,col)) > tol) res = res+1
       END DO
    END DO
    
  END FUNCTION compare_matrices

  SUBROUTINE set_matrix_to_unit(M)

    USE constants
    
    IMPLICIT NONE

    REAL(DP), INTENT(INOUT) :: M(:,:)
    INTEGER :: i

    IF (SIZE(M,1) /= SIZE(M,2)) THEN
       WRITE(ERR,"(A)") 'error(set_matrix_to_unit): non-square matrix'
       STOP
    END IF
  
    M(:,:) = 0.0_dp

    DO i=1, SIZE(M(1,:))
       
       M(i,i) = 1.0_dp
       
    END DO
    
  END SUBROUTINE set_matrix_to_unit

  SUBROUTINE write_vector (unt,a)
    REAL*8, INTENT(IN) :: a(:)
    INTEGER              :: i
    INTEGER, INTENT(IN)  :: unt
    
    DO i = 1, SIZE(a(:))
       WRITE(unt,fmt='(F20.16)',ADVANCE='NO') a(i)
       WRITE(unt,*)
    ENDDO ! i
    
  END SUBROUTINE write_vector

  SUBROUTINE write_cmatrix (unt,a)
    COMPLEX*16, INTENT(IN) :: a(:,:)
    INTEGER              :: i, k
    INTEGER, INTENT(IN)  :: unt
    
    DO i = 1, SIZE(a(:,1))
       DO k = 1, SIZE(a(1,:))
          
          WRITE(unt,fmt='(A1,F10.3,F10.3,A1)',ADVANCE='NO') '(',a(i,k),')'
          
       ENDDO ! k
       WRITE(unt,*)
    ENDDO ! i
    
  END SUBROUTINE write_cmatrix
  
  SUBROUTINE write_cvector (unt,a)
    COMPLEX*16, INTENT(IN) :: a(:)
    INTEGER              :: i
    INTEGER, INTENT(IN)  :: unt
    
    DO i = 1, SIZE(a(:))
       WRITE(unt,fmt='(F10.3,F10.3)',ADVANCE='NO') a(i)
       WRITE(unt,*)
    ENDDO ! i
    
  END SUBROUTINE write_cvector
  
  SUBROUTINE random_matrix(nr, nc, A, rmin, rmax, type)
    
    IMPLICIT NONE

    INTEGER , INTENT(IN)    :: nr, nc
    REAL*8  , INTENT(IN)    :: rmin, rmax
    CHARACTER(LEN=4) , INTENT(IN)  :: type
    REAL*8  , INTENT(INOUT) :: A(nr,nc)
    INTEGER :: rr, cc
    ! - random number
    INTEGER :: i_seed
    INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
    INTEGER, DIMENSION(1:8) :: dt_seed
    REAL :: r
    ! - done
    
    IF (nr<1 .OR. nc<1) THEN
       WRITE(ERR,*) 'error(random_matrix): invalid size', nr, nc
       STOP
    END IF

    CALL random_seed(size=i_seed)
    ALLOCATE(a_seed(1:i_seed))
    CALL random_seed(get=a_seed)
    CALL date_and_time(values=dt_seed)
    a_seed(i_seed)=dt_seed(8)
    a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
             
    CALL random_seed(put=a_seed)
    DEALLOCATE(a_seed)
    
    IF (type == 'nsym') THEN
       
       DO rr=1, nr
          DO cc=1, nc
             CALL random_number(r)
             A(rr,cc) = rmin +r*(rmax-rmin)
          END DO
       END DO

    END IF

    IF (type == 'symm' .AND. (nr==nc)) THEN
       
       DO rr=1, nr
          DO cc=rr, nc
             CALL random_number(r)
             A(rr,cc) = rmin +r*(rmax-rmin)
             A(cc,rr) = A(rr,cc)
          END DO
       END DO
       
    END IF
       
  END SUBROUTINE random_matrix
  !
  !            Routines to do mtx inversion, from Numerical
  !            Recepies, Teukolsky et al. Routines included
  !            below are MATINV, LUDCMP and LUBKSB. See chap 2
  !            of Numerical Recipes for further details
  !            Recoded in FORTRAN 90 by M. Hjorth-Jensen
  !
  SUBROUTINE matinv(a,n)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j
    REAL*8, DIMENSION(n,n), INTENT(INOUT)  :: a
    REAL*8, ALLOCATABLE :: y(:,:)
    REAL*8 :: d
    INTEGER, ALLOCATABLE :: indx(:)
    
    ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
    y=0.
    !     setup identity matrix
    DO i=1,n
       y(i,i)=1.
    ENDDO
    !     LU decompose the matrix just once
    CALL  lu_decompose(a,n,indx,d)
    
    !     Find inverse by columns
    DO j=1,n
       CALL lu_linear_equation(a,n,indx,y(:,j))
    ENDDO
    !     The original matrix a was destroyed, now we equate it with the inverse y 
    a=y
    
    DEALLOCATE ( y ); DEALLOCATE ( indx )
    
  END SUBROUTINE matinv
  
  SUBROUTINE dgeev_interface(N, A, WR, WI, VL, VR)
    
    IMPLICIT NONE
    INTEGER , INTENT(INOUT) :: N
    REAL*8  , INTENT(INOUT) :: A(N,N)
    REAL*8  , INTENT(INOUT) :: VL(N,N), VR(N,N)
    REAL*8  , INTENT(INOUT) :: WR(N),WI(N)
    
    INTEGER :: LDA, LWORK, LDVL, LDVR
    INTEGER :: INFO
    REAL*8, DIMENSION(:)  , ALLOCATABLE ::  WORK
    
    LDA   = N
    LDVL  = N
    LDVR  = N
    LWORK = -1
    ALLOCATE(WORK(1))
    
    CALL DGEEV('V', 'V', N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
    LWORK = INT(WORK(1))
    DEALLOCATE(WORK)
    ALLOCATE(WORK(LWORK))
    CALL DGEEV('V', 'V', N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
    
  END SUBROUTINE dgeev_interface
  !
  !     Given an NxN matrix A(N,N), this routine replaces it by the LU 
  !     decomposed one, where the matrix elements are stored in the same 
  !     matrix A. The array indx is  an output vector which records the row
  !     permutation effected by the partial pivoting. d is the determinant
  !
  SUBROUTINE lu_decompose(a,n,indx,d)

  IMPLICIT NONE
  INTEGER :: n, i, j, k, imax
  REAL*8 :: sum , tiny, aamax, dum, d
  REAL*8, DIMENSION(n,n) :: a
  INTEGER, DIMENSION(n) :: indx
  REAL*8, ALLOCATABLE :: vv(:)

  tiny=1.0e-20
  ALLOCATE ( vv(n) )
  D=1.
  DO i=1,n
     aamax=0.
     DO j=1,n
        IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
     ENDDO
     !     Zero is the largest element
     IF (aamax == 0.) STOP 'Singular matrix lu_decompose.'
     !     No nonzero largest element
     vv(i)=1./aamax
  ENDDO
  !     loop over columns
  DO j=1,n
     !     solves equation 2.3.12 except for i=j of Numerical Recipes
     IF (j > 1) THEN
        DO i=1,j-1
           sum=a(i,j)
           IF (i > 1)THEN
              DO k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
              ENDDO
              a(i,j)=sum
           ENDIF
        ENDDO
     ENDIF
     !    start searching for largest pivot element
     aamax=0.
     DO i=j,n
        sum=a(i,j)
        IF (j > 1)THEN
           DO k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           ENDDO
           a(i,j)=sum
        ENDIF
        dum=vv(i)*ABS(sum)
        IF (dum >= aamax) THEN
           imax=i
           aamax=dum
        ENDIF
     ENDDO
     !    interchange of rows
     IF (j /= imax)THEN
        DO k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        ENDDO
        !    change of parity for determinant
        d=-d
        vv(imax)=vv(j)
     ENDIF
     indx(j)=imax
     IF(j /= n) THEN
        IF(a(j,j) == 0.) a(j,j)=tiny
        dum=1./a(j,j)
        DO i=j+1,n
           a(i,j)=a(i,j)*dum
        ENDDO
     ENDIF
     !    set up determinant
     d=d*a(j,j)
  ENDDO
  IF(a(n,n) == 0.)  a(n,n)=tiny
  DEALLOCATE ( vv)

END SUBROUTINE lu_decompose

!     Solves set of linear equations Ax=b, A is input as an LU decompomsed
!     matrix and indx keeps track of the permutations of the rows. b is input
!     as the right-hand side vector b and returns the solution x. A, n and indx
!     are not modified by this routine. This function takes into that b can contain
!     many zeros and is therefore suitable for matrix inversion


SUBROUTINE lu_linear_equation(a,n,indx,b)

  IMPLICIT NONE
  INTEGER :: n, ii, ll, i, j
  REAL*8 :: sum 
  REAL*8, DIMENSION(n,n) :: a
  REAL*8, DIMENSION(n) :: b
  INTEGER, DIMENSION(n) :: indx

  ii=0
  !     First we solve equation 2.3.6 of numerical recipes 
  DO i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     IF (ii /= 0)THEN
        DO j=ii,i-1
           sum=sum-a(i,j)*b(j)
        ENDDO
     ELSEIF (sum /= 0.) THEN
        ii=i
     ENDIF
     b(i)=sum
  ENDDO
  !     then we solve equation 2.3.7
  DO i=n,1,-1
     sum=b(i)
     IF (i < n) THEN
        DO j=i+1,n
           sum=sum-a(i,j)*b(j)
        ENDDO
     ENDIF
     !     store a component of the solution x in the same place as b
     b(i)=sum/a(i,i)
  ENDDO

END SUBROUTINE lu_linear_equation
! 
! routine for calculating the principal square root ( and inverse ) 
! of a general matrix AA, i.e. solving the matrix equation: AA - A = 0. 
! The routine is built on the property : 
!
!            | 0  AA |      | 0    A | 
! sign(B) = (|       |) =  (|        |)
!            | 1   0 |      | A^-1 0 |
! 
! the sign of the matrix AA is calculated using Newtons iteration method,
! see Higham et. al. ref. Numerical Algorithms 15 (1997) 227-242 
! 
SUBROUTINE sqrtmat( aa, a, a_inv ,n ) 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: aa
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: a, a_inv
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: x, x0, x1, x2
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: i_mat, temp
  INTEGER :: i,j,k
  COMPLEX*16 :: d

  ALLOCATE( x(n+n,n+n), x0(n+n,n+n), x1(n+n,n+n),x2(n+n,n+n))
  ALLOCATE( i_mat(n,n),temp(n,n))
  ! setup real identity matrix only
  i_mat = (0.D0,0.D0)
  DO i = 1, n
     i_mat(i,i) = (1.d0, 0.d0)
  ENDDO

  x0 = (0.d0,0.d0)
  x1 = (0.d0,0.d0)
  x2 = (0.d0,0.d0)
  x  = (0.d0,0.d0)
  
  temp = (0.d0,0.d0)

  DO i = n+1, 2*n
     DO j = 1, n
        x0(j,i) = aa(j,i-n)
     ENDDO
  ENDDO
  DO i = 1, n
     DO j = n+1, 2*n
        x0(j,i) = i_mat(j-n,i)
     ENDDO
  ENDDO
  k = 0
  DO WHILE( MAXVAL(ABS(temp-aa)) > 1.D-14 .AND.  k < 1000 )
     x1 = x0 
     x2 = x0 
     !CALL cmplxmatinv_lapack(x2,n+n)
     CALL cmplxmatinv(x2,n+n,d)
     x = 0.5d0 * ( x1 + x2 ) 
     x0 = x 
     k = k + 1
     DO j =1, n
        DO i = 1, n
           a(i,j) = x(i,j+n)
           a_inv(i,j) = x(i+n,j) 
        ENDDO
     ENDDO
     temp = MATMUL( a,a )
  ENDDO
  DEALLOCATE(i_mat,temp); DEALLOCATE(x,x0,x1,x2)

END SUBROUTINE sqrtmat

SUBROUTINE cmplxmatinv_lapack(a,n)

  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT)  :: a
  INTEGER :: n
  INTEGER :: ipiv(n), info, lwork
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: work
  
  !ZGETRF computes an LU factorization of a general M-by-N matrix A
  CALL zgetrf(n,n,a,n,ipiv,info)
  
  lwork = -1
  ALLOCATE(work(1))
  !ZGETRI computes the inverse of a matrix using the LU factorization computed by ZGETRF.
  CALL zgetri(n,a,n,ipiv,work,lwork,info)
  
  IF (info == 0) THEN
     lwork = INT(work(1))
     DEALLOCATE(work)
  ELSE
     WRITE(ERR,"(A,I)") 'error(cmplxmatinv_lapack): lwork query failed, info = ', info
     STOP
  END IF
  
  ALLOCATE(work(lwork))
  CALL zgetri(n,a,n,ipiv,work,lwork,info)

END SUBROUTINE cmplxmatinv_lapack


!
!    F90 program library, adapted from Numerical Recipes
!    All functions have been translated to F90 from F77
!    This is the complex*16 version of the 
!    routines to do matrix inversion, from Numerical
!    Recipes, Teukolsky et al. Routines included
!    below are MATINV, LUDCMP and LUBKSB. See chap 2
!    of Numerical Recipes for further details
!
SUBROUTINE cmplxmatinv(a,n,d)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i, j
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT)  :: a
  COMPLEX*16, ALLOCATABLE :: y(:,:)
  COMPLEX*16 :: d
  INTEGER, ALLOCATABLE :: indx(:)

  ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
  y=0.
  !     setup identity matrix
  DO i=1,n
     y(i,i)=(1.d0, 0.d0) 
  ENDDO
  !     LU decompose the matrix just once
  CALL  cmplxlu_decompose(a,n,indx,d)

  !     Find inverse by columns
  DO j=1,n
     CALL cmplxlu_linear_equation(a,n,indx,y(:,j))
  ENDDO
  !     The original matrix a was destroyed, now we equate it with the inverse y 
  a=y

  DEALLOCATE ( y ); DEALLOCATE ( indx )

END SUBROUTINE cmplxmatinv
!
!     Given an NxN matrix A(N,N), this routine replaces it by the LU 
!     decomposed one, where the matrix elements are stored in the same 
!     matrix A. The array indx is  an output vector which records the row
!     permutation effected by the partial pivoting. d is the determinant
!
SUBROUTINE cmplxlu_decompose(a,n,indx,d)

  IMPLICIT NONE
  INTEGER :: n, i, j, k, imax
  COMPLEX*16 :: sum, dum, tiny, d, aamax
  COMPLEX*16, DIMENSION(n,n) :: a
  INTEGER, DIMENSION(n) :: indx
  COMPLEX*16, ALLOCATABLE :: vv(:)

  tiny= ( 1.0D-20, 0.d0 ) 
  ALLOCATE ( vv(n) )
  D=1.
  DO i=1,n
     aamax=(0.D0,0.D0)
     DO j=1,n
        IF (ABS(a(i,j)) > ABS(aamax) ) aamax=ABS(a(i,j))
     ENDDO
     !     Zero is the largest element
     IF (aamax == (0.D0,0.D0)) STOP 'Singular matrix cmplxlu_decompose.'
     !     No nonzero largest element
     vv(i)=1./aamax
  ENDDO
  !     loop over columns
  DO j=1,n
     !     solves equation 2.3.12 except for i=j of Numerical Recipes
     IF (j > 1) THEN
        DO i=1,j-1
           sum=a(i,j)
           IF (i > 1)THEN
              DO k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
              ENDDO
              a(i,j)=sum
           ENDIF
        ENDDO
     ENDIF
     !    start searching for largest pivot element
     aamax=(0.D0,0.D0)
     DO i=j,n
        sum=a(i,j)
        IF (j > 1)THEN
           DO k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           ENDDO
           a(i,j)=sum
        ENDIF
        dum=vv(i)*ABS(sum)
        IF (ABS( dum ) >= ABS( aamax) ) THEN
           imax=i
           aamax=dum
        ENDIF
     ENDDO
     !    interchange of rows
     IF (j /= imax)THEN
        DO k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        ENDDO
        !    change of parity for determinant
        d=-d
        vv(imax)=vv(j)
     ENDIF
     indx(j)=imax
     IF(j /= n) THEN
        IF(a(j,j) == 0.) a(j,j)=tiny
        dum=1./a(j,j)
        DO i=j+1,n
           a(i,j)=a(i,j)*dum
        ENDDO
     ENDIF
     !    set up determinant
     d=d*a(j,j)
  ENDDO
  IF(a(n,n) == (0.d0,0.d0) )  a(n,n)=tiny
  DEALLOCATE ( vv)

END SUBROUTINE cmplxlu_decompose
!
!     Solves set of linear equations Ax=b, A is input as an LU decompomsed
!     matrix and indx keeps track of the permutations of the rows. b is input
!     as the right-hand side vector b and returns the solution x. A, n and indx
!     are not modified by this routine. This function takes into that b can contain
!     many zeros and is therefore suitable for matrix inversion
!
SUBROUTINE cmplxlu_linear_equation(a,n,indx,b)

  IMPLICIT NONE
  INTEGER :: n, ii, ll, i, j
  COMPLEX*16 :: sum 
  COMPLEX*16, DIMENSION(n,n) :: a
  COMPLEX*16, DIMENSION(n) :: b
  INTEGER, DIMENSION(n) :: indx

  ii=0
  !     First we solve equation 2.3.6 of numerical recipes 
  DO i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     IF (ii /= 0)THEN
        DO j=ii,i-1
           sum=sum-a(i,j)*b(j)
        ENDDO
     ELSEIF (sum /= 0.) THEN
        ii=i
     ENDIF
     b(i)=sum
  ENDDO
  !     then we solve equation 2.3.7
  DO i=n,1,-1
     sum=b(i)
     IF (i < n) THEN
        DO j=i+1,n
           sum=sum-a(i,j)*b(j)
        ENDDO
     ENDIF
     !     store a component of the solution x in the same place as b
     b(i)=sum/a(i,i)
  ENDDO

END SUBROUTINE cmplxlu_linear_equation

SUBROUTINE setup_model_space(n, npp, eigvals, eigvecs, eigvals_pp, eigvecs_pp)

  USE constants
    
  IMPLICIT NONE
  
  INTEGER    , INTENT(IN)    :: n, npp
  REAL(DP)   , INTENT(IN)    :: eigvals(n), eigvecs(n,n)
  REAL(DP)   , INTENT(INOUT) :: eigvals_pp(npp), eigvecs_pp(npp,npp)
  INTEGER                    :: i, k1

  ! split the provided eigenvectors/values
  ! into model space/excluded space projections
  
  DO k1=1, npp
     ! loop over all model space coefficients of exact eigenvectors |k>
     DO i = 1, n
        IF ( i <= npp ) THEN
           eigvecs_pp(i,k1) = eigvecs(i,k1) 
           eigvals_pp(k1) = eigvals(k1) 
        ENDIF
     ENDDO
  ENDDO
  
END SUBROUTINE setup_model_space

SUBROUTINE lee_suzuki_svd(eigvals, eigvecs, np, heff)

  USE constants
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: np
  REAL(DP), DIMENSION(np,np), INTENT(IN)    :: eigvecs
  REAL(DP), DIMENSION(np)   , INTENT(IN)    :: eigvals
  REAL(DP), DIMENSION(np,np), INTENT(INOUT) :: heff
  INTEGER                                   :: i
  REAL(DP), DIMENSION(np)                   :: S
  REAL(DP), DIMENSION(np,np)                :: U, VT, E, UVT, UVTE
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: WORK
  INTEGER                                   :: LWORK, INFO
  INTEGER, DIMENSION(8*np)                  :: IWORK
  REAL(DP)                                  :: alpha, beta

  ! singular value decomposition of the eigenvectors in P-space
  !
  ! eigvecs = U*S*VT
  ALLOCATE(WORK(1))
  WORK(1) = 0.0_dp
  LWORK = -1
  CALL DGESDD('A',np,np,eigvecs,np,S,U,np,VT,np,WORK,LWORK,IWORK,INFO)
  IF (INFO == 0) THEN
     LWORK = WORK(1)
     DEALLOCATE(WORK)
     ALLOCATE(WORK(LWORK))
  ELSE
     WRITE(ERR,*) 'error(lee_suzuki_svd): unable to determine length of workspace array ', INFO
     STOP 
  END IF
  CALL DGESDD('A',np,np,eigvecs,np,S,U,np,VT,np,WORK,LWORK,IWORK,INFO)
  
  E = 0.0_dp
  !
  DO i=1,np
     E(i,i) = eigvals(i)
  END DO
  
  !
  ! Heff = U*VT*E*(U*VT)^T
  
  ! matrix-matrix
  ! UVT = U*VT
  alpha = 1.0_dp
  beta  = 0.0_dp
  CALL DGEMM('N','N',np,np,np,alpha,U,np,VT,np,beta,UVT,np)
  !
  ! matrix-matrix
  ! UVTE = UVT*E
  CALL DGEMM('N','N',np,np,np,alpha,UVT,np,E,np,beta,UVTE,np)
  !
  ! matrix-matrix
  ! heff = UVTE*(UVT)^T
  CALL DGEMM('N','T',np,np,np,alpha,UVTE,np,UVT,np,beta,heff,np)
  
  
END SUBROUTINE lee_suzuki_svd

SUBROUTINE lee_suzuki( cvec_pp, cvec_qp, ceig, np, nq, heff)
  
  USE constants
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: np, nq
  COMPLEX*16, DIMENSION(np,np), INTENT(IN) :: cvec_pp
  COMPLEX*16, DIMENSION(nq,np), INTENT(IN)  :: cvec_qp
  COMPLEX*16, DIMENSION(np), INTENT(IN) :: ceig
  COMPLEX*16, DIMENSION(np,np), INTENT(INOUT) ::  heff
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: omega2, omega_2inv, u, u_inv, cvec_pp_inv, &
       eigen_vec, vl, heff_rhs, temp1, omega, temp2
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: ceig_p, omega2_eig
  REAL(KIND=8), DIMENSION(2*np) :: rwork
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: work1
  COMPLEX*16 :: sum1, d
  INTEGER :: i_p,j_p, k, i_q
  INTEGER :: lda, ldvl, ldvr, info, lwork
  CHARACTER*1 :: jobvl, jobvr, balanc, sense
    
  !WRITE(FRUN,*) '  --  Lee-Suzuki Similiarity Transformation  --  '

  balanc = 'n';  jobvl = 'n' ;  jobvr = 'v';  sense = 'n';  lda = np
  ldvl = 1;  ldvr = np;  lwork = -1
  ALLOCATE(omega2(np,np), omega_2inv(np,np),cvec_pp_inv(np,np)) 
  ALLOCATE(eigen_vec(np,np), vl(np,np), temp1(np,np)) 
  ALLOCATE(omega(nq,np),temp2(np,nq)); ALLOCATE(ceig_p(np), omega2_eig(np))
  eigen_vec = (0.D0,0.D0) 
  temp1 = TRANSPOSE(cvec_pp) 
  ALLOCATE(work1(1))
  CALL zgeev( jobvl, jobvr, np, temp1, lda, omega2_eig, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  IF (info == 0) THEN
     lwork = INT(work1(1))
     DEALLOCATE(work1)
  ELSE
     WRITE(ERR,"(A,I)") 'error(lee_suzuki): lwork query failed, info = ', info
     STOP
  END IF
  ALLOCATE(work1(lwork))
  CALL zgeev( jobvl, jobvr, np, temp1, lda, omega2_eig, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  ! the P->Q space transformation matrix, omega 
  cvec_pp_inv = cvec_pp
  CALL cmplxmatinv(cvec_pp_inv, np, d)
  !CALL cmplxmatinv_lapack(cvec_pp_inv, np)
  DO i_p = 1, np
     DO i_q = 1, nq
        omega(i_q,i_p) = SUM( cvec_pp_inv(:,i_p)*cvec_qp(i_q,:) )
     ENDDO
  ENDDO
  ! non-hermitian effective interaction
  ! setup 2-p effective interaction in P-space
 
  heff = (0.0D0,0.0D0)

  DO i_p = 1, np
     DO j_p = 1, np 
        heff(i_p,j_p) = SUM( cvec_pp(i_p,:)*ceig(:)*cvec_pp(j_p,:) ) 
        sum1 = (0.D0,0.D0)
        DO k = 1, np
           DO i_q = 1, nq
              sum1 = sum1 + cvec_pp(i_p,k) * ceig(k) * cvec_qp(i_q,k) * omega(i_q,j_p) 
           ENDDO
        ENDDO
        heff(i_p,j_p) = heff(i_p,j_p) + sum1
     ENDDO
  ENDDO
  ! organizing the matrix (P(1 + omega * omega)P)
  omega2 = MATMUL(TRANSPOSE(cvec_pp_inv),cvec_pp_inv)
  ! calculate sqrt and inverse sqrt of matrix (P(1 + omega * omega)P)
  ALLOCATE(u(np,np), u_inv(np,np),heff_rhs(np,np)) 
  u = (0.D0,0.D0); u_inv =  (0.D0,0.D0)
  CALL sqrtmat( omega2, u, u_inv ,np ) 
  heff_rhs =  MATMUL( u, MATMUL( heff,u_inv ) ) 
  DEALLOCATE(u, u_inv) 
  heff = (0.D0,0.D0)
  heff = heff_rhs
  ! diagonalize 2p-effective shell model hamiltonian
  CALL zgeev( jobvl, jobvr, np, heff_rhs, lda, ceig_p, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  ! compare spectrum from exact and P-space diagonalization
  !WRITE(FRUN,*) 'Compare model space two-body spectrum with exact spectrum:' 
  !DO i_p = 1, np
  !   a1 = REAL( ceig_p(i_p))
  !   a2 = AIMAG( ceig_p(i_p))
  !   b1 = REAL( ceig(i_p) )
  !   b2 = AIMAG( ceig(i_p) ) 
  !   WRITE(FRUN,"(4F20.10)") A1, A2, B1, B2
  !ENDDO
  DEALLOCATE(omega2, omega_2inv, cvec_pp_inv) 
  DEALLOCATE(eigen_vec, vl, heff_rhs, temp1); DEALLOCATE(omega) 
  DEALLOCATE(temp2); DEALLOCATE(ceig_p, omega2_eig)

END SUBROUTINE lee_suzuki

END MODULE matrix_stuff
!
! for timing printout
!
MODULE time
  
  USE aux
  IMPLICIT NONE
  
CONTAINS
  
  FUNCTION make_time (time_s) RESULT (res)
    CHARACTER(LEN=50)  :: res
    REAL(DP), INTENT(IN)   :: time_s
    INTEGER            :: diff_s, diff_m, diff_h, diff_d
    INTEGER, PARAMETER :: ms = 60, hs = 3600, ds = 86400
    REAL(DP)               :: milli
    
    diff_s = time_s
    
    milli = time_s - REAL(diff_s)
    milli = 1000.0*milli

    IF (diff_s < ms) THEN
       res = TRIM(i2a(INT(diff_s)))//' s '//TRIM(i2a(INT(milli)))//' ms'
       RETURN
       
    ELSEIF (diff_s < hs) THEN
       diff_m = INT(diff_s/ms) ; diff_s = diff_s - diff_m * ms
       res = TRIM(i2a(diff_m))//' min '//TRIM(i2a(INT(diff_s)))//' s'//TRIM(i2a(INT(milli)))//' ms'
       RETURN
       
    ELSEIF (diff_s < ds) THEN
       diff_h = INT(diff_s/hs)
       diff_m = INT(diff_s - diff_h * hs) / ms
       diff_s = diff_s - diff_h * hs - diff_m * ms
       res = TRIM(i2a(diff_h))//' hr '//TRIM(i2a(diff_m))//' min '//TRIM(i2a(INT(diff_s)))//' s'//TRIM(i2a(INT(milli)))//' ms'
       RETURN
       
    ELSE
       
       diff_d = INT(diff_s/ds)
       diff_h = INT((diff_s - diff_d * ds) / hs)
       diff_m = INT(diff_s - diff_d * ds - diff_h * hs) / ms
       diff_s = diff_s - diff_d * ds - diff_h * hs - diff_m * ms
       res = TRIM(i2a(diff_d))//' d '//TRIM(i2a(diff_h))//' hr '//TRIM(i2a(diff_m))//' min '//TRIM(i2a(INT(diff_s)))//' s'//TRIM(i2a(INT(milli)))//' ms'

    ENDIF
    
  END FUNCTION make_time
    
END MODULE time
!
! Definitions of isospin and parity parameters
!
! module written by M. Kartamyshev
!
MODULE tzpi

  USE aux
  
  IMPLICIT NONE
  
  ! definition of the many-body isospin_z constants
  TYPE, PUBLIC :: tz_info
     SEQUENCE
     INTEGER :: p, n, pp, pn, nn, ppp, ppn, pnn, nnn, dtz
  END TYPE tz_info
  TYPE (tz_info), PARAMETER, PUBLIC :: isospin_z = tz_info (-1, 1, -2, 0, 2, -3, -1, 1, 3, 2)
  
  ! definition of the parity constants
  INTEGER, PARAMETER, PUBLIC :: pi_plus = 0, pi_minus = 1, pi_step = 1 
  
  ! tzpi_info TYPE is used to ease selection of tz and parity separated structures
  TYPE, PUBLIC :: tzpi_info
     SEQUENCE
     INTEGER :: tz, pi
  END TYPE tzpi_info

  ! tzpi constants
  TYPE (tzpi_info), PARAMETER, PUBLIC :: &
       pp_plus  = tzpi_info(isospin_z%pp, pi_plus),  pp_minus  = tzpi_info(isospin_z%pp, pi_minus),&
       pn_plus  = tzpi_info(isospin_z%pn, pi_plus),  pn_minus  = tzpi_info(isospin_z%pn, pi_minus),& 
       nn_plus  = tzpi_info(isospin_z%nn, pi_plus),  nn_minus  = tzpi_info(isospin_z%nn, pi_minus),&
       ppp_plus = tzpi_info(isospin_z%ppp, pi_plus), ppp_minus = tzpi_info(isospin_z%ppp, pi_minus),&
       ppn_plus = tzpi_info(isospin_z%ppn, pi_plus), ppn_minus = tzpi_info(isospin_z%ppn, pi_minus),&
       pnn_plus = tzpi_info(isospin_z%pnn, pi_plus), pnn_minus = tzpi_info(isospin_z%pnn, pi_minus),&
       nnn_plus = tzpi_info(isospin_z%nnn, pi_plus), nnn_minus = tzpi_info(isospin_z%nnn, pi_minus)
  
  INTERFACE set_parity
     MODULE PROCEDURE set_onebody_parity, set_twobody_parity, set_threebody_parity
  END INTERFACE
  
  INTERFACE OPERATOR (==)
     MODULE PROCEDURE compare_tzpi_info
  END INTERFACE
  
  PRIVATE set_onebody_parity, set_twobody_parity, set_threebody_parity

CONTAINS
  
  ! == operator for tzpi_info type variables,  a == b
  LOGICAL FUNCTION compare_tzpi_info (a, b)
    TYPE(tzpi_info), INTENT(IN) :: a, b
    compare_tzpi_info = .FALSE.
    IF (a%tz /= b%tz) RETURN
    IF (a%pi /= b%pi) RETURN
    compare_tzpi_info = .TRUE.
  END FUNCTION compare_tzpi_info
  

  ! sets values of two one-body isospin projections corresponding to a input two-body projection value
  SUBROUTINE set_twobody_tzs (tz, tz_1, tz_2)
    INTEGER, INTENT (IN) :: tz
    INTEGER, INTENT(OUT) :: tz_1, tz_2
    SELECT CASE (tz)
     CASE (isospin_z%pp)
        tz_1 = isospin_z%p; tz_2 = isospin_z%p
     CASE (isospin_z%pn)
        tz_1 = isospin_z%p; tz_2 = isospin_z%n
     CASE (isospin_z%nn)
        tz_1 = isospin_z%n; tz_2 = isospin_z%n
     CASE DEFAULT
        WRITE (*,*) 'set_twobody_tz: wrong tz value:', tz ; STOP
     END SELECT
  END SUBROUTINE set_twobody_tzs

  ! sets values of three one-body isospin projections corresponding to input three-body projection value
  SUBROUTINE set_threebody_tzs (tz, tz_1, tz_2, tz_3)
    INTEGER, INTENT (IN) :: tz
    INTEGER, INTENT(OUT) :: tz_1, tz_2, tz_3
    SELECT CASE (tz)
    CASE (isospin_z%ppp)
       tz_1 = isospin_z%p; tz_2 = isospin_z%p; tz_3 = isospin_z%p
    CASE (isospin_z%ppn)
       tz_1 = isospin_z%p; tz_2 = isospin_z%p; tz_3 = isospin_z%n
    CASE (isospin_z%pnn)
       tz_1 = isospin_z%p; tz_2 = isospin_z%n; tz_3 = isospin_z%n
    CASE (isospin_z%nnn)
       tz_1 = isospin_z%n; tz_2 = isospin_z%n; tz_3 = isospin_z%n
    CASE DEFAULT
       WRITE (*,*) 'set_threebody_tz: wrong tz value:', tz ; STOP
    END SELECT
  END SUBROUTINE set_threebody_tzs

  ! checks one-body isospin value for consistency with internal program representation
  SUBROUTINE verify_onebody_isospin_value (tz)
    INTEGER, INTENT (IN) :: tz
    SELECT CASE (tz)
    CASE (isospin_z%p)
    CASE (isospin_z%n)
    CASE DEFAULT
       WRITE(*,*) 'verify_onebody_isospin_value: wrong tz value: ', tz; STOP
    END SELECT
  END SUBROUTINE verify_onebody_isospin_value
  
  ! checks two-body isospin value for consistency with internal program representation
  SUBROUTINE verify_twobody_isospin_value (tz)
    INTEGER, INTENT (IN) :: tz
    SELECT CASE (tz)
    CASE (isospin_z%pp)
    CASE (isospin_z%pn)
    CASE (isospin_z%nn)
    CASE DEFAULT
       WRITE(*,*) 'verify_twobody_isospin_value: wrong tz value: ', tz; STOP
    END SELECT
  END SUBROUTINE verify_twobody_isospin_value

  ! checks three-body isospin value for consistency with internal program representation
  SUBROUTINE verify_threebody_isospin_value (tz)
    INTEGER, INTENT (IN) :: tz
    SELECT CASE (tz)
    CASE (isospin_z%ppp)
    CASE (isospin_z%ppn)
    CASE (isospin_z%pnn)
    CASE (isospin_z%nnn)
    CASE DEFAULT
       WRITE(*,*) 'verify_threebody_isospin_value: wrong tz value: ', tz; STOP
    END SELECT
  END SUBROUTINE verify_threebody_isospin_value

  ! verifies two-body isospin_z configuration
  SUBROUTINE verify_twobody_isospin_config (tz_1, tz_2, tz)
    INTEGER, INTENT (IN) :: tz_1, tz_2, tz
    CALL verify_onebody_isospin_value (tz_1)
    CALL verify_onebody_isospin_value (tz_2)
    CALL verify_twobody_isospin_value (tz)
    IF(tz_1 + tz_2 /= tz) THEN
       WRITE (*,*) 'verify_twobody_isospin_config: inconsistent configuration: ', tz_1, tz_2, tz; STOP
    ENDIF
  END SUBROUTINE verify_twobody_isospin_config
  
  ! verifies three-body isospin_z configuration
  SUBROUTINE verify_threebody_isospin_config (tz_1, tz_2, tz_3, tz)
    INTEGER, INTENT (IN) :: tz_1, tz_2, tz_3, tz
    CALL verify_onebody_isospin_value (tz_1)
    CALL verify_onebody_isospin_value (tz_2)
    CALL verify_onebody_isospin_value (tz_3)
    CALL verify_threebody_isospin_value (tz)
    IF(tz_1 + tz_2 + tz_3 /= tz) THEN
       WRITE (*,*) 'verify_threebody_isospin_config: inconsistent configuration: ', tz_1, tz_2, tz_3, tz; STOP
    ENDIF
  END SUBROUTINE verify_threebody_isospin_config

  ! sets values of parity for two-body configuration
  FUNCTION set_onebody_parity (l1) RESULT (res)
    INTEGER, INTENT (IN) :: l1
    INTEGER              :: res, parity
    parity = MOD(l1, 2)
    SELECT CASE (parity)
    CASE (0)
       res = pi_plus
    CASE (1)
       res = pi_minus
    CASE DEFAULT
       WRITE (*,*) 'set_onebody_parity: wrong parity value: ', parity ; STOP
    END SELECT
  END FUNCTION set_onebody_parity
  
  ! sets values of parity for two-body configuration
  ! INPUT TRUE VALUES, i.e. NOT DOUBLED
  FUNCTION set_twobody_parity (l1, l2) RESULT (res)
    INTEGER, INTENT (IN) :: l1, l2
    INTEGER              :: res, parity
    parity = MOD( (l1 + l2), 2 )
    SELECT CASE (parity)
    CASE (0)
       res = pi_plus
    CASE (1)
       res = pi_minus
    CASE DEFAULT
       WRITE (*,*) 'set_twobody_parity: wrong parity value: ', parity ; STOP
    END SELECT
  END FUNCTION set_twobody_parity

  ! sets values of parity for three-body configuration
  FUNCTION set_threebody_parity (l1, l2, l3) RESULT (res)
    INTEGER, INTENT (IN) :: l1, l2, l3
    INTEGER              :: res, parity
    parity = MOD( (l1 + l2 + l3), 2 )
    SELECT CASE (parity)
    CASE (0)
       res = pi_plus
    CASE (1)
       res = pi_minus
    CASE DEFAULT
       WRITE (*,*) 'set_threebody_parity: wrong parity value: ', parity; STOP
    END SELECT
  END FUNCTION set_threebody_parity
  
  ! checks parity variable value for consistency with internal program representation
  SUBROUTINE verify_parity_value (parity)
    INTEGER, INTENT (IN) :: parity
    SELECT CASE (parity)
    CASE (pi_plus)
    CASE (pi_minus)
    CASE DEFAULT
       WRITE (*,*) 'verify_parity_value: unknown parity value:', parity ; STOP
    END SELECT
  END SUBROUTINE verify_parity_value
  
  ! inverts parity value
  INTEGER FUNCTION invert_parity (parity)
    INTEGER, INTENT (IN) :: parity
    SELECT CASE (parity)
    CASE (pi_minus)
       invert_parity = pi_plus
    CASE (pi_plus)
       invert_parity = pi_minus
    CASE DEFAULT
       WRITE (*,*) 'invert_parity: wrong parity value: ', parity ; STOP
    END SELECT
  END FUNCTION invert_parity

  ! returns character value corresponding to a given one-body isospin projection value
  FUNCTION tz2a (tz)
    CHARACTER(LEN=1)     :: tz2a
    INTEGER, INTENT (IN) :: tz
    SELECT CASE (tz)
    CASE (isospin_z%p)       
       tz2a = 'p'
    CASE (isospin_z%n)
       tz2a = 'n'
    CASE DEFAULT
       WRITE (*,*) 'tz2a: unknown one-body isospin projection value: ', tz ; STOP
    END SELECT
  END FUNCTION tz2a

  ! returns character value corresponding to a given two-body isospin projection value
  FUNCTION tz2ab (tz)
    CHARACTER(LEN=2)     :: tz2ab
    INTEGER, INTENT (IN) :: tz
    SELECT CASE (tz)
    CASE (isospin_z%pp)
       tz2ab = 'pp'
    CASE (isospin_z%pn)
       tz2ab = 'pn'
    CASE (isospin_z%nn)
       tz2ab = 'nn'
    CASE DEFAULT
       WRITE (*,*) 'tz2ab: unknown two-body isospin projection value: ', tz ; STOP
    END SELECT
  END FUNCTION tz2ab

  ! returns character value corresponding to a given three-body isospin projection value
  FUNCTION tz2abc (tz)
    CHARACTER(LEN=3)     :: tz2abc
    INTEGER, INTENT (IN) :: tz
    
    IF (.NOT. NORMALORDER) THEN
       SELECT CASE (tz)
       CASE (isospin_z%ppp)
          tz2abc = 'ppp'
       CASE (isospin_z%ppn)
          tz2abc = 'ppn'
       CASE (isospin_z%pnn)
          ! We store the Tz=+1/2 in neutron-neutron-proton order
          ! 2012-08-30 / Andreas
          tz2abc = 'nnp'
       CASE (isospin_z%nnn)
          tz2abc = 'nnn'
       CASE DEFAULT
          WRITE (*,*) 'tz2abc: unknown three-body isospin projection value: ', tz ; STOP
       END SELECT
    ELSE IF(NORMALORDER) THEN 
       SELECT CASE (tz)
       CASE (isospin_z%ppp)
          tz2abc = 'ppp'
       CASE (isospin_z%ppn)
          tz2abc = 'pnp'
       CASE (isospin_z%pnn)
          tz2abc = 'pnn'
       CASE (isospin_z%nnn)
          tz2abc = 'nnn'
       CASE DEFAULT
          WRITE (*,*) 'tz2abc: unknown three-body isospin projection value: ', tz ; STOP
       END SELECT
    END IF
  END FUNCTION tz2abc

  ! returns character value corresponding to a given parity value
  FUNCTION pi2a (parity)
    CHARACTER(LEN=1)     :: pi2a
    INTEGER, INTENT (IN) :: parity
    SELECT CASE (parity)
    CASE (pi_plus)
       pi2a = '+'
    CASE (pi_minus)
       pi2a = '-'
    CASE DEFAULT
       WRITE (*,*) 'pi2a: unknown parity value: ', parity ; STOP
    END SELECT
  END FUNCTION pi2a

  ! returns spectroscopic notation
  FUNCTION spectroscopic_notation (nn,ll,jj)

    USE aux

    CHARACTER(LEN=5)     :: spectroscopic_notation
    INTEGER, INTENT (IN) :: nn,ll,jj

    spectroscopic_notation = TRIM(ADJUSTL(i2a(nn)))//TRIM(ADJUSTL(ll2a(ll)))//TRIM(ADJUSTL(i2a(jj)))
    
  END FUNCTION spectroscopic_notation

  FUNCTION spectroscopic_notation_z (nn,ll,jj,tz)

    USE aux

    CHARACTER(LEN=8)     :: spectroscopic_notation_z
    INTEGER, INTENT (IN) :: nn,ll,jj,tz

    spectroscopic_notation_z = TRIM(ADJUSTL(tz2a(tz)))//'('//TRIM(ADJUSTL(i2a(nn)))//TRIM(ADJUSTL(ll2a(ll)))//TRIM(ADJUSTL(i2a(jj)))//')'
    
  END FUNCTION spectroscopic_notation_z

  FUNCTION ll2a (ll)
    CHARACTER(LEN=1)     :: ll2a
    INTEGER, INTENT (IN) :: ll
    
    IF (ll>23) ll2a = i2a(ll)
    
    SELECT CASE (ll)
    CASE (0)
       ll2a = 's'
    CASE (1)
       ll2a = 'p'
    CASE (2)
       ll2a = 'd'
    CASE (3)
       ll2a = 'f'
    CASE (4)
       ll2a = 'g'
    CASE (5)
       ll2a = 'h'
    CASE (6)
       ll2a = 'i'
    CASE (7)
       ll2a = 'j'
    CASE (8)
       ll2a = 'k'
    CASE (9)
       ll2a = 'l'
    CASE (10)
       ll2a = 'm'
    CASE (11)
       ll2a = 'n'
    CASE (12)
       ll2a = 'o'
    CASE (13)
       ll2a = 'p'
    CASE (14)
       ll2a = 'q'
    CASE (15)
       ll2a = 'r'
    CASE (16)
       ll2a = 's'
    CASE (17)
       ll2a = 't'
    CASE (18)
       ll2a = 'u'
    CASE (19)
       ll2a = 'v'
    CASE (20)
       ll2a = 'w'
    CASE (21)
       ll2a = 'x'
    CASE (22)
       ll2a = 'y'
    CASE (23)
       ll2a = 'z'
    CASE DEFAULT
       ll2a = '?'
    END SELECT
  END FUNCTION ll2a

END MODULE tzpi
!
!
!
MODULE single_particle_basis
  
  USE constants
  USE aux
  USE tzpi

  IMPLICIT NONE
  
  TYPE :: splimits
     TYPE(limits_type) :: oscn, nn, ll, jj, tz, pi, e
  END TYPE splimits
  
  TYPE :: nllimits
     TYPE(limits_type) :: oscn, nn, ll
  END TYPE nllimits

  TYPE, PUBLIC :: nlorbit
     SEQUENCE
     INTEGER :: address, oscn, nn, ll
  END TYPE nlorbit

  TYPE, PUBLIC :: nlNLorbit_type
     SEQUENCE
     INTEGER :: address, oscn, nn, ll, ncm, lcm 
  END TYPE nlNLorbit_type

  TYPE, PUBLIC :: nlnlNLorbit_type
     SEQUENCE
     INTEGER :: address, oscn, oscn12, oscn3, oscncm, n12, l12, n3, l3, ncm, lcm 
  END TYPE nlnlNLorbit_type

  TYPE, PUBLIC :: sporbit
     ! nn = 0,1,2,...
     ! ll = 0,2,4,... (doubled)
     ! jj = 1,2,3,... (doubled)
     ! tz = -1, +1    (doubled)
     INTEGER :: address, oscn, nn, ll, jj, tz, pi
     REAL(DP) :: e
     LOGICAL  :: hole
  END TYPE sporbit

  TYPE, PUBLIC :: sporbit_m
     ! nn = 0,1,2,...
     ! ll = 0,2,4,... (doubled)
     ! jj = 1,2,3,... (doubled)
     ! tj = -j,..+j   (doubled)
     ! tz = -1, +1    (doubled)
     INTEGER :: address, oscn, nn, ll, jj, jz, tz, pi
     TYPE(sporbit), POINTER :: jorbit
     REAL(dp) :: e
     LOGICAL  :: hole
  END TYPE sporbit_m
  
  TYPE, PUBLIC :: spbasis_type
     INTEGER         :: amount
     TYPE(splimits)  :: limits
     TYPE(sporbit), DIMENSION(:), ALLOCATABLE :: orbit
     INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: nljtz2sporbit
  END TYPE spbasis_type
  
  TYPE, PUBLIC :: spbasis_m_type
     INTEGER         :: amount
     TYPE(splimits)  :: limits
     TYPE(sporbit_m), DIMENSION(:), ALLOCATABLE :: orbit
  END TYPE spbasis_m_type

  TYPE, PUBLIC :: nlbasis_type
     INTEGER        :: amount
     TYPE(nllimits) :: limits
     TYPE(nlorbit), DIMENSION(:), ALLOCATABLE :: orbit
  END TYPE nlbasis_type

  TYPE, PUBLIC :: nlNLbasis_type
     INTEGER        :: amount
     TYPE(nllimits) :: limits
     TYPE(nlNLorbit_type), DIMENSION(:), ALLOCATABLE :: orbit
  END TYPE nlNLbasis_type

  TYPE, PUBLIC :: nlnlNLbasis_type
     INTEGER        :: amount
     TYPE(nllimits) :: limits
     TYPE(nlnlNLorbit_type), DIMENSION(:), ALLOCATABLE :: orbit
  END TYPE nlnlNLbasis_type
  
  TYPE(spbasis_type), TARGET :: spbasis
  TYPE(spbasis_m_type), TARGET :: spbasis_m
  TYPE(nlbasis_type) ::rel_nl
  TYPE(nlbasis_type) :: rel_nl_ab ! threebody rel_nl lists, used in conjunction with howf setup
  TYPE(nlbasis_type) :: rel_nl_c  ! threebody rel_nl lists, used in conjunction with howf setup
  ! rel+cm for twobody
  ! n: relative n
  ! l: relative l
  ! N: cm n
  ! L: cm l
  TYPE(nlNLbasis_type), TARGET :: relcm_nlNL
  ! rel+cm for threebody
  ! nlnl relative n12 l12 n3 l3
  ! NL cm nl
  TYPE(nlnlNLbasis_type), TARGET :: relcm_nlnlNL
  
CONTAINS
  
  ! given the Nmax parameter read on input
  ! and accessible via constants, all the 
  ! allowed orbits in the spbasis object is 
  ! allocated. 
  SUBROUTINE construct_single_particle_basis(this)
    
    USE constants
    
    IMPLICIT NONE
    
    TYPE(spbasis_type) :: this
    INTEGER            :: n, l, tz, j, j_min, j_max, icount
    INTEGER            :: oscn, npart, max_cfg_energy
    INTEGER            :: max_ll, max_nn, nprotons, nneutrons

    ! count the number of single particle orbits given
    ! the boundaries
    
    ! n is not doubled
    ! l is doubled
    ! oscn%max is not doubled
    
    IF (nbody_basis == 'twobody'  ) npart = 2
    IF (nbody_basis == 'threebody') npart = 3
    
    max_cfg_energy = 0
    max_nn = 0
    max_ll = 0

    icount = 0
    DO n = this%limits%nn%min, this%limits%nn%max, this%limits%nn%step
       DO l = this%limits%ll%min, this%limits%ll%max, this%limits%ll%step
          
          ! for NPmax controlled space, ensure to stay within it
          IF (TRIANGULAR .AND. (2*n+l/2 > this%limits%oscn%max)        ) CYCLE
          IF (SQUARE     .AND. (npart*(2*n+l/2) > this%limits%oscn%max)) CYCLE

          j_min = l-1; j_max = l+1
          IF (j_min < 0) j_min = 1
          DO j = j_min, j_max, 2
             DO tz=this%limits%tz%min, this%limits%tz%max, this%limits%tz%step
                icount = icount + 1
             END DO
          ENDDO
       ENDDO
    ENDDO
    
    this%amount = icount
    
    IF (ALLOCATED(this%orbit)) THEN
       WRITE(ERR,"(A)") 'error(construct_single_particle_basis): this%orbit already allocated'
       STOP
    END IF
    
    ALLOCATE(this%orbit(1:this%amount))
    
    ALLOCATE(this%nljtz2sporbit(this%limits%nn%min:this%limits%nn%max,&       ! not doubled
                                this%limits%ll%min/2:this%limits%ll%max/2,&   ! not doubled
                                this%limits%ll%min-1:this%limits%ll%max+1,&   ! doubled
                                this%limits%tz%min:this%limits%tz%max))       ! doubled
                       
    this%nljtz2sporbit = irrel

    ! fill orbit array, and lookup
    icount = 0
    nprotons = 0
    nneutrons = 0
    
    DO oscn=this%limits%oscn%min, this%limits%oscn%max, this%limits%oscn%step
       DO n = this%limits%nn%min, this%limits%nn%max, this%limits%nn%step
          DO l = this%limits%ll%min, this%limits%ll%max, this%limits%ll%step
             IF (2*n+l/2 /= oscn) CYCLE
             IF (TRIANGULAR .AND. (2*n+l/2 > this%limits%oscn%max)        ) CYCLE
             IF (SQUARE     .AND. (npart*(2*n+l/2) > this%limits%oscn%max)) CYCLE
             j_min = l-1; j_max = l+1
             IF (j_min < 0) j_min = 1
             DO j = j_max, j_min, -2
                DO tz=this%limits%tz%min, this%limits%tz%max, this%limits%tz%step
                   icount = icount + 1
                   this%orbit(icount)%address = icount
                   this%orbit(icount)%oscn    = oscn
                   this%orbit(icount)%nn      = n
                   this%orbit(icount)%ll      = l
                   this%orbit(icount)%jj      = j
                   this%orbit(icount)%tz      = tz
                   this%orbit(icount)%pi      = set_parity(l/2)
                   this%orbit(icount)%e       = (oscn+1.5_dp)*hbar_omega
                   this%nljtz2sporbit(n,l/2,j,tz) = icount
                   this%orbit(icount)%hole    = .TRUE.
                   IF (tz == isospin_z%p) nprotons  = nprotons  + (j+1)
                   IF (tz == isospin_z%n) nneutrons = nneutrons + (j+1)
                   IF (nprotons  > number_of_protons  .AND. tz == isospin_z%p) this%orbit(icount)%hole    = .FALSE.
                   IF (nneutrons > number_of_neutrons .AND. tz == isospin_z%n) this%orbit(icount)%hole    = .FALSE.
                   IF (oscn > max_cfg_energy) max_cfg_energy = oscn
                   IF (l    > max_ll) max_ll = l
                   IF (n    > max_nn) max_nn = n
                END DO
             END DO
          END DO
       END DO
    END DO
    
    IF (((number_of_neutrons + number_of_protons) > (nprotons +nneutrons)) .AND. NORMALORDER) THEN
       WRITE(ERR,"(A)") 'error(construct_single_particle_basis): fermi vaccuum larger than single particle basis'
       STOP
    END IF

    IF (max_ll < this%limits%ll%max ) this%limits%ll%max = max_ll
    IF (max_nn < this%limits%nn%max ) this%limits%nn%max = max_nn

    ! do a check when using user defined lmax and/or nmax to ensure
    ! that an unecessary large NPmax was not set

    IF (SQUARE) max_cfg_energy = npart * max_cfg_energy
    
    IF (max_cfg_energy < NPmax) THEN
       WRITE(ERR,*)
       WRITE(ERR,"(A)") '(tip): your setting of nmax_lab/lmax_lab or NPmax is such that,'
       WRITE(ERR,"(A)") '       the actual model space of nbody configurations allows an'
       WRITE(ERR,"(A,I2,A)") '       energy up to ', max_cfg_energy,' while you defined the'
       WRITE(ERR,"(A,I2)") '       interaction in a space up to ', NPmax
       WRITE(ERR,"(A,I2,A)") '       Thus, set NPmax to ', max_cfg_energy,' and rerun!'
       
       ! print it anyway, just for user inspection.
       CALL print_single_particle_basis(FSP,this)
       
       STOP
    END IF

  END SUBROUTINE construct_single_particle_basis

  ! construct an m-scheme basis given a jj-scheme basis
  SUBROUTINE construct_single_particle_basis_m(this_jj, this_m)
    
    IMPLICIT NONE
    
    TYPE(spbasis_type)  , TARGET :: this_jj
    TYPE(spbasis_m_type)         :: this_m
    TYPE(sporbit), POINTER :: jjorbit
    INTEGER            :: icount, tz
    INTEGER            :: ijjorbit, im
        
    ! count the number of single particle m-orbits
    icount = 0
    DO ijjorbit=1, this_jj%amount
       jjorbit => this_jj%orbit(ijjorbit)
       icount = icount + jjorbit%jj+1
    ENDDO
    
    this_m%amount = icount
    
    IF (ALLOCATED(this_m%orbit)) THEN
       WRITE(ERR,"(A)") 'error(construct_single_particle_basis_m): this%orbit already allocated'
       STOP
    END IF
    
    ALLOCATE(this_m%orbit(1:this_m%amount))
    
    ! fill orbit array
    icount = 1
    DO tz = this_jj%limits%tz%min, this_jj%limits%tz%max, this_jj%limits%tz%step
       DO ijjorbit=1, this_jj%amount
          jjorbit => this_jj%orbit(ijjorbit)
          IF (jjorbit%tz /= tz) CYCLE
          DO im=-jjorbit%jj, jjorbit%jj, 2
             this_m%orbit(icount)%address = icount
             this_m%orbit(icount)%oscn    = jjorbit%oscn
             this_m%orbit(icount)%nn      = jjorbit%nn
             this_m%orbit(icount)%ll      = jjorbit%ll
             this_m%orbit(icount)%jj      = jjorbit%jj
             this_m%orbit(icount)%jz      = im
             this_m%orbit(icount)%tz      = jjorbit%tz
             this_m%orbit(icount)%pi      = jjorbit%pi
             this_m%orbit(icount)%e       = jjorbit%e
             this_m%orbit(icount)%hole    = jjorbit%hole
             this_m%orbit(icount)%jorbit => this_jj%orbit(ijjorbit)
             icount = icount+1
          END DO
       END DO
    END DO
    
  END SUBROUTINE construct_single_particle_basis_m
  
  SUBROUTINE zero_sporbit_m(spm) 

    TYPE(sporbit_m), INTENT(INOUT) :: spm

    spm%address = irrel
    spm%oscn    = irrel
    spm%nn = irrel
    spm%ll = irrel
    spm%jj = irrel
    spm%jz = irrel
    spm%tz = irrel
    spm%pi = irrel
    IF (ASSOCIATED(spm%jorbit)) NULLIFY(spm%jorbit)
    spm%e = -99.0_dp
        
  END SUBROUTINE zero_sporbit_m
  
  SUBROUTINE construct_nl_basis(this)
    
    IMPLICIT NONE
    
    TYPE(nlbasis_type), INTENT(INOUT) :: this
    
    CALL count_nl_states(this%amount, this%limits)
    
    IF (ALLOCATED(this%orbit)) THEN
       WRITE(ERR,"(A)") 'error(construct_nl_basis): this%orbit already allocated'
       STOP
    END IF
    
    ALLOCATE(this%orbit(1:this%amount))
    
    CALL count_nl_states(this%amount, this%limits, this)

  END SUBROUTINE construct_nl_basis

  ! given nmax and lmax, the number of 
  ! |nl> states is counted
  SUBROUTINE count_nl_states(icount, limits, this)
    
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: icount
    TYPE(nlbasis_type), INTENT(INOUT), OPTIONAL :: this
    TYPE(nllimits), INTENT(IN) :: limits
    INTEGER :: n, l, oscn

    ! count the number of single particle orbits given
    ! the boundaries above

    ! sorted in ascending order, i.e. oscn = 0,1,2,3,...
    
    ! THE REASON FOR DIVIDING WITH limits%ll%step
    ! IS THAT IT ENABLED THIS ROUTINE TO BE USED
    ! WITH BOTH DOUBLED AND NOT DOUBLED ll VALUES
    ! THE ll VALUES ARE IN GENERAL NOT DOUBLED IN
    ! THE TWOBODY SECTOR OF THE CODE, WHEREAS THEY
    ! ARE DOUBLED IN THE THREEBODY SECTOR OF THE
    ! CODE. THIS IS INCONSISTENT AND NOT VERY 
    ! ESTETHIC. SHOULD BE CHANGED!!

    icount = 0
    DO oscn = limits%oscn%min, limits%oscn%max, limits%oscn%step
       DO n = limits%nn%min, limits%nn%max, limits%nn%step
          DO l = limits%ll%min, limits%ll%max, limits%ll%step
             IF (2*n+l/limits%ll%step /= oscn) CYCLE
             icount = icount + 1
             IF (PRESENT(this)) THEN
                this%orbit(icount) = nlorbit(icount, 2*n+l/limits%ll%step, n, l)
             END IF
          ENDDO
       ENDDO
    END DO
  END SUBROUTINE count_nl_states

  SUBROUTINE construct_nlNL_basis(this)
    
    IMPLICIT NONE
    
    TYPE(nlNLbasis_type), INTENT(INOUT) :: this
    
    CALL count_nlNL_states(this%amount, this%limits)
    
    IF (ALLOCATED(this%orbit)) THEN
       WRITE(ERR,"(A)") 'error(construct_nlNL_basis): this%orbit already allocated'
       STOP
    END IF
    
    ALLOCATE(this%orbit(1:this%amount))
    
    CALL count_nlNL_states(this%amount, this%limits, this)
    
  END SUBROUTINE construct_nlNL_basis

  ! given nmax and lmax, the number of 
  ! |nlNL> states is counted
  SUBROUTINE count_nlNL_states(icount, limits, this)
    
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: icount
    TYPE(nlNLbasis_type), INTENT(INOUT), OPTIONAL :: this
    TYPE(nllimits), INTENT(IN) :: limits
    INTEGER :: nn, ll, ncm, lcm, oscn

    ! count the number of single particle orbits given
    ! the boundaries above

    ! sorted in ascending order, i.e. oscn = 0,1,2,3,...
    ! where oscn = 2n+l + 2N+L
    icount = 0
    DO oscn = limits%oscn%min, limits%oscn%max, limits%oscn%step
       DO nn = limits%nn%min, limits%nn%max, limits%nn%step
          DO ll = limits%ll%min, limits%ll%max, limits%ll%step
             DO ncm= limits%nn%min, limits%nn%max, limits%nn%step
                DO lcm = limits%ll%min, limits%ll%max, limits%ll%step
                   IF (2*nn+ll + 2*ncm+lcm /= oscn) CYCLE
                   icount = icount + 1
                   IF (PRESENT(this)) THEN
                      this%orbit(icount) = nlNLorbit_type(icount, 2*nn+ll+2*ncm+lcm, nn, ll, ncm, lcm)
                   END IF
                ENDDO
             ENDDO
          END DO
       END DO
    END DO
  END SUBROUTINE count_nlNL_states
  
  SUBROUTINE construct_nlnlNL_basis(this)
    
    IMPLICIT NONE
    
    TYPE(nlnlNLbasis_type), INTENT(INOUT) :: this
    
    CALL count_nlnlNL_states(this%amount, this%limits)
    
    IF (ALLOCATED(this%orbit)) THEN
       WRITE(ERR,"(A)") 'error(construct_nlnlNL_basis): this%orbit already allocated'
       STOP
    END IF
    
    ALLOCATE(this%orbit(1:this%amount))
    
    CALL count_nlnlNL_states(this%amount, this%limits, this)
    
  END SUBROUTINE construct_nlnlNL_basis

  ! given nmax and lmax, the number of 
  ! |nlnlNL> states is counted
  SUBROUTINE count_nlnlNL_states(icount, limits, this)
    
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: icount
    TYPE(nlnlNLbasis_type), INTENT(INOUT), OPTIONAL :: this
    TYPE(nllimits), INTENT(IN) :: limits
    INTEGER :: n12, l12, n3, l3, ncm, lcm, oscn

    ! count the number of single particle orbits given
    ! the boundaries above

    ! sorted in ascending order, i.e. oscn = 0,1,2,3,...
    ! where oscn = 2n12+l12 + 2n3+l3 + 2ncm+lcm
    icount = 0
    DO oscn = limits%oscn%min, limits%oscn%max, limits%oscn%step
       DO n12 = limits%nn%min, limits%nn%max, limits%nn%step
          DO l12 = limits%ll%min, limits%ll%max, limits%ll%step
             DO n3 = limits%nn%min, limits%nn%max, limits%nn%step
                DO l3 = limits%ll%min, limits%ll%max, limits%ll%step
                   DO ncm= limits%nn%min, limits%nn%max, limits%nn%step
                      DO lcm = limits%ll%min, limits%ll%max, limits%ll%step
                         ! NOTE, n is not doubled, l is doubled
                         IF (2*n12+l12/2 + 2*n3+l3/2 + 2*ncm+lcm/2 /= oscn) CYCLE
                         icount = icount + 1
                         IF (PRESENT(this)) THEN
                            this%orbit(icount) = nlnlNLorbit_type(icount, &
                                                                  2*n12+l12/2+2*n3+l3/2+2*ncm+lcm/2,&
                                                                  2*n12+l12/2, &
                                                                  2*n3+l3/2, &
                                                                  2*ncm+lcm/2,&
                                                                  n12, l12, n3, l3, ncm, lcm)
                         END IF
                      ENDDO
                   ENDDO
                END DO
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE count_nlnlNL_states
  
  SUBROUTINE orbit_swap(a,b)
    TYPE(sporbit), INTENT(INOUT) :: a,b
    TYPE(sporbit) :: temp_a, temp_b
    
    temp_a = a
    temp_b = b
    
    a=temp_b
    b=temp_a
    
  END SUBROUTINE orbit_swap
  
  SUBROUTINE print_single_particle_basis(unt, this)
    
    INTEGER, INTENT(IN)            :: unt
    TYPE(spbasis_type), INTENT(IN) :: this
    INTEGER :: i
    LOGICAL :: UNT_IS_OPEN
    
    INQUIRE(unt, opened=UNT_IS_OPEN)
    IF (.NOT. UNT_IS_OPEN) RETURN
    WRITE(unt,"(A)") 'P-SPACE SINGLE PARTICLE BASIS' 
    WRITE(unt,"(A14,F20.6)") 'hbar_omega  = ', hbar_omega
    IF (.NOT. NORMALORDER) WRITE(unt,"(A10,5A5,A7,A5,A13)") '#','n','l','2*j','tz','pi','','2n+l','energy'
    IF (NORMALORDER) WRITE(unt,"(A10,5A5,A7,A5,A13,A15)") '#','n','l','2*j','tz','pi','','2n+l','energy','particle/hole'
    
    DO i=1, this%amount
       WRITE(unt,"(I10,3I5,2A5,A7,I5,F13.2)",ADVANCE='NO') i, this%orbit(i)%nn, &
                                              this%orbit(i)%ll/2 , &
                                              this%orbit(i)%jj   , &
                                              tz2a(this%orbit(i)%tz) ,&
                                              pi2a(this%orbit(i)%pi) ,&
                                              spectroscopic_notation(this%orbit(i)%nn , &
                                                                     this%orbit(i)%ll/2 , &
                                                                     this%orbit(i)%jj), &
                                              this%orbit(i)%oscn, &
                                              this%orbit(i)%e
       IF (NORMALORDER) THEN
          IF (this%orbit(i)%hole) THEN
             WRITE(unt,"(A15)") 'holes'
          ELSE
             WRITE(unt,"(A15)") 'particles'
          END IF
       ELSE
          WRITE(unt,*)
       END IF
       
    END DO
    
  END SUBROUTINE print_single_particle_basis

  SUBROUTINE print_single_particle_basis_m(unt, this)
    
    INTEGER, INTENT(IN)            :: unt
    TYPE(spbasis_m_type), INTENT(IN) :: this
    INTEGER :: i
    LOGICAL :: UNT_IS_OPEN
    
    INQUIRE(unt, opened=UNT_IS_OPEN)
    IF (.NOT. UNT_IS_OPEN) RETURN
    WRITE(unt,"(A)") 'P-SPACE SINGLE PARTICLE M-SCHEME BASIS' 
    IF (.NOT. NORMALORDER) WRITE(unt,"(A10,6A5,A7,A5,A13)") '#','n','l','2*j','2*jz','tz','pi','','2n+l','energy'
    IF (NORMALORDER) WRITE(unt,"(A10,6A5,A7,A5,A13,A15)") '#','n','l','2*j','2*jz','tz','pi','','2n+l','energy','particle/hole'
    DO i=1, this%amount
       IF (model_space == 'square') THEN
          IF (nbody_basis == 'twobody' .AND. (this%orbit(i)%oscn > this%limits%oscn%max/2)) CYCLE
          IF (nbody_basis == 'threebody' .AND. (this%orbit(i)%oscn > this%limits%oscn%max/3)) CYCLE
       END IF
       WRITE(unt,"(I10,4I5,2A5,A7,I5,F13.2)", ADVANCE='NO') i, &
                                              this%orbit(i)%nn, &
                                              this%orbit(i)%ll/2 , &
                                              this%orbit(i)%jj   , &
                                              this%orbit(i)%jz   , &
                                              tz2a(this%orbit(i)%tz) ,&
                                              pi2a(this%orbit(i)%pi) ,&
                                              spectroscopic_notation(this%orbit(i)%nn , &
                                                                     this%orbit(i)%ll/2 , &
                                                                     this%orbit(i)%jj), &
                                              this%orbit(i)%oscn, &
                                              this%orbit(i)%e

       IF (NORMALORDER) THEN
          IF (this%orbit(i)%hole) THEN
             WRITE(unt,"(A15)") 'holes'
          ELSE
             WRITE(unt,"(A15)") 'particles'
          END IF
       ELSE
          WRITE(unt,*)
       END IF

    END DO
    
  END SUBROUTINE print_single_particle_basis_m

  FUNCTION single_particle_basis_memory(spbasis) RESULT(sizeof_spbasis_type_char)
    
    IMPLICIT NONE

    TYPE(spbasis_type), INTENT(IN) :: spbasis
    INTEGER, PARAMETER :: sizeof_sporbit_type = 7*sizeof_int+sizeof_double
    INTEGER, PARAMETER :: sizeof_splimits     = 6*sizeof_limits_type
    INTEGER            :: sizeof_spbasis_type
    CHARACTER(LEN=100) :: sizeof_spbasis_type_char

    sizeof_spbasis_type = spbasis%amount*sizeof_sporbit_type + &
                          sizeof_splimits + sizeof_int ! Bytes

    sizeof_spbasis_type_char = 'sizeof_sp_basis_type: '//&
         TRIM(ADJUSTL(r2a(REAL(sizeof_spbasis_type)/2.0**10)))// ' KBytes'

  END FUNCTION single_particle_basis_memory

  SUBROUTINE print_nl_basis(unt, this)
    
    INTEGER, INTENT(IN)            :: unt
    TYPE(nlbasis_type), INTENT(IN) :: this
    INTEGER :: i
    LOGICAL :: UNT_IS_OPEN
    
    INQUIRE(unt, opened=UNT_IS_OPEN)
    IF (.NOT. UNT_IS_OPEN) RETURN
    WRITE(unt,"(A)") 'NL BASIS'
    WRITE(unt,"(A6,3A5)") '#','oscn','n','l'
    
    DO i=1, this%amount
       WRITE(unt,"(I6,3I5)") i, this%orbit(i)%oscn, &
                                 this%orbit(i)%nn  , &
                                 this%orbit(i)%ll 
    END DO

  END SUBROUTINE print_nl_basis
  
  SUBROUTINE print_nlNL_basis(unt, this)
    
    INTEGER, INTENT(IN)            :: unt
    TYPE(nlNLbasis_type), INTENT(IN) :: this
    INTEGER :: i
    LOGICAL :: UNT_IS_OPEN
    
    INQUIRE(unt, opened=UNT_IS_OPEN)
    IF (.NOT. UNT_IS_OPEN) RETURN
    WRITE(unt,"(A)") 'NL(rel)NL(cm) BASIS'
    WRITE(unt,"(A10,5A5)") '#','oscn','n','l','ncm','lcm'
    
    DO i=1, this%amount
       WRITE(unt,"(I10,5I5)") i, this%orbit(i)%oscn, &
                                 this%orbit(i)%nn  , &
                                 this%orbit(i)%ll  , &
                                 this%orbit(i)%ncm  , &
                                 this%orbit(i)%lcm
    END DO

  END SUBROUTINE print_nlNL_basis
  
  SUBROUTINE print_nlnlNL_basis(unt, this)
    
    INTEGER, INTENT(IN)            :: unt
    TYPE(nlnlNLbasis_type), INTENT(IN) :: this
    INTEGER :: i
    LOGICAL :: UNT_IS_OPEN
    
    INQUIRE(unt, opened=UNT_IS_OPEN)
    IF (.NOT. UNT_IS_OPEN) RETURN
    WRITE(unt,"(A)") 'nab lab, nc lc, ncm lcm BASIS'
    WRITE(unt,"(A10,10A8)") '#','oscn','oscnab','oscnc','oscncm','nab','2*lab','nc','2*lc','ncm','2*lcm'
    
    DO i=1, this%amount
       WRITE(unt,"(I10,10I8)") i, this%orbit(i)%oscn  , &
                                  this%orbit(i)%oscn12, &
                                  this%orbit(i)%oscn3 , &
                                  this%orbit(i)%oscncm, &
                                  this%orbit(i)%n12   , &
                                  this%orbit(i)%l12   , &
                                  this%orbit(i)%n3    , &
                                  this%orbit(i)%l3    , &
                                  this%orbit(i)%ncm   , &
                                  this%orbit(i)%lcm 
       
    END DO
    
  END SUBROUTINE print_nlnlNL_basis
  
  FUNCTION nlbasis_memory(nlbasis) RESULT(sizeof_nlbasis_type_char)
    
    IMPLICIT NONE
    
    TYPE(nlbasis_type), INTENT(IN) :: nlbasis
    INTEGER, PARAMETER :: sizeof_nlorbit_type = 4*sizeof_int
    INTEGER, PARAMETER :: sizeof_nllimits     = 3*sizeof_limits_type
    INTEGER            :: sizeof_nlbasis_type
    CHARACTER(LEN=100) :: sizeof_nlbasis_type_char

    sizeof_nlbasis_type = nlbasis%amount*sizeof_nlorbit_type + &
                          sizeof_nllimits + sizeof_int ! Bytes

    sizeof_nlbasis_type_char = 'sizeof_nlbasis_type: '//&
         TRIM(ADJUSTL(r2a(REAL(sizeof_nlbasis_type)/2.0**10)))// ' KBytes'

  END FUNCTION nlbasis_memory
  
  FUNCTION nlNLbasis_memory(nlNLbasis) RESULT(sizeof_nlNLbasis_type_char)
    
    IMPLICIT NONE
    
    TYPE(nlNLbasis_type), INTENT(IN) :: nlNLbasis
    INTEGER, PARAMETER :: sizeof_nlNLorbit_type = 6*sizeof_int
    INTEGER, PARAMETER :: sizeof_nllimits     = 3*sizeof_limits_type
    INTEGER            :: sizeof_nlNLbasis_type
    CHARACTER(LEN=100) :: sizeof_nlNLbasis_type_char

    sizeof_nlNLbasis_type = nlNLbasis%amount*sizeof_nlNLorbit_type + &
                          sizeof_nllimits + sizeof_int ! Bytes

    sizeof_nlNLbasis_type_char = 'sizeof_nlNLbasis_type: '//&
         TRIM(ADJUSTL(r2a(REAL(sizeof_nlNLbasis_type)/2.0**10)))// ' KBytes'

  END FUNCTION nlNLbasis_memory
     
END MODULE single_particle_basis

MODULE harmonic_oscillator

  USE constants
  USE gauss_legendre_mesh
    
  IMPLICIT NONE
  
  TYPE, PUBLIC :: howf_values
     REAL(DP) :: integral, norm
     REAL(DP), DIMENSION(:), ALLOCATABLE :: value
  END TYPE howf_values
  
  TYPE, PUBLIC :: howf_evaluated
     INTEGER                                        :: amount
     TYPE(howf_values), DIMENSION(:,:), ALLOCATABLE :: nl_number
     TYPE(gauleg_mesh)                              :: mesh
     REAL(DP)                                       :: b
  END TYPE howf_evaluated
  
  ! use this type for storing 'assembled' wave functions
  TYPE, PUBLIC :: wfunc_type
     !                   x,L  
     REAL(DP), DIMENSION(:,:), ALLOCATABLE :: psi
     TYPE(gauleg_mesh) :: mesh
     REAL(DP)          :: norm
     LOGICAL           :: TENSOR
  END TYPE wfunc_type

  INTERFACE harmonic_oscillator_energy
     MODULE PROCEDURE onebody_harmonic_oscillator_energy, &
                      twobody_harmonic_oscillator_energy, &
                      threebody_harmonic_oscillator_energy
  END INTERFACE

CONTAINS
  
  SUBROUTINE set_howf_mesh_and_b(this, b, Nmesh, min, max)

    IMPLICIT NONE
    TYPE(howf_evaluated), INTENT(INOUT):: this 
    INTEGER, INTENT(IN) :: Nmesh
    REAL(DP), INTENT(IN) :: b, min, max

    this%mesh%info = mesh_info(Nmesh, min, max)
    CALL setup_gauleg_mesh(this%mesh)
    
    this%b = b
    
  END SUBROUTINE set_howf_mesh_and_b
  
  SUBROUTINE allocate_radial_wavefunction(this, nlbasis)
    
    USE single_particle_basis

    IMPLICIT NONE
    TYPE(howf_evaluated), INTENT(INOUT):: this    ! resulting harmonic oscillator wf's
    TYPE(nlbasis_type)  , INTENT(IN)   :: nlbasis ! list of harmonic oscillator quantum numbers |nl>

    IF (ALLOCATED(this%nl_number)) THEN
       WRITE(ERR,"(A)") 'error(allocated_radial_wavefunction): this%nl_number already allocated'
       STOP
    END IF
    
    ALLOCATE(this%nl_number(0:nlbasis%limits%nn%max,0:nlbasis%limits%ll%max))
    
  END SUBROUTINE allocate_radial_wavefunction
  
  SUBROUTINE setup_radial_wavefunction(this, nlbasis)
    
    USE special_functions
    USE single_particle_basis
    
    IMPLICIT NONE

    TYPE(howf_evaluated), INTENT(INOUT):: this     ! resulting harmonic oscillator wf's
    TYPE(nlbasis_type)  , INTENT(IN)   :: nlbasis ! list of harmonic oscillator quantum numbers |nl>
    INTEGER  :: nn, ll, ip, iorbit
    
    ! allocate the arrays for the radial points, for each allowed n/l 

    DO iorbit=1, nlbasis%amount
       
       nn = nlbasis%orbit(iorbit)%nn          
       ll = nlbasis%orbit(iorbit)%ll          
       IF(ALLOCATED(this%nl_number(nn, ll)%value)) CYCLE
       ALLOCATE (this%nl_number(nn, ll)%value( 1:this%mesh%info%amount ))
       
       this%nl_number(nn, ll)%value(:) = 0.0_dp
       
       DO ip = 1, this%mesh%info%amount
          this%nl_number(nn, ll)%value(ip) = &
               rnl_laguerre(nn, ll, this%mesh%pnt(ip)%x, this%b)

       ENDDO ! ip        
       
       ! EVALUATE INTEGRALS AND NORMS OF H.O. FUNCTIONS
       
       this%nl_number(nn, ll)%integral = &
            SUM(this%nl_number(nn, ll)%value(:)*this%mesh%pnt(:)%xxw)
       
       this%nl_number(nn, ll)%norm = &
            SUM(this%nl_number(nn, ll)%value(:)**2*this%mesh%pnt(:)%xxw)
       
       IF (ABS(this%nl_number(nn, ll)%norm-1.0_dp) > 1.d-1) THEN
          WRITE(ERR,"(A,I,A,I,A,F)") 'warning(setup_radial_wavefunction): for n', nn, ' l', ll, ' norm = ', this%nl_number(nn,ll)%norm
       END IF

       
    ENDDO ! confs
    
  END SUBROUTINE setup_radial_wavefunction
  
  SUBROUTINE setup_radial_wavefunction_with_scaled_argument(this, nlbasis, scale)
    
    USE special_functions
    USE single_particle_basis
    
    IMPLICIT NONE

    TYPE(howf_evaluated), INTENT(INOUT):: this     ! resulting harmonic oscillator wf's
    TYPE(nlbasis_type)  , INTENT(IN)   :: nlbasis ! list of harmonic oscillator quantum numbers |nl>
    REAL(DP)            , INTENT(IN)   :: scale   ! Rnl(scale*r)
    INTEGER  :: nn, ll, ip, iorbit
    
    ! allocate the arrays for the radial points, for each allowed n/l 

    DO iorbit=1, nlbasis%amount
       
       nn = nlbasis%orbit(iorbit)%nn          
       ll = nlbasis%orbit(iorbit)%ll          
       IF(ALLOCATED(this%nl_number(nn, ll)%value)) CYCLE
       ALLOCATE (this%nl_number(nn, ll)%value( 1:this%mesh%info%amount ))
       
       this%nl_number(nn, ll)%value(:) = 0.0_dp
       
       DO ip = 1, this%mesh%info%amount
          this%nl_number(nn, ll)%value(ip) = &
               rnl_laguerre(nn, ll, this%mesh%pnt(ip)%x*scale, this%b)
       ENDDO ! ip        
       
       ! EVALUATE INTEGRALS AND NORMS OF H.O. FUNCTIONS
       
       this%nl_number(nn, ll)%integral = &
            SUM(this%nl_number(nn, ll)%value(:)*this%mesh%pnt(:)%xxw)
       
       this%nl_number(nn, ll)%norm = &
            SUM(this%nl_number(nn, ll)%value(:)**2*this%mesh%pnt(:)%xxw)
       
       !IF (ABS(this%nl_number(nn, ll)%norm-1.0_dp) > 1.d-1) THEN
       !   WRITE(ERR,"(A,I,A,I,A,F)") 'warning(setup_radial_wavefunction): for n', nn, ' l', ll, ' norm = ', this%nl_number(nn,ll)%norm
       !END IF

       
    ENDDO ! confs
    
  END SUBROUTINE setup_radial_wavefunction_with_scaled_argument

  SUBROUTINE print_radial_wavefunctions_norm(unt, a ,b, c, sp)
    
    USE single_particle_basis

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unt
    TYPE(howf_evaluated), INTENT(IN) :: a,b,c
    TYPE(nlbasis_type), INTENT(IN) :: sp
    INTEGER :: wf,nn,ll
    LOGICAL :: UNT_IS_OPEN
    
    INQUIRE(unt, opened=UNT_IS_OPEN) 
    IF (.NOT. UNT_IS_OPEN) RETURN
    
    WRITE(unt,"(A)") 'WAVEFUNCTION NORMS'
    WRITE(unt,"(2A6,3A15)") 'n','l','pp','pn','nn'
    
    DO wf=1, sp%amount
       nn=sp%orbit(wf)%nn
       ll=sp%orbit(wf)%ll
       WRITE(unt,"(2I6,3F15.6)") nn,ll,a%nl_number(nn, ll)%norm,b%nl_number(nn, ll)%norm,c%nl_number(nn, ll)%norm 
    END DO
    
  END SUBROUTINE print_radial_wavefunctions_norm

  SUBROUTINE print_radial_wavefunction_norm(unt, a,sp)
    
    USE single_particle_basis

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unt
    TYPE(howf_evaluated), INTENT(IN) :: a
    TYPE(nlbasis_type), INTENT(IN) :: sp
    INTEGER :: wf,nn,ll
    LOGICAL :: UNT_IS_OPEN
    
    INQUIRE(unt, opened=UNT_IS_OPEN) 
    IF (.NOT. UNT_IS_OPEN) RETURN
    
    WRITE(unt,"(A)") 'WAVEFUNCTION NORM'
    WRITE(unt,"(2A6,A15)") 'n','l'
    
    DO wf=1, sp%amount
       nn=sp%orbit(wf)%nn
       ll=sp%orbit(wf)%ll
       WRITE(unt,"(2I6,F15.6)") nn,ll,a%nl_number(nn, ll)%norm
    END DO
    
  END SUBROUTINE print_radial_wavefunction_norm

  SUBROUTINE print_radial_wavefunction(unt, this, nn, ll)
    
    USE constants
    
    INTEGER , INTENT(IN) :: unt, nn, ll
    TYPE(howf_evaluated), INTENT(IN) :: this
    INTEGER :: ii
    
    DO ii=1, this%mesh%info%amount
       WRITE(unt, "(2E16.6)") this%mesh%pnt(ii)%x, this%nl_number(nn, ll)%value(ii)
    END DO
    
  END SUBROUTINE print_radial_wavefunction
  
  FUNCTION radial_wavefunction_memory(howf) RESULT(sizeof_howf_evaluated_char)
    
    IMPLICIT NONE

    TYPE(howf_evaluated), INTENT(IN) :: howf
    INTEGER :: mesh_amount        
    INTEGER :: n_amount           
    INTEGER :: l_amount           
    INTEGER :: sizeof_howf_values 
    INTEGER :: sizeof_gauleg_mesh 
    INTEGER :: sizeof_howf_evaluated
    CHARACTER(LEN=100) :: sizeof_howf_evaluated_char

    mesh_amount        = howf%mesh%info%amount                   
    n_amount           = SIZE(howf%nl_number(:,1))               
    l_amount           = SIZE(howf%nl_number(1,:))               
    sizeof_howf_values = (2+mesh_amount)*sizeof_double           
    sizeof_gauleg_mesh = (2+mesh_amount)*sizeof_double+sizeof_int
    
    sizeof_howf_evaluated = sizeof_int + (n_amount*l_amount)*sizeof_howf_values +&
         sizeof_gauleg_mesh + sizeof_double! Bytes

    sizeof_howf_evaluated_char = 'sizeof_howf_evaluated: '//&
         TRIM(ADJUSTL(r2a(REAL(sizeof_howf_evaluated)/2.0**20)))// ' MBytes'

  END FUNCTION radial_wavefunction_memory
    
  FUNCTION onebody_harmonic_oscillator_energy(n, l) RESULT(res)
    
    IMPLICIT NONE
    INTEGER  :: n, l ! NOT DOUBLED
    REAL(DP) :: res

    res = T_fac*(2*n + l + 1.5_dp)*hbar_omega

  END FUNCTION onebody_harmonic_oscillator_energy

  FUNCTION twobody_harmonic_oscillator_energy(n1, l1, n2, l2) RESULT(res)
    
    IMPLICIT NONE
    INTEGER  :: n1, l1, n2, l2 ! NOT DOUBLED
    REAL(DP) :: res

    res = T_fac*(2*n1 + l1 + 2*n2 + l2 + 3.0_dp)*hbar_omega

  END FUNCTION twobody_harmonic_oscillator_energy

  FUNCTION threebody_harmonic_oscillator_energy(n1, l1, n2, l2, n3, l3) RESULT(res)
    
    IMPLICIT NONE
    INTEGER  :: n1, l1, n2, l2, n3, l3 ! NOT DOUBLED
    REAL(DP) :: res

    res = T_fac*(2*n1 + l1 + 2*n2 + l2 + 2*n3 + l3 + 4.5_dp)*hbar_omega

  END FUNCTION threebody_harmonic_oscillator_energy
    
  FUNCTION harmonic_oscillator_recoil_energy(bn, bl, kn, kl) RESULT(res)
    
    IMPLICIT NONE
    INTEGER :: bn, bl, kn, kl
    REAL(DP) :: res, recoil
    
    recoil = 0.0_dp
    IF (bl == kl) THEN
       IF (bn == kn-1) THEN
          recoil = -SQRT((bn+1.0_dp)*(bn + bl + 1.5_dp))
       END IF
       IF (bn == kn) THEN
          recoil = +2.0_dp*bn + bl + 1.5_dp
       END IF
       IF (bn == kn+1) THEN
          recoil = -SQRT(bn*(bn + bl + 0.5_dp))
       END IF
    END IF
    
    res =  T_fac*hbar_omega*recoil/REAL(mass_nucleus, kind=dp)
    
  END FUNCTION harmonic_oscillator_recoil_energy

  FUNCTION harmonic_oscillator_kinetic_energy(bn, bl, kn, kl) RESULT(res)
    
    IMPLICIT NONE
    INTEGER :: bn, bl, kn, kl
    REAL(DP) :: res, kinetic
    
    kinetic = 0.0_dp
    IF (bl == kl) THEN
       IF (bn == kn-1) THEN
          kinetic = +SQRT((bn+1.0_dp)*(bn + bl + 1.5_dp))
       END IF
       IF (bn == kn) THEN
          kinetic = +2.0_dp*bn + bl + 1.5_dp
       END IF
       IF (bn == kn+1) THEN
          kinetic = +SQRT(bn*(bn + bl + 0.5_dp))
       END IF
    END IF
    
    ! KAI MOD EKIN
    res =  T_fac*0.5_dp*hbar_omega*kinetic

  END FUNCTION harmonic_oscillator_kinetic_energy

  ! needs to be verified
  FUNCTION harmonic_oscillator_potential_energy(bn, bl, kn, kl) RESULT(res)
    
    IMPLICIT NONE
    INTEGER :: bn, bl, kn, kl
    REAL(DP) :: res, potential
    
    potential = 0.0_dp
    IF (bl == kl) THEN
       IF (bn == kn-1) THEN
          potential = -SQRT((bn+1.0_dp)*(bn + bl + 1.5_dp))
       END IF
       IF (bn == kn) THEN
          potential = +2.0_dp*bn + bl + 1.5_dp
       END IF
       IF (bn == kn+1) THEN
          potential = -SQRT(bn*(bn + bl + 0.5_dp))
       END IF
    END IF
    
    res =  0.5_dp*hbar_omega*potential

  END FUNCTION harmonic_oscillator_potential_energy
  
  FUNCTION ho_overlap(bra, ket, mesh) RESULT(res)
    TYPE(howf_values), INTENT(IN) :: bra, ket
    TYPE(gauleg_mesh), INTENT(IN) :: mesh
    REAL(DP) :: res
    
    res = SUM (bra%value(:)* ket%value(:)* mesh%pnt(:)%xxw)
    
  END FUNCTION ho_overlap

  SUBROUTINE determine_HO_mesh(tol, n, Np, b, cutoff) 
    
    USE constants
    USE special_functions
    USE gauss_legendre_mesh

    IMPLICIT NONE
    REAL(DP), INTENT(IN)                         :: tol 
    INTEGER, INTENT(IN)                          :: n
    INTEGER, INTENT(INOUT)                       :: Np
    REAL(DP), INTENT(IN)                         :: b, cutoff
    TYPE(gauleg_mesh)                            :: mesh
    INTEGER                                      :: ir, step_Np
    REAL(DP)                                     :: r, eps
    REAL(DP), ALLOCATABLE, DIMENSION(:)          :: wf
    
    step_Np = 2
    
    Np = 10
    ALLOCATE(wf(1:Np))
    wf(:) = 0.0_dp
    
    ! checking when n = Ncut/2 l = 0 is converged within tol
    eps = 1.0_dp
    DO WHILE(eps>tol)
       
       Np = Np + step_Np
       DEALLOCATE(wf) ; ALLOCATE(wf(1:Np))
       wf(:) = 0.0_dp
       
       mesh%info = mesh_info(Np, zero, cutoff)
       CALL setup_gauleg_mesh(mesh)
       
       DO ir = 1, mesh%info%amount
          r = mesh%pnt(ir)%x
          wf(ir) = rnl_laguerre(n, 1, r, b) ! <- use ell = 1
       END DO ! ip        
       eps = ABS(1.0_dp - SUM(wf(:)**2*mesh%pnt(:)%xxw))
      
       IF (Np>1000) THEN
          WRITE(ERR,*) 'error(determine_HO_mesh): Cant find nof HO mesh points. '
          WRITE(ERR,*) '                          Try increasing radial cutoff'
          STOP
          WRITE(SCR,*) eps, Np
       END IF
    END DO
    WRITE(FRUN,"(A,I4,A)") '                    -  HO wf norm convergence requires ', Np, ' quadrature points'
    
  END SUBROUTINE determine_HO_mesh

    SUBROUTINE deallocate_wfunc(psi)
    
    IMPLICIT NONE
    
    TYPE(wfunc_type), INTENT(INOUT) :: psi
    
    IF (ALLOCATED(psi%psi)) DEALLOCATE(psi%psi)
    IF (ALLOCATED(psi%mesh%pnt)) DEALLOCATE(psi%mesh%pnt)
    psi%norm = 0.0_dp
    psi%TENSOR = .FALSE.
    
  END SUBROUTINE deallocate_wfunc

  SUBROUTINE print_wfunc(unt, psi)
        
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: unt
    TYPE(wfunc_type), INTENT(IN) :: psi
    INTEGER :: i
    
    WRITE(unt,"(A,F)") '# norm = ', psi%norm

    IF (ALLOCATED(psi%psi) ) THEN
       
       IF (SIZE(psi%psi(1,:)) == 2) THEN

          WRITE(unt,"(3A20)") '#x','psi[1]','psi[2]'
          DO i=1, psi%mesh%info%amount
             WRITE(unt,"(3F20.9)") psi%mesh%pnt(i)%x, psi%psi(i,1), psi%psi(i,2)
                                   
          END DO
          
       END IF

       IF (SIZE(psi%psi(1,:)) == 1) THEN

          WRITE(unt,"(2A20)") '#x','psi_x[1]'
          DO i=1, psi%mesh%info%amount
             WRITE(unt,"(2F20.9)") psi%mesh%pnt(i)%x, psi%psi(i,1)
                                   
          END DO
          
       END IF

    END IF

        
  END SUBROUTINE print_wfunc

  FUNCTION compute_wfunc_norm(psi) RESULT(norm)
    
    IMPLICIT NONE

    TYPE(wfunc_type), INTENT(IN) :: psi
    REAL(DP) :: norm
    
    norm = 0.0_dp

    IF (psi%TENSOR) THEN
       
       norm = SUM(psi%mesh%pnt(:)%w*(psi%psi(:,1)**2 + psi%psi(:,2)**2))
       
    END IF

    IF (.NOT. psi%TENSOR) THEN
       
       norm = SUM(psi%mesh%pnt(:)%w*(psi%psi(:,1)**2))
       
    END IF

  END FUNCTION compute_wfunc_norm
  
  SUBROUTINE compute_twobody_matter_radius(psi, rd)
    
    IMPLICIT NONE
    
    TYPE(wfunc_type), INTENT(IN)    :: psi
    REAL(DP)        , INTENT(INOUT) :: rd
    
    rd = 0.0_dp
    
    IF (psi%TENSOR) THEN
       
       rd = 0.5_dp*SQRT(SUM(psi%mesh%pnt(:)%xxw*(psi%psi(:,1)**2 + psi%psi(:,2)**2)))
       
    END IF
    
    IF (.NOT. psi%TENSOR) THEN
       
       rd = 0.5_dp*SQRT(SUM(psi%mesh%pnt(:)%xxw*(psi%psi(:,1)**2)))
       
    END IF
    
  END SUBROUTINE compute_twobody_matter_radius

END MODULE harmonic_oscillator
!
! contain methods to read a given input file
!
MODULE inifile
  IMPLICIT NONE
  PUBLIC
  INTEGER, PARAMETER :: ini_max_name_len = 4096
  INTEGER, PARAMETER :: ini_max_string_len =4096
  LOGICAL :: ini_fail_on_not_found = .FALSE.
  LOGICAL :: ini_echo_read = .FALSE.
  TYPE tnamevalue
     !no known way to make character string pointers..
     CHARACTER(ini_max_name_len)  :: name
     CHARACTER(ini_max_string_len):: value
  END TYPE tnamevalue

  TYPE tnamevalue_pointer
     TYPE(tnamevalue), POINTER :: p
  END TYPE tnamevalue_pointer

  TYPE tnamevaluelist
     INTEGER count
     INTEGER delta
     INTEGER capacity
     TYPE(tnamevalue_pointer), DIMENSION(:), POINTER :: items
  END TYPE tnamevaluelist

  TYPE tinifile
     LOGICAL slashcomments
     TYPE (tnamevaluelist) :: l, readvalues
  END TYPE tinifile

  TYPE(tinifile) :: defini

CONTAINS

  SUBROUTINE tnamevaluelist_init(l)
    TYPE (tnamevaluelist) :: l

    l%count = 0
    l%capacity = 0
    l%delta = 4096
    NULLIFY(l%items)

  END SUBROUTINE tnamevaluelist_init

  SUBROUTINE tnamevaluelist_clear(l)
    TYPE (tnamevaluelist) :: l
    INTEGER i, status

    DO i=l%count,1,-1
       DEALLOCATE (l%items(i)%p, stat = status)
    END DO
    DEALLOCATE (l%items, stat = status)
    CALL tnamevaluelist_init(l)

  END SUBROUTINE tnamevaluelist_clear

  SUBROUTINE tnamevaluelist_valueof(l, aname, avalue)
    TYPE (tnamevaluelist) :: l
    CHARACTER(len=*), INTENT(in) :: aname
    CHARACTER(len=*) :: avalue
    INTEGER i

    DO i=1, l%count
       IF (l%items(i)%p%name == aname) THEN
          avalue = l%items(i)%p%value 
          RETURN
       END IF
    END DO
    avalue = ''

  END SUBROUTINE tnamevaluelist_valueof

  SUBROUTINE tnamevaluelist_add(l, aname, avalue)
    TYPE (tnamevaluelist) :: l
    CHARACTER(len=*), INTENT(in) :: aname, avalue

    IF (l%count == l%capacity) CALL tnamevaluelist_setcapacity(l, l%capacity + l%delta)
    l%count = l%count + 1
    ALLOCATE(l%items(l%count)%p)
    l%items(l%count)%p%name = aname
    l%items(l%count)%p%value = avalue

  END SUBROUTINE tnamevaluelist_add

  SUBROUTINE tnamevaluelist_setcapacity(l, c)
    TYPE (tnamevaluelist) :: l
    INTEGER c
    TYPE(tnamevalue_pointer), DIMENSION(:), POINTER :: tmpitems

    IF (l%count > 0) THEN
       IF (c < l%count) STOP 'tnamevaluelist_setcapacity: smaller than count'
       ALLOCATE(tmpitems(l%count))
       tmpitems = l%items(1:l%count)
       DEALLOCATE(l%items)
       ALLOCATE(l%items(c))
       l%items(1:l%count) = tmpitems
       DEALLOCATE(tmpitems)
    ELSE
       ALLOCATE(l%items(c))
    END IF
    l%capacity = c

  END SUBROUTINE tnamevaluelist_setcapacity

  SUBROUTINE tnamevaluelist_delete(l, i)
    TYPE (tnamevaluelist) :: l
    INTEGER, INTENT(in) :: i

    DEALLOCATE(l%items(i)%p)
    IF (l%count > 1) l%items(i:l%count-1) = l%items(i+1:l%count)
    l%count = l%count -1

  END SUBROUTINE tnamevaluelist_delete

  SUBROUTINE ini_namevalue_add(ini,ainline)
    TYPE(tinifile) :: ini
    CHARACTER (len=*), INTENT(in) :: ainline
    INTEGER eqpos, slashpos, lastpos
    CHARACTER (len=LEN(ainline)) :: aname, s, inline

    inline=TRIM(ADJUSTL(ainline))
    eqpos = SCAN(inline,'=')
    IF (eqpos/=0 .AND. inline(1:1)/='#' .AND. inline(1:7) /= 'comment' ) THEN
       aname = TRIM(inline(1:eqpos-1))
       s = ADJUSTL(inline(eqpos+1:)) 
       IF (ini%slashcomments) THEN
          slashpos=SCAN(s,'/')
          IF (slashpos /= 0) THEN
             s  = s(1:slashpos-1)
          END IF
       END IF
       lastpos=LEN_TRIM(s)
       IF (lastpos>1) THEN
          IF (s(1:1)=='''' .AND. s(lastpos:lastpos)=='''') THEN
             s = s(2:lastpos-1)
          END IF
       END IF
       CALL tnamevaluelist_add(ini%l, aname, s)
    END IF

  END SUBROUTINE ini_namevalue_add


  SUBROUTINE ini_open(filename, unit_id,  error, slash_comments)
    CHARACTER (len=*), INTENT(in) :: filename
    INTEGER, INTENT(in) :: unit_id
    LOGICAL, OPTIONAL, INTENT(out) :: error
    LOGICAL, OPTIONAL, INTENT(in) :: slash_comments
    LOGICAL aerror

    CALL tnamevaluelist_init(defini%l)
    CALL tnamevaluelist_init(defini%readvalues)
    IF (PRESENT(slash_comments)) THEN
       CALL ini_open_file(defini,filename,unit_id,aerror,slash_comments)
    ELSE
       CALL ini_open_file(defini,filename,unit_id,aerror)
    END IF
    IF (PRESENT(error)) THEN
       error = aerror
    ELSE
       IF (aerror) THEN
          WRITE (*,*) 'ini_open: error opening file ' // TRIM(filename)
          STOP
       END IF
    END IF

  END SUBROUTINE ini_open


  SUBROUTINE ini_open_file(ini, filename, unit_id,  error, slash_comments)
    TYPE(tinifile) :: ini
    CHARACTER (len=*), INTENT(in) :: filename
    INTEGER, INTENT(in) :: unit_id
    LOGICAL, INTENT(out) :: error
    LOGICAL, OPTIONAL, INTENT(in) :: slash_comments
    CHARACTER (len=4096) :: inline

    CALL tnamevaluelist_init(ini%l)
    CALL tnamevaluelist_init(ini%readvalues)

    IF (PRESENT(slash_comments)) THEN
       ini%slashcomments = slash_comments
    ELSE
       ini%slashcomments = .FALSE.
    END IF
    OPEN(unit=unit_id,file=filename,form='formatted',status='old', err=500)
    DO 
       READ (unit_id,'(a)',END=400) inline
       IF (inline == 'end') EXIT;
       IF (inline /= '') CALL ini_namevalue_add(ini,inline) 
    END DO
400 CLOSE(unit_id)
    error=.FALSE.
    RETURN
500 error=.TRUE.

  END SUBROUTINE ini_open_file

  SUBROUTINE ini_open_fromlines(ini, lines, numlines, slash_comments)
    TYPE(tinifile) :: ini
    INTEGER, INTENT(in) :: numlines
    CHARACTER (len=*), DIMENSION(numlines), INTENT(in) :: lines
    LOGICAL, INTENT(in) :: slash_comments
    INTEGER i

    CALL tnamevaluelist_init(ini%l)
    ini%slashcomments = slash_comments
    DO i=1,numlines
       CALL ini_namevalue_add(ini,lines(i))
    END DO

  END  SUBROUTINE ini_open_fromlines

  SUBROUTINE ini_close

    CALL ini_close_file(defini)

  END SUBROUTINE ini_close


  SUBROUTINE ini_close_file(ini)
    TYPE(tinifile) :: ini

    CALL tnamevaluelist_clear(ini%l)
    CALL tnamevaluelist_clear(ini%readvalues)

  END  SUBROUTINE ini_close_file


  FUNCTION ini_read_string(key, notfoundfail) RESULT(avalue)
    CHARACTER (len=*), INTENT(in) :: key
    LOGICAL, OPTIONAL, INTENT(in) :: notfoundfail
    CHARACTER(len=ini_max_string_len) :: avalue

    IF (PRESENT(notfoundfail)) THEN
       avalue = ini_read_string_file(defini, key, notfoundfail)
    ELSE
       avalue = ini_read_string_file(defini, key)
    END IF

  END FUNCTION ini_read_string


  FUNCTION ini_read_string_file(ini, key, notfoundfail) RESULT(avalue)
    TYPE(tinifile) :: ini
    CHARACTER (len=*), INTENT(in) :: key
    LOGICAL, OPTIONAL, INTENT(in) :: notfoundfail
    CHARACTER(len=ini_max_string_len) :: avalue

    CALL tnamevaluelist_valueof(ini%l, key, avalue)
    IF (avalue/='') THEN
       CALL  tnamevaluelist_add(ini%readvalues, key, avalue)
       IF (ini_echo_read) WRITE (*,*) TRIM(key)//' = ',TRIM(avalue)
       RETURN
    END IF
    IF (ini_fail_on_not_found) THEN
       WRITE(*,*) 'key not found : '//key
       STOP
    END IF
    IF (PRESENT(notfoundfail)) THEN
       IF (notfoundfail) THEN
          WRITE(*,*) 'key not found : '//key
          STOP
       END IF
    END IF

  END FUNCTION ini_read_string_file


  FUNCTION ini_read_int(key, default)
    INTEGER, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    INTEGER ini_read_int

    IF (PRESENT(default)) THEN
       ini_read_int = ini_read_int_file(defini, key, default)
    ELSE
       ini_read_int = ini_read_int_file(defini, key)
    END IF

  END FUNCTION ini_read_int


  FUNCTION ini_read_int_file(ini, key, default)
    TYPE(tinifile) :: ini
    INTEGER ini_read_int_file
    INTEGER, OPTIONAL, INTENT(in) :: default
    CHARACTER  (len=*), INTENT(in) :: key
    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini, key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_int_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       IF (VERIFY(TRIM(s),'-+0123456789') /= 0) GOTO 10
       READ (s,*, err = 10) ini_read_int_file
    END IF
    RETURN
10  WRITE (*,*) 'error reading integer for key: '//key
    STOP

  END FUNCTION ini_read_int_file

  FUNCTION ini_read_double(key, default)
    DOUBLE PRECISION, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    DOUBLE PRECISION ini_read_double

    IF (PRESENT(default)) THEN
       ini_read_double = ini_read_double_file(defini, key, default)
    ELSE
       ini_read_double = ini_read_double_file(defini, key)
    END IF

  END FUNCTION ini_read_double


  FUNCTION ini_read_double_file(ini,key, default)
    TYPE(tinifile) :: ini
    DOUBLE PRECISION ini_read_double_file 
    DOUBLE PRECISION, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini,key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_double_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       READ (s,*, err=10) ini_read_double_file
    END IF
    RETURN
10  WRITE (*,*) 'error reading double for key: '//key
    STOP

  END FUNCTION ini_read_double_file


  FUNCTION ini_read_real(key, default)
    REAL, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    REAL ini_read_real

    IF (PRESENT(default)) THEN
       ini_read_real = ini_read_real_file(defini, key, default)
    ELSE
       ini_read_real = ini_read_real_file(defini, key)
    END IF

  END FUNCTION ini_read_real

  FUNCTION ini_read_real_file(ini,key, default)
    TYPE(tinifile) :: ini
    REAL ini_read_real_file 
    REAL, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini,key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_real_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       READ (s,*, err=10) ini_read_real_file
    END IF
    RETURN
10  WRITE (*,*) 'error reading double for key: '//key
    STOP

  END FUNCTION ini_read_real_file

  FUNCTION ini_read_logical(key, default)
    LOGICAL, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    LOGICAL ini_read_logical

    IF (PRESENT(default)) THEN
       ini_read_logical = ini_read_logical_file(defini, key, default)
    ELSE
       ini_read_logical = ini_read_logical_file(defini, key)
    END IF

  END FUNCTION ini_read_logical

  FUNCTION ini_read_logical_file(ini, key, default)
    TYPE(tinifile) :: ini
    LOGICAL ini_read_logical_file
    LOGICAL, OPTIONAL, INTENT(in) :: default
    CHARACTER  (len=*), INTENT(in) :: key

    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini,key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_logical_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       IF (VERIFY(TRIM(s),'10tf') /= 0) GOTO 10  
       READ (s,*, err = 10) ini_read_logical_file
    END IF

    RETURN

10  WRITE (*,*) 'error reading logical for key: '//key
    STOP
  END FUNCTION ini_read_logical_file


  SUBROUTINE ini_savereadvalues(afile,unit_id)
    CHARACTER(len=*)  :: afile
    INTEGER, INTENT(in) :: unit_id

    CALL ini_savereadvalues_file(defini, afile, unit_id)

  END SUBROUTINE ini_savereadvalues

  SUBROUTINE ini_savereadvaluestofile(unit_id)
    INTEGER, INTENT(in) :: unit_id
    
    CALL ini_savereadvalues_existing_file(defini, unit_id)

  END SUBROUTINE ini_savereadvaluestofile
  
  SUBROUTINE ini_savereadvalues_file(ini, afile, unit_id)
    TYPE(tinifile) :: ini
    CHARACTER(len=*), INTENT(in) :: afile
    INTEGER, INTENT(in) :: unit_id
    INTEGER i

    OPEN(unit=unit_id,file=afile,form='formatted',status='replace', err=500)

    DO i=1, ini%readvalues%count
       WRITE (unit_id,'(a)') TRIM(ini%readvalues%items(i)%p%name) // ' = ' &
            //TRIM(ini%readvalues%items(i)%p%value)

    END DO

    CLOSE(unit_id)
    RETURN

500 WRITE(*,*) 'ini_savereadvalues_file: error creating '//TRIM(afile)

  END SUBROUTINE ini_savereadvalues_file

  SUBROUTINE ini_savereadvalues_existing_file(ini, unit_id)
    TYPE(tinifile) :: ini
    INTEGER, INTENT(in) :: unit_id
    INTEGER i

    DO i=1, ini%readvalues%count
       WRITE (unit_id,'(a)') TRIM(ini%readvalues%items(i)%p%name) // ' = ' &
            //TRIM(ini%readvalues%items(i)%p%value)

    END DO

    RETURN

  END SUBROUTINE ini_savereadvalues_existing_file

END MODULE inifile
