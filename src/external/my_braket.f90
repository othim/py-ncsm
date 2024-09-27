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

