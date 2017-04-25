module func  
!use nrtype
use adaptint
use param
contains

function comoving_distance(z1,z2)
	use adaptint
	implicit none

	double precision,intent(in)::z1,z2
	double precision::comoving_distance,tH,dH,res,res2
	double precision::RELACC,ABSACC,acc
	integer::MAXRUL,IPARM,ifail,N
	double precision,dimension(390)::alpha
	relacc=1.0E-9
	absacc=0.0d0
	maxrul=9
	iparm=0
	call D01ARF(z1,z2,drdz,RELACC,ABSACC,MAXRUL,IPARM,ACC,res2,N,&
		     &ALPHA,IFAIL)
	comoving_distance=res2/megaparsec

end function comoving_distance

function drdz(z) 
	implicit none
	double precision,intent(in)::z
	double precision::drdz,c_light
  	c_light=3.d+10
	drdz=- c_light/hubble_constant_dp(z) 
end function drdz



!
!!function rtoz returns the Z value corresponding to the comoving distance, R, in MPC, traveled by a photon emitted at zsource
FUNCTION rtoz (R, zsource)
IMPLICIT NONE
double precision,INTENT(IN)::R, zsource
double precision::rtoz, Rguess, zmax, zmin, zguess
REAL(8):: epsilon = 1e-6
REAL(8),PARAMETER::TINY= 1e-30
integer::i

  if (R < TINY) then
     rtoz= zsource
  end if

  Rguess = 0.0
  zmax = zsource
  zmin = 0.0

  DO i=1,1000
    zguess = (zmax+zmin)/2.0
    Rguess = comoving_distance(zsource,zguess)

    if ( abs((Rguess-R)/R) > epsilon) then
      if (Rguess > R) then
	 zmin = zguess
      else 
	 zmax = zguess
      end if
    else 
	exit
    end if
  END DO
  

  rtoz= zguess

END FUNCTION rtoz



function hubble_constant_dp(z) 
implicit none
double precision,intent(in)::z
double precision::hubble_constant_dp
hubble_constant_dp=Ho*SQRT(omega_k*(1.0+z)**2+omega_m*(1.0+z)**3+&
         &omega_r*(1.0+z)**4+omega_l)
end function hubble_constant_dp





end module func
