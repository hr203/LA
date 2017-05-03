
!!This will optimize the distribution of jobs into many processes.
!
SUBROUTINE para_range(njob, nprocs, irank, myoffset, myonset)
	use param
	implicit none
	integer:: njob, nprocs, irank, myoffset, myonset,chunksize, i, nj, m

	if(njob<nprocs) then
		if(irank<njob) then
			myoffset=irank+1
			myonset=irank+1
		else
			myoffset=0
			myonset=0
		end if
	else
		myoffset=1
		nj=njob
		DO i=1, nprocs
			m=int(nj/(nprocs-i+1))
			myonset=myoffset+m-1
			if(irank==(i-1)) then
				exit
			else
			 	myoffset=myonset+1
				nj=nj-m
			end if
		END DO
	end if

	write(*,*) 'JOB DISTRIBUTION, RANK, offset, onset', irank, myoffset, myonset

END SUBROUTINE para_range


!This write the three dimensional 4bit array to file in binary form
!
SUBROUTINE write_binary_4bit(arr, ifile, n1, n2, n3)
	IMPLICIT NONE
 	character(*)::ifile
	integer, intent(in)::n1, n2, n3
	real(4), dimension(n1, n2, n3)::arr
	integer::fstat

	open(unit=17,file=trim(ifile),status='replace',iostat=fstat,access='stream')
	write(17) arr
	close(17)

	write(*,*) 'New file created :', ifile

END SUBROUTINE write_binary_4bit

!This write 2D slice.
!!
SUBROUTINE write_slice(arr, ifile, n1, n2, n3, n)
IMPLICIT NONE
 	character(*)::ifile
	integer, intent(in)::n1, n2, n3, n
	real(4), dimension(n1, n2, n3)::arr
	integer::fstat,i,j

	open(unit=17,file=trim(ifile),status='replace',iostat=fstat,form='formatted')

	DO i=1, n1
	DO j=1, n2
	write(17, '(17f20.4)') real(i), real(j), arr(i,j,n)
	END DO
	END DO

	 close(17)


END SUBROUTINE write_slice



!!This determine time between two redshift z1, z2
!!z1, z2 are redshfit and delt is time gap in Myr
!!
SUBROUTINE timegap(z1,z2, delt)
	use adaptint
	use param
	IMPLICIT NONE
	real(8)::z1,z2, delt

	double precision::RELACC,ABSACC,acc, res
	integer::MAXRUL,IPARM,ifail,N
	double precision,dimension(390)::alpha
	interface
	FUNCTION fn_timegap(z)
	IMPLICIT NONE
		real(8),intent(in)::z
	real(8)::fn_timegap

	END FUNCTION fn_timegap
	end interface
	relacc=1.0E-9
	absacc=0.0d0
	maxrul=9
	iparm=0

	call D01ARF(z1,z2,fn_timegap,RELACC,ABSACC,MAXRUL,IPARM,ACC,res,N,&
		     &ALPHA,IFAIL)
	delt=res/megayear

END SUBROUTINE timegap

FUNCTION fn_timegap(z)
	IMPLICIT NONE
	real(8),intent(in)::z
	real(8)::fn_timegap
	interface
	function hubble_constant(z) 
	use param
	implicit none
	real(8),intent(in)::z
	real(8)::hubble_constant
	end function hubble_constant
	end interface

	fn_timegap=1.d0/(1.d0+z)/hubble_constant(z) 


END FUNCTION fn_timegap

!!! return hubble costant at redshift z in sec-1
function hubble_constant(z) 
	use param
	implicit none
	real(8),intent(in)::z
	real(8)::hubble_constant
	hubble_constant=Ho*SQRT(omega_k*(1.d0+z)**2+omega_m*(1.d0+z)**3+&
         &omega_r*(1.d0+z)**4+omega_l)
end function hubble_constant

