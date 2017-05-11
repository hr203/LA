!!!!read the halofinds redshifts for the simulation
!
subroutine read_redshift
use param
	implicit none
	integer::fstat,i
	real(4)::arg


	open(11,file=checkpoints,status='old',iostat=fstat)
	if(fstat /= 0) then
		write(*,*) 'File missing!!ERROR!!', checkpoints
		stop
	end if

	DO num_checkpoints=1,max_input
	    read(unit=11,err=51,end=41,fmt='(f20.10)') arg
		z_checkpoint(num_checkpoints)=dble(arg)
	END DO
	41  num_checkpoints=num_checkpoints-1
	51  close(11)


!	write(*,*) 'halofinds to recompose:'
	!  do i=1,num_checkpoints
	!    write(*,'(f5.1)') z_checkpoint(i)
	!  enddo



end subroutine read_redshift

!!calculate the time between the simulation snapshots..
!!
subroutine cal_age_timestep
	use param
	implicit none
	real(8)::z1, z2, delt
	integer::i

	z2=z_checkpoint(1)

	DO i=1,num_checkpoints-1
	z1=z_checkpoint(i+1)
	call timegap(z1,z2,delt)
	age_checkpoint(i)=delt
	END DO
	age_checkpoint(num_checkpoints)=age_checkpoint(num_checkpoints-1) + 10.0

!	write(*,*) 'Age to recompose:'
!	  do i=1,num_checkpoints
!	    write(*,'(f5.1)') age_checkpoint(i)
!	  enddo

end subroutine cal_age_timestep


!! 

SUBROUTINE readcube
	use param
	implicit none
  	character(180) ::ofile5,ifile1
	integer::q,q1,count,i,j,k,fs,n2p,count1,num, fstat
	real(8)::rad,radp
	integer::spdim1
	integer,parameter::max_in=10000


	ifile1=cellpath//'cells432.dat'
	open(unit=31,file=ifile1, status='old', form='formatted', iostat=fstat)
	if(fstat /= 0) then
		write(*,*) 'File missing!!ERROR!', ifile1
		stop
	end if
	count=0
	count1=0

	DO i=1,max_in
		read(31,*,end=100) n2p,radp
		count=count + n2p
		count1=count1 + 1
	END DO
	100  print*, 'max num cells',count, count1


	spdim1=count1
	spdim=count
	close(31)


	count=0
	ifile1=cellpath//'cube432.dat'
	open(unit=31,file=ifile1, access='stream', status='old', iostat=fstat)
	if(fstat /= 0) then
		write(*,*) 'File missing!!ERROR!', ifile1
		stop
	end if
	read(31) cube
	

					
	close(31)
	write(*,*) 'done cube reading, maxval cube',maxval(cube), cube(1,:),  cube(2,:)
	write(*,*) 'spdim=',spdim

END SUBROUTINE readcube

!!To check periodicity in the box
!!
subroutine check_boundary(x, n, i, j, k)
	implicit none

	integer, dimension(3), intent(in)::x
	integer, intent(in)::n
	integer, intent(out)::i,j,k

	i=x(1)
	j=x(2)
	k=x(3)

	if(i<1) i=i+n
	if(j<1) j=j+n
	if(k<1) k=k+n

	if(i>n) i=i-n
	if(j>n) j=j-n
	if(k>n) k=k-n

end subroutine check_boundary


