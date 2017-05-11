subroutine count_num_halo(myid, z, nhalo)
	use param
	use param_halo
	implicit none
	integer::myid
	real(8)::z

  	character(7) :: z_s
  	character(250) :: ifile
	real(8) ::halo_max,halo_min,M,Mass, Meff, mass_gal, M1, z_ns
	integer::nhalo,i,j,k,cnt,nhalo_max, ll, xx, yy, cnt1, fstat, nhalo_ns

	integer, dimension(n_cell, n_cell, n_cell)::halo_trace_ns
	real(8), dimension(n_cell, n_cell, n_cell)::halo_mass_new_ns

	halo_trace=0 !! This tracks the grid with halo..1 mean halo present

!!choose the halo file : consider feedback, large/small source

	
	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)

!	ifile=halolist_path//z_s(1:len_trim(z_s))//"-coarsest_wsubgrid_sources.dat"
	ifile=halolist_path//z_s(1:len_trim(z_s))//"-halos.dat"
	write(*,*) 'halo file, z', ifile, z
	open(unit=13,file=ifile, status='old')
	read(13, *) nhalo
	close(13)

end subroutine count_num_halo


!!This subroutine reads the list of halo at redshift z

!!This track age
!!The halo file should be given in grid positions, sum all haloes >Mcut at every grid points and place.
!!Sould not have many similar sources at the same grid points, otherwise age calculation will be wrong.
!!
subroutine read_halofile_track_age(myid, z)
	use param
	use param_halo
	implicit none
	integer::myid
	real(8)::z

  	character(7) :: z_s
  	character(250) :: ifile
	real(8) ::halo_max,halo_min,M,Mass, Meff, mass_gal, M1, z_ns, M2
	integer::nhalo,i,j,k,cnt,nhalo_max, ll, xx, yy, cnt1, fstat, nhalo_ns

	integer, dimension(n_cell, n_cell, n_cell)::halo_trace_ns
	real(8), dimension(n_cell, n_cell, n_cell)::halo_mass_new_ns
integer, parameter::nmax=100000000
	halo_trace=0 !! This tracks the grid with halo..1 mean halo present

!!choose the halo file : consider feedback, large/small source

	
	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)

        !ifile=halolist_path//z_s(1:len_trim(z_s))//"-coarsest_wsubgrid_sources.dat"
	ifile=halolist_path//z_s(1:len_trim(z_s))//"-halos.dat"
	write(*,*) 'halo file, z', ifile, z
	open(unit=13,file=ifile, status='old')
	read(13, *) nhalo_max

	write(*,*) 'Max number of halo', nhalo_max


	nhalo=0
	halo_max=0.0
	halo_min=1.0d+20
	halo_mass_new=0.0		!!This is the new array of halo mass


!!Read the haloes :

	DO ll=1, nmax!nhalo_max ! halo loop
		read(13,*,end=110)  i,j,k,M1,M, M2
		mass=M*M_grid_s  ! mass of the halo in solar mass



	if(mass >mmin_cut) then
		halo_trace(i,j,k)=1 !To identify the cells with source	
		mass_gal=Mass
			
		nhalo=nhalo+1
		if(mass_gal .lt. halo_min) halo_min=mass_gal
		if(mass_gal .gt. halo_max) halo_max=mass_gal


		halo_mass_new(i,j,k)=halo_mass_new(i,j,k)+real(mass_gal)

	end if
				
	END DO !l
	110 print*,'Halo considered',nhalo


	close(13)
	num_halo=nhalo

	write(*,*) 'mim max mass halo', halo_min, halo_max

	DO i=1, n_cell
	DO j=1, n_cell
	DO k=1, n_cell

	if(halo_trace(i,j,k)==1) then
			halo_age_new(i,j,k)=halo_age_new(i,j,k)+1
	end if
	END DO
	END DO
	END DO


!!Check if halo is moving to different grid or not.. estimate the age of the halo..
	cnt=0
	cnt1=0
	DO k=1, n_cell
	DO j=1, n_cell
	Do i=1, n_cell
		xx=halo_age_old(i,j,k)
		yy=halo_age_new(i,j,k)
		if((xx .ne. 0) .and. (xx==yy) ) then
			call halo_move(i,j,k, fstat)
			cnt1=cnt1+1
			if(fstat==0) cnt=cnt+1
		end if
	END DO
	END DO
	END DO

!!This loop can not be merge with previous loops, as previous array correct for halo movement.
!!This one 

	DO k=1, n_cell
	DO j=1, n_cell
	Do i=1, n_cell
		if(halo_trace(i,j,k)==1 .and. halo_age_new(i,j,k)>1) then
			if(halo_mass_old(i,j,k)==0.0) halo_mass_old(i,j,k)=halo_mass_new(i,j,k)
		end if
	END DO
	END DO
	END DO

	write(*,*) 'Number of halo moved, rejected, frac of rejection', cnt1, cnt, real(cnt)/real(cnt1)


!!Now we will calculate the effective mass of the halo over time, this is to account the mass evlution of the halo which is in general exponential.
!!However, we will not taken into account the dennisty contrast evolution around the source over the time and will work only with the current density distribution at that redshift.

	dat_overlap=0.0	!!This array will carry the mass, X(:), etc.. 
	cnt=0


	DO k=1,n_cell
	DO j=1,n_Cell
	DO i=1,n_cell
		xx=halo_trace(i,j,k)
		if(xx ==1 ) then
			cnt=cnt+1
			Meff=(halo_mass_new(i,j,k)+halo_mass_old(i,j,k)*(real(halo_age_new(i,j,k))-1.))/real(halo_age_new(i,j,k))

			halo_mass_new(i,j,k)=Meff
			dat_overlap(1,cnt)=Meff 
			dat_overlap(2,cnt)=real(i)+0.5	
			dat_overlap(3,cnt)=real(j)+0.5	
			dat_overlap(4,cnt)=real(k)+0.5
			dat_overlap(5,cnt)=real(age_checkpoint(halo_age_new(i,j,k)))
		end if
	END DO
	END DO
	END DO

	WRITE(*,*) 'max halo age',maxval(dat_overlap(5,1:cnt))
	print*,'Read_halofile exit. nhalo',cnt


!!Save the mass and age arrays for the next time step.

	halo_age_old=halo_age_new
	halo_mass_old=halo_mass_new

end subroutine read_halofile_track_age

!!In case halo move to new grid
!!
SUBROUTINE halo_move(i,j,k, fstat)
	use param
	use param_halo
	implicit none
	integer, intent(out)::fstat
	integer,intent(in)::i,j,k
	integer::r1,r2,r3,n
	integer, dimension(3)::r

!!Posibility of halo moving to new grid
	fstat=0
	DO n=2,100
		r(1)=i + cube(n,1)
		r(2)=j + cube(n,2)
		r(3)=k + cube(n,3)
		call check_boundary(r, n_cell, r1, r2, r3)

		if((halo_age_new(r1,r2,r3)==1) .and. (halo_trace(r1,r2,r3)==1) ) then !! possibility of halo move to new cell
			halo_age_new(r1,r2,r3)=halo_age_old(i,j,k) + 1	!!This halo will not be considered for next halo move, as the age have been changed to more than 1
			halo_age_new(i,j,k)=0
			halo_age_old(i,j,k)=0


			halo_mass_old(r1, r2, r3)=halo_mass_old(i,j,k)
			halo_mass_old(i,j,k)=0.0

			fstat=1
			exit
		end if
	END DO

!!Possibility of halo merge to existing halo..

if(fstat==0) then
	DO n=2,100
		r(1)=i + cube(n,1)
		r(2)=j + cube(n,2)
		r(3)=k + cube(n,3)
		call check_boundary(r, n_cell, r1, r2, r3)

		if( (halo_trace(r1,r2,r3)==1) ) then !! possibility of halo merge
		if(halo_age_new(r1,r2,r3) .lt. halo_age_new(i,j,k)) then
			halo_age_new(r1,r2,r3)=halo_age_new(i,j,k) + 1
		end if
			halo_age_new(i,j,k)=0
			halo_age_old(i,j,k)=0
	

			halo_mass_old(i,j,k)=0.0
			fstat=1
			exit
		end if
	END DO

end if

!!This part is only to reduce the error..consider if the halo mass becme less than previous step

if(fstat==0) then
	DO n=2,100
		r(1)=i + cube(n,1)
		r(2)=j + cube(n,2)
		r(3)=k + cube(n,3)
		call check_boundary(r, n_cell, r1, r2, r3)

		if( (halo_trace(r1,r2,r3)==1) ) then !! possibility of halo mass reduced..
		if(halo_age_new(r1,r2,r3) .lt. halo_age_new(i,j,k)) then
			halo_age_new(r1,r2,r3)=halo_age_new(i,j,k) + 1

			halo_age_new(i,j,k)=0
			halo_age_old(i,j,k)=0
	

			halo_mass_old(i,j,k)=0.0
			fstat=1
			exit
		end if
		end if
	END DO

end if

!!This part is if the halo vanishes.. then the ionized region will be there and reduce as recombination proceed
!!No age increase, we add an artificial source to keep the ionizeed bubble.
if(fstat==0) then
!	halo_age_new(i,j,k)=halo_age_old(i,j,k)
!	halo_mass_new(i,j,k)=halo_mass_old(i,j,k)	!!removing this fake introduction of halo.. halo may fade away
!	halo_trace(i,j,k)=1
	fstat=1
end if

if(fstat==0) then
 write(*,*) '!!!!WARNING!!!.. halo age is not calculated'
!write(*,*) halo_mass_old(i,j,k), halo_mass_new(i,j,k), maxval(halo_mass_new(i-5:i+5,j-5:j+5,k-5:k+5)), minval(halo_mass_new(i-3:i+3,j-3:j+3,k-3:k+3)), i,j,k
stop
end if


END SUBROUTINE halo_move

subroutine z_relevant(z,check_min)
     use param
     implicit none
     integer::counter
     real(8),intent(in)::z,check_min
     real(8)::z_max

     z_max = (z+1)/0.75 - 1
     write(*,*) "Redshift range: ", z_max, " to ",z
     do counter=1,num_checkpoints
        if (z_checkpoint(counter)>=z .and. z_checkpoint(counter)<=z_max) then
            comm_redshifts(counter)=z_checkpoint(counter)
        else
            comm_redshifts(counter)=0.0
        endif
     enddo

end subroutine

