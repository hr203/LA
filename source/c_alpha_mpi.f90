!! make -j8 -f Make_c_alpha
!! mpirun -np 4 ./c_alpha

!!This program calculate the x_alpha maps for a given source distribution
!!Assumption 1/r^2 fall of the Lyalpha photns number density
!!
program mpi_sim
	use param
	use param_halo
	implicit none 
     	include 'mpif.h'
      	integer  ierr, offset, onset, i, j, k, tag1, &
     	&         tag2, tag3, tag4, source, check,check2,count,check_min,counter
      	real(8) ::   mysum, mysum1, total_sum,total_sum1,af
!        real(8),dimension(:),allocatable::comm_redshifts2
	real(4) :: time
      	integer  status(MPI_STATUS_SIZE),fstat
  	character(180) :: ifile,ofile
  	character(7) :: z_s
	logical,parameter::restart=.false.		!! To restart at some particular slice
	INTEGER,PARAMETER::restart_id=15
	real(4), dimension(n_cell, n_cell, n_cell):: dummy_arr


!C ***** Initializations *****
      	call MPI_INIT(ierr)
      	call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
      	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

	tag1=1
	tag2=2
	tag3=3
	tag4=4


!C***** Master task only ****** read
      if (rank .eq. 0) then
		call readcube
		call read_redshift
		call cal_age_timestep

		halo_age_new=0
		halo_mass_old=0.0
		write(*,*) '# of cell =', n_cell
		write(*,*) '# of checkpoints = ',num_checkpoints
	
		open(unit=16, file='x_alpha.dat')


		call MPI_BCAST(num_checkpoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(spdim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(cube,max_cell_read*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	else
		call MPI_BCAST(num_checkpoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(spdim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(cube,max_cell_read*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	end if

	if(restart) then
		write(*,*) 'Program is restarted : redshift',z_checkpoint(restart_id)
		check_min=restart_id
	else
		check_min=1
	end if

DO check=check_min,num_checkpoints
        if (num_checkpoints>max_input) write(*,*) "WARNING comm_redshifts not large enough"
	comm_redshift=z_checkpoint(check)
        call z_relevant(comm_redshift,comm_redshifts,check_min)
	IF(rank ==0) THEN
		write(*,*) ''
		write(*,*) 'Processing redshift = ',comm_redshift
                write(*,*) "Relevant redshifts here are:",pack(comm_redshifts,comm_redshifts/=0.0)

		call count_num_halo(rank, comm_redshift, num_halo)
		call MPI_BCAST(comm_redshift,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(num_halo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(num_halo>0) then
		allocate(dat_overlap(5,num_halo))

		call read_halofile_track_age(rank, comm_redshift)

		call cpu_time(time)
		write(*,*) 'Halofile reading : done : time ',time





		DO i=1,numtasks-1
			call para_range(num_halo, numtasks, i, offset, onset)
		     	call MPI_SEND(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, ierr) 
		     	call MPI_SEND(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, ierr)
	if((offset >0) .and. (onset>=offset)) then
		     	call MPI_SEND(dat_overlap(1,offset),5*(onset-offset+1), MPI_REAL, i, tag4, MPI_COMM_WORLD, ierr) 
		     	call MPI_SEND(matrix_alpha, n_cell*n_cell*n_cell, MPI_REAL, i, tag3, MPI_COMM_WORLD, ierr) 
	end if
		end do 


		call para_range(num_halo, numtasks, rank, offset, onset)
		call visual_calpha(offset,onset,rank)

		dummy_arr=matrix_alpha
		call cpu_time(time)
		write(*,*) 'RANK0 : done :time ', time

		do i=1, numtasks-1
		     	call MPI_RECV(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, status, ierr)
		     	call MPI_RECV(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, status, ierr)
	if((offset >0) .and. (onset>=offset)) then
		     	call MPI_RECV(dat_overlap(1,offset),5*(onset-offset+1), MPI_REAL, i, tag4, MPI_COMM_WORLD, status, ierr)
		     	call MPI_RECV(matrix_alpha, n_cell*n_cell*n_cell, MPI_REAL, i, tag3, MPI_COMM_WORLD, status, ierr)
			dummy_arr=dummy_arr+matrix_alpha
	end if
			write(*,*) 'Recv,rank=',i
		end do 



	
		matrix_alpha=dummy_arr
		write(*,*) 'max,min alpha',maxval(matrix_alpha(:,:,:)),minval(matrix_alpha)
		count=0
		DO i=1,n_cell
		DO j=1,n_cell
		DO k=1,n_cell
			if(matrix_alpha(i,j,k)<10.0) then
			count=count+1
			end if
		END DO
		END DO
		END DO

		af=1.0-dble(count)/dble(n_cell)**3.0
		write(*,*) 'coupled frac=',af

		write(z_s,'(f7.3)') comm_redshift
		z_s=adjustl(z_s)
		ofile=output_path//trim(z_s)//"matrix_alphastar.dat"
		call write_binary_4bit(matrix_alpha, ofile, n_cell, n_cell, n_cell)

		!ofile=output_path//trim(z_s)//"matrix_alphastar_slice.dat"
		!call write_slice(matrix_alpha, ofile, n_cell, n_cell, n_cell, n_cell/2)

		call cpu_time(time)
		write(*,*) 'Finish : time, z', time, comm_redshift

		write(16,*) comm_redshift, af
		deallocate(dat_overlap)
		if(af> 0.99) then
			write(*,*) 'coupling completed'
			close(16)
			call mpi_abort(mpi_comm_world, ierr)
		end if
end if
	ELSE

		call MPI_BCAST(comm_redshift,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(num_halo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(num_halo>0) then
		allocate(dat_overlap(5,num_halo))
end if

	END IF 

	if(rank .ne. 0) then
		call MPI_RECV(offset, 1, MPI_INTEGER, 0, tag1, MPI_COMM_WORLD, status, ierr)
		call MPI_RECV(onset, 1, MPI_INTEGER, 0, tag2, MPI_COMM_WORLD, status, ierr)
	if((offset >0) .and. (onset>=offset)) then
		call MPI_RECV(dat_overlap(1,offset),5*(onset-offset+1), MPI_REAL, 0, tag4, MPI_COMM_WORLD, status, ierr)
		call MPI_RECV(matrix_alpha, n_cell*n_cell*n_cell, MPI_REAL, 0, tag3, MPI_COMM_WORLD, status, ierr)

         
		call visual_calpha(offset,onset,rank)
	end if


		call MPI_SEND(offset, 1, MPI_INTEGER, 0, tag1, MPI_COMM_WORLD, ierr) 
		call MPI_SEND(onset, 1, MPI_INTEGER, 0, tag2, MPI_COMM_WORLD, ierr)
	if((offset >0) .and. (onset>=offset)) then
		call MPI_SEND(dat_overlap(1,offset),5*(onset-offset+1), MPI_REAL, 0, tag4, MPI_COMM_WORLD, ierr)
		call MPI_SEND(matrix_alpha, n_cell*n_cell*n_cell, MPI_REAL, 0, tag3, MPI_COMM_WORLD, ierr)

	end if
		deallocate(dat_overlap)
	end if
!end loop here
end do

      call MPI_FINALIZE(ierr)


end program


