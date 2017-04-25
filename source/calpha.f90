!!This include the SED for calculating the number of Lyamn alpha photons
!!SED is divided into two files. "lamda_mod.dat" contains all the lambda bins of the SED. Must contain the wavelengths below Lyman-alpha.
!! "lum_mod.dat" contains the Luminosity.
!! 'lamda' in angstrom, 'lum' in ev/s/lamda.. lamda should be in increasing order
!! One needs to replace this with prefered SED.. keeping the requirements ..
!!
SUBROUTINE read_SED(lamda, lum) 
	use param
	IMPLICIT NONE
	real(8),dimension(n_sed), intent(out)::lamda,lum
  	character(180) :: ifile,ofile
	integer::i
	real(8)::xz

	ifile=pspath//"lamda_mod.dat"
	open(unit=21,file=ifile)
	DO i=1,n_sed
		READ(21,*) xz
		lamda(i)=xz
	END DO
	close(21)

	ifile=pspath//"lum_mod.dat"
	open(unit=21,file=ifile)
	DO i=1,n_sed
		READ(21,*) xz
		lum(i)=xz
	END DO
	close(21)

END SUBROUTINE read_SED

!! This subroutine calculates the Lyman-alpha coupling coefficient.
!! This takes the halos subset (myoffset,myonset) from rank 'myid'.
! 
subroutine visual_calpha(myoffset,myonset,myid)
	use param
	use func
	IMPLICIT NONE
	INTEGER ::myoffset,myonset,myid
	real(8),dimension(n_sed)::lamda,lum
	real(8),dimension(n_cell)::arr
	integer,dimension(3)::x1,r
	real(8)::xz,sim_unit,max_dist,res,vol1,r_dist
	integer::i,j,k,no1,no2,pq,num,q, r1, r2, r3
  	character(180) :: ifile,ofile


	if((myoffset >0) .and. (myonset>=myoffset)) then

	call zeffectlyman(comm_redshift,arr)
	write(*,*) 'visual_calpha called.. nalpha=',maxval(arr)
write(*,*) 'rank, off on set', myid, myoffset, myonset
write(*,*) MAXVAL(DAT_OVERLAP(1,myoffset:myonset))
	sim_unit=box/n_cell/hlittle/(1.0+comm_redshift) !mpc physical

	matrix_alpha=0.0

	DO i=myoffset,myonset  !halo
		x1(1:3)=int(dat_overlap(2:4,i)) + 1

!write(*,*) i
		no1=1
		no2=0
		max_dist=dble(dat_overlap(5,i))*Megayear*c_light/Megaparsec/sim_unit	!!Distance in sim unit travel within the age

		if(max_dist > dble(n_cell)) max_dist = dble(n_cell)


		DO j=1,int(max_dist)
			r_dist=dble(j)
			res=arr(j)*dble(dat_overlap(1,i))
			if(res < 0.01) exit
			vol1=4.0/3.0*pi*r_dist**3.0
			no2=int(vol1)
			if(no2>spdim) exit

			IF(no2 >= no1) then
				DO q=no1,no2
					r(:)=x1(:)+cube(q,:)
					call check_boundary(r, n_cell, r1, r2, r3)
					matrix_alpha(r1,r2,r3)=matrix_alpha(r1,r2,r3) + real(res)
				END DO
				no1=no2+1	
			END IF
		END DO
	END DO !halo
	END IF

	write(*,*) 'rank max,min alpha',myid, maxval(matrix_alpha(:,:,:)),minval(matrix_alpha)


end subroutine visual_calpha

!! This subruitne will call the SED per solar mass and calculate the number of Lyman alpha photon 
!! This uses 1/r/r relation to calculate x_alpha
!!
SUBROUTINE zeffectlyman(z,alpha)
	use param
	use func
	IMPLICIT NONE
	real(8),intent(in)::z
	real(8),dimension(n_cell),intent(out)::alpha
	real(8),dimension(n_sed)::lamda,lum
	real(8)::r,zr,lr,l1,l0,e1,e0,en,nu_al,del_l,rp
	integer::i,j,k
	real(8)::xz,res,del_nu, nlya_emission

	call read_SED(lamda, lum)

	del_nu=c_light/angstrom_cm/lamda_alpha/lamda_alpha !!frequency difference for del lamda = 1 A
	nu_al=E_lyalpha/hplanck_ev !Lyman alpha freq in Hz


	DO i=1,n_cell	!!Do we need to increase the size.. photon escaping?? One can try using a increased value
		r=dble(i)*box/hlittle/n_cell !comoving distace at grid points in Mpc
		zr=rtoz(r,z)
		lr=lamda_alpha*(1.0+zr)/(1.0+z)	!!wavelength in the continuum SEd which redshifted to Lyman alpha at distance r
		if(lr<lamda_hi) then
			alpha(i) = 0.0	!!above EHI, photons will be absorbed for ionization

		else 
			DO j=1,n_sed
				if(lamda(j)>=lr) then
					l1=lamda(j)
					l0=lamda(j-1)
					e1=lum(j)
					e0=lum(j-1)
					exit
				end if
			END DO

			en=e0+ (lr-l0)*(e1-e0)/(l1-l0) !linear interpolation !energy per A per sec in erg
			del_l=(1.0+zr)/(1.0+z)! wave length differ for which  at r the wavelength differ will be 1 A!c_light/nu_al/nu_al*(1+z)/(1+zr)*1e8 !del lamda in A for unit frequency at r
			rp=r/(1.0+zr)*megaparsec
			alpha(i)=en*del_l/4.0/pi/rp/rp/del_nu*erg_ev/E_lyalpha
			alpha(i)=1.66e11/(1.0+zr)*alpha(i)*f_esc_lya  ! ly coeffi for unit mass

		end if
	END DO


END SUBROUTINE zeffectlyman

