MODULE param
	!USE const


!!Constants..
	real(8), parameter:: hplanck_erg = 6.6260695729E-27
	real(8), parameter:: hplanck_ev = 4.13566751691E-15
	real(8), parameter:: erg_ev = 6.2415E+11
	real(8), parameter:: E_HI = 13.60	
	real(8), parameter:: E_lyalpha = 10.20		
	real(8), parameter:: E_HeI = 24.590		
	real(8), parameter:: E_HeII = 54.420	
	real(8), parameter:: E_max = 10000.00	
	real(8), parameter:: Megayear = 3.1536E+13
	real(8), parameter:: Megaparsec = 3.08568025E+24                 !megaparsec in cm
	real(8), parameter:: kiloparsec = 3.08568025E+21			!kiloparsec in cm
	real(8), parameter:: angstrom_cm = 1.d-8
	real(8), parameter:: pi = 3.14159265359
	real(8), parameter:: C_light = 3.0e+10 
	real(8), parameter:: lamda_alpha = 1216.0	!!In A
	real(8), parameter::lamda_hi = 912.0	!!In A
  	real(8), parameter :: grav_const = 6.672d-8 !cm^3/g/s^2
  	real(8), parameter :: Gkgs = grav_const/1.0d3
  	real(8), parameter :: M_sun = 1.989d33
	real(8), parameter :: m_proton = 1.6726d-24	!!gram
	real(8),parameter::Kboltz_ev=8.62911325d-5 		

!! COSMOLOGICAL PARAMETERS .. suman da paper : WMAP5

	real(8), parameter :: omega_l = 0.73
	real(8), parameter :: omega_m = 1.0 - omega_l
	real(8), parameter :: omega_k = 1.0 - omega_l - omega_m
	real(8), parameter :: omega_r = 5.0E-5
	real(8), parameter :: omega_b = 0.044	!!not 0.046
	real(8), parameter :: omega_ch = 0.7
	real(8), parameter :: bias = 1
	real(8), parameter :: power_index = 2.0 
	real(8), parameter :: hlittle = 0.7
	real(8), parameter :: Ho = (hlittle*3.2407e-18) ! s^-1 at z=0 
	real(8), parameter :: hydrogen_frac = 0.752
	real(8), parameter :: f_he = 1.0 - hydrogen_frac
	real(8), parameter :: z_dec = 151.0
	real(8), parameter :: tcmb = 2.73

!!!SIMULATION PARAMETERS


	real(8), parameter :: box = 244.0 !!box length in cMpc/h unit
	integer, parameter :: nc = 8000
	integer, parameter :: n_cell = 250	!!grid number in simulation box

!!RELATED TO CUBEP3M SIMULATION


	real(8), parameter :: h0kgs = hlittle*3.2407e-18 
	real(8), parameter :: rho_ckgs = (3.0*h0kgs**2)/(8.0*pi*Gkgs)
	REAL(8), parameter::den_unit0=omega_m*rho_ckgs*(dble(n_cell)/dble(nc))**3.0	!Density 


  	real(8), parameter :: Vol_m3 = (box/hlittle*(Megaparsec/100.0))**3.0
  	real(8), parameter :: M_grid = (omega_m*rho_ckgs*Vol_m3)/dble(nc**3.0)
  	real(8), parameter :: M_grid_s = M_grid/(M_sun/1.0d3) !mass unit in M_solar, particle mass=8*M_grid_s


!!*******************************************************************YOU HAVE TO DEFINE ***********************************


	integer, parameter :: max_input=300	!!Number of redshift bins
	real(8), dimension(max_input) :: z_checkpoint,age_checkpoint,comm_redshifts

	!integer,parameter :: max_halos=80000000
	!real(4), dimension(12, max_halos):: dat_overlap	!!This carry all the information about the haloes
	real(4), dimension(:, :), allocatable:: dat_overlap	!!This carry all the information about the haloes

	integer, parameter:: max_cell_read=42210291	!!This is the dimension of array cube, which helps to create the sphere around the sources.
	integer, dimension(max_cell_read, 3)::cube	!!This array will be used for creating spherical bubbles


!! SIMULATION PARAMETERS

 	character(*), parameter :: comm_path    = '../' 
 	character(*), parameter :: halolist_path = comm_path!'halo_files/'	
 	character(*), parameter :: output_path    = comm_path//'output/' 
 	character(*), parameter :: pspath = comm_path//'input/'	!!put the SED files here
 	character(*), parameter :: checkpoints =comm_path//'red.dat' !!This file contains the redshifts to be considered for the simulation.. follow decending order
 	character(*), parameter :: cellpath=comm_path//'cellpath/'	!!This contains the cubes to generate sphere around a point



!!Information for selecting the 1D prfile


	real(8), parameter :: f_esc_lya = 0.01	!!This should be 1.. 0.1 to incorporate f_star!!
	integer, parameter :: n_sed=141		!!dim of the arry to read the SED

	real(8), parameter :: mmin_cut = 5d9	!!Minimum mass of halo to be considered in the simulation

!! Large Arrays used in the simulation 

	real(4), dimension(n_cell, n_cell, n_cell):: matrix_alpha




	real(8):: comm_redshift, gap_timestep
	integer :: spdim, num_halo, rank, numtasks, num_checkpoints, redshift_loop_index




END MODULE param

module param_halo
use param

	real(4), dimension(n_cell, n_cell, n_cell):: halo_mass_old, halo_mass_new
	integer, dimension(n_cell, n_cell, n_cell):: halo_trace, halo_age_old, halo_age_new
end module param_halo
