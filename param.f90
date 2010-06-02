module param
  use types,only:rprec
  $if ($MPI)
  use mpi
  $endif
  implicit none

  save

  private rprec  !--this is dumb.
  public
  
!---------------------------------------------------
! MPI PARAMETERS
!---------------------------------------------------
  $if ($MPI)
  $define $MPI_LOGICAL .true.
  $define $NPROC 4
  $else
  $define $MPI_LOGICAL .false.
  $define $NPROC 1
  $endif

  logical, parameter :: USE_MPI = $MPI_LOGICAL

  $undefine $MPI_LOGICAL

  $if ($MPI)
  integer :: status(MPI_STATUS_SIZE)
  $endif

  character (*), parameter :: path = './'

  !--this stuff must be defined, even if not using MPI
  character (8) :: chcoord  !--holds character representation of coord
  integer, parameter :: nproc = $NPROC  !--this must be 1 if no MPI
  integer :: ierr
  integer :: comm
  integer :: up, down
  integer :: global_rank
  integer :: MPI_RPREC, MPI_CPREC
  integer :: rank = -1   !--init to bogus (so its defined, even if no MPI)
  integer :: coord = -1  !--same here
  integer :: rank_of_coord(0:nproc-1), coord_of_rank(0:nproc-1)
  !--end mpi stuff
  
!---------------------------------------------------
! COMPUTATIONAL DOMAIN PARAMETERS
!---------------------------------------------------  
! characteristic length is H=z_i and characteristic velocity is u_star  
!   L_x, L_y, L_z, dx, dy, dz are non-dim. using H

  integer,parameter:: nx=32,ny=32,nz=(33-1)/nproc + 1
  integer, parameter :: nz_tot = (nz - 1) * nproc + 1
  integer,parameter:: nx2=3*nx/2,ny2=3*ny/2
  integer,parameter:: lh=nx/2+1,ld=2*lh,lh_big=nx2/2+1,ld_big=2*lh_big

  !this value is dimensional [m]:
  real(rprec),parameter::z_i=1._rprec   !dimensions in meters, height of BL
    
  !these values should be non-dimensionalized by z_i: 
  !set as multiple of BL height (z_i) then non-dimensionalized by z_i
  real(rprec),parameter::L_x=4.*z_i/z_i           
    !real(rprec),parameter::L_y=1.*z_i/z_i          
    !real(rprec),parameter::L_z=1./nproc * z_i/z_i
  real(rprec),parameter::L_y=(ny - 1.)/(nx - 1.)*L_x               ! ensure dy=dx
  real(rprec),parameter::L_z=(nz_tot - 1./2.)/(nx - 1.)/nproc*L_x  ! ensure dz = dx

  !these values are also non-dimensionalized by z_i:
    real(rprec),parameter::dz=nproc*L_z/(nz_tot-1./2.)
    real(rprec),parameter::dx=L_x/(nx-1),dy=L_y/(ny-1) ! Need to fix this wrt to autowrapping

  integer, parameter :: iBOGUS = -1234567890  !--NOT a new Apple product
  real (rprec), parameter :: BOGUS = -1234567890._rprec
  real(rprec),parameter::pi=3.1415926535897932384626433_rprec
  
!---------------------------------------------------
! MODEL PARAMETERS
!---------------------------------------------------   
  !Model type: 1->Smagorinsky; 2->Dynamic; 3->Scale dependent
  !            4->Lagrangian scale-sim   5-> Lagragian scale-dep
  !Models type: 1->static prandtl, 2->Dynamic
  integer,parameter::model=1,models=1,nnn=2
  !Cs is the Smagorinsky Constant
  !Co and nnn are used in the mason model for smagorisky coeff
  real(kind=rprec),parameter::Co=0.16_rprec

  !Test filter type: 1->cut off 2->Gaussian 3->Top-hat
  integer,parameter::ifilter=2

  ! u_star=0.45 m/s if coriolis_forcing=.FALSE. and =ug if coriolis_forcing=.TRUE.
  real(rprec),parameter::u_star=0.45_rprec,Pr=.4_rprec

  !--Coriolis stuff
  ! coriol=non-dim coriolis parameter,
  ! ug=horiz geostrophic vel, vg=transverse geostrophic vel
  logical,parameter::coriolis_forcing=.false.
  real(rprec),parameter::coriol=9.125E-05*z_i/u_star,      &
       ug=u_star/u_star,vg=0._rprec/u_star
	   
  real(rprec),parameter::vonk=0.4_rprec 
  integer,parameter::c_count=10000,p_count=10000
  integer, parameter :: cs_count = 5  !--tsteps between dynamic Cs updates	   
  
  ! nu_molec is dimensional m^2/s
  real(rprec),parameter::nu_molec=1.14e-5_rprec  
	   
  logical,parameter::use_bldg=.false.
  logical,parameter::molec=.false.,sgs=.true.,dns_bc=.false.  
  
!---------------------------------------------------
! TIMESTEP PARAMETERS
!---------------------------------------------------   
  integer, parameter :: nsteps = 100
 
  real (rprec), parameter :: dt = 2.e-5_rprec      !dt=2.e-4 usually works for 64^3
  real (rprec), parameter :: dt_dim = dt*z_i/u_star     !dimensional time step in seconds                                 
  
  integer :: jt                 ! global time-step counter
  integer :: jt_total           !--used for cumulative time (see io module)
  real(rprec) :: total_time, total_time_dim

  ! time advance parameters (AB2)
  real (rprec), parameter :: tadv1 = 1.5_rprec, tadv2 = 1._rprec - tadv1
  
!---------------------------------------------------
! BOUNDARY/INITIAL CONDITION PARAMETERS
!---------------------------------------------------  
  !--initu = true to read from a file; false to create with random noise
  logical, parameter :: initu = .false.
  !--initlag = true to initialize cs, FLM & FMM; false to read from vel.out
  logical, parameter :: inilag = .true.

  ! ubc: upper boundary condition: ubc=0 stress free lid, ubc=1 sponge
  integer,parameter::ubc=0

  !'wall', 'stress free'
  character (*), parameter :: lbc_mom = 'wall'

  ! prescribed inflow: constant or read from file
  ! read from file is not working properly
  logical,parameter::inflow=.false.
  logical, parameter :: use_fringe_forcing = .false.  
  
  ! position of right end of buffer region, as a fraction of L_x
  real (rprec), parameter :: buff_end = 1._rprec
  ! length of buffer region as a fraction of L_x
  real (rprec), parameter :: buff_len = 0.25_rprec
  
  real (rprec), parameter :: face_avg = 1.0_rprec

  logical, parameter :: read_inflow_file = .false.
  logical, parameter :: write_inflow_file = .false.

  ! records at position jx_s
  integer, parameter :: jt_start_write = 6

  ! forcing along top and bottom bdrys
  ! if inflow is true and force_top_bot is true, then the top & bottom
  ! velocities are forced to the inflow velocity
  logical, parameter :: force_top_bot = .false.

  logical, parameter :: use_mean_p_force = .true.
  real (rprec), parameter :: mean_p_force = 1._rprec * 1./(nproc*L_z)
  
!---------------------------------------------------
! DATA OUTPUT PARAMETERS
!--------------------------------------------------- 
  !records time-averaged data to files ./output/*_avg.dat
  logical, parameter :: tavg_calc = .true.
  integer, parameter :: tavg_nstart = 1, tavg_nend = nsteps

!  Turns instantaneous velocity recording on or off
  logical, parameter :: point_calc = .true.
  integer, parameter :: point_nstart = 1, point_nend = nsteps, point_nskip = 10
  integer, parameter :: point_nloc = 2
  real(rprec), save, dimension(3,point_nloc) :: point_loc = (/ &
      (/ L_x/2., L_y/2., 2._rprec /), &
      (/ 3._rprec, 2._rprec, 2._rprec /) &
      /)

  !  domain instantaneous output
  logical, parameter :: domain_calc = .true.
  integer, parameter :: domain_nstart = 100, domain_nend = nsteps, domain_nskip = 100
  
  !  x-plane instantaneous output
  logical, parameter :: xplane_calc   = .true.
  integer, parameter :: xplane_nstart = 100, xplane_nend = nsteps, xplane_nskip  = 100
  integer, parameter :: xplane_nloc   = 2
  real(rprec), save, dimension(xplane_nloc) :: xplane_loc = (/ 1.0, 3.0 /)

  !  y-plane instantaneous output
  logical, parameter :: yplane_calc   = .true.
  integer, parameter :: yplane_nstart = 100, yplane_nend = nsteps, yplane_nskip  = 100
  integer, parameter :: yplane_nloc   = 2
  real(rprec), save, dimension(yplane_nloc) :: yplane_loc = (/ 1.0, 3.0 /)  

  !  z-plane instantaneous output
  logical, parameter :: zplane_calc   = .true.
  integer, parameter :: zplane_nstart = 100, zplane_nend = nsteps, zplane_nskip  = 100
  integer, parameter :: zplane_nloc   = 7
  real(rprec), save, dimension(zplane_nloc) :: zplane_loc = (/ 0.733347, 1.550644, 1.959293, &
                                                            2.163617, 2.265780, 2.316861, &
                                                            2.342401 /)
 
  !------xxxxxxxxx--SCALARS_PARAMETERS--xxxxxxxxx---------------
  ! S_FLAG=1 for Theta and q, =0 for no scalars
  !logical,parameter::S_FLAG=.TRUE.,coupling_flag=.FALSE.,mo_flag=.TRUE.
  logical,parameter::S_FLAG=.false.
  !integer,parameter::DYN_init=2, SCAL_init=5, no_days=1
  integer,parameter::DYN_init=100, SCAL_init=5, no_days=1
  !integer,parameter::DYN_init=1, SCAL_init=1, no_days=1
  integer,parameter::patch_flag=1, remote_flag=0, time_start=0
  ! initu=.TRUE. & initsc=.FALSE read velocity fields from a binary file
  ! initu=.TRUE. & initsc=.TRUE. read velocity & scalar fields from a binary file
  ! initu=.FALSE. & S_FLAG=.TRUE. initialize velocity & scalar fields using ran
  ! initu=.FALSE. & S_FLAG=.FALSE. initialize velocity fields using ran
  logical,parameter::initsc=.false.
  ! lbc=0: prescribed surface temperature, lbc=1 prescribed surface flux
  ! (wT=0.06 Km/s)
  integer,parameter :: lbc=0
  ! Added a new parameter - sflux_flag for passive scalars with bldngs
  logical,parameter :: sflux_flag=.false.
  logical,parameter :: wt_evolution_flag=.FALSE.
  logical,parameter :: test_phase=.FALSE., vec_map=.FALSE., smag_sc=.FALSE.
  logical,parameter :: check_dt=.TRUE.
  integer,parameter :: stencil_pts=4
  logical,parameter :: coarse_grain_flag=.FALSE.
  !inversion strength (K/m)
  real(kind=rprec),parameter::g=9.81_rprec, inv_strength=0._rprec
  real(kind=rprec),parameter::theta_top=300._rprec,T_scale=300._rprec&
       ,wt_s=20._rprec,T_init=300._rprec
  real(kind=rprec),parameter::cap_thick=80._rprec, z_decay=1._rprec


end module param
