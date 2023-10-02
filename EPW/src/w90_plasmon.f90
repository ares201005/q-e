  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  module w90_plasmon
  !-----------------------------------------------------------------------

!! new notes: (Mar 2020)
!!    etf(:,:),               &!  interpolated eigenvalues (nbnd, nkqf)
!!    etf_k(:,:),             &!  Saved interpolated KS eigenenergies for later used in q-parallelization (nbnd, nkqf)

  USE kinds,      ONLY : dp, i4b
  USE epwcom,     ONLY : nbndsub, rlx_approx,efermi_read, fermi_energy
  USE phcom,      ONLY : nmodes
  USE elph2,      ONLY : etf, etf_k, etf_ks,               &
                         xqf, xkf, wkf,                    &
                         nkqf, nkf, nqf, nkqtotf, nqtotf,  &
                         ibndmin, ibndmax,dmef,            & 
                         epsi, efnew
  
  USE constants_epw, ONLY : one, two, zero, czero,         &
                            pi, twopi, hbar, hbarJ,        &
                            electron_SI, kb, ryd2ev,       &
                            ci, cone, czero 
!  USE mp,            ONLY : mp_barrier, mp_sum
!  USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool 
!  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : stdout, ionode, ionode_id
  
  implicit none
  
  real(kind=dp), parameter :: lightSpeed=2.99792458e8_dp !m/s
  real(kind=dp), parameter :: hbarc=hbarJ*lightSpeed/electron_SI*1.e10_dp !ev*A

  real(kind=dp), allocatable :: eig_k(:,:)
  real(kind=dp), allocatable :: eig_kall(:,:)
  real(kind=dp), allocatable :: xk(:,:)
  real(kind=dp), allocatable :: xk_all(:,:)

  !
!!  real(kind=dp),   allocatable :: qpoints(:,:)    ! plasmon wave vectors
!!  real(kind=dp),   allocatable :: eig(:)          ! eigen energy for certain k
!!  real(kind=dp),   allocatable :: del_eig(:,:,:)  ! derivative of E_{nk}
  real(kind=dp),   allocatable :: jhot_h_k(:,:)   ! hot hole distribution      (E,k)
  real(kind=dp),   allocatable :: jhot_e_k(:,:)   ! hot electron distribution  (E,k)
  real(kind=dp),   allocatable :: jhot_h(:)       ! hot hole distribution, (E)
  real(kind=dp),   allocatable :: jhot_e(:)       ! hot electron distribution, (E)
  real(kind=dp),   allocatable :: dhot_h(:)       ! hot hole density distribution
  real(kind=dp),   allocatable :: dhot_e(:)       ! hot electron density distribution
  real(kind=dp),   allocatable :: dos(:)          ! density of states
  real(kind=dp),   allocatable :: hot_freq_list(:)!
  real(kind=dp), allocatable :: occ(:)          ! occupancy
  real(kind=dp), allocatable :: fnk(:,:)      ! f_nk per pool
  real(kind=dp), allocatable :: fnk_all(:,:)      ! f_nk, total
  real(kind=dp), allocatable :: fnk0(:,:)     ! f_nk at t0
  real(kind=dp), allocatable :: dfnk(:,:)     ! d f_nk / dt
  real(kind=dp), allocatable :: fnk_e(:,:)      ! f_nk 
  real(kind=dp), allocatable :: fnk_h(:,:)      ! f_nk 
  real(kind=dp), allocatable :: dfnk_e(:,:)     ! d f_nk / dt
  real(kind=dp), allocatable :: dfnk_h(:,:)     ! d f_nk / dt
  real(kind=dp), allocatable :: trap(:,:)     ! trapping matrix
!!  complex(kind=dp),allocatable :: Ceps(:,:)       ! e-plasmon coupling matrix
!!  complex(kind=dp),allocatable :: Cepn(:,:)       ! e-phonon  coupling matrix
!!  complex(kind=dp),allocatable :: eSelf(:,:)     ! e-e scattering collision integral
!!  !
  integer,         allocatable :: eeklist(:,:)
  !
  real(kind=dp)    :: relaxtime,hot_tmax
  real(kind=dp)    :: qvalue, q0       ! plasmon wavenumber
  real(kind=dp)    :: Eplas            ! plasmon energy
  real(kind=dp)    :: telec            ! tempearture, in Kelvin
  real(kind=dp)    :: eta_smr,eta_smr2 ! adapative smear constant
  real(kind=dp)    :: dt, tmax         ! time step and tmax
  real(kind=dp)    :: trap_rate
  REAL(kind=DP) :: ef0
  !! Fermi energy level
!  complex(kind=dp) :: epsilonw         ! complex optical dielectric constant

  integer          :: nksqtotf            ! nkqtotf/2
  integer          :: nksqf            ! nkqf/2
  integer          :: num_wann         !
  integer          :: hot_nfreq
  integer          :: ntrap
  integer          :: nee
  !
  logical          :: tDynamics        ! do dynamic calculation or not
  logical          :: transit_dist     ! distribution of transitin rate (plasmon decay rate)
  logical          :: lepi             ! consider electron-phonon interaction or not
  logical          :: ltrap            ! consider trapping or not
  !

  contains


  !-----------------------------------------
  !
  !-----------------------------------------
    subroutine plasmon_setup
!!    use w90_parameters, only : num_wann, berry_kmesh
    implicit none
    !===========================================================!
    !         set up of plasmon module                          !
    !===========================================================!
    real(kind=dp) :: hot_freq_min
    real(kind=dp) :: hot_freq_max
    real(kind=dp) :: hot_freq_step=0.01_dp
    integer :: i
  
    telec = 300.0_dp      ! room temperature
    Eplas = 3.0_dp
    q0    = Eplas/hbarc   ! unit A
  
    hot_freq_min = - Eplas-1.0_dp
    hot_freq_max =   Eplas+1.0_dp
  
    hot_nfreq=nint((hot_freq_max-hot_freq_min)/hot_freq_step)
    allocate( hot_freq_list(hot_nfreq) )
    do i = 1, hot_nfreq
      hot_freq_list(i)=hot_freq_min+(i-1)*hot_freq_step
    enddo 
 
    num_wann = ibndmax-ibndmin+1

!!    allocate( qpoints(3,nqs) )
!!    allocate( eig(num_wann) )
    
    tDynamics      = .true.
    transit_dist   = .true.
    lepi           = .true.
    ltrap          = .false.
!!  
!!    ntrap = 9
!!    !trap_rate = 1.E-2_dp
!!   
!!    allocate(Ceps(num_wann,num_wann))

    allocate(occ(num_wann) )

    relaxtime = 1.E-3_dp
    hot_tmax  = 1.e3_dp 

    write(6,*) 'nkf      =', nkf
    write(6,*) 'nqf      =', nqf
    write(6,*) 'nkqf     =', nkqf
    write(6,*) 'nkqtotf  =', nkqtotf 
    write(6,*) 'nqtotf   =', nqtotf 
  !
  ! Fermi level and corresponding DOS
  !
  IF ( efermi_read ) THEN
    !
    ef0 = fermi_energy
    !
  ELSE
    !
    ef0 = efnew
    !
  ENDIF
  !

    if(tDynamics) then
      ! determine k mesh and number of k points
      allocate(fnk0(num_wann,nkf) )
      allocate( fnk(num_wann,nkf) )
      allocate( fnk_all(num_wann,nkf) )
      allocate(dfnk(num_wann,nkf) )
      allocate( fnk_e(num_wann,nkf) )
      allocate( fnk_h(num_wann,nkf) )
      allocate(dfnk_e(num_wann,nkf) )
      allocate(dfnk_h(num_wann,nkf) )
      allocate( trap(num_wann,nkf) )
  
      fnk    = zero
      fnk0   = zero
      fnk_e  = zero
      fnk_h  = zero
      dfnk   = zero
      dfnk_e = zero
      dfnk_h = zero

      tmax=hot_tmax
      dt=0.02_dp
    endif

    write(stdout,'(A,L)')    "  relaxation time approximation: ", rlx_approx
    write(stdout,'(A,f9.4)') "  relax time                     ", relaxtime
    write(stdout,'(A,f9.4)') "  trapping rate                  ", trap_rate
    write(stdout,'(A,f9.4)') "  maximum propagation time:      ", tmax
    write(stdout,'(A,f9.4)') "  propagation time step:         ", dt
    call flush(stdout)
 
    end subroutine plasmon_setup
  
  
    !===========================================================!
    !                   PUBLIC PROCEDURES                       ! 
    !===========================================================!
  
    subroutine plasmon_main (iqq, iq, totq)
    !============================================================!
    !                                                            !
    !============================================================!
!!      USE io_global,     ONLY : stdout, ionode

!!    use w90_constants,     only : pi,cone,czero,ci
!!    use w90_utility,       only : w0gauss,utility_rotate
!!    use w90_comms,         only : on_root
!!    use w90_io,            only : stdout,io_file_unit,io_error,seedname
!!    use w90_postw90_common,only : nrpts, irvec, ndegen, &
!!                                  get_occ_telec,kmesh_spacing, &
!!                                  fourier_R_to_k
!!    use w90_wan_ham,       only : get_D_h,get_eig_deleig
!!    use w90_parameters,    only : num_wann,nfermi,fermi_energy_list,  &
!!                                  kubo_adpt_smr_fac,kubo_adpt_smr_max,&
!!                                  kubo_smr_index
!!    use w90_eph

  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool 
  USE mp_world,      ONLY : mpime

    implicit none
    INTEGER, INTENT(IN) :: iqq, iq
    !! Q-point index from selecq.fmt window
    INTEGER, INTENT(IN) :: totq
    !! Total number of q-points from the selecq.fmt grid. 


    REAL(kind=DP), ALLOCATABLE :: xkf_all(:,:)
    REAL(kind=DP), ALLOCATABLE :: etf_all(:,:)

!!    !
!!    real(kind=dp)     :: omega, Ept
!!    real(kind=dp)     :: rdum_r(3),rdum_i(3)
!!    real(kind=dp)     :: kweight, kpt(3),qpt(2),db1,db2,db3,Delta_k,vdum(3)
!!    real(kind=dp)     :: joint_level_spacing,arg
!!    real(kind=dp)     :: rfac1, rfac2, occ_prod,delta,delta1,delta2
    integer           :: loop_xyz,loop_x,loop_y,loop_z,q_xy
    integer           :: m,n,i,j,k,idir,ir,ifreq
    integer           :: ik
!!    integer           :: file_unit1,file_unit2
!!    character(len=24) :: file_name1,file_name2
!!    !
!!    complex(kind=dp)  :: ztmp
!!    !
!!    complex(kind=dp), allocatable :: HH(:,:)
!!    complex(kind=dp), allocatable :: delHH(:,:,:)
!!    complex(kind=dp), allocatable :: UU(:,:)
!!    complex(kind=dp), allocatable :: D_h(:,:,:)
!!    complex(kind=dp), allocatable :: vv_r(:,:,:,:)
!!    complex(kind=dp), allocatable :: vv_k(:,:,:)
!!  
    !
    if (ionode) call printHeader(stdout)

    if(ionode) then
       write(stdout,'(1x,a)')&
            '------------------------------------------'
       write(stdout,'(1x,a)')&
            '! Properties calculated in module plasmon!'
       write(stdout,'(1x,a)')&
            '------------------------------------------'
       write(stdout,'(/,3x,a)') '* initial distribution'
       write(stdout,'(/,3x,a)') '* dynamics/relaxation'
    endif
  
    !
    call plasmon_setup
    nksqtotf = nkqtotf/2
    nksqf = nkqf/2

  ! The k points are distributed among pools: here we collect them
    allocate(xk_all(3,nksqtotf), eig_kall(nbndsub, nksqtotf))
    allocate(xk(3,nksqf), eig_k(nbndsub, nksqf))
    xk_all  = zero
    xk      = zero
    eig_k   = zero
    eig_kall= zero

!if(iqq == totq) then
    ALLOCATE ( xkf_all      ( 3,       nkqtotf ), &
               etf_all      ( nbndsub, nkqtotf ) )

    xkf_all(:,:) = zero
    etf_all(:,:) = zero
    !

    write(6,*) 'nkf, pool_id', nkf, my_pool_id, npool
    write(6,*) 'xkf'
    do i = 1, nkqf
      write(6,*) i, xkf(:,i)
    enddo


    !do ik = 1, nkf
    !   ikk = 2 * ik - 1
    !   xkf_all(:,ik+lower_band-1) = xkf(:,ikk)
    !enddo

    write(stdout,*) 'size (xkf, xkf_all)', size(xkf), size(xkf_all)
    write(stdout,*) 'pool gathering!', allocated(xkf), allocated(etf), nkqtotf,nkqf
    call flush(stdout)

    

#if defined(__MPI)
    !CALL mp_barrier(inter_pool_comm)
    !
    ! note that poolgather2 works with the doubled grid (k and k+q)
    !
    write(stdout,*) 'pool gathering-0!'
    CALL poolgather2 ( 3,       nkqtotf, nkqf, xkf,    xkf_all  )
    !CALL poolgather_yz ( 3,       nkqtotf, nkqf, xkf,    xkf_all  )
    write(stdout,*) 'pool gathering-1!'
    call flush(stdout)

    CALL poolgather2 ( nbndsub, nkqtotf, nkqf, etf,    etf_all  )
    !!CALL poolgather_yz ( nbndsub, nkqtotf, nkqf, etf,    etf_all  )
    write(stdout,*) 'pool gathering-2!'
    !call flush(stdout)

    !write(stdout,*) 'mp_barrier!'
    !CALL mp_barrier(inter_pool_comm)
    !write(stdout,*) 'mp_barrier,done!'
    call flush(stdout)
    !
#else
    !
    xkf_all = xkf
    etf_all = etf
    !
#endif


    write(6,*) 'xf_all'
    write(stdout,*) 'nbndsub=',nbndsub
    write(stdout,*) 'ibndmin=',ibndmin
    write(stdout,*) 'ibndmax=',ibndmax
    write(stdout,*) 'fermi=  ',ef0,ryd2ev
    do i = 1, nksqtotf
      xk_all(:,i) = xkf_all(:,2*i-1)
      eig_kall(:,i) = etf_all(:,2*i-1)
      write(stdout,'(I5,3f9.4)') i, xk_all(:,i)
      do j = ibndmin, ibndmax
         write(stdout,'(I5,8e15.6)') j,ryd2ev*(eig_kall(j,i)-ef0)
      enddo
    enddo

    write(6,*) 'xk and eig_k of this pool'
    do i = 1, nksqf
      xk(:,i) = xkf(:,2*i-1)
      eig_k(:,i) = etf(:,2*i-1)
      write(stdout,'(I5,3f9.4)') i, xk(:,i)
      do j = ibndmin, ibndmax
         write(stdout,'(I5,8e15.6)') j,ryd2ev*(eig_k(j,i)-ef0)
      enddo
    enddo

    ! convert eig to ev
    eig_k = eig_k * ryd2ev
    eig_kall = eig_kall * ryd2ev

!endif
 
    
!!    call get_qpoints
    !
!!    db1=1.0_dp/real(kmesh(1),dp)
!!    db2=1.0_dp/real(kmesh(2),dp)
!!    db3=1.0_dp/real(kmesh(3),dp)
!!    kweight=db1*db2*db3
!!    !
!!    allocate( vv_r(num_wann,num_wann,nrpts,3))
!!    allocate( vv_k(num_wann,num_wann,3))
!!    vv_r=0.d0
!!    vv_k=czero
!!  
!!    ! read vv_r from file
!!    file_name1=trim(seedname)//'-p.dat.ZY'
!!    file_unit1=io_file_unit()
!!    open(file_unit1,file=file_name1, status='old')
!!    read(file_unit1,*) m
!!    read(file_unit1,*) n
!!    if(m/=num_wann) then
!!      write(stdout,*) "Error:mismathc of num_wann"
!!      stop
!!    else if(n/=nrpts) then
!!      write(stdout,*) "Error:mismathc of nrpts"
!!      stop
!!    endif
!!    read(file_unit1,'(15I5)') (ndegen(i),i=1,nrpts)
!!    do ir=1,nrpts
!!      do n=1,num_wann
!!      do m=1,num_wann
!!        read(file_unit1,*) irvec(1:3,ir), i,j, (rdum_r(idir),rdum_i(idir),idir=1,3)
!!        do idir=1,3
!!          vv_r(i,j,ir,idir)=dcmplx(rdum_r(idir),rdum_i(idir))
!!        enddo
!!      enddo
!!      enddo
!!    enddo 
!!    close(file_unit1) 
!!    !
!!    file_unit1=io_file_unit()
!!    file_name1="dielectric.dat"
!!    open(file_unit1,file=file_name1,STATUS='unknown')
!!    omega =0.1_dp
!!    do while(omega<=9.0_dp)
!!       call dielectric(1,'Au',omega,epsilonw)
!!       write(file_unit1,'(f9.4,2X,2es15.6)') omega, dble(epsilonw), dimag(epsilonw)
!!   
!!       omega = omega + 0.1_dp
!!    enddo
!!    close(file_unit1)
!!  
!!  
!!    file_unit2=io_file_unit()
!!    file_name2="dispersion.dat"
!!    open(file_unit2,file=file_name2,STATUS='unknown')
!!    omega =0.1_dp
!!    do while(omega<=7.0_dp)
!!       call dielectric(1,'Au',omega,epsilonw)
!!       ztmp = omega*sqrt(epsilonw/(epsilonw+cone))
!!       qvalue = dble(ztmp)
!!       write(file_unit2,'(2es15.6)') qvalue, omega
!!       omega = omega + 0.1_dp
!!    enddo
!!    close(file_unit2)
!!  
!!    ! determine qmesh: momemtum vector of plasmon
!!    
!!    !--------------------------------------------------------------
!!    !initial distribution I: transition rate distribution (Gamma) !
!!    !--------------------------------------------------------------
!!    
!!    allocate(HH(num_wann,num_wann))
!!    allocate(delHH(num_wann,num_wann,3))
!!    allocate(UU(num_wann,num_wann))
!!    allocate(D_h(num_wann,num_wann,3))
!!    allocate(del_eig(num_wann,3,nkf))
!!  
    if(transit_dist) then
      allocate( jhot_h_k(hot_nfreq, nkf) )
      allocate( jhot_e_k(hot_nfreq, nkf) )
      allocate( jhot_h(hot_nfreq) )
      allocate( jhot_e(hot_nfreq) )
      jhot_h_k= zero
      jhot_e_k= zero
      jhot_h  = zero
      jhot_e  = zero
    endif
  
    if(tDynamics) then
      allocate( dhot_h(hot_nfreq) )
      allocate( dhot_e(hot_nfreq) )
      allocate( dos(hot_nfreq) )
      dhot_h = zero
      dhot_e = zero
      dos    = zero
    endif
  
    do ik=1,nkf

!!      loop_x= loop_xyz/(kmesh(2)*kmesh(3))
!!      loop_y=(loop_xyz-loop_x*(kmesh(2)*kmesh(3)))/kmesh(3)
!!      loop_z= loop_xyz-loop_x*(kmesh(2)*kmesh(3))-loop_y*kmesh(3)
!!      kpt(1)=loop_x*db1
!!      kpt(2)=loop_y*db2
!!      kpt(3)=loop_z*db3
!!  
!!      xkf(1:3,loop_xyz+1)=kpt(1:3)
!!  
!!      !write(3000,'(I5,3f15.7)') loop_xyz+1, xkf(1:3,loop_xyz+1)
!!  
!!      !write(stdout,*) db1,db2,db3
!!      !write(stdout,*) loop_x,loop_y,loop_z
!!      !write(stdout,*) kpt
!!      
!!      call get_eig_deleig(kpt,eig,del_eig(1:num_wann,1:3,loop_xyz+1),HH,delHH,UU)
!!  
!!      ! get the trapping matrix in the Hamiltonian gauge
!!      ! O^H(k)= U^\dag(k) O^w(k) U(k) where O^w is in wannier representation
!!      if(ltrap) then
!!      trap(:,:,loop_xyz+1)=czero
!!      do i=1,num_wann
!!        do j=1,num_wann
!!          ztmp = czero
!!          if(eig(i)<fermi_energy_list(1).or.eig(j)<fermi_energy_list(1)) cycle
!!          do m=1,ntrap ! the first ntrap orbital is connected to trapping site
!!            ztmp = ztmp + dconjg(UU(m,j))*trap_rate*UU(m,i)
!!          enddo
!!          trap(j,i,loop_xyz+1) = ztmp
!!        enddo
!!      enddo
!!      endif
!!     
!!      !write(stdout,*) HH
!!      !write(stdout,*) UU
!!      !write(stdout,*) eig
!!  
!!      Delta_k=kmesh_spacing(kmesh)
!!      if(tDynamics) eig_k(:,loop_xyz+1)=eig(:)
!!      
!!      call get_occ_telec(eig,occ,efnew,telec)
!!      call get_D_h(delHH,UU,eig,D_h)
!!  
!!      ! (1) transform the momentum matrix from reald space to k space
!!      ! (2) transform the momentum matrxi to the eigen space of HH_k
!!      call fourier_R_to_k(kpt,vv_r(:,:,:,1), vv_k(:,:,1),0)
!!      call fourier_R_to_k(kpt,vv_r(:,:,:,2), vv_k(:,:,2),0)
!!      call fourier_R_to_k(kpt,vv_r(:,:,:,3), vv_k(:,:,3),0)
!!  
!!      do i = 1, 3     
!!        vv_k(:,:,i)=utility_rotate(vv_k(:,:,i),UU,num_wann)
!!      enddo
!!  
!!      do q_xy=0,product(qmesh)-1
!!        ! get e-plasmon coupling matrxi M(k,q), not finished!
!!        Ceps=czero
!!        ! make a subroutine for calculating e-plasmon coupling
!!        ! matrix for given plasmon energy(momentum) and electron
!!        ! momentum matrix element in wannier presentation. 
!!        
!!        Ceps(:,:) = vv_k(:,:,1)*0.5_dp
!!        !Ceps = cone
!!  
!!        do m=1,num_wann
!!          if(tDynamics) then 
!!             rho(m,m,loop_xyz+1)=occ(m)*cone
!!             rho0(m,m,loop_xyz+1)=occ(m)*cone
!!             rho_e(m,m,loop_xyz+1)=occ(m)*cone
!!             rho_h(m,m,loop_xyz+1)=(1.0_dp-occ(m))*cone
!!          endif
!!          do n=1,num_wann
!!             if(m==n) cycle
!!             ! Eq.(35) YWVS07 
!!             vdum(:)=del_eig(m,:,loop_xyz+1)-del_eig(n,:,loop_xyz+1)
!!             !if(m==n) vdum(:)=del_eig(m,:,loop_xyz+1)
!!             joint_level_spacing=sqrt(dot_product(vdum(:),vdum(:)))*Delta_k
!!             eta_smr=min(joint_level_spacing*kubo_adpt_smr_fac,kubo_adpt_smr_max)
!!             if(eta_smr<1.d-3) eta_smr=1.d-3
!!  
!!             vdum(:)=del_eig(m,:,loop_xyz+1)
!!             joint_level_spacing=sqrt(dot_product(vdum(:),vdum(:)))*Delta_k
!!             eta_smr2=min(joint_level_spacing*kubo_adpt_smr_fac,kubo_adpt_smr_max)
!!             if(eta_smr2<1.d-3) eta_smr2=1.d-3
!!  
!!  
!!             rfac1=real(Ceps(n,m)*Ceps(m,n),dp)
!!             occ_prod=occ(n)*(1.0_dp-occ(m))
!!  
!!             if(transit_dist) then
!!                !transit_rate
!!               do ifreq=1,hot_nfreq
!!                 Ept=hot_freq_list(ifreq)
!!  
!!                 arg=(eig(m)-eig(n)-real(Eplas,dp))/eta_smr
!!                 delta=w0gauss(arg,kubo_smr_index)/eta_smr
!!  
!!                 arg=(eig(m)-fermi_energy_list(1)-real(Ept,dp))/eta_smr
!!                 delta1=w0gauss(arg,kubo_smr_index)/eta_smr
!!  
!!                 arg=(eig(n)-fermi_energy_list(1)-real(Ept,dp))/eta_smr
!!                 delta2=w0gauss(arg,kubo_smr_index)/eta_smr
!!  
!!                 if(eig(m)>=fermi_energy_list(1)) &
!!                   jhot_e_k(ifreq,loop_xyz+1)=jhot_e_k(ifreq,loop_xyz+1)+occ_prod*rfac1*delta*delta1
!!                 if(eig(n)<=fermi_energy_list(1)) &
!!                   jhot_h_k(ifreq,loop_xyz+1)=jhot_h_k(ifreq,loop_xyz+1)+occ_prod*rfac1*delta*delta2
!!                 !
!!               enddo
!!             endif
!!  
!!             if(tDynamics) then
!!               ! initial distribution for propagation
!!              
!!               if(eig(m)>fermi_energy_list(1)) then
!!                 arg=(eig(m)-eig(n)-real(Eplas,dp))/eta_smr
!!               else
!!                 arg=(eig(m)-eig(n)+real(Eplas,dp))/eta_smr
!!               endif
!!  
!!               !lorentzian shape
!!               !delta=1.0_dp/(1.0_dp+arg*arg)/pi
!!               !delta=delta/eta_smr*rfac1
!!  
!!               !gaussian shape
!!               delta=w0gauss(arg,kubo_smr_index)/eta_smr*rfac1
!!  
!!               if(eig(m)>fermi_energy_list(1)) then
!!                 drho(m,m,loop_xyz+1)  =drho(m,m,loop_xyz+1)  +occ(n)*delta*cone
!!                 drho_e(m,m,loop_xyz+1)=drho_e(m,m,loop_xyz+1)+occ(n)*delta*cone
!!               else
!!                 drho(m,m,loop_xyz+1)  =drho(m,m,loop_xyz+1)  -(1.0_dp-occ(n))*delta*cone
!!                 drho_h(m,m,loop_xyz+1)=drho_h(m,m,loop_xyz+1)+(1.0_dp-occ(n))*delta*cone
!!               endif
!!               !
!!               do ifreq=1,hot_nfreq
!!                 Ept=hot_freq_list(ifreq)
!!                 arg=(eig(m)-fermi_energy_list(1)-real(Ept,dp))/eta_smr2
!!                 delta2=w0gauss(arg,kubo_smr_index)/eta_smr2
!!                 if(eig(m)>fermi_energy_list(1)) then
!!                   dhot_e(ifreq)=dhot_e(ifreq)+occ(n)*delta*delta2*kweight
!!                 else
!!                   dhot_h(ifreq)=dhot_h(ifreq)+(1.0_dp-occ(n))*delta*delta2*kweight
!!                 endif
!!               enddo
!!               !
!!             endif
!!          enddo
!!        enddo
!!  
!!      enddo  ! q_xy
!!      ! sum over different k
!!      if(transit_dist) then
!!        jhot_h(:)=jhot_h(:)+jhot_h_k(:,loop_xyz+1)*kweight
!!        jhot_e(:)=jhot_e(:)+jhot_e_k(:,loop_xyz+1)*kweight
!!      endif
    enddo  ! loop_xyz
!! 



    if(tDynamics) then
      allocate(eeklist(8,1000000))
      write(stdout,*) 'calculate ee klist'
      write(stdout,*) num_wann, nkqf, size(etf), dabs(maxval(etf))
      call flush(stdout)
      call ee_klist(nkqf/2,nksqtotf, num_wann, xk, eig_k, xk_all, eig_kall, eeklist)
    endif

!!    if(transit_dist) then
!!       !
!!       ! hot-carrier distribution
!!       !
!!       call printfile(zero,'-jhot.dat.new',0,hot_nfreq,hot_freq_list,jhot_h,jhot_e,2)
!!    endif

    if(tDynamics) then
      !
      ! hot-carrier density distribution
      !
!!      !call sumk(nkf,hot_nfreq,Delta_k,kweight,eig_k,del_eig,hot_freq_list,fermi_energy_list(1),drho_e,dhot_e)
!!      !call sumk(nkf,hot_nfreq,Delta_k,kweight,eig_k,del_eig,hot_freq_list,fermi_energy_list(1),drho_h,dhot_h)
!!  
!!      !call sumk2(nkf,hot_nfreq,Delta_k,kweight,eig_k,del_eig,hot_freq_list,fermi_energy_list(1),drho,dhot_e,dhot_h)
!!      call printfile(zero,'-dhot.dat',0,hot_nfreq,hot_freq_list,dhot_h,dhot_e,2)
!!  
      !
      ! equilibrium electron density distribution
      !
!!      call sumk(nkf,hot_nfreq,Delta_k,kweight,eig_k,del_eig,hot_freq_list,fermi_energy_list(1),rho,dhot_h)
!!      call printfile(zero,'-fd.dat',0,hot_nfreq,hot_freq_list,dhot_h,dhot_e,1)
    endif
  
!!    if(tDynamics) then 
!!      fnk = fnk + dfnk
!!     fnk_e = fnk_e + dfnk_e
!!      fnk_h = fnk_h + dfnk_h
!!    endif
!!  
!!    if(lepi) then
!!      !---------------------------------------------------
!!      ! calculate phonon modes and e-ph coupling matrix  !
!!      ! or get the self-energy from wannier projected    !
!!      ! real-space electron-phonon coupling matrix       !
!!      ! g(R_e,R_p)                                       !
!!      !---------------------------------------------------
!!      call eph_setup(nkf,nkf,zero)
!!  
!!      call self_elec(nkf,nkf,nmode,num_wann,zero,fermi_energy_list(1),wqf,eig_k,freq,epmk)
!!  
!!      !compare eigenvalues
!!      open(120,file="bands-compare.dat",status="replace")
!!      do i=1,nkf
!!        do m=1,num_wann
!!          write(120,"(2f15.5)") eig_k(m,i)-fermi_energy_list(1),eigk2(m,i)
!!        enddo
!!      enddo
!!      !
!!    else
!!      if(tDynamics) then
!!         allocate( relax_time(num_wann,nkf) )
!!         relax_time = relaxtime 
!!      endif
!!    endif
!!  
!!    ! test e-e collision
!!    if(tDynamics) then
!!      write(stdout,*) 'calculate ee scattering is included!'
!!      allocate(eSelf(num_wann,nkf))
!!      !call eeself(num_wann,nkf,xkf,eig_k,rho,kweight,eSelf)
!!      !stop
!!    endif
!!  
!!    if(tDynamics) then
!!      !--------------
!!      !  dynamics   !
!!      !--------------
!!      call propagation
!!    endif
!!  
!!    !
!!    deallocate(HH,delHH,UU,del_eig)
!!    deallocate(eig, qpoints,Ceps)
!!    if(tDynamics) then
!!      deallocate(eig_k,occ,rho_e,rho_h,drho_e,drho_h)
!!      deallocate(dhot_e,dhot_h,dos)
!!    endif
    !
    if(ionode) call printEnding(stdout)
    !
    end subroutine plasmon_main
!!  
!!  
!!    !===========================================================!
!!    !                   PRIVATE PROCEDURES                      ! 
!!    !===========================================================!
!!  
!!    subroutine get_qpoints
!!    !==============================================================!
!!    ! momemtum vector of plasmon, obey the disperison relation     !
!!    !                                                              !
!!    ! k^2=k^2_x + k^2_y;                                           !
!!    ! and k = w/c sqrt{epsilon(w)/(epsilon(w)+1)}                  !
!!    !==============================================================!
!!    use w90_constants, only : dp, cone
!!    use w90_io,        only : stdout
!!    implicit none
!!    !
!!    integer :: i,j
!!    complex(kind=dp) :: ztmp
!!  
!!    !
!!    ! dielectric constant epsilon(w)
!!    !
!!    call dielectric(1,'Au',Eplas,epsilonw)
!!    write(stdout,'(A,f9.4,A,2es15.5)') &
!!          '  dielectric constant at frequency', &
!!          Eplas, " eV is ",dble(epsilonw), dimag(epsilonw)
!!       
!!    !
!!    ! qvalue = (w/c) * sqrt{epsilon(w)/(epsilon(w)+1)}
!!    !
!!    ztmp = Eplas*sqrt(epsilonw/(epsilonw+cone))
!!    qvalue = dble(ztmp)
!!  
!!    do i=1,nqs
!!  
!!    enddo
!!  
!!    end subroutine get_qpoints
!!  
!!  
!!    !==============================================================!
!!    !                                                              !
!!    !==============================================================!
!!    subroutine sumk(nkf,nfreq,Delta_k,kweight,eig_k,del_eig,freq_list,efermi,rhok,rhoe)
!!    use w90_constants,  only : dp
!!    use w90_parameters, only : num_wann,kubo_adpt_smr_max,kubo_adpt_smr_fac, &
!!                               kubo_smr_index
!!    use w90_utility,    only : w0gauss
!!    implicit none
!!    integer,         intent(in) :: nkf,nfreq
!!    real(kind=dp),   intent(in) :: Delta_k,kweight
!!    real(kind=dp),   intent(in) :: eig_k(num_wann,nkf)
!!    real(kind=dp),   intent(in) :: del_eig(num_wann,3,nkf)
!!    real(kind=dp),   intent(in) :: freq_list(nfreq)
!!    real(kind=dp),   intent(in) :: efermi
!!    complex(kind=dp),intent(in) :: rhok(num_wann,num_wann,nkf)
!!    real(kind=dp),   intent(out):: rhoe(nfreq)
!!    !
!!    real(kind=dp) :: joint_level_spacing
!!    real(kind=dp) :: vdum(3)
!!    real(kind=dp) :: arg,Ept,delta,eta_smr
!!    integer       :: m,n,k,ifreq
!!    !
!!  
!!    rhoe=zero
!!    do k=1,nkf
!!      do m=1,num_wann
!!        vdum(:)=del_eig(m,:,k)
!!        joint_level_spacing=sqrt(dot_product(vdum(:),vdum(:)))*Delta_k
!!        eta_smr=min(joint_level_spacing*kubo_adpt_smr_fac,kubo_adpt_smr_max)
!!        if(eta_smr<1.d-3) eta_smr=1.d-3
!!        !
!!        do ifreq=1,nfreq
!!          Ept=freq_list(ifreq)
!!          arg=(eig_k(m,k)-efermi-real(Ept,dp))/eta_smr
!!          delta=w0gauss(arg,kubo_smr_index)/eta_smr
!!          rhoe(ifreq)=rhoe(ifreq)+real(rhok(m,m,k),dp)*delta*kweight
!!        enddo
!!        !
!!      enddo
!!    enddo
!!  
!!    end subroutine sumk 
!!    !
!!  
!!    subroutine sumelec(nkf,kweight,efermi,eig_k,rhok,ntot)
!!    use w90_constants,  only : dp
!!    use w90_parameters, only : num_wann
!!    !
!!    implicit none
!!    integer :: nkf
!!    real(kind=dp),   intent(in) :: kweight
!!    real(kind=dp),   intent(in) :: efermi
!!    real(kind=dp),   intent(in) :: eig_k(num_wann,nkf)
!!    complex(kind=dp),intent(in) :: rhok(num_wann,num_wann,nkf)
!!    real(kind=dp),   intent(out):: ntot(2)
!!    !
!!    integer :: m,k
!!    !
!!    ntot=zero
!!    do k=1,nkf
!!      do m=1,num_wann
!!        if(eig_k(m,k)>=efermi) ntot(1) = ntot(1) + dble(rhok(m,m,k)) * kweight
!!        if(eig_k(m,k)<=efermi) ntot(2) = ntot(2) + dble(1.0_dp-rhok(m,m,k)) * kweight
!!      enddo
!!    enddo
!!    !
!!    end subroutine sumelec
!!  
!!    !==============================================================!
!!    !                                                              !
!!    !                                                              !
!!    !==============================================================!
!!    subroutine sumk2(nkf,nfreq,Delta_k,kweight,eig_k,del_eig,freq_list,efermi,rhok,rhoe,rhoh)
!!    use w90_constants,  only : dp
!!    use w90_parameters, only : num_wann,kubo_adpt_smr_max,kubo_adpt_smr_fac, &
!!                               kubo_smr_index
!!    use w90_utility,    only : w0gauss
!!    implicit none
!!    integer,         intent(in) :: nkf,nfreq
!!    real(kind=dp),   intent(in) :: Delta_k,kweight
!!    real(kind=dp),   intent(in) :: eig_k(num_wann,nkf)
!!    real(kind=dp),   intent(in) :: del_eig(num_wann,3,nkf)
!!    real(kind=dp),   intent(in) :: freq_list(nfreq),efermi
!!    complex(kind=dp),intent(in) :: rhok(num_wann,num_wann,nkf)
!!    real(kind=dp),   intent(out):: rhoe(nfreq)
!!    real(kind=dp),   intent(out):: rhoh(nfreq)
!!    !
!!    real(kind=dp) :: joint_level_spacing
!!    real(kind=dp) :: vdum(3)
!!    real(kind=dp) :: arg,Ept,delta,eta_smr
!!    integer       :: m,n,k,ifreq
!!    !
!!  
!!    rhoe=zero
!!    rhoh=zero
!!    do k=1,nkf
!!      do m=1,num_wann
!!        vdum(:)=del_eig(m,:,k)
!!        joint_level_spacing=sqrt(dot_product(vdum(:),vdum(:)))*Delta_k
!!        eta_smr=min(joint_level_spacing*kubo_adpt_smr_fac,kubo_adpt_smr_max)
!!        if(eta_smr<1.d-3) eta_smr=1.d-3
!!        !
!!        do ifreq=1,nfreq
!!          Ept=freq_list(ifreq)
!!          arg=(eig_k(m,k)-efermi-real(Ept,dp))/eta_smr
!!          delta=w0gauss(arg,kubo_smr_index)/eta_smr
!!          if(eig_k(m,k)>efermi) then
!!            rhoe(ifreq)=rhoe(ifreq)+real(rhok(m,m,k),dp)*delta*kweight
!!          else
!!            rhoh(ifreq)=rhoh(ifreq)-real(rhok(m,m,k),dp)*delta*kweight
!!          endif
!!        enddo
!!        !
!!      enddo
!!    enddo
!!  
!!    end subroutine sumk2
!!    !
!!   
!!  
!!    !==============================================================!
!!    ! initial condition for propagation                            !
!!    ! the electron density is calculated from:                     !
!!    !                                                              !
!!    ! rho^<_nk=  -i* int G^r_nk \Sigma^<_{nk} G^a_{nk} dE /2Pi     !
!!    ! rho^>_nk=  -i* int G^r_nk \Sigma^>_{nk} G^a_{nk} dE /2Pi     !
!!    !                                                              !
!!    !                                                              !
!!    ! \Sigma^<_{nk}= M[NG^<(E-w) + (N+1) G^<(E+w)]M                !
!!    ! \Sigma^>_{nk}= M[NG^>(E+w) + (N+1) G^>(E-w)]M                !
!!    !                                                              !
!!    !==============================================================!
!!    subroutine initialize
!!    use w90_constants,  only : dp, ci,czero,cone
!!    use w90_parameters, only : num_wann
!!    implicit none
!!    
!!    !
!!    ! can be merged with the transition rate subroutine
!!    ! 
!!  
!!  
!!    end subroutine initialize
!!  
!!  
!!    subroutine propagation
!!    use w90_constants,     only : dp, ci,czero,cone
!!    use w90_parameters,    only : num_wann,fermi_energy_list
!!    use w90_io,            only : stdout,io_file_unit,io_error,seedname
!!    use w90_postw90_common,only : kmesh_spacing
!!    implicit none
!!    !
!!    complex(kind=dp), allocatable :: rho_d(:,:,:)
!!    complex(kind=dp), allocatable :: ztmp0(:,:,:)
!!    complex(kind=dp), allocatable :: ztmp1(:,:,:)
!!    !
!!    real(kind=dp) :: tt, dt2,dt3,dt6
!!    real(kind=dp) :: kweight, Delta_k
!!    real(kind=dp) :: ntot(2), n_exci, n_cool,n_trap
!!    integer           :: file_unit,ifreq
!!    integer           :: i,k
!!    character(len=24) :: file_name
!!    !
!!    allocate( ztmp0(num_wann,num_wann,nkf) )
!!    allocate( ztmp1(num_wann,num_wann,nkf) )
!!    allocate( rho_d(num_wann,num_wann,nkf) )
!!   
!!    kweight=1.0_dp/(real(kmesh(1),dp)*real(kmesh(2),dp)*real(kmesh(3),dp))
!!    Delta_k=kmesh_spacing(kmesh) 
!!    tt=0.0_dp
!!    dt2=dt/2.0_dp
!!    dt3=dt/3.0_dp
!!    dt6=dt/6.0_dp
!!      
!!    file_name=trim(seedname)//'-dhot-td.dat'
!!    write(stdout,'(/,3x,a)') '* '//file_name
!!    file_unit=io_file_unit()
!!    open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
!!    close(file_unit)
!!  
!!    write(stdout,'(2x,A)') "======================"
!!    write(stdout,'(2x,A)') "| propagation begins |"
!!    write(stdout,'(2x,A)') "======================"
!!  
!!    call sumelec(nkf,kweight,fermi_energy_list(1),eig_k,rho,ntot)
!!    write(stdout,'(2x,A,2f15.7)') "total number of electron is:", ntot
!!    n_cool=ntot(2)
!!    n_exci=ntot(1)
!!  
!!    trap_rho = 0.d0
!!  
!!    do 
!!      if(mod(nint(tt/dt),10)==0) then
!!        write(stdout,'(2x, A,f12.4)') "current time:", tt
!!      endif
!!  
!!      !----------------------------------
!!      ! statistics                      !
!!      ! get average occ at given energy !
!!      !----------------------------------
!!  
!!      !if( mod(nint(tt/1.e-2_dp),nint(200.0_dp/1.e-2_dp)/100)==0 ) then
!!      if( mod(nint(tt/1.e-2_dp),500) ==0 ) then
!!        !call sumk2(nkf,hot_nfreq,Delta_k,kweight,eig_k,del_eig,hot_freq_list,fermi_energy_list(1),drho,dhot_e,dhot_h)
!!        call sumk(nkf,hot_nfreq,Delta_k,kweight,eig_k,del_eig,hot_freq_list,fermi_energy_list(1),rho,dhot_h)
!!        write(stdout,'(2x,A,f9.3)') 'Time=',tt
!!        call printfile(tt,file_name,file_unit,hot_nfreq,hot_freq_list,dhot_h,dhot_e,3)
!!  
!!        call sumelec(nkf,kweight,fermi_energy_list(1),eig_k,rho,ntot)
!!        write(stdout,'(2x,A,f15.7,2f15.7)') "total number of electron at time ", tt, ntot
!!        write(stdout,'(2x,A,f15.7,2f15.7)') "number of trapped electron       ", tt, trap_rho
!!      endif
!!  
!!  
!!      if(rlx_approx.and.(.not.ltrap)) then
!!        call derivtr(tt,rho_d)
!!        rho=rho_d
!!      else
!!        call derivtr(tt,rho_d)
!!        ztmp0 = rho
!!        rho   = ztmp0 + dt2*rho_d
!!        ztmp1 = ztmp0 + dt6*rho_d
!!        trap_rho(1)=trap_rho(1)+trap_rho(2)*dt6*kweight
!!  
!!        call derivtr(tt+dt2,rho_d)
!!        rho   = ztmp0 + dt2*rho_d
!!        ztmp1 = ztmp1 + dt3*rho_d
!!        trap_rho(1)=trap_rho(1)+trap_rho(2)*dt3*kweight
!!        
!!        call derivtr(tt+dt2,rho_d)
!!        rho   = ztmp0 + dt *rho_d
!!        ztmp1 = ztmp1 + dt3*rho_d
!!        trap_rho(1)=trap_rho(1)+trap_rho(2)*dt3*kweight
!!     
!!        call derivtr(tt+dt,rho_d)
!!        rho = ztmp1 + dt6*rho_d
!!        trap_rho(1)=trap_rho(1)+trap_rho(2)*dt6*kweight
!!      endif
!!  
!!  
!!      !debug
!!      if( abs(tt-tmax)<1.d-6 ) then
!!        write(stdout,'(/,A,/)') 'occupation number'
!!        do k=1,nkf
!!          do i=1,num_wann
!!            !write(stdout,'(2I5,3e15.7)') i,k, eig_k(i,k)-fermi_energy_list(1), rho(i,i,k)
!!          enddo
!!        enddo
!!      endif
!!      !
!!  
!!      if(tt>=tmax) exit
!!      tt=tt+dt
!!    enddo
!!  
!!    n_cool=n_cool-ntot(2)
!!    n_trap=n_exci-n_cool
!!   
!!    write(stdout,'(2x,A,f15.7)') "number of electron excited:    ", n_exci
!!    write(stdout,'(2x,A,f15.7)') "number of electron trapped:    ", n_trap
!!    write(stdout,'(2x,A,f15.7)') "number of electron trapped:    ", trap_rho(1)
!!    write(stdout,'(2x,A,f15.7)') "number of electron cooled down:", n_cool
!!    write(stdout,'(2x,A,f15.7)') "tapping efficiency:            ", n_trap/n_exci
!!  
!!    write(stdout,'(2x,A)') "======================"
!!    write(stdout,'(2x,A)') "!  propagation ends  !"
!!    write(stdout,'(2x,A)') "======================"
!!    
!!    deallocate( rho_d,ztmp0,ztmp1 )
!!    
!!    end subroutine propagation
!!  
    !==============================================================!
    ! derivative of DM in eigen spaces                             !
    ! initial condition is calculated from integration of          !
    ! Green's function                                             !
    !                                                              !
    ! EOM is                                                       !
    ! i\partial_t \rho = [H(k),\rho(k)] - Q[\rho(k)]               !
    !                                                              !
    ! coherent part                                                !
    !                                                              !
    ! [H(k),\rho(k)], H(k) is diagonal,                            !
    ! ==> [H(k),\rho(k)]_{ji}= (E_j-E_i)\rho_{ji}                  !
    !                                                              !
    !==============================================================!
    subroutine derivtr(tt,fnk_d)
!!    use w90_parameters,    only : num_wann,nfermi,fermi_energy_list
!!    use w90_postw90_common,only : get_occ_telec
!!    use w90_eph
    implicit none
    real(kind=dp),    intent(in)  :: tt
    real(kind=dp),    intent(out) :: fnk_d(num_wann,nkf)  ! derivative of density matrix
    !
    real(kind=dp),   allocatable :: phi(:,:)
    complex(kind=dp), allocatable :: zmat0(:,:)
    !
    integer          :: i, j,k,kq
    real(kind=dp)    :: Ediff, kweight
    complex(kind=dp) :: ztmp(3)
  
    fnk_d = zero
  
    !----------------------------------------------------
    ! coherent part, [h,\rho]=E_j\rho_{ji}-\rho_{ji}E_i |
    !                        =(E_j-E_i)\rho_{ji}        |
    ! = 0 for BTE                                       !
    !----------------------------------------------------
  
  
    !----------------------------------------------------
    ! incoherent part                                   |
    !----------------------------------------------------
    allocate(phi(num_wann,nkf))
!!    allocate(zmat0(num_wann,num_wann))
!!  
    if(lepi) then
      !----------------------------------------------------
      ! calculate dissipation matrix by EOMs method       |
      ! sum over all phonon modes                         |
      !----------------------------------------------------
!!      call dissipative(nkf,nkf,xkf,eig_k,fnk,phi)
  
      do k=1,nkf
        fnk_d(:,k)=fnk_d(:,k)+phi(:,k)
      enddo
    endif
  
    if(rlx_approx) then
      if(.not.ltrap) then
        ! analytic solution, rho(t) = rho_eq + drho(-t/tau)
        ! rho_eq = occ
!!        do k = 1, nkf
!!          ! get occ 
!!          call get_occ_telec(eig_k(:,k),occ,fermi_energy_list(1),telec)
!!          do i = 1, num_wann
!!             rho_d(i,i,k) = drho(i,i,k)*exp(-tt*relax_time(i,k)*2.d0/hbar) + occ(i)
!!            !rho_d(i,i,k) = rho_d(i,i,k) - ci*(rho(i,i,k)-occ(i))/relax_time(i,k)
!!          enddo
!!        enddo
        return
      !else
      !  drho=rho 
      !  do k=1,nkf
      !    call get_occ_telec(eig_k(:,k),occ,fermi_energy_list(1),telec)
      !    do i = 1, num_wann
      !      drho(i,i,k) = drho(i,i,k) - occ(i)*cone
      !      rho_d(i,i,k) = rho_d(i,i,k) - ci*drho(i,i,k)*relax_time(i,k)*2.d0
      !    enddo 
      !  enddo
      endif
    endif
  
!!    kweight=1.0_dp/real(kmesh(1),dp)/real(kmesh(2),dp)/real(kmesh(3),dp)
!!  
!!    call eeself(num_wann,nkf,xkf,eig_k,rho,kweight,eSelf)
!!    do k=1,nkf
!!      do i=1,num_wann
!!        rho_d(i,i,k) = rho_d(i,i,k) + ci * eSelf(i,k)
!!      enddo
!!    enddo
!!  
!!    if(.true.) then
!!      call dissipative2(nkf,eig_k,1.0_dp,rho,phi)
!!      do k=1,nkf
!!        rho_d(:,:,k) = rho_d(:,:,k) + (phi(:,:,k)-dconjg(transpose(phi(:,:,k)))) * relax_time(1,k)
!!      enddo
!!    endif
!!  
!!    if(ltrap) then
!!      if(tt<0.02_dp) write(stdout,*) 'test, trap_rate=', trap_rate
!!      do k=1,nkf
!!        zmat0 = rho(:,:,k) - rho0(:,:,k)
!!        !call zgemm('n','n',num_wann,num_wann,num_wann,cone,trap(:,:,k),num_wann,drho(:,:,k),num_wann,czero,phi(:,:,k),num_wann)
!!        call zgemm('n','n',num_wann,num_wann,num_wann,cone,trap(:,:,k),num_wann,zmat0,num_wann,czero,phi(:,:,k),num_wann)
!!        rho_d(:,:,k)=rho_d(:,:,k)-ci*(phi(:,:,k)+dconjg(transpose(phi(:,:,k))))
!!      enddo
!!    endif
!!  
!!    trap_rho(2)=0.d0
!!    do k=1,nkf
!!      rho_d(:,:,k)=-ci*rho_d(:,:,k)/hbar
!!      do i=1,num_wann
!!        trap_rho(2)=trap_rho(2)+2.d0*dble(phi(i,i,k))/hbar
!!      enddo
!!    enddo
  
    deallocate(phi,zmat0)
    ! 
    end subroutine derivtr

  
!!    !-------------------------------
!!    ! print distributions          !
!!    !-------------------------------
!!    subroutine printfile(tt,suffix,fu,np,list,var1,var2,iop)
!!    use w90_io,         only : stdout,io_file_unit,io_error,seedname
!!    use w90_parameters, only : fermi_energy_list
!!    implicit none
!!    real(kind=dp),    intent(in) :: tt
!!    character(len=*), intent(in) :: suffix
!!    integer,          intent(in) :: fu,np
!!    real(kind=dp),    intent(in) :: list(np),var1(np),var2(np)
!!    integer,          intent(in) :: iop
!!    
!!    character(len=24) :: file_name
!!    integer           :: file_unit,i
!!  
!!    if(iop==3) then
!!      file_unit=fu
!!      file_name=suffix
!!      open(file_unit,FILE=file_name,FORM='FORMATTED',position="append")
!!    else
!!      file_name=trim(seedname)//suffix
!!      file_unit=io_file_unit()
!!      open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
!!      write(file_unit,'(A,E16.8)') 'E-fermi:',fermi_energy_list(1)
!!    endif
!!  
!!    write(stdout,'(/,3x,a)') '* '//file_name
!!  
!!    do i=1,np
!!      if(iop==1) then
!!        write(file_unit,'(2E16.8)') list(i),var1(i)
!!      else if(iop==2) then
!!        write(file_unit,'(3E16.8)') list(i),var1(i),var2(i)
!!      else if(iop==3) then
!!        write(file_unit,'(f9.3,2x,3E16.8)') tt,list(i),var1(i),var2(i)
!!      endif
!!    enddo
!!  
!!    if(iop==3) write(file_unit,*) 
!!    close(file_unit)
!!  
!!    end subroutine printfile
!!  
    !------------------------------------------------------------------------
    ! |M(k1,k2,k3,k4)|^2 = e^4/[q^4 \epslon^2_0 \epsilon^2_r] 
    !
    ! S^pm_e(k_1) = 2pi/hbar * \sum_{k1,k2,k3} M(k1,k2,k3,k4) x 
    ! (0.5 \pm 0.5 - f(k2)) x [0.5 \mp 0.5 - f(k3)] x [0.5 \mp 0.5 - f(k4)]
    !------------------------------------------------------------------------
!!    subroutine eeself(num_wann,nkf,xkf,eig_k,rho,kweight,selfp,selfm)
!!    implicit none
!!    integer, intent(in) :: num_wann,nkf
!!    real*8,  intent(in) :: xkf(3,nkf)
!!    real*8,  intent(in) :: eig_k(num_wann,nkf)
!!    real*8,  intent(in) :: rho(num_wann,num_wann,nkf)
!!    real*8,  intent(in) :: kweight
!!    real*8,  intent(out):: selfp(num_wann,nkf)
!!    real*8,  intent(out):: selfm(num_wann,nkf)
!!    !
!!    integer :: i, j, m, n
!!    integer :: k1, k2, k3, k4
!!    integer :: ist, jst, kst, lst
!!    real*8  :: ees(2), dtmp, dk(3), dE, diele, mk4
!!    real*8  :: traceSelf
!!    
!!    selfp = zero
!!    selfm = zero
!!    diele = 1.d0
!!  
!!    traceSelf = 0.d0
!!
!!    do i = 1, nee
!!       k1 = eeklist(1,nee) 
!!       k2 = eeklist(2,nee) 
!!       k3 = eeklist(3,nee) 
!!       k4 = eeklist(4,nee) 
!!           
!!       ist = eeklist(5,nee) 
!!       jst = eeklist(6,nee) 
!!       kst = eeklist(7,nee) 
!!       lst = eeklist(8,nee)
!!
!!       dk = xkf(:,k1) - xkf(:,k3)
!!       mk4 = dot_product(dk, dk)
!!       mk4 = mk4 * mk4
!!       if (mk4 < 1.e-4_dp) cycle
!!
!!       !call dielectric(dk,dE,diele))
!!       mk4 = one / mk4 / diele
!!
!!       selfp(ist,k1) = selfp(ist,k1) + mk4 * (one - fnk(jst,k2)) * fnk(kst,k3) * fnk(lst,k4)
!!       selfm(ist,k1) = selfm(ist,k1) - mk4 * fnk(jst,k2) * (one - fnk(kst,k3)) * (one - fnk(lst,k4))
!!    enddo
!!
!!    end subroutine eeself
!!  
!!  
!!    subroutine eecore(iband,jband,k1,k3,nkf,num_wann,xkf,eig_k,kweight,ees)
!!    implicit none
!!    integer,  intent(in) :: iband, jband, k1,k3, nkf
!!    integer,  intent(in) :: num_wann
!!    real*8,   intent(in) :: xkf(3,nkf)
!!    real*8,   intent(in) :: eig_k(num_wann,nkf)
!!    real*8,   intent(in) :: kweight
!!    real*8,   intent(out):: ees(2)
!!    !
!!    integer :: k2,k4,m,n
!!    real*8  :: crit
!!    real*8  :: dk(3), dE, vec(3), dtmp
!!   
!!    ees = 0.d0
!!  
!!    dk = xkf(:,k1) - xkf(:,k3) ! k1 - k3
!!    dE = eig_k(iband,k1) - eig_k(jband,k3)
!!  
!!    crit = 1.d-4
!!  
!!    do k2=1, nkf
!!      do k4=1,nkf
!!        ! k2 - k4
!!        vec = xkf(:,k2) - xkf(:,k4) + dk
!!        dtmp = sqrt(dot_product(vec,vec))
!!        if(dtmp<crit) then
!!          do m=1,num_wann
!!            do n=1,num_wann
!!              dtmp = eig_k(m,k2) - eig_k(n,k4) + dE
!!              if(dabs(dtmp)<crit) then
!!                ees(1) = ees(1) + (1.d0-fnk(m,k2)) * fnk(n,k4)
!!                ees(2) = ees(2) + (1.d0-fnk(n,k4)) * fnk(m,k2)
!!              endif
!!            enddo
!!          enddo
!!        endif
!!        !
!!      enddo
!!    enddo
!!  
!!    ees = ees * kweight * kweight 
!!  
!!    end subroutine eecore
!!  
!!    !-------------------------------------
!!    !                                    !
!!    !-------------------------------------

    subroutine ee_klist(nkf,nkftot, num_wann,xk,eig_k,xk_all,eig_kall,eeklist)
    implicit none
    integer,       intent(in) :: nkf, nkftot, num_wann
    real(kind=dp), intent(in) :: xk(3,nkf)
    real(kind=dp), intent(in) :: eig_k(num_wann,nkf)
    real(kind=dp), intent(in) :: xk_all(3,nkftot)
    real(kind=dp), intent(in) :: eig_kall(num_wann,nkftot)
    integer,    intent(inout):: eeklist(8,*)
  
    real*8  :: vec(3),dist,critk, crite
    integer :: k1,k2,k3,k4
    integer :: ist, jst, kst, lst, ijst, klst
  
    nee=0
  
    critk = 1.d-6
    crite = 1.d-6
 
!    write(6,*) 'nkf, pool_id', nkf, my_pool_id, npool
    do k1 = 1,nkf
    do k2 = 1, nkftot
      if(k2==k1) cycle
      do k3 = 1, nkftot
      do k4 = 1, nkftot
        if(k3==k4) cycle
        vec = xk(:,k1) + xk_all(:,k2) - xk_all(:,k3) - xk_all(:,k4)
        dist=sqrt(dot_product(vec,vec))
  
        if(dist<critk) then
          write(stdout,'(4(I5,2x),e15.6)') k1,k2,k3,k4, dist
          call flush(stdout)
        endif

        if(dist<critk) then
          ijst = 0
          do ist=ibndmin,ibndmax
          do jst=ibndmin,ibndmax
            if (jst > ist) cycle
            ijst = ijst + 1

            klst = 0
            do kst=ibndmin,ibndmax
            do lst=ibndmin,ibndmax
               if (lst > kst ) cycle
               klst = klst + 1
               if (klst > ijst) cycle
               
               if(ist==jst.and.ist==kst.and.ist==lst) cycle
               if(ist==kst.and.k1==k3.and.jst==lst.and.k2==k4) cycle

               if( dabs(eig_k(ist,k1)-eig_kall(jst,k2)) > 3.0_dp) cycle
               if( dabs(eig_k(ist,k1)-eig_kall(kst,k3)) > 5.0_dp) cycle
               if( dabs(eig_kall(jst,k2)-eig_kall(lst,k4)) > 5.0_dp) cycle

               !if(dabs(eig_k(ist,k1)-eig_k(kst,k3))<crite) cycle
               dist=eig_k(ist,k1)+eig_kall(jst,k2) - eig_kall(kst,k3) - eig_kall(lst,k4)
               if(abs(dist)<crite) then
                 nee = nee + 1
                 write(stdout,'(I6,1x,8I4,3e15.6)') nee, k1,k2,k3,k4,ist,jst,kst,lst,eig_k(ist,k1)-eig_kall(kst,k3),&
                                                    eig_kall(ist,k1),eig_kall(lst,k4)
                 eeklist(1,nee) = k1
                 eeklist(2,nee) = k2
                 eeklist(3,nee) = k3
                 eeklist(4,nee) = k4
  
                 eeklist(5,nee) = ist
                 eeklist(6,nee) = jst
                 eeklist(7,nee) = kst
                 eeklist(8,nee) = lst
               endif
            enddo
            enddo
          enddo
          enddo
        endif  ! check k1+k2-k3-k4=0
      enddo
      enddo
    enddo
    enddo  
 
    write(6,'(A,I7,A)') 'totally', nee, ' iterms in the ee klist!'
  
    end subroutine ee_klist
  
  

  !-----------------------------------------------------------------------
  SUBROUTINE relaxbte ( )
    implicit none
  !-----------------------------------------------------------------------
  !!  compute hot-elecrtron relaxation dynamics by employing the ab-initio
  !!  boltzman equation, where electron-electron and electron-phonon scattering
  !!  are considered. 
  !! 
  !!  YZ, Mar. 2020, at Los Alamos National Laboratory
  !! 
  !-----------------------------------------------------------------------

 
  END SUBROUTINE relaxbte



!------------------------------------------------
  
    subroutine dielectric(model, metal, omega, epsilonw)
    !==============================================================!
    !                                                              !
    ! Reference: Applied Optics 37, 5271 (1998).                   !
    !                                                              !
    ! A. Lorentz-Drude Model                                       !
    !                                                              !
    ! epsilon(w) = epsilon_f(w) + epsilon_b(w)                     !
    !                                                              !
    ! This separates explicitly the intraband effects (ususally    !
    ! referred to as free electron effects) from interband effects !
    ! (ususally referred to as bound electron effects).            !
    !                                                              !
    !                             \Omega^2_p                       !
    ! Where, epsilon_f(w) = 1 - --------------,                    !
    !                           w ( w - i G_0 )                    !
    !                                                              !
    ! where \Omega_p = \sqrt{f_0}*omegap                           !
    !                                                              !
    !                               f_i w^2_p                      !
    ! epsilon_b(w) = \sum_j  ------------------------              !
    !                        (w^2_j - w^2) + i w G_j               !
    !                                                              !
    ! B. Brendel-Bormann Model                                     !
    !                                                              !
    !                     \Omega^2_p                               !
    ! epsilon_f(w) = 1 - ------------- + \sum_j x_j(w),            !
    !                     w(w - i G_0)                             !
    ! where                                                        !
    !                                                              !
    ! x_j(w) =                                                     !
    !                                                              !
    !==============================================================!
    implicit none
    !
    integer,           intent(in) :: model     ! choice of dielectric models
    character*2,       intent(in) :: metal     ! choice of metal
    real(kind=dp),     intent(in) :: omega     ! frequency 
    complex(kind=dp),  intent(out):: epsilonw  ! dielectric at given frequency
    !
    integer, parameter :: nLB = 5
    !
    real(kind=dp)  :: dtmp
    real(kind=dp)  :: wp, f0, G0, omegaj(nLB), fj(nLB), Gj(nLB)
    integer :: i, j, k
    
    epsilonw = ci
    
    ! parameters
    !
    if ( metal=='Ag' ) then
     wp =  9.01_dp 
     f0 = 0.845_dp
     G0 = 0.048_dp
     data fj(1:nLB)     /0.065_dp, 0.124_dp, 0.011_dp, 0.840_dp, 5.646_dp/
     data Gj(1:nLB)     /3.886_dp, 0.452_dp, 0.065_dp, 0.916_dp, 2.419_dp/
     data omegaj(1:nLB) /0.816_dp, 4.481_dp, 8.185_dp, 9.083_dp, 20.29_dp/
    else if( metal == 'Au' ) then
     wp = 9.03_dp
     f0 = 0.760_dp
     G0 = 0.053_dp
     data fj(1:nLB)     /0.024_dp,0.010_dp,0.071_dp,0.601_dp,4.384_dp/
     data Gj(1:nLB)     /0.241_dp,0.345_dp,0.870_dp,2.494_dp,2.214_dp/
     data omegaj(1:nLB) /0.415_dp,0.830_dp,2.969_dp,4.304_dp,12.32_dp/
    else
      write(stdout,*) 'warning: metal is not fund!'
      stop
    endif
    !
    if(model==1) then
      !write(stdout,*) 'test1'
      epsilonw = cone - f0 * wp * wp / omega / (omega - ci * G0) 
      !do i = 1, nLB
      !  epsilonw = epsilonw + fj(i) * wp * wp / (omegaj(i)*omegaj(i) - omega * omega + ci * omega * Gj(i) )
      !enddo
    else
    
    endif
    !
    end subroutine dielectric


    !-------------------------------
    !         print Header         !
    !-------------------------------
    subroutine printHeader(funit)
    implicit none
    !
    integer, intent(in) :: funit
  
    character*10 :: date,time
    character(len=*), parameter :: Header = "(&
            '                                                                                '/&
            '********************************************************************************'/&
            '**                                                                            **'/&
            '**                               PHD-K                                        **'/&
            '**                         Copyright (c)  Yu Zhang, PhD,                      **'/&
            '**                           Email: zhyhku@gmail.com                          **'/&
            '**                                                                            **'/&
            '**      Theoretical Division, Los Alamos National Laboratory, Los Alamos      **'/&
            '**                                                                            **'/&
            '**   This propgram is used to calculate the hot-carrier generation from       **'/&
            '**   generation plasmon decay and its injection to other materials via        **'/&
            '**   the interface. The electron transport is calculated within the NEGF      **'/&
            '**   formalism. The electronic structure employs tight-binding, DFTB or       **'/&
            '**   ab-initio tight-binding (localized wannier function as basis).           **'/&
            '**                                                                            **'/&
            '**   modify the description above!!!!!                                        **'/&
            '**                                                                            **'/&
            '********************************************************************************')"
    !
    ! Print Title and Time
    !
    write(funit,Header)
    call date_and_time(date,time)
    write(funit,1000) date(7:8),date(5:6),date(1:4),time(1:2),time(3:4),time(5:6)
    !
    1000 format( ' Calculation Started on ',' ',A2,"/",A2,"/",A4,"  ",A2,":",A2,":",A2,/)
    end subroutine printHeader
    !
    subroutine printEnding(funit)
    implicit none
    integer, intent(in) :: funit
    !
    character*10 :: date,time
    character(len=*), parameter :: Header = "(&
            '********************************************************************************'/&
            '**                                                                            **'/&
            '**                             Program Ended                                  **'/&
            '**                                                                            **'/&
            '********************************************************************************')"
    !
    call date_and_time(date,time)
    write(funit,1000) date(7:8),date(5:6),date(1:4),time(1:2),time(3:4),time(5:6)
    write(funit,Header)
    !
    1000 format( /,' Calculation Finished on ',' ',A2,"/",A2,"/",A4,"  ",A2,":",A2,":",A2)
    end subroutine printEnding

end module w90_plasmon
