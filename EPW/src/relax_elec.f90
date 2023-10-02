  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            

  !-----------------------------------------------------------------------
  SUBROUTINE relax_elec_new ( )
  !-----------------------------------------------------------------------
  !! 
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, shortrange, &
                            fsthick, ngaussw, degaussw, &
                            wmin_specfun,wmax_specfun,nw_specfun, &
                            eps_acustic, efermi_read, fermi_energy,&
                            restart, restart_freq, &
                            nomega, &
                            rlx_dt, rlx_tmax,Epump,Th_cond, mobility
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, &
                            epf17, wkf, nkf, wf, wqf, xkf, nkqtotf,&
                            esigmar_all, esigmai_all, a_all,  &
                            he_ij,he_all2, he_all, edosef, &
                            edos_all, vdos_all,jdos  !ZY
  USE transportcom,  ONLY : lower_bnd
  USE control_flags, ONLY : iverbosity
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, kelvin2eV, two, zero, hbar, pi
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : inter_pool_comm
  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : stdout, ionode, ionode_id
  !
  implicit none
  !
  ! Local variables 
  !
  REAL(kind=DP) :: Te
  !! electron temperature
  REAL(kind=DP) :: eptemp
  ! temperature for the electronic Fermi occupations in the e-p calculation 
  REAL(kind=DP) :: Tenv ! environment temperature
  Real(kind=dp) :: beta,Tl, Tl1, Tl0, dTl
  !! lattice temperature
  !REAL(kind=DP) :: mobility
  !!! electron mobility
  !
  INTEGER :: n
  !! Integer for the degenerate average over eigenstates
  INTEGER :: iw
  !! Counter on the frequency
  !
  Real(kind=dp) :: omega 
  !! omega is the volume of the primitive cell in a.u. 
  !! in restart mode, omega is not read from the pw save, so, we read it from
  !! volume.dat file!!! remember to prepare this file 
  REAL(kind=DP) :: tt,dt,dt2,dt3,dt6
  REAL(kind=DP) :: tmax
  REAL(kind=DP) :: ww
  REAL(kind=DP) :: dw 
  !! Frequency intervals
  REAL(kind=DP) :: tmp
  !! Temporary variable to store real part of Sigma for the degenerate average
  REAL(kind=DP) :: tmp2
  !! Temporary variable to store imag part of Sigma for the degenerate average
  REAL(kind=DP) :: tmp3
  !! Temporary variable to store Z for the degenerate average
  REAL(kind=DP) :: diff1,diff2
  real(kind=dp) :: Nhot(2)
  !
  real(kind=DP), external :: efermig, dos_ef, wgauss
  real(kind=DP), external :: w0gauss

  !!
  real(kind=dp), allocatable :: fe0(:)  ! f(e) at t=0 (including the excitation)
  real(kind=dp), allocatable :: feq(:)  ! f(e) at equilibrium (temperature given by eptemp)
  real(kind=dp), allocatable :: fe(:)   ! f(e) 
  real(kind=dp), allocatable :: dfe(:)   ! d f(e) / dt
  real(kind=dp), allocatable :: tau(:) ! lifetime
  logical,       allocatable :: gettau(:) 
  real(kind=dp), allocatable :: dvec0(:) 
  real(kind=dp), allocatable :: dvec1(:) 
  !

  if (ionode) call printHeader(stdout)

  allocate(feq(nw_specfun),fe0(nw_specfun),fe(nw_specfun),dfe(nw_specfun))
  allocate(dvec1(nw_specfun),dvec0(nw_specfun))
  allocate(gettau(nw_specfun),tau(nw_specfun))

  Te = 300.0_dp
  Te = Te * kelvin2eV / ryd2ev
  Tl = Te
  Tenv = Te
  eptemp = Te

  open(20,file='volume.dat',status='old')
  read(20,*) omega
  write(6,*) 'unit-cell volume=', omega
  close(20)

  !Th_cond = 1.0_dp ! from input
  !Th_cond = 5x10^14 W/m^3/K
  ! = 5x10^14 / (1.602x10^{-19}) ev / 10^15 fs / ((1/0.0529)^3x10^27 au^3) K
  ! = 5/1.502 10^-9 /(1/0.0529)^3 ev/au^3/K
  ! since Tl absorbed Kb, we need to divide kB here
  ! Th_cond = 5/1.602 *(0.0529)^3 10^-9 / (8.617333262145 x 10^-5)  /au^3
  !         = 5/1.602 *(0.0520)^3 / 8.617333262145 * 10^-4
  !
  !Th_cond = 5.0_dp/1.602 * 0.0529_dp**3 / 8.617333262145_dp * 1.e-4

  ! initial distribution
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1)

  DO iw = 1, nw_specfun
    !
    ww = wmin_specfun + dble (iw-1) * dw
    feq(iw) = wgauss(-ww/eptemp, -99) 
    fe(iw) = wgauss(-ww/eptemp, -99) 
    !
  ENDDO


  ! rescale the maximum value to 1.0
  tmp = max( maxval(jdos(:,1)/edos_all),maxval(jdos(:,2)/edos_all) ) 
  jdos = jdos/tmp 
  IF (mpime.eq.ionode_id) THEN
    ! Write to file
    write(stdout,*) 'total hot e/h'
    write(stdout,*) sum(jdos(:,1)), sum(jdos(:,2))
    write(stdout,*) sum(jdos(:,1)/edos_all), sum(jdos(:,2)/edos_all)
    open(666,file='jdos.dat',status='replace')
    open(667,file='vdos.dat',status='replace')
    write(666,'(A9,5A15)') 'energy','jdos-e', 'jdos-h','dos','f(t=0)', 'he(w)'
  ENDIF

  !--------------------------
  ! hot-carrier excitation
  !--------------------------
  do iw = 1, nw_specfun
     ww = wmin_specfun + dble(iw-1) * dw
     ! hot
     dfe(iw) = jdos(iw,1) / edos_all(iw) * Epump - & ! electron
               jdos(iw,2) / edos_all(iw) * Epump     ! hole
     fe(iw) = fe(iw) + dfe(iw)
     IF (mpime.eq.ionode_id) THEN
        write(666,'(f9.4, 5e16.7)') ww, jdos(iw,1), jdos(iw,2),edos_all(iw),& !dfe(iw),
                                   fe(iw), he_all(iw)
     ENDIF
  enddo

  IF (mpime.eq.ionode_id) THEN
     do iw = 1, nomega
       ww =  dble(iw-1) * dw
       write(667,'(f9.4, 1e16.7)') ww, vdos_all(iw)
     enddo
  ENDIF
    
  IF (mpime.eq.ionode_id) close(666)
  IF (mpime.eq.ionode_id) close(667)

  !DO iw = 1, nw_specfun
  !  ww = wmin_specfun + dble (iw-1) * dw
  !  write(stdout,'(I5,2e15.6)') iw, ww, fe(iw)
  !ENDDO


  !
  tmp = 2.0 * pi / (hbar * 1.0e15 / ryd2ev)
  write(stdout,*) 'scaling factor for he_all=', tmp

  !he_all = he_all * omega
  he_all = he_all * 2.0 * pi / (hbar * 1.0e15 /ryd2ev)
  he_all2= he_all2* 2.0 * pi / (hbar * 1.0e15 /ryd2ev)
  he_ij = he_ij * pi / (hbar * 1.0e15 /ryd2ev)

  ! mobility (eV^-1) 
  mobility = mobility / ryd2ev
  !mobility = 0.021_dp / ryd2ev

  !--------------------
  ! relaxation
  !--------------------

  tmax = rlx_tmax
  dt   = rlx_dt
  write(stdout,*) 'time-step=', dt
  tt = 0.0_dp

  dt2 = dt / 2.0_dp
  dt3 = dt / 3.0_dp
  dt6 = dt / 6.0_dp

  tau    = 0.e0_dp
  gettau = .false.
  fe0    = fe ! salve fe at t=0

  IF (mpime.eq.ionode_id) THEN
     open(999,file='occ.dat',status='replace')
  endif
  !tmax = 1.e4_dp
  !
  do
    !
    if ( mod(nint(tt/dt),nint(1.0_dp/dt)) == 0 ) then
      ! fit fe with effective Te
      call getTe(nw_specfun,dw,wmin_specfun,fe,Te)
      Nhot = 0.0_dp
      do iw = 1, nw_specfun
         ww = wmin_specfun + dble(iw-1) * dw
         if (ww>0.0_dp) Nhot(1) = Nhot(1) + (fe(iw) - feq(iw)) * edos_all(iw) * dw
         if (ww<0.0_dp) Nhot(2) = Nhot(2) - (fe(iw) - feq(iw)) * edos_all(iw) * dw

      enddo

      write(stdout,'(A,f9.2,1x,4e16.7,A)') 'current time=', tt, Nhot,Tl*ryd2ev/kelvin2eV, Te*ryd2ev/kelvin2eV,' K'

      IF (mpime.eq.ionode_id) THEN
        write(999,'(A,f9.2,1x,4e16.7,A)') 'current time=', tt, Nhot,Tl*ryd2ev/kelvin2eV, Te*ryd2ev/kelvin2eV,' K'
        DO iw = 1, nw_specfun
          ww = wmin_specfun + dble (iw-1) * dw
          write(999,'(I5,f10.4,e20.8)') iw, ww, fe(iw) !,dfe(iw)
        ENDDO
      endif
      !
    endif

    !---------------------------
    ! calculate the lifetime
    !---------------------------
    do iw = 1, nw_specfun
      ww = wmin_specfun + dble(iw-1) * dw
      diff1 = abs( fe(iw) - feq(iw) ) 
      diff2 = abs( fe0(iw) - feq(iw) ) * exp(-1.0_dp)
      ! stop updating the lifetime once diff1 < diff2
      if (diff1 > diff2) then 
        tau(iw) = tt
      else 
        if (.not.gettau(iw) ) then
          gettau(iw) = .true.
          write(stdout,'(A,f9.3,A,e15.6)') 'lifetime of energy ', ww*ryd2ev, " is ", tau(iw)
        endif
      endif
    enddo

    call dfedt(tt,nw_specfun,dw,wmin_specfun,mobility,tl,Tenv, fe,dfe,dTl)
    dvec0 = fe
    fe    = dvec0 + dt2 * dfe
    dvec1 = dvec0 + dt6 * dfe
    !
    Tl0   = Tl
    Tl    = Tl0   + dt2 * dTl
    Tl1   = Tl0   + dt6 * dTl

    call dfedt(tt+dt2,nw_specfun,dw,wmin_specfun,mobility,tl,Tenv, fe,dfe,dTl)
    fe    = dvec0 + dt2 * dfe
    dvec1 = dvec1 + dt3 * dfe
    !
    Tl    = Tl0   + dt2 * dTl
    Tl1   = Tl1   + dt3 * dTl

    call dfedt(tt+dt2,nw_specfun,dw,wmin_specfun,mobility,tl,Tenv, fe,dfe,dTl)
    fe    = dvec0 + dt  * dfe
    dvec1 = dvec1 + dt3 * dfe
    !
    Tl    = Tl0   + dt  * dTl
    Tl1   = Tl1   + dt3 * dTl
   
    call dfedt(tt+dt,nw_specfun,dw,wmin_specfun,mobility,tl,Tenv, fe,dfe,dTl)
    fe   = dvec1 + dt6 * dfe
    !
    Tl   = Tl1   + dt6 * dTl

    if(tt>=tmax) exit
    tt=tt+dt

    !exit
    !
  enddo

  !
  if(ionode) then
  DO iw = 1, nw_specfun
    ww = wmin_specfun + dble (iw-1) * dw
    write(stdout,'(I5,4e15.6)') iw, ww, fe(iw),fe(iw)-feq(iw),tau(iw)
  ENDDO
  endif


  if(ionode) call printEnding(stdout)
  !
  deallocate(fe,dfe,fe0,feq,gettau,tau,dvec0,dvec1)
  RETURN
  !
  END SUBROUTINE relax_elec_new
  

  subroutine dfedt(tt,np, dw, ebot, mobility, tl, Tenv, fe,dfe,dTl)
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, kelvin2eV, two, zero, pi, hbar, ci
  USE epwcom,        ONLY : ngaussw, degaussw, eps_acustic,nomega, alpha_heat,Th_cond
  uSE elph2,         ONLY : vdos_all,he_ij,he_all,he_all2,edosef,edos_all
  implicit none
  real(kind=dp), intent(in)  :: tt
  integer,       intent(in)  :: np
  real(kind=dp), intent(in)  :: dw, ebot
  real(kind=dp), intent(in)  :: mobility,tl, Tenv
  real(kind=dp), intent(in)  :: fe(np)
  real(kind=dp), intent(out) :: dfe(np),dTl
  !
  integer       :: i, j, k,l, ijklmax
  integer       :: dij
  real(kind=dp) :: ei,ej,ek,el,Ediff
  real(kind=dp) :: delta0, delta, tmp0, tmp,nocc
  real(kind=dp) :: hbar_ryfs
  real(kind=dp) :: dfe0(np),dfe1(np)
  real(kind=dp) :: demax  ! maximal |ei+ej-ek-el| allowed in the calculaiton
                          ! larger |ei+ej-ek-el| is, smaller contribution to the
                          ! sum
  real(kind=dp) :: Eng2l  ! energy transferred from electron to lattice
  real(kind=dp) :: Eng2env! energy transferred from lattice to enviroment
  real(kind=dp) :: heat_capacity
  REAL(kind=DP), EXTERNAL :: w0gauss
  real(kind=DP), external :: wgauss
  !
  logical       :: lapprox
  !
  dfe = 0.0_dp
  demax = 5.0*degaussw
  ijklmax = nint(demax/dw)

  !    (reduced hbar)  to (ev*fs)  to (ryd * fs)
  hbar_ryfs = hbar  * 1.0e15 / ryd2ev
  dij = nint(0.5_dp/ryd2ev/dw)
  if (dij < 5) dij = 5

  ! e_i = ebot + (i-1) * dw
  ! so, e_i + e_j - e_k - e_l = 
  !   (i+j-k-l) * dw
  ! |(i+j-k-l) * dw| <= demax -->
  ! |i+j-k-l| <= (demax/dw) = ijlkmax
  ! ==>
  ! i+j = k + l \pm ijklmax
  ! or 
  ! k + l = i + j \mp ijklmax
  ! i.e., |k| <= |i+j+ijklmax|
  ! l = i+j + ijklmax - k
  !

  dfe  = 0.0_dp

  !write(stdout,*) 'test-delta function'
  !tmp = 0.0_dp
  !do i = 1, np
  !  ei  = ebot + dble(i-1) * dw

  !  delta = w0gauss(ei/degaussw, 0)/degaussw
  !  tmp = tmp + delta * dw
  !  write(stdout,'(A,I4,3e15.6)') 'i, ei, delta, tmp', i, ei, delta, tmp
  !enddo

  dfe0 = 0.0_dp
  dfe1 = 0.0_dp
  !
  ! a falst algorithim is need to accelerate the four indelx integral
  ! 1) F_ee
  ! paralle, to be done!!!
  !write(stdout,*) 'test-zy: e-e'
  do i = 1, np
    ei = ebot + dble (i-1) * dw

    do j = 1, np
      ej = ebot + dble (j-1) * dw

      !(k+l) = i+j, so k < i+j
      do k = 1, min(np,i+j+ijklmax)
        ek = ebot + dble (k-1) * dw

        do l = max(1,i+j-ijklmax-k), min(np,i+j+ijklmax-k)
          el = ebot + dble (l-1) * dw
          ! get delta(e1+e2-e3-e4)
          Ediff = ei + ej - ek - el
          !delta0= w0gauss(Ediff/degaussw, 0)/degaussw
          !delta0= aimag(1.0_dp / (Ediff - ci * 1.0) )
          delta0= aimag(1.0_dp / (Ediff - ci * degaussw) )
          delta = delta0 * edos_all(j) * edos_all(k) * edos_all(l) / (edosef**3.0)

          dfe0(i) = dfe0(i) + delta*(fe(k)*fe(l) * (one-fe(i))*(one-fe(j)) - &
                                     fe(i)*fe(j) * (one-fe(k))*(one-fe(l))) 
          !if (mod(i,20)==0.and.mod(j,50)==0) write(stdout,'(4I5,4e20.9)') l,k,j,i,Ediff,delta0, delta, dfe0(i)
        enddo
      enddo
      !if (mod(i,10)==0.and.mod(j,10)==0) write(stdout,'(2I5,2e20.9)') j,i,dfe0(i)
    enddo
    
  enddo
  
  ! mobility : 1/E
  ! dw*dw*dw : E^3
  ! delta    : 1/E
  ! hbar     : 1/(E*Time)
  ! 2*De/hbar

  dfe = dfe0 * dw * dw *dw * 2.0 * mobility / hbar_ryfs

  !tmp = 0.0_dp
  !tmp0= 0.0_dp
  !do i = 1, np
  !  tmp = tmp + dfe(i)
  !  tmp0= tmp0+ dfe(i) * edos(i)
  !  write(stdout,'(I5,3e20.9)') i, edos(i), dfe0(i), dfe(i)
  !enddo
  !write(stdout,*) 'trace of dfe, (e-e)=', tmp,tmp0

  !dfe = 0.0_dp

  !------------
  ! e-phonon
  !------------

  lapprox = .true.
  lapprox = .false.

  !write(stdout,*) 'test e-ph'
  Eng2l = 0.0_dp
  do i = 1, np
    ei = ebot + dble (i-1) * dw

    !-------------------------------------- 
    ! without approximation
    !-------------------------------------- 

    if(.not.lapprox) then
      !do j = max(1,i-dij), min(np,i+dij)
      do j = 1, nomega
        !ej = ei + dble (j-1) * dw
        !ej = ei - dble (j-1) * dw
        Ediff = dble(j-1) * dw
        !Ediff = dble(j-i) * dw
      
        if (Ediff < dw) cycle
        if (Ediff < eps_acustic) cycle
      
        nocc = wgauss(-Ediff/Tl,-99)      ! 1/(exp(x)+1)
        ! one - two * occ = [exp(x) - 1] / [exp(x) + 1]
        nocc = nocc / (one - two * nocc)

        !if(i==1) write(stdout,*) 'Ediff,nocc', Ediff, nocc
      
        ! 1). ei > ej
        !if (ei>ej) then
        !if ( i>j ) then
        if ( (i-j+1) >= 1) then
          ! ej + w --> ei, .i.e, \delta(ei-ej-w)
          tmp = he_all(i-j+1) !* edos_all(i-j+1)/edosef
          !tmp = he_all(i) * he_all(i-j+1) * dw / edosef !edos_all(i-j+1)/edosef * dw
          tmp = he_ij(i,i-j+1) * dw / edosef !edos_all(i-j+1)/edosef * dw
          tmp = he_ij(i,i-j+1) * dw
          dfe1(i) = dfe1(i) - (fe(i) * (one - fe(i-j+1)) * (nocc+1.0_dp) - &
                               fe(i-j+1) * (one - fe(i)) * nocc) * tmp !* vdos_all(j) 
          !dfe1(i) = dfe1(i) - (fe(i) * (one - fe(j)) * (nocc+1.0_dp) - &
          !                     fe(j) * (one - fe(i)) * nocc) * he_ij(j,i) * dw 
                               !fe(j) * (one - fe(i)) * nocc) * he_all2(i)
        endif
        !else

        !if ( i<j ) then
        if ( (i+j-1) <= np) then
          ! ei + w --> ej, .i.e. \delta(ei-ej+w)
          !tmp = he_all(i+j-1) !* edos_all(i+j-1)/edosef
          !tmp = he_all(i) * he_all(i+j-1) * dw / edosef !edos_all(i+j-1)/edosef * dw
          tmp = he_ij(i+j-1,i) * dw / edosef !edos_all(i+j-1)/edosef * dw
          tmp = he_ij(i+j-1,i) * dw 
          dfe1(i) = dfe1(i) + (fe(i+j-1) * (one - fe(i)) * (nocc+1.0_dp) - &
                               fe(i) * (one - fe(i+j-1)) * nocc) * tmp !* vdos_all(j) 
          !dfe1(i) = dfe1(i) + (fe(j) * (one - fe(i)) * (nocc+1.0_dp) - &
          !                     fe(i) * (one - fe(j)) * nocc) * he_ij(i,j) * dw 
        endif
      enddo
      !\int e dos(e) * df/dt_{eph} de
      Eng2l = Eng2l - ei * dfe1(i) * dw
      !Eng2l = Eng2l + dfe1(i) * dw !* Ediff

    else

    !-------------------------------------- 
    ! approximaiton with Talor expansion
    !-------------------------------------- 
    ! = 1/g(e) * d/de[H(e)(f(e)(1-fe)+ KbTl df/de]
    ! f(e) = 1/(exp(beta*(e-mu)) + 1.0)
    ! df/de = [exp(eta*(e-mu))+1]^{-2} * exp(beta*(e-mu)) * beta

      ! finite-difference calculation of d f/de
      if (i==1) then
        tmp = (fe(i+1) - fe(i)) / dw
      else if (i==np) then
        tmp = (fe(i) - fe(i-1)) / dw
      else
        tmp = (fe(i+1) - fe(i-1)) / dw / 2.0
      endif
      dfe1(i) = (fe(i) * (one - fe(i)) + Tl * tmp) * he_all(i) 

      Eng2l = Eng2l + dfe1(i) * dw
      !write(stdout,'(I5,3e20.8)') i, tmp, dfe1(i)
    endif
  enddo

  ! e-e part is correct, confirmed, May 11
  !dfe = 0.0_dp

  ! C_l(T_l) dT_l/dt = Eng2l
  ! calculate heat capacity
  call heatcapacity(nomega, dw, 0.0_dp, Tl, vdos_all, heat_capacity)

  Eng2env = Th_cond * (Tl - Tenv)
  
  dTl = alpha_heat * (Eng2l - Eng2env) / heat_capacity
  if ( mod(nint(tt/0.1),nint(1.0_dp/0.1)) == 0 ) then
    write(stdout,*) 
    write(stdout,'(A,7e15.6)') ' test dTl', Tl, Eng2l, Eng2env, heat_capacity, dTl
  endif
  !dTl = 0.0_dp

  ! finite-difference
  ! d/dw : 1/E
  ! edos : 1/L^3/E ==> d/dw / edos = L^3 
  ! H(e) : 1/Time/L^3
  !
  !write(stdout,*) 'test dfe'
  tmp = 0.0_dp
  tmp0= 0.0_dp
  do i = 1, np
    !-------------------------------------- 
    ! approximaiton with Taylor expansion
    !-------------------------------------- 
    if (lapprox) then

      if (i == 1) then
         dfe(i) = dfe(i) + (dfe1(i+1) - dfe1(i)) / dw / edos_all(i)
      else if (i==np) then
         dfe(i) = dfe(i) + (dfe1(i) - dfe1(i-1)) / dw / edos_all(i)
      else
         dfe(i) = dfe(i) + (dfe1(i+1) - dfe1(i-1)) / dw / 2.0_dp / edos_all(i)
      endif

    else

      !-------------------------------------- 
      ! without approximation
      !-------------------------------------- 
      dfe(i) = dfe(i) + dfe1(i) / edos_all(i)

    endif


    tmp = tmp + dfe(i)
    tmp0= tmp0+ dfe(i) * edos_all(i)
    !write(stdout,'(I5,3e20.9)') i, edos(i), dfe1(i), dfe(i)
  enddo
  if ( mod(nint(tt/0.1),nint(1.0_dp/0.1)) == 0 ) then
    write(stdout,*) 'trace of dfe=', tmp,tmp0
  endif

  end subroutine dfedt

  !-----------------------------------------------------------------------
  ! calculate heat capacity
  !-----------------------------------------------------------------------
  subroutine heatcapacity(np, dw, ebot, Tl, dos, heat_capacity)
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, kelvin2eV, two, zero, pi, ci, eps6, eps8
  USE epwcom,        ONLY : ngaussw, degaussw, eps_acustic
  implicit none
  integer,       intent(in)  :: np
  real(kind=dp), intent(in)  :: dw, ebot
  real(kind=dp), intent(in)  :: Tl, dos(np)
  real(kind=dp), intent(out) :: heat_capacity

  integer       :: i
  real(kind=dp) :: wq, x, dnq, occ
  real(kind=DP), external :: wgauss

  heat_capacity = 0.0_dp
  do i = 1, np
    wq = ebot + (i-1) * dw
    if (wq > eps_acustic) then
       
       x = wq/Tl
       occ = wgauss(-x,-99)          ! 1/(exp(-x)+1)
       ! one - two * occ = [exp(-x) - 1] / [exp(-x) + 1]
       ! one - occ = exp(-x) /[exp(-x) + 1] 
       occ = occ / (one - two * occ) ! 1/(exp(-x)-1)

       if ( x .lt. -100.0) then
       else if (x .gt. 100.0) then
       else
       endif

       !
       ! docc/de = - [exp(x) -1] ^{-2} * exp(x) * dx/T (x=w/k_bT)
       !         = occ * occ *exp(x) * wq / Tl**2

       dnq = 1.0_dp ! derivative of Bose-Einstein function
       dnq = occ * occ * exp(x) * wq / Tl / Tl !kb is canceled when kB * Eng2/heat_capacity 

       heat_capacity = heat_capacity + wq * dos(i) * dnq
    endif
  enddo
  !
  end subroutine heatcapacity

  !-----------------------------------------------------------------------
  subroutine getTe(np, dw, ebot, fe, Te)
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, kelvin2eV, two, zero, pi, ci, eps6, eps8
  implicit none
  
  integer,       intent(in)    :: np
  real(kind=dp), intent(in)    :: dw, ebot
  real(kind=dp), intent(in)    :: fe(np)
  real(kind=dp), intent(inout) :: Te
  !
  real(kind=dp) :: beta(1)
  real(kind=dp), allocatable :: x(:)
  !
  integer :: i, j

  allocate(x(np))

  beta = 1.0_dp / Te

  do i = 1, np
    x(i) = ebot + dble(i-1) * dw
  enddo

  call fit_fermi(np,x,fe,beta)

  write(stdout,*) 'fitted temperature=', Te, 1.0_dp /beta(1)
  Te = 1.0_dp /beta(1)
  !Te = Te

  deallocate(x)

  end subroutine getTe
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  SUBROUTINE relax_elec ( xkf_all, etf_all)
  !-----------------------------------------------------------------------
  !! 
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : linewidth_elself
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, shortrange, &
                            fsthick, eptemp, ngaussw, degaussw, &
                            eps_acustic, efermi_read, fermi_energy,&
                            restart, restart_freq
  USE pwcom,         ONLY : ef !, nelec, isk
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, &
                            nkf, epf17, wf, wqf, xkf, nkqtotf, &
                            sigmar_all, sigmai_all, sigmai_mode, zi_all, efnew
  USE transportcom,  ONLY : lower_bnd
  USE control_flags, ONLY : iverbosity
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6, eps8
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : inter_pool_comm
  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : stdout, ionode, ionode_id
  !
  implicit none
  !
  REAL(kind=DP), intent(in) :: xkf_all(3,nkqtotf)
  REAL(kind=DP), intent(in) :: etf_all(nbndsub,nkqtotf)
  !
  ! Local variables 
  !
  INTEGER :: n
  !! Integer for the degenerate average over eigenstates
  INTEGER :: ik
  !! Counter on the k-point index 
  INTEGER :: ikk
  !! k-point index
  INTEGER :: ikq
  !! q-point index 
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: jbnd
  !! Counter on bands
  INTEGER :: imode
  !! Counter on mode
  !! Number of states on the Fermi surface
  INTEGER :: nksqtotf
  !! Total number of k+q points 
  integer :: nee
  ! 
  REAL(kind=DP) :: tmp
  !! Temporary variable to store real part of Sigma for the degenerate average
  REAL(kind=DP) :: tmp2
  !! Temporary variable to store imag part of Sigma for the degenerate average
  REAL(kind=DP) :: tmp3
  !! Temporary variable to store Z for the degenerate average
  REAL(kind=DP) :: ekk2
  !! Temporary variable to the eigenenergies for the degenerate average
  REAL(kind=DP) :: sigmar_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the real-part of Sigma 
  REAL(kind=DP) :: sigmai_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the imag-part of Sigma 
  REAL(kind=DP) :: zi_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the Z
  REAL(kind=DP) :: g2
  !! Electron-phonon matrix elements squared in Ry^2
  REAL(kind=DP) :: ekk
  !! Eigen energy on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ekq
  !! Eigen energy of k+q on the fine grid relative to the Fermi level
  REAL(kind=DP) :: wq
  !! Phonon frequency on the fine grid
  REAL(kind=DP) :: ef0
  !! Fermi energy level
  REAL(kind=DP) :: wgq
  !! Bose occupation factor $n_{q\nu}(T)$
  REAL(kind=DP) :: wgkq
  !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
  REAL(kind=DP) :: weight
  !! Self-energy factor 
  !!$$ N_q \Re( \frac{f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q\nu} - i\delta }) $$ 
  !!$$ + N_q \Re( \frac{1- f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q\nu} - i\delta }) $$ 
  REAL(kind=DP) :: w0g1
  !! Dirac delta for the imaginary part of $\Sigma$
  REAL(kind=DP) :: w0g2
  !! Dirac delta for the imaginary part of $\Sigma$
  REAL(kind=DP) :: inv_wq
  !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
  REAL(kind=DP) :: inv_eptemp0
  !! Inverse of temperature define for efficiency reasons
  REAL(kind=DP) :: g2_tmp
  !! If the phonon frequency is too small discart g
  REAL(kind=DP) :: inv_degaussw
  !! Inverse of the smearing for efficiency reasons
  !REAL(kind=DP), external :: efermig
  !! Function to compute the Fermi energy 
  REAL(kind=DP), external :: dos_ef
  !! Function to compute the Density of States at the Fermi level
  REAL(kind=DP), external :: wgauss
  !! Fermi-Dirac distribution function (when -99)
  REAL(kind=DP), external :: w0gauss
  !! This function computes the derivative of the Fermi-Dirac function
  !! It is therefore an approximation for a delta function
  !  
  ! SP: Define the inverse so that we can efficiently multiply instead of
  ! dividing

  !
  real(kind=dp), allocatable :: eig_k(:,:)
  real(kind=dp), allocatable :: xk_all(:,:)
  integer,       allocatable :: eeklist(:,:)

  ! 
  inv_eptemp0 = 1.0/eptemp
  inv_degaussw = 1.0/degaussw
  !
  if (ionode) call printHeader(stdout)
  call flush(stdout)

  nksqtotf = nkqtotf/2 ! odd-even for k,k+q
  
  allocate(xk_all(3,nksqtotf), eig_k(nbndsub, nksqtotf))
  xk_all = zero
  eig_k  = zero

  write(6,*) 'xf_all'
  !do ik = 1, nksqtotf
  !  xk_all(:,ik) = xkf_all(:,2*ik-1)
  !  eig_k(:,ik) = etf_all(:,2*ik-1)
  !  write(6,*) ik, xk_all(:,ik)
  !enddo

  if(.true.) then
    allocate(eeklist(8,1000000))
    write(stdout,*) 'calculate ee klist'
    write(stdout,*) nbndsub, nkqf, size(etf), dabs(maxval(etf))
    !call flush(stdout)
    !call ee_klist(nkqf/2,nksqtotf, nbndsub, nee, xk_all, eig_k, eeklist)
    !stop
  endif

  
  deallocate(xk_all, eig_k,eeklist)

  if(ionode) call printEnding(stdout)

  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  102 FORMAT(5x,'E( ',i3,' )=',f9.4,' eV   Re[Sigma]=',f15.6,' meV Im[Sigma]=',f15.6,' meV     Z=',f15.6,' lam=',f15.6)
  !
  RETURN
  !
  END SUBROUTINE relax_elec
  
  subroutine ee_klist(nkf,nkftot,num_wann,nee,xkf,eig_k,eeklist)
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  implicit none
  integer,       intent(in) :: nkf, nkftot, num_wann
  integer,      intent(out) :: nee  
  real(kind=dp), intent(in) :: xkf(3,nkftot)
  real(kind=dp), intent(in) :: eig_k(num_wann,nkftot)
  integer,    intent(inout) :: eeklist(8,*)
  
  real*8  :: vec(3),dist,crit
  integer :: k1,k2,k3,k4
  integer :: ist, jst, kst, lst, ijst, klst
  
  nee=0
  
  crit = 1.d-6
 
!  write(6,*) 'nkf, pool_id', nkf, my_pool_id, npool
  do k1 = 1,nkf
  do k2 = 1, nkftot
    if(k2==k1) cycle
    do k3 = 1, nkftot
    do k4 = 1, nkftot
      if(k3==k4) cycle
      vec = xkf(:,k1) + xkf(:,k2) - xkf(:,k3) - xkf(:,k4)
      dist=sqrt(dot_product(vec,vec))
  
      if(dist<crit) then
        !write(6,'(12f9.4)') xkf(:,k1), xkf(:,k2), xkf(:,k3), xkf(:,k4)
        write(stdout,'(4(I5,2x),e15.6)') k1,k2,k3,k4, dist
        call flush(stdout)
      endif

      if(dist<crit) then
        ijst = 0
        do ist=1,num_wann
        do jst=1,num_wann
          if (jst > ist) cycle
          ijst = ijst + 1

          klst = 0
          do kst=1,num_wann
          do lst=1,num_wann
             if (lst > kst ) cycle
             klst = klst + 1
             if (klst > ijst) cycle
 
             !if(dabs(eig_k(ist,k1)-eig_k(kst,k3))<crit) cycle
             dist=eig_k(ist,k1)+eig_k(jst,k2) - eig_k(kst,k3) - eig_k(lst,k4)
             if(abs(dist)<crit) then
               nee = nee + 1
               write(stdout,'(9I5,2e15.6)') nee, k1,k2,k3,k4,ist,jst,kst,lst,eig_k(ist,k1)-eig_k(kst,k3), dist
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
 
  write(6,'(A,I4,A)') 'totally', nee, ' iterms in the ee klist!'
  
  end subroutine ee_klist
 

!
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


