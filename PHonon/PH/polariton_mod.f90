!-------------------------------------------------------------------------------
!  Phonon-cavity coupling for dynmat.x (\Gamma-point phonon-polaritons)
!   minimal post-DFPT module (CBOA / linear response)
!  Author: Yu Zhang @ LANL
!-------------------------------------------------------------------------------
module polariton_mod
  use kinds,     only : DP
  use constants, only : RY_TO_THZ, RY_TO_CMM1, BOHR_RADIUS_ANGS
  implicit none
contains


subroutine phonon_polariton_nonanal(nat, nat_blk, itau_blk, epsil, q, zeu, omega, dyn, &
    ncav, omega_in, polvec, lambda_in, wpol, evec_pol, phot_frac)

  !-----------------------------------------------------------------------
  !     add the nonanalytical term with coupling to photon
  !
  use kinds, only: dp
  use constants, only: pi, fpi, e2
 implicit none
 integer, intent(in) :: nat, nat_blk, itau_blk(nat)
 !  nat: number of atoms in the cell (in the supercell in the case
 !       of a dyn.mat. constructed in the mass approximation)
 !  nat_blk: number of atoms in the original cell (the same as nat if
 !       we are not using the mass approximation to build a supercell)
 !  itau_blk(na): atom in the original cell corresponding to
 !                atom na in the supercell
 !
 complex(DP), intent(inout) :: dyn(3,3,nat,nat) ! dynamical matrix
 real(DP), intent(in) :: q(3),  &! polarization vector
      &       epsil(3,3),     &! dielectric constant tensor
      &       zeu(3,3,nat_blk),   &! effective charges tensor
      &       omega            ! unit cell volume
 !
 integer, intent(in) :: ncav
 real(DP), intent(in) :: omega_in(ncav)                     ! photon freq (user units)
 real(DP), intent(in) :: polvec(3,ncav)                     ! unit polarization vectors
 real(DP), intent(in) :: lambda_in(ncav)                    ! optional direct lambda (Gaussian a.u.), <=0 -> build from Vmode
 real(DP), intent(out) :: wpol(3*nat + ncav)                ! polariton freqs (Ry)
 real(DP), intent(out) :: evec_pol(3*nat+ncav, 3*nat+ncav)  ! polariton eigenvectors (dimensionless)
 real(DP), intent(out) :: phot_frac(3*nat+ncav)             ! photon fraction per mode
 !
 ! local variables
 !
 real(DP) zag(3),zbg(3),  &! eff. charges  times g-vector
      &       qeq              !  <q| epsil | q>
 integer na,nb,nc,           &! counters on atoms
      &  na_blk,nb_blk,      &! as above for the original cell
      &  i,j                  ! counters on cartesian coordinates
 !
 qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+    &
        q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+    &
        q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
 !
 write(6, *) "q in phonon_polariton_nonanal is:", q(1), q(2), q(3)

 !
 do na = 1,nat
    na_blk = itau_blk(na)
    do nc = 1, ncav

       do i=1,3
          !
          zag(i) = q(1)*zeu(1,i,na_blk) +  q(2)*zeu(2,i,na_blk) + &
                   q(3)*zeu(3,i,na_blk)
          ! for cavity
          zbg(i) = 0.0 ! TBA.

       enddo
       !
       ! phonon-photon interaction term
       do i = 1,3
          do j = 1,3
          end do
       end do
       !
    end do
 end do
 ! photon term
 do i = 1, ncav

 enddo
 !
 return

end subroutine phonon_polariton_nonanal

subroutine build_polaritons( &
   nat, nmodes, amass_atom, vol_bohr3, w2, z_in, zstar, epsinf, &
   ncav, omega_in, omega_units, polvec, lambda_in, vmode_ang3, eps_ext, &
   nout, wpol, evec_pol, phot_frac)
  !---------------- arguments ----------------
  integer, intent(in) :: nat, nmodes, ncav
  real(DP), intent(in) :: amass_atom(nat)                       ! amu
  real(DP), intent(in) :: vol_bohr3                             ! cell volume (bohr^3)
  real(DP), intent(in) :: w2(nmodes)                            ! phonon squared freqs (Ry)
  real(DP), intent(in) :: z_in(3*nat,nmodes)                    ! phonon eigenvectors (mass-weighted style)
  real(DP), intent(in) :: zstar(nat,3,3)                        ! Born Z*
  real(DP), intent(in) :: epsinf(3,3)                           ! epsilon_infty (ion-clamped)
  real(DP), intent(in) :: omega_in(ncav)                        ! photon freq (user units)
  character(len=*), intent(in) :: omega_units                   ! 'Ry','THz','cm-1','meV','eV'
  real(DP), intent(in) :: polvec(3,ncav)                        ! unit polarization vectors
  real(DP), intent(in) :: lambda_in(ncav)                       ! optional direct lambda (Gaussian a.u.), <=0 -> build from Vmode
  real(DP), intent(in) :: vmode_ang3(ncav)                      ! effective mode volumes [] (optional)
  real(DP), intent(in) :: eps_ext                               ! background dielectric of cavity
  integer, intent(out) :: nout
  real(DP), intent(out) :: wpol(nmodes+ncav)                    ! polariton freqs (Ry)
  real(DP), intent(out) :: evec_pol(nmodes+ncav,nmodes+ncav)    ! polariton eigenvectors (dimensionless)
  real(DP), intent(out) :: phot_frac(nmodes+ncav)               ! photon fraction per mode
  !---------------- locals ----------------
  integer :: i, j, iat, a, b, nu, mu, al, be, ndim, info, lwork
  real(DP) :: pi, ang2bohr, vol_ang3, convw, normp
  real(DP) :: qeq, I3(3,3), chi(3,3), lam(3,ncav)
  real(DP) :: wph(ncav), S(nmodes,ncav), X(ncav,ncav), IX(ncav,ncav)
  real(DP), allocatable :: K(:,:), work(:), eval(:)
  real(DP) :: dnu(3,nmodes), mhalf
  !
  pi = dacos(-1.0_DP)
  !
  ! ---- unit conversions (to Rydberg atomic units) ----
  select case (trim(adjustl(omega_units)))
  case ('Ry');   convw = 1.0_DP
  case ('THz');  convw = 1.0_DP/RY_TO_THZ
  case ('cm-1'); convw = 1.0_DP/RY_TO_CMM1
  case ('meV');  convw = 1.0_DP/13.605693009D0/1.0D-3
  case ('eV');   convw = 1.0_DP/13.605693009D0
  case default;  convw = 1.0_DP
  end select
  do al=1,ncav
     wph(al) = omega_in(al)*convw
  enddo
  !
  ! ---- (Gaussian a.u.) ----
  I3 = 0.0_DP; do a=1,3; I3(a,a)=1.0_DP; enddo
  !
  chi = (epsinf - I3) * (vol_bohr3/(4.0_DP*pi))
  !
  !qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+    &
  !       q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+    &
  !       q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
  !
  ! ---- cavity couplings \lambda ----
  ang2bohr = 1.0_DP/BOHR_RADIUS_ANGS
  do al=1,ncav
     normp = max(1.0D-14, sqrt(polvec(1,al)**2+polvec(2,al)**2+polvec(3,al)**2))
     if (lambda_in(al) > 0.0_DP) then
        lam(:,al) = polvec(:,al)/normp * lambda_in(al)
     else if (vmode_ang3(al) > 0.0_DP) then
        ! \lambda =\vec{e} * sqrt(4\pi/ (Veff * \varepsilon_{ext})), with Veff in bohr^3
        vol_ang3 = max(1.0D-24, vmode_ang3(al))
        lam(:,al) = polvec(:,al)/normp * sqrt( 4.0_DP*pi / ( (vol_ang3*(ang2bohr**3)) * max(eps_ext,1.0D-12) ) )
     else
        ! if both lamda and voluem are zero, set lam zero
        lam(:,al) = 0.0_DP
     endif
  enddo
  WRITE(6, *) 'DEBUG-YZ: lam = ', lam
  !
  ! ---- mode effective dipoles (using Z* and eigenvectors) ----
  dnu = 0.0_DP
  do nu=1,nmodes
     i=0
     do iat=1, nat
        do a=1,3
           i = i+1
           mhalf = 1.0_DP/sqrt(amass_atom(iat))
           do b=1,3
              dnu(b,nu) = dnu(b,nu) + zstar(iat, a, b) * z_in(i,nu) * mhalf
           enddo
        enddo
     enddo
  enddo
  WRITE(6,*) "DEBUG-YZ: dnu = ", dnu
  !
  ! ---- X = \lambda^T \chi \lambda
  X = 0.0_DP
  do al=1,ncav
     do be=1,ncav
        do a=1,3; do b=1,3
           X(al,be) = X(al,be) + lam(a,al)*chi(a,b)*lam(b,be)
        enddo; enddo
     enddo
  enddo
  !
  IX = X
  do al=1,ncav
    IX(al,al)=IX(al,al) + 1.0_DP
 enddo
  call inv_spd(IX, ncav, info)   ! (I+X)^{-1} (SPD expected)
  WRITE(6, *) "DEBUG-YZ: X = ", X
  WRITE(6, *) "DEBUG-YZ: IX = ", IX
  !
  if (info /= 0) stop 'polariton_mod: inversion failed'
  !
  S = 0.0_DP
  do nu=1,nmodes
     do al=1,ncav
        S(nu,al) = lam(1,al)*dnu(1,nu) + lam(2,al)*dnu(2,nu) + lam(3,al)*dnu(3,nu)
     enddo
  enddo
  WRITE(6,*) "S = ", S

  !
  ! ---- Assemble symmetric K ----
  ndim = nmodes + ncav
  allocate(K(ndim,ndim)); K=0.0_DP
  WRITE(6, *) "S = ", S
  ! K_QQ
  do nu=1,nmodes
     K(nu,nu) = w2(nu)
     do mu=1,nmodes
        do al=1,ncav
           do be=1,ncav
              K(nu,mu) = K(nu,mu) + S(nu,al)*IX(al,be)*S(mu,be)
           enddo
        enddo
     enddo
     WRITE(6, '(A, I5, e15.7)') "KQQ(i,i) = ", nu, K(nu, nu)
  enddo
  !WRITE(6, *) "KQQ = ", K(1:nmodes, 1:nmodes)
  !
  ! TODO, currently we only considered the |0, 1> photon excitation space
  ! make it general to arbitray truncation
  !
  ! K_Qq (=K_qQ^T)
  do nu=1,nmodes
     do al=1,ncav
        do be=1,ncav
           K(nu, nmodes+al) = K(nu, nmodes+al) - S(nu,be)*IX(be,al)*wph(al)
           K(nmodes+al, nu) = K(nu, nmodes+al)
        enddo
     WRITE(6, '(A, 2I5, e15.7)') "KQq(i,i) = ", nu, al, K(nu, nmodes+al)
     enddo
  enddo
  ! WRITE(6, *) "KQq = ", K(1:nmodes, nmodes+1:nmodes+ncav)
  ! K_qq
  do al=1,ncav
     do be=1,ncav
        K(nmodes+al, nmodes+be) = wph(al)*IX(al,be)*wph(be)
     enddo
     WRITE(6, '(A, I5, e15.7)') "Kqq(i,i) = ", al, K(nmodes+al, nmodes+al)
  enddo
  !WRITE(6, *) "Kqq = ", K(nmodes+1:nmodes+1, nmodes+1:nmodes+1)
  !
  ! ---- Diagonalize K (symmetric) ----
  allocate(eval(ndim))
  lwork = max(1, 3*ndim)
  allocate(work(lwork))
  call dsyev('V','U',ndim,K,ndim,eval,work,lwork,info)
  if (info /= 0) stop 'polariton_mod: diagonalization failed'
  !
  nout = ndim
  write(6,*) " Eigenvalue of dynamical matrix K"
  do i=1,ndim
     wpol(i) = sqrt(max(0.0_DP, eval(i)))
     write(6, '(I5, e15.7)') i, eval(i)
  enddo
  evec_pol(:,:) = K(:,:)     ! eigenvectors returned in K
  !
  phot_frac = 0.0_DP
  do i=1,ndim
     do al=1,ncav
        phot_frac(i) = phot_frac(i) + evec_pol(nmodes+al,i)**2
     enddo
  enddo
  !
  deallocate(K, eval, work)
end subroutine build_polaritons

!---- small SPD inverse using Cholesky (LAPACK) ----
subroutine inv_spd(A, n, info)
  use kinds, only : DP
  implicit none
  integer, intent(in) :: n
  real(DP), intent(inout) :: A(n,n)
  integer, intent(out) :: info
  call dpotrf('U', n, A, n, info)
  if (info /= 0) return
  call dpotri('U', n, A, n, info)
  if (info /= 0) return
  ! fill lower triangle
  A = 0.5D0*(A + transpose(A))
end subroutine inv_spd


subroutine construct_polariton_dynmat(nat, nmodes, ncav, amass_atom, w2, zstar, z_in, chi, lam, K)
  implicit none
  integer,  intent(in) :: nat
  integer,  intent(in) :: nmodes
  integer,  intent(in) :: ncav
  real(DP), intent(in) :: amass_atom(nat)                       ! amu
  real(DP), intent(in) :: w2(nmodes)                            ! phonon squared freqs (Ry)
  real(DP), intent(in) :: zstar(nat,3,3)                        ! Born Z*
  real(DP), intent(in) :: z_in(3*nat,nmodes)                    ! phonon eigenvectors (mass-weighted style)
  real(DP), intent(in) :: chi(3, 3)
  real(DP), intent(in) :: lam(3,ncav)
  real(DP), intent(inout) :: K(nmodes + ncav, nmodes + ncav)

  !
  integer  :: i, j, iat, a, b, nu, mu, al, be, ndim, info, lwork
  real(DP) :: dnu(3,nmodes), mhalf
  real(DP) :: wph(ncav), S(nmodes,ncav), X(ncav,ncav), IX(ncav,ncav)

  ndim = nmodes + ncav
  ! ---- mode effective dipoles (using Z* and eigenvectors) ----
  dnu = 0.0_DP
  do nu=1,nmodes
     i=0
     do iat=1, nat
        do a=1,3
           i = i+1
           mhalf = 1.0_DP/sqrt(amass_atom(iat))
           do b=1,3
              dnu(b,nu) = dnu(b,nu) + zstar(iat, a, b) * z_in(i,nu) * mhalf
           enddo
        enddo
     enddo
  enddo
  !
  ! ---- X = \lambda^T \chi \lambda  and  S = \lambda\cdot d ----
  X = 0.0_DP
  do al=1,ncav
     do be=1,ncav
        do a=1,3; do b=1,3
           X(al,be) = X(al,be) + lam(a,al)*chi(a,b)*lam(b,be)
        enddo; enddo
     enddo
  enddo
  IX = X
  do al=1,ncav
    IX(al,al) = IX(al,al) + 1.0_DP
  enddo
  call inv_spd(IX, ncav, info)   ! (I+X)^{-1} (SPD expected)
  if (info /= 0) stop 'polariton_mod: inversion failed'
  !
  S = 0.0_DP
  do nu=1,nmodes
     do al=1,ncav
        S(nu,al) = lam(1,al)*dnu(1,nu) + lam(2,al)*dnu(2,nu) + lam(3,al)*dnu(3,nu)
     enddo
  enddo
  !
  ! ---- Assemble symmetric K ----
  WRITE(6, *) "S = ", S
  ! K_QQ
  do nu=1,nmodes
     K(nu,nu) = w2(nu)
     do mu=1,nmodes
        do al=1,ncav
           do be=1,ncav
              K(nu,mu) = K(nu,mu) + S(nu,al)*IX(al,be)*S(mu,be)
           enddo
        enddo
     enddo
  enddo
  WRITE(6, *) "KQQ = ", K(1:nmodes, 1:nmodes)
  !
  ! TODO, currently we only considered the |0, 1> photon excitation space
  ! make it general to arbitray truncation
  !
  ! K_Qq (=K_qQ^T)
  do nu=1,nmodes
     do al=1,ncav
        do be=1,ncav
           K(nu, nmodes+al) = K(nu, nmodes+al) - S(nu,be)*IX(be,al)*wph(al)
           K(nmodes+al, nu) = K(nu, nmodes+al)
        enddo
     enddo
  enddo
  ! K_qq
  do al=1,ncav
     do be=1,ncav
        K(nmodes+al, nmodes+be) = wph(al)*IX(al,be)*wph(be)
     enddo
  enddo
  WRITE(6, *) "Kqq = ", K(nmodes+1:nmodes+1, nmodes+1:nmodes+1)
end subroutine construct_polariton_dynmat


subroutine build_twisted_polaritons( &
   nat, nmodes, amass_atom, vol_bohr3, w2, z_in, zstar, epsinf, &
   twist_angle, thickness, distance, &
   ncav, omega_in, omega_units, polvec, lambda_in, vmode_ang3, eps_ext, &
   nout, wpol, evec_pol, phot_frac)
  !---------------- arguments ----------------
  integer, intent(in) :: nat, nmodes, ncav
  real(DP), intent(in) :: amass_atom(nat)                       ! amu
  real(DP), intent(in) :: vol_bohr3                             ! cell volume (bohr^3)
  real(DP), intent(in) :: w2(nmodes)                            ! phonon squared freqs (Ry)
  real(DP), intent(in) :: z_in(3*nat,nmodes)                    ! phonon eigenvectors (mass-weighted style)
  real(DP), intent(in) :: zstar(nat,3,3)                        ! Born Z*
  real(DP), intent(in) :: epsinf(3,3)                           ! epsilon_infty (ion-clamped)
  real(DP), intent(in) :: twist_angle                           ! twist angle
  real(DP), intent(in) :: thickness                             ! thickness of the slab ( = 1, for by-layer)
  real(DP), intent(in) :: distance                              ! distance between the two layers (slabs)
  real(DP), intent(in) :: omega_in(ncav)                        ! photon freq (user units)
  character(len=*), intent(in) :: omega_units                   ! 'Ry','THz','cm-1','meV','eV'
  real(DP), intent(in) :: polvec(3,ncav)                        ! unit polarization vectors
  real(DP), intent(in) :: lambda_in(ncav)                       ! optional direct lambda (Gaussian a.u.), <=0 -> build from Vmode
  real(DP), intent(in) :: vmode_ang3(ncav)                      ! effective mode volumes [Å^3] (optional)
  real(DP), intent(in) :: eps_ext                               ! background dielectric of cavity
  integer, intent(out) :: nout
  real(DP), intent(out) :: wpol(nmodes+ncav)                    ! polariton freqs (Ry)
  real(DP), intent(out) :: evec_pol(nmodes+ncav,nmodes+ncav)    ! polariton eigenvectors (dimensionless)
  real(DP), intent(out) :: phot_frac(nmodes+ncav)               ! photon fraction per mode
  !---------------- locals ----------------
  integer :: i, j, iat, a, b, nu, mu, al, be, ndim, info, lwork
  real(DP) :: pi, ang2bohr, vol_ang3, convw, normp
  real(DP) :: I3(3,3), chi(3,3), lam(3,ncav)
  real(DP) :: wph(ncav), S(nmodes,ncav), X(ncav,ncav), IX(ncav,ncav)
  real(DP), allocatable :: K0(:, :), K(:,:), work(:), eval(:)
  real(DP) :: dnu(3,nmodes), mhalf
  !
  pi = dacos(-1.0_DP)
  !
  ! ---- unit conversions (to Rydberg atomic units) ----
  select case (trim(adjustl(omega_units)))
  case ('Ry');   convw = 1.0_DP
  case ('THz');  convw = 1.0_DP/RY_TO_THZ
  case ('cm-1'); convw = 1.0_DP/RY_TO_CMM1
  case ('meV');  convw = 1.0_DP/13.605693009D0/1.0D-3
  case ('eV');   convw = 1.0_DP/13.605693009D0
  case default;  convw = 1.0_DP
  end select
  do al=1,ncav
     wph(al) = omega_in(al)*convw
  enddo

  ! compute the phonon-polariton of twisted bi-layer (or by-slab)
  ! TBA
  I3 = 0.0_DP; do a=1,3; I3(a,a)=1.0_DP; enddo
  chi = (epsinf - I3) * (vol_bohr3/(4.0_DP*pi))
  !
  ! ---- cavity couplings λ ----
  ang2bohr = 1.0_DP/BOHR_RADIUS_ANGS
  do al=1,ncav
     normp = max(1.0D-14, sqrt(polvec(1,al)**2+polvec(2,al)**2+polvec(3,al)**2))
     if (lambda_in(al) > 0.0_DP) then
        lam(:,al) = polvec(:,al)/normp * lambda_in(al)
     else if (vmode_ang3(al) > 0.0_DP) then
        ! \lambda = \ver{e} \cdot sqrt{4pi} / (Veff \varepsilon_{\mathrm{ext}} with Veff in bohr^3
        vol_ang3 = max(1.0D-24, vmode_ang3(al))
        lam(:,al) = polvec(:,al)/normp * sqrt( 4.0_DP*pi / ( (vol_ang3*(ang2bohr**3)) * max(eps_ext,1.0D-12) ) )
     else
        ! if both lamda and voluem are zero, set lam zero
        lam(:,al) = 0.0_DP
     endif
  enddo
  ! WRITE(6, *) 'DEBUG-YZ: lam = ', lam
  !
  ndim = 2 * nmodes + ncav
  allocate(K0(nmodes+ncav, nmodes+ncav))
  K0 = 0.0_DP
  call construct_polariton_dynmat(nat, nmodes, ncav, amass_atom, w2, zstar, z_in, chi, lam, K0)
  !
  allocate(K(2*nmodes + ncav, 2*nmodes + ncav))

  ! 1) first layer/slab
  K(1:nmodes, 1:nmodes) = K0(1:nmodes, 1:nmodes)

  ! 2) rotate second K for the second layer/slab
  ! K(nmodes + 1: 2*nmodes, nmodes+1:2*nmodes) =

  ! 3) add the coupling between two-layers/slabs
  ! K(1:nmodes, nmodes+1:2*nmodes) = xx

  ! 4) coupling between each layer/slab and cavity photon
  ! K(1:nmodes, 2*nmodes + ncav) = xx
  ! K(2*nmodes + ncav, 1:nmodes) = xx
  ! K(1+nmodes:2*nmodes, 2*nmodes + ncav) = xx
  ! K(2*nmodes + ncav, 1+nmodes:2*nmodes) = xx


  ! ---- Diagonalize K (symmetric) ----
  allocate(eval(ndim))
  lwork = max(1, 3*ndim)
  allocate(work(lwork))
  call dsyev('V','U',ndim,K,ndim,eval,work,lwork,info)
  if (info /= 0) stop 'polariton_mod: diagonalization failed'
  !
  nout = ndim
  do i=1,ndim
     wpol(i) = sqrt(max(0.0_DP, eval(i)))
  enddo
  evec_pol(:,:) = K(:,:)     ! eigenvectors returned in K
  !
  phot_frac = 0.0_DP
  do i=1,ndim
     do al=1,ncav
        phot_frac(i) = phot_frac(i) + evec_pol(nmodes+al,i)**2
     enddo
  enddo
  !
  deallocate(K, eval, work)


end subroutine build_twisted_polaritons

end module polariton_mod
