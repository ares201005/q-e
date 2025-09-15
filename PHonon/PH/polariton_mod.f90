!-------------------------------------------------------------------------------
!  Phonon–cavity coupling for dynmat.x (Γ-point phonon-polaritons)
!    —  minimal post-DFPT module (CBOA / linear response)
!  Author: Yu Zhang @ LANL
!-------------------------------------------------------------------------------
module polariton_mod
  use kinds,     only : DP
  use constants, only : RY_TO_THZ, RY_TO_CMM1, BOHR_RADIUS_ANGS
  implicit none
contains

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
  real(DP), intent(in) :: lambda_in(ncav)                       ! optional direct λ (Gaussian a.u.), <=0 -> build from Vmode
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
  ! ---- χ from ε∞ : χ = (ε∞ - I) * Ω_cell / (4π)  (Gaussian a.u.) ----
  I3 = 0.0_DP; do a=1,3; I3(a,a)=1.0_DP; enddo
  chi = (epsinf - I3) * (vol_bohr3/(4.0_DP*pi))
  !
  ! ---- cavity couplings λ ----
  ang2bohr = 1.0_DP/BOHR_RADIUS_ANGS
  do al=1,ncav
     normp = max(1.0D-14, sqrt(polvec(1,al)**2+polvec(2,al)**2+polvec(3,al)**2))
     if (lambda_in(al) > 0.0_DP) then
        lam(:,al) = polvec(:,al)/normp * lambda_in(al)
     else
        ! λ = ê * sqrt(4π / (Veff * ε_ext)), with Veff in bohr^3
        vol_ang3 = max(1.0D-24, vmode_ang3(al))
        lam(:,al) = polvec(:,al)/normp * sqrt( 4.0_DP*pi / ( (vol_ang3*(ang2bohr**3)) * max(eps_ext,1.0D-12) ) )
     endif
  enddo
  write(6, *) 'DEBUG-YZ: lam = ', lam
  !
  ! ---- mode effective dipoles dν (using Z* and eigenvectors) ----
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
  ! ---- X = λ^T χ λ  and  S = λ·d ----
  X = 0.0_DP
  do al=1,ncav
     do be=1,ncav
        do a=1,3; do b=1,3
           X(al,be) = X(al,be) + lam(a,al)*chi(a,b)*lam(b,be)
        enddo; enddo
     enddo
  enddo
  IX = X
  do al=1,ncav; IX(al,al)=IX(al,al)+1.0_DP; enddo
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
  ndim = nmodes + ncav
  allocate(K(ndim,ndim)); K=0.0_DP
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
  !
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

end module polariton_mod
