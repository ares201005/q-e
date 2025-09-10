!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
program dynmat
  !--------------------------------------------------------------------
  !! This program:
  !
  !! * reads a dynamical matrix file produced by the phonon code;
  !! * adds the nonanalytical part (if Z* and epsilon are read from file),
  !!   applies the chosen Acoustic Sum Rule (if q=0);
  !! * diagonalise the dynamical matrix; 
  !! * calculates IR and Raman cross sections (if Z* and Raman tensors
  !!   are read from file, respectively);
  !! * writes the results to files, both for inspection and for plotting.
  !
  !! Input data (namelist "input"):
  !
  !! * \(\text{fildyn} [character]: input file containing the dynamical matrix
  !!   (default: fildyn='matdyn')
  !! * \(q(3)) - [real]: calculate LO modes (add nonanalytic terms) along
  !!   the direction q (cartesian axis, default: q=(0,0,0) )
  !! * \(\text{amass}(\text{nt})\) - [real]: mass for atom type nt, amu
  !!   (default: amass is read from file fildyn)
  !! * \(\text{asr}\) - [character]: indicates the type of Acoustic Sum Rule imposed:
  !!    * 'no': no Acoustic Sum Rules imposed (default)
  !!    * 'simple':  previous implementation of the asr used
  !!      (3 translational asr imposed by correction of
  !!      the diagonal elements of the dynamical matrix)
  !!    * 'crystal': 3 translational asr imposed by optimized
  !!      correction of the dyn. matrix (projection).
  !!    * 'one-dim': 3 translational asr + 1 rotational asr
  !!      imposed by optimized correction of the dyn. mat. (the
  !!      rotation axis is the direction of periodicity; it
  !!      will work only if this axis considered is one of
  !!      the cartesian axis).
  !!    * 'zero-dim': 3 translational asr + 3 rotational asr
  !!      imposed by optimized correction of the dyn. mat.
  !!      Note that in certain cases, not all the rotational asr
  !!      can be applied (e.g. if there are only 2 atoms in a
  !!      molecule or if all the atoms are aligned, etc.).
  !!      In these cases the supplementary asr are cancelled
  !!      during the orthonormalization procedure (see below).
  !!      Finally, in all cases except 'no' a simple correction
  !!      on the effective charges is performed (same as in the
  !!      previous implementation).
  !! * \(\text{axis}\) - [integer]: indicates the rotation axis for a 1D system
  !!   (1=Ox, 2=Oy, 3=Oz ; default =3)
  !! * \(\text{lperm}\) - [logical]: TRUE to calculate Gamma-point mode contributions to
  !!   dielectric permittivity tensor (default: lperm=.false.)
  !! * \(\text{lplasma}\) - [logical]: TRUE to calculate Gamma-point mode effective plasma 
  !!   frequencies, automatically triggers lperm = TRUE
  !!   (default: lplasma=.false.)
  !! * \(\text{filout} - [character]: output file containing phonon frequencies and normalized
  !!   phonon displacements (i.e. eigenvectors divided by the
  !!   square root of the mass and then normalized; they are
  !!   not orthogonal). Default: filout='dynmat.out'
  !! * \(\text{fileig}\) - [character]: output file containing phonon frequencies and eigenvectors
  !!   of the dynamical matrix (they are orthogonal). Default: fileig=' '
  !! * \(\text{filmol}\) - [character]: as above, in a format suitable for 'molden'
  !!   (default: filmol='dynmat.mold')
  !! * \(\text{filxsf}\) - [character]: as above, in axsf format suitable for xcrysden
  !!   (default: filxsf='dynmat.axsf')
  !! * \(\text{loto_2d}\) - [logical]: set to TRUE to activate two-dimensional treatment of
  !!   LO-TO splitting.
  !
  USE kinds,       ONLY : DP
  USE mp,          ONLY : mp_bcast
  USE mp_global,   ONLY : mp_startup, mp_global_end
  USE mp_world,    ONLY : world_comm
  USE io_global,   ONLY : ionode, ionode_id, stdout
  USE environment, ONLY : environment_start, environment_end
  USE io_dyn_mat,  ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_dyn_mat, read_dyn_mat_tail
  USE constants,   ONLY : amu_ry
  USE dynamical,  ONLY : dyn, m_loc, ityp, tau, zstar, dchi_dtau
  USE rigid,       ONLY : dyndiag, nonanal, remove_dyn_interaction
  ! USE cell_base,  only : omega        ! unit-cell volume [bohr^3] filled by readers
  ! USE polariton_mod, only : build_polaritons
  !
  implicit none
  !
  integer, parameter :: ntypx = 10
  character(len=256):: fildyn, filout, filmol, filxsf, fileig
  character(len=3) :: atm(ntypx)
  character(len=10) :: asr
  logical :: lread, gamma, loto_2d
  complex(DP), allocatable :: z(:,:)
  real(DP) :: amass(ntypx), amass_(ntypx), eps0(3,3), a0, omega, &
       at(3,3), bg(3,3), q(3), q_(3)
  real(DP), allocatable :: w2(:)
  integer :: nat, na, nt, ntyp, iout, axis, nspin_mag, ios
  real(DP) :: celldm(6)
  logical :: xmldyn, lrigid, lraman, lperm, lplasma, remove_interaction_blocks
  logical, external :: has_xml
  integer :: ibrav, nqs
  integer, allocatable :: itau(:)
  !
  ! YZ: variables for phonon polaritons
  ! ---- polariton additions (by Yu Zhang) ----
  logical :: lcavity=.false., print_ir=.true.
  integer :: ncav=0, ios_pol=0
  character(len=8) :: cav_omega_units='cm-1'
  real(DP) :: eps_ext=1.0_DP
  real(DP), allocatable :: cav_omega(:), cav_lambda(:), cav_vmode(:), cav_pol(:,:)
  !
  namelist /POLARITON/ lcavity, ncav, cav_omega, cav_omega_units, cav_pol, &
                       cav_lambda, cav_vmode, eps_ext, print_ir
  !
  namelist /input/ amass, asr, axis, fildyn, filout, filmol, filxsf, &
                   fileig, lperm, lplasma, q, loto_2d, remove_interaction_blocks
  !
  ! code is parallel-compatible but not parallel
  !
  CALL mp_startup()
  CALL environment_start('DYNMAT')
  !
  IF (ionode) CALL input_from_file ( )
  !
  asr  = 'no'
  axis = 3
  fildyn='matdyn'
  filout='dynmat.out'
  filmol='dynmat.mold'
  filxsf='dynmat.axsf'
  fileig=' '
  amass(:)=0.0d0
  q(:)=0.0d0
  lperm=.false.
  lplasma=.false.
  loto_2d=.false.
  remove_interaction_blocks = .false.
  !
  IF (ionode) read (5,input, iostat=ios)
  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('dynmat', 'reading input namelist', ABS(ios))
  !
  IF (ionode) read (5, POLARITON, iostat=ios_pol)
  CALL mp_bcast(ios_pol, ionode_id, world_comm)
  CALL errore('dynmat', 'reading polariton namelist', ABS(ios_pol))
  !
  CALL mp_bcast(asr,ionode_id, world_comm)
  CALL mp_bcast(axis,ionode_id, world_comm)
  CALL mp_bcast(amass,ionode_id, world_comm)
  CALL mp_bcast(fildyn,ionode_id, world_comm)
  CALL mp_bcast(filout,ionode_id, world_comm)
  CALL mp_bcast(filmol,ionode_id, world_comm)
  CALL mp_bcast(fileig,ionode_id, world_comm)
  CALL mp_bcast(filxsf,ionode_id, world_comm)
  CALL mp_bcast(q,ionode_id, world_comm)
  CALL mp_bcast(remove_interaction_blocks, ionode_id, world_comm)
  !
  IF (ionode) inquire(file=fildyn,exist=lread)
  CALL mp_bcast(lread, ionode_id, world_comm)
  IF (lread) THEN
     IF (ionode) WRITE(6,'(/5x,a,a)') 'Reading Dynamical Matrix from file '&
                                     , TRIM(fildyn)
  ELSE
     CALL errore('dynmat', 'File '//TRIM(fildyn)//' not found', 1)
  END IF
  !
  ntyp = ntypx ! avoids spurious out-of-bound errors
  xmldyn=has_xml(fildyn)
  IF (xmldyn) THEN
     CALL read_dyn_mat_param(fildyn,ntyp,nat)
     ALLOCATE (m_loc(3,nat))
     ALLOCATE (tau(3,nat))
     ALLOCATE (ityp(nat))
     ALLOCATE (zstar(3,3,nat))
     ALLOCATE (dchi_dtau(3,3,3,nat) )
     CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, atm, amass_, tau, ityp, &
             m_loc, nqs, lrigid, eps0, zstar, lraman, dchi_dtau)
     IF (nqs /= 1) CALL errore('dynmat','only q=0 matrix allowed',1)
     a0=celldm(1) ! define alat
     ALLOCATE (dyn(3,3,nat,nat) )
     CALL read_dyn_mat(nat,1,q_,dyn(:,:,:,:))
     CALL read_dyn_mat_tail(nat)
     IF(asr.ne.'no') THEN
         CALL set_asr ( asr, axis, nat, tau, dyn, zstar )
     END IF
     IF (ionode) THEN
        DO nt=1, ntyp
           IF (amass(nt) <= 0.0d0) amass(nt)=amass_(nt)
        END DO
     END IF
  ELSE
     IF (ionode) THEN
        CALL readmat2 ( fildyn, asr, axis, nat, ntyp, atm, a0, &
                        at, omega, amass_, eps0, q_ )
        DO nt=1, ntyp
           IF (amass(nt) <= 0.0d0) amass(nt)=amass_(nt)/amu_ry
        END DO
     END IF
  ENDIF
  IF (remove_interaction_blocks)  CALL remove_dyn_interaction(dyn, nat) 
  !
  IF (ionode) THEN
     !
     ! from now on, execute on a single processor
     !
     gamma = ( abs( q_(1)**2+q_(2)**2+q_(3)**2 ) < 1.0d-8 )
     !
     IF (gamma .and. .not.loto_2d) THEN
        ALLOCATE (itau(nat))
        DO na=1,nat
           itau(na)=na
        END DO
        CALL nonanal ( nat, nat, itau, eps0, q, zstar, omega, dyn )
        DEALLOCATE (itau)
     END IF
     !
     ALLOCATE ( z(3*nat,3*nat), w2(3*nat) )
     CALL dyndiag(nat,ntyp,amass,ityp,dyn,w2,z)
     !
     IF (filout.eq.' ') then
        iout=6
     ELSE
        iout=4
        OPEN (unit=iout,file=filout,status='unknown',form='formatted')
     END IF
     CALL writemodes(nat,q_,w2,z,iout)

     !
     !! if (lcavity .and. ios_pol==0) then
     !!    integer :: nmodes, i, nout, na
     !!    real(DP), allocatable :: amass_atom(:), zreal(:,:), wpol(:), evec_pol(:,:), phot_frac(:)
     !!    nmodes = 3*nat
     !!    allocate(amass_atom(nat))
     !!    do na=1,nat
     !!       amass_atom(na) = amass(ityp(na))
     !!    enddo
     !!    allocate(zreal(nmodes,nmodes))
     !!    do i=1,nmodes
     !!       zreal(:,i) = real(z(:,i), kind=DP)   ! take real part (Γ phonons)
     !!    enddo
     !!    if (.not.allocated(cav_pol)) then
     !!       allocate(cav_pol(3,max(1,ncav))); cav_pol=0.0_DP
     !!       if (ncav>=1) cav_pol(:,1) = (/0.0_DP,0.0_DP,1.0_DP/)  ! default z-pol
     !!    endif
     !!    if (.not.allocated(cav_omega)) then
     !!       allocate(cav_omega(max(1,ncav))); cav_omega=0.0_DP
     !!    endif
     !!    if (.not.allocated(cav_lambda)) then
     !!       allocate(cav_lambda(max(1,ncav))); cav_lambda=0.0_DP
     !!    endif
     !!    if (.not.allocated(cav_vmode)) then
     !!       allocate(cav_vmode(max(1,ncav))); cav_vmode=0.0_DP
     !!    endif
     !!    allocate(wpol(nmodes+ncav), evec_pol(nmodes+ncav,nmodes+ncav), phot_frac(nmodes+ncav))
     !!    call build_polaritons( nat, nmodes, amass_atom, omega, w2, zreal, zstar, eps0, &
     !!         ncav, cav_omega, cav_omega_units, cav_pol, cav_lambda, cav_vmode, eps_ext, &
     !!         nout, wpol, evec_pol, phot_frac )
     !!    write(iout,'(/,a)') ' ===== CAVITY–POLARITON SUMMARY (Γ) ====='
     !!    write(iout,'(a)')    '  #    freq(cm-1)    photon_frac'
     !!    do i=1,nout
     !!       write(iout,'(i3,2x,f12.4,3x,f7.3)') i, wpol(i)*RY_TO_CMM1, phot_frac(i)
     !!    enddo
     !!    write(iout,'(a,/,a)') '  (frequencies converted from Ry to cm^-1 using RY_TO_CMM1)', &
     !!                          '  Note: modes beyond 3*nat are mostly photonic.'
     !!    deallocate(amass_atom, zreal, wpol, evec_pol, phot_frac)
     !! endif
     !

     IF(iout .ne. 6) close(unit=iout)
     IF (fileig .ne. ' ') THEN
       OPEN (unit=15,file=TRIM(fileig),status='unknown',form='formatted')
       CALL write_eigenvectors (nat,ntyp,amass,ityp,q_,w2,z,15)
       CLOSE (unit=15)
     ENDIF
     CALL writemolden (filmol, gamma, nat, atm, a0, tau, ityp, w2, z)
     CALL writexsf (filxsf, gamma, nat, atm, a0, at, tau, ityp, z)
     IF (gamma) THEN 
        CALL RamanIR (nat, omega, w2, z, zstar, eps0, dchi_dtau)
        IF (lperm .OR. lplasma) THEN
            CALL polar_mode_permittivity(nat,eps0,z,zstar,w2,omega, &
                                         lplasma)
            IF ( ABS( q(1)**2+q(2)**2+q(3)**2 ) > 1.0d-8 ) &
               WRITE(6,'(5x,a)') 'BEWARE: phonon contribution to &
               & permittivity computed with TO-LO splitting'
        ENDIF
     ENDIF
  ENDIF
  !
  IF (xmldyn) THEN
     DEALLOCATE (m_loc)
     DEALLOCATE (tau)
     DEALLOCATE (ityp)
     DEALLOCATE (zstar)
     DEALLOCATE (dchi_dtau)
     DEALLOCATE (dyn)
  ENDIF
  !
  CALL environment_end('DYNMAT')
  !
  CALL mp_global_end()
  !
end program dynmat
