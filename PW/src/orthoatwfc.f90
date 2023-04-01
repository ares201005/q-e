!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE orthoUwfc(save_wfcatom)
  !-----------------------------------------------------------------------
  !
  ! This routine saves to buffer "iunhub" atomic wavefunctions having an
  ! associated Hubbard U term * S, for DFT+U(+V) calculations. Same for 
  ! "iunhub2" but without S (this is then used to computed Hubbard forces 
  ! and stresses). Atomic wavefunctions
  ! are orthogonalized if desired, depending upon the value of "Hubbard_projectors"
  ! "swfcatom" must NOT be allocated on input.
  !
  ! If save_wfcatom == .TRUE., also write atomic wavefunctions before
  ! applying S to buffer.
  !
  USE kinds,      ONLY : DP
  USE buffers,    ONLY : get_buffer, save_buffer
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunhub, iunhub_noS, nwordwfcU
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc, swfcatom
  USE klist,      ONLY : nks, xk, ngk, igk_k
  USE ldaU,       ONLY : Hubbard_projectors, wfcU, nwfcU, copy_U_wfc
  USE wvfct,      ONLY : npwx
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, &
                         bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only, use_gpu
  USE noncollin_module, ONLY : noncolin, npol
  USE mp_bands,         ONLY : use_bgrp_in_hpsi
  USE becmod_gpum,      ONLY : becp_d
  USE becmod_subs_gpum, ONLY : using_becp_auto, using_becp_d_auto, calbec_gpu
  USE uspp_init,        ONLY : init_us_2
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: save_wfcatom
  !! If .TRUE., write atomic wavefunction before applying S to buffer
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol, npw
  ! ik: the k point under consideration
  ! ibnd: counter on bands
  LOGICAL :: orthogonalize_wfc, normalize_only, save_flag
  COMPLEX(DP) , ALLOCATABLE :: wfcatom (:,:)
  !
  IF ( Hubbard_projectors == "pseudo" ) THEN
     WRITE( stdout,*) 'Beta functions used for Hubbard projectors'
     RETURN
  ELSE IF (Hubbard_projectors=="wf") THEN
     !
     ! Read Wannier functions from file (produced by pmw.x).
     !
     WRITE( stdout,*) 'Hubbard projectors are read from file produced by pmw.x'
     DO ik = 1, nks
        CALL get_buffer (wfcU, nwordwfcU, iunhub, ik)
     END DO
     RETURN
     !
  ELSE IF (Hubbard_projectors=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     normalize_only = .FALSE.
     WRITE( stdout,'(/5x,a,/)') 'Atomic wfc used for Hubbard projectors are NOT orthogonalized'
  ELSE IF (Hubbard_projectors=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.
     WRITE( stdout,'(/5x,a,/)') 'Atomic wfc used for Hubbard projectors are orthogonalized'
     IF (gamma_only) CALL errore('orthoUwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE IF (Hubbard_projectors=="norm-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .TRUE.
     WRITE( stdout,'(/5x,a,/)') 'Atomic wfc used for Hubbard projectors are normalized but NOT orthogonalized'
     IF (gamma_only) CALL errore('orthoUwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE
     WRITE(stdout,'(/5x,"Hubbard_projectors = ",a)') Hubbard_projectors
     CALL errore ("orthoUwfc"," This type of Hubbard projectors is not valid",1)
  END IF
  !
  ALLOCATE ( wfcatom(npwx*npol, natomwfc), swfcatom(npwx*npol, natomwfc) )
  !$acc enter data create(wfcatom, swfcatom)
  !
  save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi=.false.
  !
  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type (nkb,natomwfc, becp)
  CALL using_becp_auto(2)
  !
  DO ik = 1, nks
     !
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown (ik, wfcatom)
       !$acc update device(wfcatom)
     ELSE
       IF(use_gpu) THEN
         !$acc host_data use_device(wfcatom)
         CALL atomic_wfc_gpu( ik, wfcatom )
         !$acc end host_data
       ELSE
         CALL atomic_wfc (ik, wfcatom)
       END IF
     ENDIF
     npw = ngk (ik)
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb, use_gpu)
     if(use_gpu) then 
       CALL using_becp_d_auto(2)
       !$acc host_data use_device(vkb, wfcatom)
       CALL calbec_gpu( npw, vkb, wfcatom, becp_d )
       !$acc end host_data

       !$acc host_data use_device(wfcatom, swfcatom)
       CALL s_psi_gpu( npwx, npw, natomwfc, wfcatom, swfcatom )
       !$acc end host_data
     else
       CALL calbec (npw, vkb, wfcatom, becp)
       CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     end if 
     !
     IF (orthogonalize_wfc) CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, .FALSE. )
     !
     ! copy S * atomic wavefunctions with Hubbard U term only in wfcU
     ! (this is used during the self-consistent solution of Kohn-Sham equations)
     ! save to unit iunhub
     !
     !$acc update host(swfcatom)
     CALL copy_U_wfc (swfcatom, noncolin)
     IF ( nks > 1 ) CALL save_buffer (wfcU, nwordwfcU, iunhub, ik)
     !
     ! If save_wfcatom=.TRUE. copy the orthonormalized wfcatom to wfcU and save
     ! to unit iunhubnoS
     !
     IF (save_wfcatom.and..not.use_gpu) THEN
        IF (orthogonalize_wfc) CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, .TRUE. )
        CALL copy_U_wfc (wfcatom, noncolin)
        CALL save_buffer (wfcU, nwordwfcU, iunhub_noS, ik)
     ENDIF
     !
  ENDDO
  !$acc exit data delete(wfcatom, swfcatom)
  DEALLOCATE (wfcatom, swfcatom)
  CALL deallocate_bec_type ( becp )
  CALL using_becp_auto(2)
  !
  use_bgrp_in_hpsi = save_flag
  !
  RETURN
  !
END SUBROUTINE orthoUwfc
!
!-----------------------------------------------------------------------
SUBROUTINE orthoUwfc_k (ik, lflag)
  !-----------------------------------------------------------------------
  !
  ! For a given k point "ik", this routine computes (ortho-)atomic wavefunctions 
  ! having an associated Hubbard U term * S, for DFT+U(+V) calculations. 
  ! Also without S (this is then used to computed Hubbard forces and stresses). 
  ! wfcatom and swfcatom must be allocated on input.
  ! Beta functions vkb must be already computed before.
  !
  ! lflag=.TRUE.  : wfcU = O^{-1/2}  \phi (w/o ultrasoft S)
  ! lflag=.FALSE. : wfcU = O^{-1/2} S\phi (w/  ultrasoft S)
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : iunhub, nwordwfcU
  USE ions_base,        ONLY : nat
  USE basis,            ONLY : natomwfc, wfcatom, swfcatom
  USE klist,            ONLY : nks, xk, ngk, igk_k
  USE ldaU,             ONLY : Hubbard_projectors, wfcU, nwfcU, copy_U_wfc
  USE wvfct,            ONLY : npwx
  USE uspp,             ONLY : nkb, vkb
  USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, &
                               bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin 
  USE becmod_subs_gpum, ONLY : using_becp_auto
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik ! the k point under consideration
  LOGICAL, INTENT(IN) :: lflag
  !
  INTEGER :: ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
             l, lm, ltot, ntot, ipol, npw
  LOGICAL :: orthogonalize_wfc, normalize_only, save_flag
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)

  IF ( Hubbard_projectors == "pseudo" ) THEN
     CALL errore ("orthoUwfc_k","Hubbard_projectors=pseudo is not supported",1)
  ELSE IF (Hubbard_projectors=="wf") THEN
     CALL errore ("orthoUwfc_k","Hubbard_projectors=wf is not supported",1)
  ELSE IF (Hubbard_projectors=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     normalize_only = .FALSE.
  ELSE IF (Hubbard_projectors=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.    
     IF (gamma_only) CALL errore('orthoUwfc_k', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE IF (Hubbard_projectors=="norm-atomic") THEN
     CALL errore ("orthoUwfc_k","Hubbard_projectors=norm-atomic is not supported",1)
  ELSE
     WRITE(stdout,'(/5x,"Hubbard_projectors = ",a)') Hubbard_projectors
     CALL errore ("orthoUwfc_k"," this Hubbard_projectors type is not valid",1)
  END IF
  !
  ! Compute atomic wfc at this k (phi)
  IF (noncolin) THEN
     CALL atomic_wfc_nc_updown (ik, wfcatom)
  ELSE
     CALL atomic_wfc (ik, wfcatom)
  ENDIF
  !
  IF (Hubbard_projectors=="ortho-atomic") THEN
     ALLOCATE(aux(npwx,natomwfc))
     ! Copy atomic wfcs (phi)
     aux(:,:) = wfcatom(:,:)
  ENDIF
  !
  ! Number of plane waves at this k point
  npw = ngk(ik)
  !
  IF (orthogonalize_wfc .OR. .NOT.lflag) THEN
     ! Allocate the array becp = <beta|wfcatom>
     CALL allocate_bec_type (nkb,natomwfc, becp)
     CALL using_becp_auto(2)
     CALL calbec (npw, vkb, wfcatom, becp)
     ! Calculate swfcatom = S * phi
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     CALL deallocate_bec_type (becp)
     CALL using_becp_auto(2)
  ENDIF
  !
  ! Compute the overlap matrix
  ! lflag=.FALSE. : On the output wfcatom are unchanged, swfcatom = O^{-1/2} S\phi.
  ! lflag=.TRUE.  : On the output wfcatom = O^{-1/2} \phi (no ultrasoft S), swfcatom are unchanged.
  IF (orthogonalize_wfc) THEN
     !$acc data copy(wfcatom, swfcatom)
     CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, lflag )
     !$acc end data
  END IF
  !
  IF (lflag) THEN
     ! Copy (ortho-)atomic wavefunctions with Hubbard U term only
     ! in wfcU (no ultrasoft S): wfcatom = O^{-1/2} \phi.
     CALL copy_U_wfc (wfcatom, noncolin)
  ELSE
     ! Copy (ortho-)atomic wavefunctions with Hubbard U term only
     ! in wfcU (with ultrasoft S): swfcatom = O^{-1/2} S\phi.
     CALL copy_U_wfc (swfcatom, noncolin)
  ENDIF
  !
  IF (Hubbard_projectors=="ortho-atomic") THEN
     ! Copy atomic wfcs
     wfcatom(:,:) = aux(:,:)
     DEALLOCATE(aux)
  ENDIF
  !
  RETURN
  !   
END SUBROUTINE orthoUwfc_k
!
!-----------------------------------------------------------------------
SUBROUTINE orthoatwfc (orthogonalize_wfc)
  !-----------------------------------------------------------------------
  !
  ! This routine calculates atomic wavefunctions, orthogonalizes them
  ! if "orthogonalize_wfc" is .true., saves them into buffer "iunsat".
  ! "swfcatom" must be allocated on input.
  ! Useful for options "wannier" and "one_atom_occupations"
  !
  USE kinds,            ONLY : DP
  USE buffers,          ONLY : save_buffer
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : iunsat, nwordatwfc
  USE ions_base,        ONLY : nat
  USE basis,            ONLY : natomwfc, swfcatom
  USE klist,            ONLY : nks, xk, ngk, igk_k
  USE wvfct,            ONLY : npwx
  USE uspp,             ONLY : nkb, vkb
  USE becmod_gpum,      ONLY : becp_d
  USE becmod_subs_gpum, ONLY : using_becp_auto, using_becp_d_auto, calbec_gpu
  USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, &
                               bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only, use_gpu
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp_init,        ONLY : init_us_2
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in) :: orthogonalize_wfc
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol, npw
  ! ik: the k point under consideration
  ! ibnd: counter on bands
  LOGICAL :: normalize_only = .FALSE.
  COMPLEX(DP) , ALLOCATABLE :: wfcatom (:,:)
  
  normalize_only=.FALSE.
  ALLOCATE (wfcatom( npwx*npol, natomwfc))
  !$acc enter data create(wfcatom, swfcatom)

  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type (nkb,natomwfc, becp) 
  
  DO ik = 1, nks
     
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown (ik, wfcatom)
       !$acc update device(wfcatom)
     ELSE
       IF(use_gpu) THEN 
         !$acc host_data use_device(wfcatom)
         CALL atomic_wfc_gpu( ik, wfcatom )
         !$acc end host_data
       ELSE
         CALL atomic_wfc (ik, wfcatom)
       END IF
     ENDIF
     npw = ngk (ik)
     !
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb, use_gpu)
     !
     IF(use_gpu) THEN 
       CALL using_becp_auto(2)
       CALL using_becp_d_auto(2)
       !$acc host_data use_device(vkb, wfcatom)
       CALL calbec_gpu( npw, vkb, wfcatom, becp_d )
       !$acc end host_data
       !
       !$acc host_data use_device(wfcatom, swfcatom)
       CALL s_psi_gpu( npwx, npw, natomwfc, wfcatom, swfcatom )
       !$acc end host_data
     ELSE
       CALL calbec (npw, vkb, wfcatom, becp) 
       CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     END IF

     IF (orthogonalize_wfc) CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, .FALSE. )
     !
     ! write S * atomic wfc to unit iunsat
     !
     !$acc update host(swfcatom)
     CALL save_buffer (swfcatom, nwordatwfc, iunsat, ik)
     !
  ENDDO
  !$acc exit data delete(wfcatom, swfcatom)
  DEALLOCATE (wfcatom)
  CALL deallocate_bec_type ( becp )
  !
  RETURN
     
END SUBROUTINE orthoatwfc
!
!-----------------------------------------------------------------------
SUBROUTINE ortho_swfc ( npw, normalize_only, m, wfc, swfc, lflag )
  !-----------------------------------------------------------------------
  !
  ! On input : 
  ! wfc (npwx*npol,m) =  \phi = a set of "m" (atomic) wavefcts
  ! swfc(npwx*npol,m) = S\phi 
  ! normalize_only    = only normalize, do not orthonormalize
  !
  ! On output this routine will compute the overlap matrix O: 
  ! O_ij = <wfc_i|S|wfc_j> = <wfc_i|swfc_j>
  ! If lflag=.FALSE. : wfc are unchanged,   swfc = O^{-1/2} S\phi.
  ! If lflag=.TRUE.  : wfc = O^{-1/2} \phi, swfc are unchanged.
  !
  USE kinds,            ONLY : DP
  USE wvfct,            ONLY : npwx
  USE mp_bands,         ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
  USE mp,               ONLY : mp_sum
  USE noncollin_module, ONLY : noncolin, npol
  USE force_mod,        ONLY : eigenval, eigenvect, overlap_inv
  USE control_flags,    ONLY : use_gpu
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: m, npw
  LOGICAL, INTENT(IN) :: normalize_only
  COMPLEX(dp), INTENT(INOUT) :: wfc (npwx*npol,m)
  COMPLEX(dp), INTENT(INOUT) :: swfc(npwx*npol,m)
  LOGICAL, INTENT(IN) :: lflag
  !
  ! ... local variables
  !
  COMPLEX(DP) :: temp 
  COMPLEX(DP) , ALLOCATABLE ::  work (:,:), overlap (:,:)
  REAL(DP) , ALLOCATABLE :: e (:)
  !$acc declare device_resident(work, overlap, e)
#if defined(__CUDA)
  COMPLEX(DP) , ALLOCATABLE ::  s(:,:)
  !$acc declare device_resident(s)
#endif
  INTEGER :: i, j, k, ipol

  ALLOCATE (overlap( m , m))    
  ALLOCATE (work   ( m , m))    
  ALLOCATE (e      ( m))    
#if defined(__CUDA)
  ALLOCATE(s(m,m))
#endif
  ! 
  !$acc kernels
  overlap(:,:) = (0.d0,0.d0)
  work(:,:) = (0.d0,0.d0)
  !$acc end kernels
  !
  ! calculate overlap matrix
  !
  IF (noncolin) THEN
     !$acc host_data use_device(wfc, swfc, overlap)
     CALL MYZGEMM ('c', 'n', m, m, npwx*npol, (1.d0, 0.d0), wfc, &
          npwx*npol, swfc, npwx*npol, (0.d0,0.d0), overlap, m)
     !$acc end host_data
  ELSE
     !$acc host_data use_device(wfc, swfc, overlap)
     CALL MYZGEMM ('c', 'n', m, m, npw, (1.d0, 0.d0), wfc, &
          npwx, swfc, npwx, (0.d0, 0.d0), overlap, m)
     !$acc end host_data
  END IF
  !
  !$acc host_data use_device(overlap)
  CALL mp_sum(  overlap, intra_bgrp_comm )
  !$acc end host_data
  !
  IF ( normalize_only ) THEN
     !$acc parallel
     !$acc loop gang
     DO i = 1, m
        !$acc loop vector
        DO j = i+1, m
           overlap(i,j) = CMPLX(0.d0,0.d0, kind=dp)
           overlap(j,i) = CMPLX(0.d0,0.d0, kind=dp)
        ENDDO
     ENDDO
     !$acc end parallel
  END IF
  !
  ! find O^(-1/2) (actually, its transpose)
  !
!civn: ZHEEV not available in cuBLAS/cuSOLVER?
  IF(use_gpu) THEN
    !
#if defined(__CUDA)
    ! s_d = CMPLX(0.d0,0.d0, kind=dp)  ! fused below
    !$acc kernels
    s(:,:) = CMPLX(0.d0,0.d0, kind=dp) 
    DO i = 1, m
       s(i,i) = CMPLX(1.d0,0.d0, kind=dp)
    ENDDO
    !$acc end kernels
    ! THIS SHOULD BE A SIMPLE CDIAGH (NOT GENERALIZED!) DRIVER NEEDED IN LAXLIB
    !$acc host_data use_device(overlap, s, e, work)
    CALL laxlib_cdiaghg_gpu( m, m, overlap, s, m, e, work, me_bgrp, &
                             root_bgrp, intra_bgrp_comm )
    !$acc end host_data
#endif
    !
  ELSE
    CALL cdiagh (m, overlap, m, e, work)
  END IF 
  !
  !$acc parallel loop collapse(2)
  DO i = 1, m
    DO j = 1, m
      s(i,j) = work(i,j) * (1.d0/SQRT(e(j)))
    END DO
  END DO 
  !$acc host_data use_device(s, work, overlap)
  Call MYZGEMM( 'n', 'c', m, m, m, (1.d0, 0.d0), s, m, work, m, (0.d0,0.d0), overlap, m)
  !$acc end host_data
  !
  IF (lflag) THEN
     !
     ! Save quantities which are needed for 
     ! calculations of Hubbard forces and stress
     !$acc kernels copyout(eigenval, eigenvect, overlap_inv)
     eigenval(:) = e(:)
     eigenvect(:,:) = work(:,:)
     overlap_inv(:,:) = overlap(:,:)
     !$acc end kernels
     !
  END IF 
  !
  DEALLOCATE( work )
  !
  ALLOCATE( work(m, npwx*npol ) )
  !$acc kernels
  work(:,:) = (0.d0,0.d0)
  !$acc end kernels
  !
  IF (lflag) THEN
     !
     ! Transform atomic orbitals WITHOUT the ultrasoft S operator 
     ! O^(-1/2) \psi (note the transposition):
     ! \phi_I = \sum_J O^{-1/2}_JI \phi_J
     !
     IF(noncolin) THEN 
       !$acc host_data use_device(overlap, wfc, work)
       CALL MYZGEMM('n', 't', m, npwx*npol, m, (1.d0,0.d0), overlap, m, wfc, npwx*npol, (0.d0,0.d0), work, m )
       !$acc end host_data
       !$acc parallel loop collapse(2) 
       DO i = 1, npwx*npol
         DO j = 1, m
           wfc(i,j) = work(j,i)
         END DO 
       END DO
     ELSE
       !$acc host_data use_device(overlap, wfc, work)
       CALL MYZGEMM('n', 't', m, npw, m, (1.d0,0.d0), overlap, m, wfc, npwx*npol, (0.d0,0.d0), work, m )
       !$acc end host_data
       !$acc parallel loop collapse(2)
       DO i = 1, npw
         DO j = 1, m
           wfc(i,j) = work(j,i)
         END DO 
       END DO
     END IF
    
     !
     !
  ELSE
     !
     ! Transform atomic orbitals WITH the ultrasoft S operator 
     ! O^(-1/2) \Spsi (note the transposition):
     ! \Sphi_I = \sum_J O^{-1/2}_JI \Sphi_J
     ! FIXME: can be done in a faster way by using wfc as work space 
     !
     IF(noncolin) THEN 
       !$acc host_data use_device(overlap, swfc, work)
       CALL MYZGEMM('n', 't', m, npwx*npol, m, (1.d0,0.d0), overlap, m, swfc, npwx*npol, (0.d0,0.d0), work, m )
       !$acc end host_data 
       !$acc parallel loop collapse(2)
       DO i = 1, npwx*npol
         DO j = 1, m
           swfc(i,j) = work(j,i)
         END DO 
       END DO
     ELSE
       !$acc host_data use_device(overlap, swfc, work)
       CALL MYZGEMM('n', 't', m, npw, m, (1.d0,0.d0), overlap, m, swfc, npwx*npol, (0.d0,0.d0), work, m )
       !$acc end host_data
       !$acc parallel loop collapse(2)
       DO i = 1, npw
         DO j = 1, m
           swfc(i,j) = work(j,i)
         END DO 
       END DO
     END IF
     !
  ENDIF
  !
#if defined(__CUDA)
  DEALLOCATE(s)
#endif
  DEALLOCATE (overlap)
  DEALLOCATE (work)
  DEALLOCATE (e)
  !
  RETURN
  !      
END SUBROUTINE ortho_swfc
!
!-----------------------------------------------------------------------
SUBROUTINE calculate_doverlap_inv (m, e, work, doverlap, doverlap_inv)
  !---------------------------------------------------------------------
  !! This routine computes the derivative of O^{-1/2}, i.e.
  !! [d((O^{-1/2}))]_IJ, where O_IJ is the overlap matrix. 
  !! Note, on the input this routine requires dO (not transposed).
  !! The solution is written in a closed form by solving the Lyapunov
  !! equation (a particular case of the Sylvester equation).
  !! See Eq. (32) in PRB 105, 199901(E) (2022).
  !! See Eq. (32) in PRB 102, 235159 (2020).
  !! Written by I. Timrov (June 2020)
  !
  USE kinds,       ONLY : DP
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER, INTENT(IN)      :: m
  !! The total number of atomic functions 
  REAL(DP), INTENT(IN)     :: e(m)
  !! The eigenvalues of the overlap matrix
  COMPLEX(DP), INTENT(IN)  :: work(m,m)
  !! The eigenvectors of the overlap matrix
  COMPLEX(DP), INTENT(IN)  :: doverlap(m,m)
  !! The derivative of the overlap matrix O_IJ (not transposed)  
  COMPLEX(DP), INTENT(OUT) :: doverlap_inv(m,m)
  !! The derivative of transposed O^{-1/2}
  !
  ! Local variables
  INTEGER :: m1, m2, m3, m4
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  !$acc declare device_resident(aux)
  !! eigenvectors of the overlap matrix
  !! auxiliary array
  !
  ALLOCATE (aux(m,m))
  !
  ! Compute (work^H) * doverlap * work 
  ! and put the result back in doverlap
  !
  ! Compute aux = doverlap * work
  !$acc host_data use_device(doverlap, work, aux)
  CALL MYZGEMM('N','N', m, m, m, (1.d0,0.d0), doverlap, &
              m, work, m, (0.d0,0.d0), aux, m)
  !$acc end host_data
  ! Compute (work^H) * aux
  !$acc host_data use_device(work, aux, doverlap)
  CALL MYZGEMM('C','N', m, m, m, (1.d0,0.d0), work, &
              m, aux, m, (0.d0,0.d0), doverlap, m)
  !$acc end host_data
  !
  !$acc parallel loop collapse(2)
  DO m1 = 1, m
     DO m2 = 1, m
        aux(m1,m2) = doverlap(m1,m2) / &
                    (e(m1)*DSQRT(e(m2))+e(m2)*DSQRT(e(m1)))
     ENDDO
  ENDDO
  !
  ! Compute work * aux * (work^H)
  !
  ! Compute doverlap = aux * (work^H)
  !$acc host_data use_device(aux, work, doverlap)
  CALL MYZGEMM('N','C', m, m, m, (1.d0,0.d0), aux, &
              m, work, m, (0.d0,0.d0), doverlap, m)
  !$acc end host_data
  ! Compute doverlap_inv = work * doverlap
  !$acc host_data use_device(work, doverlap, doverlap_inv)
  CALL MYZGEMM('N','N', m, m, m, (-1.d0,0.d0), work, &
              m, doverlap, m, (0.d0,0.d0), doverlap_inv, m)
  !$acc end host_data
  !
  DEALLOCATE (aux)
  !
  RETURN
  !
END SUBROUTINE calculate_doverlap_inv
