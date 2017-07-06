      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_find_neigh_old
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine find the neighbours of a sub domain.
      !                 It does that by comparing the coordinates of the
      !                 bounding box of the sub domains using an N-square
      !                 algorithm ! Moreover, since the routine is called before
      !                 any mapping has been performed, ALL the subs (nsubs)
      !                 must be searched !
      !
      !  Input        : min_phys(:)  (F) the min. extent of the physical domain
      !                 max_phys(:)  (F) the max. extent of the physical domain
      !                 nsubs        (I) the total number of (real) sub domains
      !                 bcdef(:)     (I) boundary condition definition
      !
      !  Input/output : min_sub(:,:) (F) the min. extent of the sub domain
      !                 max_sub(:,:) (F) the max. extent of the sub domain
      !
      !  Output       : nneigh(:)    (I) nneigh(isub) returns the total number
      !                                  of neighbouring subs of isub
      !               : ineigh(:,:)  (I) points to the nneigh(:) subs of isub
      !                 info         (I) return status
      !
      !  Remarks      : This routine offers plenty of performance improvements
      !                 If it becomes too expensive we could do some ranking
      !                 of the sub domains to search the neighbours more
      !                 efficiently.
      !
      !                 The side effect of this routine is that the lists
      !                 min_sub(:) and max_sub(:) are extended to include
      !                 ghost sub domains. This is not used beyond this routine.
      !
      !                 The lists that are built contain unique entries.
      !                 This means that a sub can occur at most once in the
      !                 neighbor list of a particular other sub. This needs
      !                 to be taken into account when doing periodic
      !                 ghosts.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_find_neigh_old.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:02  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/10/01 16:08:59  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.1  2004/07/26 08:55:30  ivos
      !  Renamed.
      !
      !  Revision 1.8  2004/07/26 07:42:37  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.7  2004/04/08 11:45:10  ivos
      !  bugfix: the same sub is added at most once to any list. This
      !  prevents the same mesh block from being sent as ghost twice.
      !
      !  Revision 1.6  2004/04/01 15:23:42  walther
      !  Added a comment in the header regarding the efficiency or lack thereof.
      !
      !  Revision 1.5  2004/01/26 17:23:46  walther
      !  Updated header and error checking.
      !
      !  Revision 1.4  2004/01/23 17:24:14  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.3  2003/12/12 15:54:00  ivos
      !  Removed topoid from argument list as it is not needed.
      !
      !  Revision 1.2  2003/12/09 13:42:00  hiebers
      !  added draft of 2d version, needs validation
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_find_neigh_old_s(min_phys,max_phys,bcdef, &
     &           min_sub,max_sub,nsubs,nneigh,ineigh,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_find_neigh_old_d(min_phys,max_phys,bcdef, &
     &           min_sub,max_sub,nsubs,nneigh,ineigh,info)
#endif
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_copy_image_subs
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys,max_phys
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: bcdef
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      INTEGER , DIMENSION(:,:), POINTER       :: ineigh
      INTEGER , DIMENSION(:  ), POINTER       :: nneigh
      INTEGER                 , INTENT(IN   ) :: nsubs
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim) :: len_phys
      REAL(MK):: mean_npbx,rmean_npbx,max_npbx,var_npbx
      REAL(MK):: mx1,mx2,mx3
      REAL(MK):: mn1,mn2,mn3
      INTEGER , DIMENSION(:), POINTER :: subid
      INTEGER , DIMENSION(ppm_dim) :: ldc
      INTEGER               :: nsubsplus
      INTEGER               :: i,j,ii,jj,k,iopt
      INTEGER               :: istat,isize
      REAL(MK)              :: t0
      LOGICAL               :: isin
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_find_neigh_old',t0,info)

      !-------------------------------------------------------------------------
      !  Allocate memory for subdomain IDs
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldc(1) = nsubs
      NULLIFY(subid)
      CALL ppm_alloc(subid,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_find_neigh_old',     &
     &        'allocation of subid failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the ID
      !-------------------------------------------------------------------------
      DO i=1,nsubs
         subid(i) = i
      ENDDO

      !-------------------------------------------------------------------------
      !  Add ghost domains for the periodic system
      !-------------------------------------------------------------------------
      CALL ppm_copy_image_subs(min_phys,max_phys,bcdef, &
     &   min_sub,max_sub,nsubs,subid,nsubsplus,info)
      IF (info.NE.0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Allocate memory for the neighbours of the subdomains
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ldc(1) = nsubsplus
      CALL ppm_alloc(nneigh,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_find_neigh_old',     &
     &        'allocation of subid failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the number of neighbours to zero
      !-------------------------------------------------------------------------
      DO i=1,nsubsplus
         nneigh(i) = 0
      ENDDO

      !-------------------------------------------------------------------------
      !  And allocate the pointers to the neighbours
      !-------------------------------------------------------------------------
      ldc(1) = 26 ! guessing
      ldc(2) = nsubsplus
      CALL ppm_alloc(ineigh,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_find_neigh_old',     &
     &        'allocation of subid failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the neighbor sub indices to undefined
      !-------------------------------------------------------------------------
      DO i=1,nsubsplus
          DO j=1,ldc(1)
              ineigh(j,i) = ppm_param_undefined
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Save the current leading dimension of ineigh
      !-------------------------------------------------------------------------
      isize = SIZE(ineigh,1)

      !-------------------------------------------------------------------------
      !  Searching the neighbours
      !-------------------------------------------------------------------------
      IF (ppm_dim.EQ.2) THEN
         !----------------------------------------------------------------------
         !  two dimensions
         !----------------------------------------------------------------------
         DO i=1,nsubs
            DO j=1,nsubsplus
               jj = subid(j)
               IF (i.LT.jj) THEN
                  mx1 = MIN(max_sub(1,i),max_sub(1,j))
                  mx2 = MIN(max_sub(2,i),max_sub(2,j))

                  mn1 = MAX(min_sub(1,i),min_sub(1,j))
                  mn2 = MAX(min_sub(2,i),min_sub(2,j))
                  !-------------------------------------------------------------
                  !
                  !-------------------------------------------------------------
                  IF (mx1.GE.min_sub(1,i).AND.mx1.LE.max_sub(1,i).AND. &
     &                mx1.GE.min_sub(1,j).AND.mx1.LE.max_sub(1,j).AND. &
     &                mx2.GE.min_sub(2,i).AND.mx2.LE.max_sub(2,i).AND. &
     &                mx2.GE.min_sub(2,j).AND.mx2.LE.max_sub(2,j)) THEN

                     !----------------------------------------------------------
                     !  Check that this sub is not already in the list
                     !----------------------------------------------------------
                     isin = .FALSE.
                     DO k=1,nneigh(i)
                        IF (ineigh(k,i) .EQ. jj) THEN
                           isin = .TRUE.
                           EXIT
                        ENDIF
                     ENDDO
                     IF (.NOT. isin) THEN
                        !-------------------------------------------------------
                        !  Check if we have enough memory allocated
                        !-------------------------------------------------------
                        IF (nneigh(i)+1.GT.isize.OR.nneigh(jj)+1.GT.isize) THEN
                           !----------------------------------------------------
                           !  if not, reallocate and update the isize variable
                           !----------------------------------------------------
                           iopt   = ppm_param_alloc_grow_preserve
                           ldc(1) = isize + 100
                           ldc(2) = nsubsplus
                           CALL ppm_alloc(ineigh,ldc,iopt,info)
                           IF (info.NE.0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,      &
     &                            'ppm_find_neigh_old',       &
     &                            'allocation of ineigh failed',__LINE__,info)
                              GOTO 9999
                           ENDIF
                           !----------------------------------------------------
                           !  Save the size
                           !----------------------------------------------------
                           isize  = ldc(1)
                        ENDIF

                        !-------------------------------------------------------
                        !  Store the neighbour set
                        !-------------------------------------------------------
                        nneigh(i)             = nneigh(i) + 1
                        ineigh(nneigh(i),i)   = jj
                        nneigh(jj)            = nneigh(jj) + 1
                        ineigh(nneigh(jj),jj) = i
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ELSE
         !----------------------------------------------------------------------
         !  three dimensions
         !----------------------------------------------------------------------
         DO i=1,nsubs
            DO j=1,nsubsplus
               jj = subid(j)
               IF (i.LT.jj) THEN
                  mx1 = MIN(max_sub(1,i),max_sub(1,j))
                  mx2 = MIN(max_sub(2,i),max_sub(2,j))
                  mx3 = MIN(max_sub(3,i),max_sub(3,j))

                  mn1 = MAX(min_sub(1,i),min_sub(1,j))
                  mn2 = MAX(min_sub(2,i),min_sub(2,j))
                  mn3 = MAX(min_sub(3,i),min_sub(3,j))
                  !-------------------------------------------------------------
                  !
                  !-------------------------------------------------------------
                  IF (mx1.GE.min_sub(1,i).AND.mx1.LE.max_sub(1,i).AND. &
     &                mx1.GE.min_sub(1,j).AND.mx1.LE.max_sub(1,j).AND. &
     &                mx2.GE.min_sub(2,i).AND.mx2.LE.max_sub(2,i).AND. &
     &                mx2.GE.min_sub(2,j).AND.mx2.LE.max_sub(2,j).AND. &
     &                mx3.GE.min_sub(3,i).AND.mx3.LE.max_sub(3,i).AND. &
     &                mx3.GE.min_sub(3,j).AND.mx3.LE.max_sub(3,j)) THEN

                     !----------------------------------------------------------
                     !  Check that this sub is not already in the list
                     !----------------------------------------------------------
                     isin = .FALSE.
                     DO k=1,nneigh(i)
                        IF (ineigh(k,i) .EQ. jj) THEN
                           isin = .TRUE.
                           EXIT
                        ENDIF
                     ENDDO
                     IF (.NOT. isin) THEN
                        !-------------------------------------------------------
                        !  Chech if we have enough memory allocated
                        !-------------------------------------------------------
                        IF (nneigh(i)+1.GT.isize.OR.nneigh(jj)+1.GT.isize) THEN
                           !----------------------------------------------------
                           !  if not, reallocate and update the isize variable
                           !----------------------------------------------------
                           iopt   = ppm_param_alloc_grow_preserve
                           ldc(1) = isize + 100
                           ldc(2) = nsubsplus
                           CALL ppm_alloc(ineigh,ldc,iopt,info)
                           IF (info.NE.0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,     &
     &                            'ppm_find_neigh_old',      &
     &                            'allocation of ineigh failed',__LINE__,info)
                              GOTO 9999
                           ENDIF
                           !----------------------------------------------------
                           !  Save the size
                           !----------------------------------------------------
                           isize  = ldc(1)
                        ENDIF

                        !-------------------------------------------------------
                        !  Store the neighbour set
                        !-------------------------------------------------------
                        nneigh(i)             = nneigh(i) + 1
                        ineigh(nneigh(i),i)   = jj
                        nneigh(jj)            = nneigh(jj) + 1
                        ineigh(nneigh(jj),jj) = i
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Free the memory again
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(subid,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_find_neigh_old',     &
     &        'deallocation of subid failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_find_neigh_old',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_find_neigh_old_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_find_neigh_old_d
#endif
