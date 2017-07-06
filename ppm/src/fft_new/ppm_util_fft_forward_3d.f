      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_util_fft_forward_3d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Parallel forward FFT for 3D arrays. This does the 
      !                 complete transform. Call ppm_util_fft_init before
      !                 using this routine!
      !                 It takes the field quantity in DATA_fv which is 
      !                 assumed to be periodic outside the domain. 
      !                 This routine comes in two versions: 
      !                    ppm_util_fft_forward_3d is the standard version
      !                    ppm_util_fft_forward_pencil_3d is the pencil version
      !
      !                 The standard version performs the first FFT on an
      !                 xy-slab topology, followed by one on z pencils. The
      !                 most efficient input topology thus is an xy-slab
      !                 topology.
      !                 The pencil version performs three FFTs: the first
      !                 on x-pencils, the second on y-pencils, and the third
      !                 on z-pencils. It is provided for use with external
      !                 FFT libraries that do not provide 2d transforms.
      !                 The most efficient input topology is an x-pencil
      !                 topology.
      !
      !  Input        : lda_fv      (I) leading dimension (vector case only)
      !                 mesh_id_user(I) mesh ID of the current data field mesh
      !                 ghostsize(3)(I) ghostsize in all directions
      !
      !  Input/output : DATA_fv(:,:,:,:[,:]) (F) field data. FFT of the data on
      !                                          output.
      !
      !  Output       : info       (I) return status. 0 on success.
      !                   
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
<<<<<<< ppm_util_fft_forward_3d.f
      !  $Log: ppm_util_fft_forward_3d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
=======
      !  $Log: ppm_util_fft_forward_3d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.16  2006/09/04 18:34:59  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.15  2005/01/07 12:34:51  ivos
      !  fixed typo: isign_fft_fft to isign_fft.
      !
      !  Revision 1.14  2004/11/03 11:07:36  hiebers
      !  Removed __NEED4SPEED, exchanged __SXF90 by __MATHKEISAN
      ! (ffts are library specific not compiler specific)
      !
      !  Revision 1.13  2004/10/01 16:09:13  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.12  2004/08/27 14:03:35  michaebe
      !  Included calls to the advanced fftw interface.
      !   These are still only used
      !  if __NEED4SPEED is defined. Should be changed in the near future.
      !
      !  Revision 1.11  2004/08/26 10:47:25  hiebers
      !  inserted correction of the last entries
      !
      !  Revision 1.10  2004/08/25 07:54:15  hiebers
      !  fixed memory leak
      !
      !  Revision 1.9  2004/08/23 12:30:02  hiebers
      !  adjusted field format to ppm convention
      !
      !  Revision 1.8  2004/07/26 08:12:52  hiebers
      !  added SX MathKeisan Interface (NEC)
      !
      !  Revision 1.7  2004/07/26 07:42:33  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.6  2004/02/25 11:51:55  hiebers
      !  bug fix in array size
      !
      !  Revision 1.5  2004/02/19 15:43:03  walther
      !  Changed DOUBLE COMPLEX to COMPLEX(MK).
      !
      !  Revision 1.4  2004/02/11 16:31:20  ivos
      !  Changed include of fftw.f from a cpp include to a f90 INCLUDE.
      !
      !  Revision 1.3  2004/02/11 15:29:50  ivos
      !  Some heavy cosmetics: header formatted, includes merged, continuation
      !  characters moved to proper column, indentation fixed, blank lines and
      !  spaces in argument lists removed, too long lines in log wrapped.
      !  Bugfix: arguments are now checked BEFORE they are assigned to Nx_in...
      !
      !  Revision 1.2  2004/02/11 10:13:23  hiebers
      !  changed arguments, included test on info , included ppm_define.h, 
      !  shortened lines to 80 characters, excluded module_mesh, 
      !  included fftw3.f
      !
>>>>>>> 1.16
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __CASE == __SLAB

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_forward_3d_ss(DATA_fv,mesh_id_user, &
     &  ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_forward_3d_sd(DATA_fv,mesh_id_user, &
     &  ghostsize,info)
#endif

#endif


#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_forward_3d_vs(DATA_fv,lda_fv,mesh_id_user, &
        ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_forward_3d_vd(DATA_fv,lda_fv,mesh_id_user, & 
        ghostsize,info)
#endif
#endif

#else

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_forward_pencil_3d_ss(DATA_fv,mesh_id_user, &
     &    ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_forward_pencil_3d_sd(DATA_fv,mesh_id_user, &
     &    ghostsize,info)
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_forward_pencil_3d_vs(DATA_fv,lda_fv, &
     &    mesh_id_user,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_forward_pencil_3d_vd(DATA_fv,lda_fv, &
     &    mesh_id_user,ghostsize,info)
#endif
#endif

#endif

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_data_util_fft
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_util_fft_step_fd
      USE ppm_module_util_fft_map

      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
      REAL(MK), DIMENSION(:,:,:,:),    POINTER     :: DATA_fv
#elif __DIM == __VFIELD
      REAL(MK), DIMENSION(:,:,:,:,:),  POINTER     :: DATA_fv
#endif
      INTEGER                    , INTENT(IN)      :: mesh_id_user
#if   __DIM == __VFIELD
      INTEGER                    , INTENT(IN)      :: lda_fv
#endif
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                   :: t0
      INTEGER, DIMENSION(3)                      :: lda
      INTEGER                                    :: i,j,k,l,m 
      INTEGER                                    :: ires,jres,lres
      REAL(MK)                                   :: rN
#if   __DIM == __SFIELD
      INTEGER, PARAMETER                         :: lda_fv = 1
#endif
#if   __DIM == __SFIELD
      INTEGER , DIMENSION(4  )                   :: lda_DATA_fv
      COMPLEX(MK), DIMENSION(:,:,:,:), POINTER   :: DATA_fv_com
      INTEGER , DIMENSION(4  )                   :: lda_DATA_fv_com
#elif __DIM == __VFIELD
      INTEGER , DIMENSION(5  )                   :: lda_DATA_fv
      COMPLEX(MK), DIMENSION(:,:,:,:,:), POINTER :: DATA_fv_com
      INTEGER , DIMENSION(5  )                   :: lda_DATA_fv_com
#endif
      REAL(MK), DIMENSION(:,:,:), POINTER        :: data_in
      COMPLEX(MK), DIMENSION(:,:,:), POINTER     :: data_com
      COMPLEX(MK), DIMENSION(:,:,:), POINTER     :: FFT_x,FFT_xy,FFT_xyz
      REAL(MK), DIMENSION(:,:,:), POINTER        :: Reslt
      INTEGER , DIMENSION(:,:)  , POINTER        :: ndata,ndata_trans
      INTEGER , DIMENSION(:,:)  , POINTER        :: ndata_ypen,ndata_zpen
      INTEGER                                    :: yhmax,zhmax,n
      INTEGER                                    :: mesh_id_internal
      INTEGER                                    :: mesh_id_xpen,mesh_id_ypen
      INTEGER                                    :: mesh_id_xpen_complex
      INTEGER                                    :: mesh_id_zpen
      INTEGER                                    :: iopt,topo_id_user
      INTEGER                                    :: topo_id_internal
      INTEGER                                    :: topo_id_xpen,topo_id_ypen
      INTEGER                                    :: topo_id_zpen
      INTEGER                                    :: idom,nsublist
      INTEGER, DIMENSION(3)                      :: Nm,Nm_com
      CHARACTER(LEN=ppm_char)                    :: mesg
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_fft_forward_3d',t0,info)

      !-------------------------------------------------------------------------
      ! Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_util_fft_forward_3d',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ppm_field_topoid .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_forward_3d',  &
     &            'No field topology has been defined so far',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

#if  !(defined(__FFTW) | defined(__MATHKEISAN))
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
#ifndef __FFTW
      CALL ppm_error(ppm_err_nofftw,'ppm_util_fft_forward_3d',  &
     &    'PPM was compiled without fftw support',__LINE__,info)
#endif
#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_util_fft_forward_3d',  &
     &    'PPM was compiled without MATHKEISAN support',__LINE__,info)
#endif
      GOTO 9999   
#else

#if __CASE == __SLAB
      !-------------------------------------------------------------------------
      !  Map the data onto the xy-slab topology if needed 
      !-------------------------------------------------------------------------
      IF (.NOT.Its_xyslab_topo) THEN
#if   __DIM == __SFIELD
          CALL ppm_util_fft_map(DATA_fv,topo_ids_tmp,mesh_ids_tmp,   &
     &        ghostsize, info)
#elif __DIM == __VFIELD
          CALL ppm_util_fft_map(DATA_fv,lda_fv,topo_ids_tmp,mesh_ids_tmp,   &
     &        ghostsize, info)
#endif
          topo_ids_tmp(1) = topo_id_xpen
          topo_ids_tmp(2) = topo_id_zpen

          topo_id_internal = ppm_field_topoid
          mesh_id_internal = ppm_meshid(topo_id_internal)%internal(mesh_id_xpen)
      ENDIF
#else
      !-------------------------------------------------------------------------
      !  Map the data onto the x-pencil topology if needed
      !-------------------------------------------------------------------------
      IF (.NOT.Its_xpencil_topo) THEN
#if   __DIM == __SFIELD
          CALL ppm_util_fft_map(DATA_fv,topo_ids_tmp,mesh_ids_tmp,  &
     &        ghostsize,info)
#elif __DIM == __VFIELD
          CALL ppm_util_fft_map(DATA_fv,lda_fv,topo_ids_tmp,mesh_ids_tmp,  &
     &        ghostsize, info)
#endif
          topo_ids_tmp(1) = topo_id_xpen
          topo_ids_tmp(2) = topo_id_ypen

<<<<<<< ppm_util_fft_forward_3d.f
          topo_id_internal = ppm_field_topoid
          mesh_id_internal = ppm_meshid(topo_id_internal)%internal(mesh_id_xpen)
      ENDIF
#endif
=======
>>>>>>> 1.16

      !-------------------------------------------------------------------------
      !  Allocate complex array
      !-------------------------------------------------------------------------
      yhmax = 0
      zhmax = 0
      DO i=1,ppm_nsublist(topo_id_internal)
          idom = ppm_isublist(i,topo_id_internal)
          IF (ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes(2,idom) &
     &        .GT.yhmax) THEN
              yhmax = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%   &
     &            nnodes(2,idom)
          ENDIF
          IF (ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes(3,idom) &
     &        .GT.zhmax) THEN
              zhmax = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%   &
     &            nnodes(3,idom)
          ENDIF
      ENDDO

#if   __DIM == __SFIELD
      lda_DATA_fv_com(1) = Nm_com(1)
      lda_DATA_fv_com(2) = yhmax   
      lda_DATA_fv_com(3) = zhmax   
      lda_DATA_fv_com(4) = nsublist
#elif __DIM == __VFIELD
      lda_DATA_fv_com(1) = lda_fv
      lda_DATA_fv_com(2) = Nm_com(1)
      lda_DATA_fv_com(3) = yhmax   
      lda_DATA_fv_com(4) = zhmax   
      lda_DATA_fv_com(5) = nsublist
#endif

      iopt = ppm_param_alloc_fit
      CALL ppm_alloc(DATA_fv_com,lda_DATA_fv_com,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_forward_3d',     &
     &        'data array DATA_FV_COM',__LINE__,info)
          GOTO 9999
      ENDIF     

      !-------------------------------------------------------------------------
      !  FFT - Transformation in x-direction (PENCIL CASE) or on xy slabs
      !        (SLAB CASE).
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      lda(1) = 3
      lda(2) = ppm_nsubs(topo_id_internal)
      CALL ppm_alloc(ndata,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_forward_3d',     &
     &        'ndata array',__LINE__,info)
          GOTO 9999
      ENDIF     
      ndata = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes

      !-------------------------------------------------------------------------
      !  Treat each subdomain separately
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsublist(topo_id_internal)
          idom = ppm_isublist(k,topo_id_internal)
          !---------------------------------------------------------------------
          !  SIMONE: do we really need to alloc inside the loop ?
          !---------------------------------------------------------------------
          CALL ppm_alloc(data_in,ndata(:,idom),iopt,info)
          IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_util_fft_forward_3d',     &
     &           'data_in array',__LINE__,info)
              GOTO 9999
          ENDIF

#if   __DIM == __SFIELD
          !---------------------------------------------------------------------
          !  SCALAR CASE: Copy in, do the 1d transform, copy out
          !---------------------------------------------------------------------
          DO l=1,ndata(3,idom)
              DO j=1,ndata(2,idom)
                  DO i=1,ndata(1,idom)
                      data_in(i,j,l) =  DATA_fv(i,j,l,k)
                  ENDDO
              ENDDO
          ENDDO
          lda = ndata(:,idom)

#if __CASE == __SLAB
          CALL ppm_util_fft_slab_fd(data_in,lda,FFT_x,info) 
#else
          CALL ppm_util_fft_step_fd(data_in,lda,FFT_x,info) 
#endif

          DO l=1,lda(3)
              DO j=1,lda(2)
                  DO i=1,lda(1)
                      DATA_fv_com(i,j,l,k) = FFT_x(i,j,l)
                  ENDDO
              ENDDO
          ENDDO

          !---------------------------------------------------------------------
          !  SIMONE: do we really need to dealloc inside the loop ?
          !---------------------------------------------------------------------
          iopt = ppm_param_dealloc
          CALL ppm_alloc(FFT_x,ndata(:,idom),iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',  &
     &            'FFT_x',__LINE__,info)
          ENDIF

#elif __DIM == __VFIELD
          !---------------------------------------------------------------------
          !  VECTOR CASE: Copy in, do the 1d transform, copy out
          !  Do separate transforms for each component.
          !---------------------------------------------------------------------
          DO n=1,lda_fv  
              DO l=1,ndata(3,idom)
                  DO j=1,ndata(2,idom)
                      DO i=1,ndata(1,idom)
                          data_in(i,j,l) =  DATA_fv(n,i,j,l,k)
                      ENDDO
                  ENDDO
              ENDDO
              lda = ndata(:,idom)

#if __CASE == __SLAB
              CALL ppm_util_fft_slab_fd(data_in,lda,FFT_x,info) 
#else
              CALL ppm_util_fft_step_fd(data_in,lda,FFT_x,info) 
#endif

              DO l=1,lda(3)
                  DO j=1,lda(2)
                      DO i=1,lda(1)
                          DATA_fv_com(n,i,j,l,k) = FFT_x(i,j,l)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO

          !---------------------------------------------------------------------
          !  SIMONE: do we really need to dealloc inside the loop ?
          !---------------------------------------------------------------------
          iopt = ppm_param_dealloc
          CALL ppm_alloc(data_in,ndata(:,idom),iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &            'data_in',__LINE__,info)
          ENDIF
          CALL ppm_alloc(FFT_x,ndata(:,idom),iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &            'data_in',__LINE__,info)
          ENDIF

#endif
      ENDDO           ! loop over subdomains

#if __CASE == __SLAB
      !-------------------------------------------------------------------------
      !  Skip the x-y transpose and the y pencil transform in the slab case
      !-------------------------------------------------------------------------
      GOTO 1000
#endif

      !-------------------------------------------------------------------------
      !  Transpose x-direction and y-direction
      !-------------------------------------------------------------------------
      topo_ids_tmp(1) = topo_id_xpen
      topo_ids_tmp(2) = topo_id_ypen

      mesh_ids_tmp(1) = mesh_id_xpen_complex
      mesh_ids_tmp(2) = mesh_id_ypen

      !-------------------------------------------------------------------------
      !  Map the data onto the y pencil topology
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
      CALL ppm_util_fft_map(DATA_fv_com,topo_ids_tmp,mesh_ids_tmp,    &
     &    ghostsize,info)
#elif __DIM == __VFIELD
      CALL ppm_util_fft_map(DATA_fv_com,lda_fv,topo_ids_tmp,mesh_ids_tmp,    &
     &    ghostsize, info)
#endif

      DO k=1,ppm_nsublist(ppm_field_topoid)
          idom   = ppm_isublist(k,ppm_field_topoid)     
          lda(1) = 3
          lda(2) = idom
          iopt   = ppm_param_alloc_fit
          CALL ppm_alloc(ndata_trans,lda,iopt,info)
          
          ndata_trans(1,idom) = ndata_ypen(2,idom)
          ndata_trans(2,idom) = ndata_ypen(1,idom)
          ndata_trans(3,idom) = ndata_ypen(3,idom)

#if  __DIM == __VFIELD
          DO n=1,lda_fv
#endif
              !-----------------------------------------------------------------
              !  SIMONE: alloc within the loop necessary ?
              !-----------------------------------------------------------------
              CALL ppm_alloc(data_com,ndata_trans(:,idom),iopt,info)

              DO l=1,ndata_ypen(3,idom)
                  DO j=1,ndata_ypen(2,idom)
                      DO i=1,ndata_ypen(1,idom)
#if  __DIM == __SFIELD
                          data_com(j,i,l)= DATA_fv_com(i,j,l,k)
#elif __DIM == __VFIELD
                          data_com(j,i,l)= DATA_fv_com(n,i,j,l,k)
#endif
                      ENDDO
                  ENDDO
              ENDDO
           
              !-----------------------------------------------------------------
              !  FFT - Transformation in y-direction
              !-----------------------------------------------------------------

              CALL ppm_util_fft_step_fd(data_com,ndata_trans(:,idom),  &
     &            FFT_xy,info) 

              !-----------------------------------------------------------------
              !  Transpose back x-direction and y-direction
              !-----------------------------------------------------------------
              DO l=1, ndata_ypen(3,idom)
                  DO i=1, ndata_ypen(1,idom)
                      DO j=1, ndata_ypen(2,idom)

#if  __DIM == __SFIELD
                          DATA_fv_com(i,j,l,k) = FFT_xy(j,i,l)
#elif __DIM == __VFIELD
                          DATA_fv_com(n,i,j,l,k) = FFT_xy(j,i,l)
#endif
                      ENDDO
                  ENDDO
              ENDDO

              !-----------------------------------------------------------------
              !  SIMONE: alloc within the loop necessary ?
              !-----------------------------------------------------------------
              iopt = ppm_param_dealloc
              CALL ppm_alloc(data_com,ndata_trans(:,idom),iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &                'data_com',__LINE__,info)
              ENDIF
              CALL ppm_alloc(FFT_xy, ndata_trans(:,idom),iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &                'FFT_xy',__LINE__,info)
              ENDIF

#if  __DIM == __VFIELD
          ENDDO                      ! loop over lda
#endif
      ENDDO               ! loop over subdomains

1000  CONTINUE         ! SLAB CASE continues

      !-------------------------------------------------------------------------
      !  Transpose x-direction and z-direction
      !-------------------------------------------------------------------------
      topo_ids_tmp(1) = topo_id_ypen 
      topo_ids_tmp(2) = topo_id_zpen 

      mesh_ids_tmp(1) = mesh_id_ypen
      mesh_ids_tmp(2) = mesh_id_zpen
 
#if __CASE == __SLAB
      topo_ids_tmp(1) = topo_id_xpen 
      mesh_ids_tmp(1) = mesh_id_xpen_complex
#endif

      !-------------------------------------------------------------------------
      !  Map data onto z pencils
      !-------------------------------------------------------------------------
#if __DIM == __SFIELD
      CALL ppm_util_fft_map(DATA_fv_com,topo_ids_tmp,mesh_ids_tmp,   &
     &    ghostsize,info)
#elif __DIM == __VFIELD
      CALL ppm_util_fft_map(DATA_fv_com,lda_fv,topo_ids_tmp,   &
     &    mesh_ids_tmp,ghostsize,info)
#endif

      DO k=1,ppm_nsublist(ppm_field_topoid)
          idom   = ppm_isublist(k,ppm_field_topoid)     
          lda(1) = 3
          lda(2) = idom
          iopt   = ppm_param_alloc_fit

          !---------------------------------------------------------------------
          !  SIMONE: alloc within the loop needed ?
          !---------------------------------------------------------------------
          CALL ppm_alloc(ndata_trans,lda,iopt,info)
          IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_util_fft_forward_3d',     &
     &            'ndata_trans',__LINE__,info)
             GOTO 9999
          ENDIF

          ndata_trans(1,idom) = ndata_zpen(3,idom)
          ndata_trans(2,idom) = ndata_zpen(2,idom)
          ndata_trans(3,idom) = ndata_zpen(1,idom)

          CALL ppm_alloc(data_com,ndata_trans(:,idom),iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_fft_forward_3d',     &
     &            'data_com',__LINE__,info)
              GOTO 9999
          ENDIF

#if __DIM == __VFIELD
          DO n=1,lda_fv
#endif

              DO l=1, ndata_zpen(3,idom)
                  DO j=1, ndata_zpen(2,idom)
                      DO i=1, ndata_zpen(1,idom)

#if __DIM == __SFIELD
                          data_com(l,j,i) = DATA_fv_com(i,j,l,k)
#elif __DIM == __VFIELD
                          data_com(l,j,i) = DATA_fv_com(n,i,j,l,k)
#endif

                      ENDDO
                  ENDDO
              ENDDO
          
              !-----------------------------------------------------------------
              !  FFT - Transformation in z-direction
              !-----------------------------------------------------------------
#if __CASE == __SLAB
              CALL ppm_util_fft_z_fd(data_com,ndata_trans(:,idom),FFT_xyz,info) 
#else
              CALL ppm_util_fft_z_fd(data_com,ndata_trans(:,idom),FFT_xyz,info) 
#endif
#if __DIM == __VFIELD
          ENDDO            ! loop over lda
#endif
      ENDDO          ! loop over subdomains

      !-------------------------------------------------------------------------
      !  Deallocate memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(Reslt,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &        'RESLT',__LINE__,info)
      ENDIF
      CALL ppm_alloc(data_in,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &        'DATA_IN',__LINE__,info)
      ENDIF
      CALL ppm_alloc(data_com,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &        'DATA_COM',__LINE__,info)
      ENDIF
      CALL ppm_alloc(FFT_x,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &        'FFT_X',__LINE__,info)
      ENDIF
      CALL ppm_alloc(FFT_xy,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &        'FFT_XY',__LINE__,info)
      ENDIF
      CALL ppm_alloc(FFT_xyz,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &        'FFT_XYZ',__LINE__,info)
      ENDIF
      CALL ppm_alloc(DATA_fv_com,lda_DATA_fv_com,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_forward_3d',   &
     &        'DATA_FV_COM',__LINE__,info)
      ENDIF

#endif
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_fft_forward_3d',t0,info)
      RETURN

#if   __CASE == __SLAB
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_forward_3d_ss
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_forward_3d_sd
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_forward_3d_vs
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_forward_3d_vd
#endif
#endif

#else

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_forward_pencil_3d_ss
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_forward_pencil_3d_sd
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_forward_pencil_3d_vs
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_forward_pencil_3d_vd
#endif
#endif

#endif

