      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_util_fft_backward_2d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Parallel backward FFT for 2D arrays. This does the 
      !                 complete transform. Call ppm_util_fft_init before
      !                 using this routine!
      !                 It takes the field quantity omega in DATA_fv which is 
      !                 assumed to be periodic outside the domain. The most 
      !                 efficient initial topology for the transform is 
      !                 an x-pencil topology. After performing an FFT on the  
      !                 x-pencils, the data are mapped onto y-pencils where 
      !                 another a FFT is performed.
      !
      !  Input        : lda_fv      (I) leading dimension (vector case only)
      !                 mesh_id_user(I) mesh ID of the current data field mesh
      !                 topo_ids(2) (I) Temporary topologies for the FFTs as
      !                                 created and returned by
      !                                 ppm_util_fft_init.
      !                                 topo_ids(1) is an x-pencil topology, 
      !                                 topo_ids(2) is a y-pencil topology.
      !                 mesh_ids(3) (I) Temporary meshes for the FFTs as 
      !                                 created and returned by
      !                                 ppm_util_fft_init.
      !                                 mesh_ids(1): x-pencils, real
      !                                 mesh_ids(2): x-pencils, complex
      !                                 mesh_ids(3): y-pencils, complex
      !                 ghostsize(2)(I) ghostsize in all directions
      !
      !  Input/output : DATA_fv(:,:,:[,:]) (F) field data. FFT of the data on
      !                                        output.
      !
      !  Output       : info       (I) return status. 0 on success.
      !                   
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
<<<<<<< ppm_util_fft_backward_2d.f
      !  $Log: ppm_util_fft_backward_2d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
=======
      !  $Log: ppm_util_fft_backward_2d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:58  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.12  2005/02/16 10:25:03  hiebers
      !  use fftw_many_dft to speed up ffts
      !
      !  Revision 1.11  2004/11/03 11:10:30  hiebers
      !   exchanged __SXF90 by __MATHKEISAN (ffts are library specific
      !  not compiler specific)
      !
      !  Revision 1.10  2004/10/01 16:09:13  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.9  2004/08/25 07:54:34  hiebers
      !  fixed memory leak
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
      !  Revision 1.5  2004/02/19 15:43:02  walther
      !  Changed DOUBLE COMPLEX to COMPLEX(MK).
      !
      !  Revision 1.4  2004/02/11 16:31:19  ivos
      !  Changed include of fftw.f from a cpp include to a f90 INCLUDE.
      !
      !  Revision 1.3  2004/02/11 15:29:49  ivos
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
>>>>>>> 1.13
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_backward_2d_ss(DATA_fv, mesh_id_user,&
         topo_ids, mesh_ids, ghostsize,  info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_backward_2d_sd(DATA_fv, mesh_id_user,& 
         topo_ids, mesh_ids,  ghostsize,info)
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_backward_2d_vs(DATA_fv, lda_fv, mesh_id_user, &
     &  topo_ids, mesh_ids, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_backward_2d_vd(DATA_fv, lda_fv, mesh_id_user, &
     &  topo_ids, mesh_ids, ghostsize, info)
#endif
#endif

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE   ppm_module_data
      USE   ppm_module_data_mesh
      USE   ppm_module_data_util_fft
      USE   ppm_module_substart
      USE   ppm_module_substop
      USE   ppm_module_write
      USE   ppm_module_error
      USE   ppm_module_alloc
      USE   ppm_module_util_fft_step_fd
  
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                    , INTENT(IN)      :: mesh_id_user
#if   __DIM == __SFIELD
      REAL(MK), DIMENSION(:,:,:),  POINTER         :: DATA_fv
#elif   __DIM == __VFIELD
      REAL(MK), DIMENSION(:,:,:,:),  POINTER       :: DATA_fv
      INTEGER                    , INTENT(IN)      :: lda_fv
#endif
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: mesh_ids 
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
<<<<<<< ppm_util_fft_backward_2d.f
      REAL(MK)                                :: t0
      INTEGER, DIMENSION(2)                   :: lda
      INTEGER                                 :: k,i,j,n
      REAL(MK)                                :: rN
#if   __DIM == __SFIELD
      INTEGER, PARAMETER                       :: lda_fv = 1
#endif
#if   __DIM == __SFIELD
      COMPLEX(MK), DIMENSION(:,:,:), POINTER   :: DATA_fv_com
      INTEGER , DIMENSION(3  )                 :: lda_DATA_fv_com
#elif   __DIM == __VFIELD
      COMPLEX(MK), DIMENSION(:,:,:,:), POINTER :: DATA_fv_com
      INTEGER , DIMENSION(4  )                 :: lda_DATA_fv_com
#endif
      REAL(MK), DIMENSION(:, :),       POINTER :: data_in 
      COMPLEX(MK), DIMENSION(:,:),     POINTER :: data_com
      COMPLEX(MK), DIMENSION(:,:),     POINTER :: FFT_x, FFT_xy
      REAL(MK), DIMENSION(:,:),        POINTER :: Result
      REAL(MK), DIMENSION(2,1)                 :: xp
      INTEGER                                  :: Npart
      INTEGER                                  :: decomp, assign
      REAL(MK), DIMENSION(2  )                 :: min_phys, max_phys
      REAL(MK), DIMENSION(2  )                 :: length
      REAL(MK), DIMENSION(2  )                 :: length_phys
      INTEGER , DIMENSION(4  )                 :: bcdef 
      INTEGER                                  :: nsubs,topo_id, mesh_id
      INTEGER                                  :: mesh_id_internal
      INTEGER                                  :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                  :: mesh_id_xpen_complex
      REAL(MK), DIMENSION(:,:), POINTER        :: min_sub,max_sub
      REAL(MK), DIMENSION(:  ), POINTER        :: cost
      INTEGER , DIMENSION(:,:), POINTER        :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER        :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER        :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER        :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:  ), POINTER        :: sub2proc
      INTEGER , DIMENSION(:  ), POINTER        :: isublist
      INTEGER                                  :: nsublist, idom
      INTEGER                                  :: dim, yhmax
      INTEGER                                  :: iopt
      INTEGER                                  :: topo_id_user, topo_id_internal
      INTEGER                                  :: topo_id_xpen, topo_id_ypen
      INTEGER, DIMENSION(2)                    :: Nm, Nm_com, Nm_poisson
      CHARACTER(LEN=ppm_char)                  :: mesg
      
=======
      ! input data
      COMPLEX(MK), DIMENSION(:,:)       , INTENT(IN   ) :: data_in
      ! size of array
      INTEGER, DIMENSION(:)         , INTENT(INOUT) :: lda
      ! output data, inverse fast fourier transformed
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:)      , POINTER       :: data_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX| __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:)   , POINTER       :: data_out 
#endif
      INTEGER                       , INTENT(  OUT) :: info

>>>>>>> 1.13
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
<<<<<<< ppm_util_fft_backward_2d.f
      CALL substart('ppm_util_fft_backward_2d',t0,info)
=======
      ! timer
      REAL(MK)                                :: t0
      ! counters
      INTEGER                                 :: i,j,iopt
      ! size of the data_in 
      INTEGER                                 :: Nx_in, Ny_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out

#ifdef __FFTW
      ! FFTW Plan
      INTEGER*8                        :: Plan      
      INTEGER                          :: mbistride, mbrank, mbidist, mbiembed 
      INTEGER                          :: mboembed, mbhowmany, mbodist
#endif
>>>>>>> 1.13

<<<<<<< ppm_util_fft_backward_2d.f
=======
#ifdef __MATHKEISAN
      ! MATHKEISAN variables for MathKeisan FFTs
      INTEGER                                 :: isign_fft,isys
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION

#elif __KIND == __SINGLE_PRECISION_COMPLEX| __KIND == __DOUBLE_PRECISION_COMPLEX
      ! parameters for cfft (default=1)
      INTEGER                                 :: incx, incy
#endif

      ! scale_fft of the transformation
      REAL(MK)                                :: scale_fft
      ! working storage
      REAL(MK), DIMENSION(:),POINTER          :: table, work
      ! the size of the working storage
      INTEGER, DIMENSION(1)                   :: lda_table, lda_work

#endif

>>>>>>> 1.13
      !-------------------------------------------------------------------------
      ! Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_util_fft_backward_2d',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ppm_field_topoid .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_backward_2d',  &
     &            'No field topology has been defined so far',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if FFTW-Library is available of if NEC Library
      !-------------------------------------------------------------------------
#if  !(defined(__FFTW) | defined(__MATHKEISAN))
      info = ppm_error_error

#ifndef __FFTW
              CALL ppm_error(ppm_err_nofftw,'ppm_util_fft_backward_2d',  &
     &            'fdsolver needs FFTW Library  ' ,__LINE__,info)
#endif

#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_util_fft_backward_2d',  &
     &    'PPM was compiled without MATHKEISAN support',__LINE__,info)
#endif
      GOTO 9999
#else

      !-------------------------------------------------------------------------
      ! Allocate complex array
      !-------------------------------------------------------------------------
      yhmax = 0
      DO i=1,ppm_nsublist(topo_id_internal)
          idom = ppm_isublist(i,topo_id_internal)
          IF (ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes(2,idom)  &
    &           .GT.yhmax) THEN
              yhmax = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes(2,idom)
          ENDIF
      ENDDO

#if   __DIM == __SFIELD
      lda_DATA_fv_com(1)= Nm_com(1)
      lda_DATA_fv_com(2)= yhmax   
      lda_DATA_fv_com(3)= nsublist
#elif __DIM == __VFIELD
      lda_DATA_fv_com(1)= lda_fv
      lda_DATA_fv_com(2)= Nm_com(1)
      lda_DATA_fv_com(3)= yhmax   
      lda_DATA_fv_com(4)= nsublist
#endif   
    
      iopt = ppm_param_alloc_fit
      CALL ppm_alloc(DATA_fv_com, lda_DATA_fv_com, iopt,info)

      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_backward_2d',     &
     &        'data array',__LINE__,info)
          GOTO 9999
<<<<<<< ppm_util_fft_backward_2d.f
      ENDIF     
      
=======
      ENDIF
     

      !-------------------------------------------------------------------------
      !  NEC version - Use MathKeisan Library 1.5
      !-------------------------------------------------------------------------




#ifdef __MATHKEISAN

>>>>>>> 1.13
      !-------------------------------------------------------------------------
      !  FFT - Transformation in x-direction
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      lda(1)=2
      lda(2)=ppm_nsubs(topo_id_internal)
      CALL ppm_alloc(ndata,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_backward_2d',     &
     &        'ndata array',__LINE__,info)
          GOTO 9999
      ENDIF     
      ndata = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes

      DO k=1,ppm_nsublist(topo_id_internal)
          idom = ppm_isublist(k,topo_id_internal)
          CALL  ppm_alloc(data_in, ndata(:,idom), ppm_param_alloc_fit, info)

          IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_util_fft_backward_2d',     &
     &        'data_in array',__LINE__,info)
          GOTO 9999
          ENDIF     
  
#if  __DIM == __VFIELD       
      DO n=1,lda_fv
#endif

          DO i=1, ndata(1,idom)
            DO j=1, ndata(2,idom)

#if  __DIM == __SFIELD
              data_in(i,j) =  DATA_fv(i,j,k)
#elif __DIM == __VFIELD
              data_in(i,j) =  DATA_fv(n,i,j,k)
#endif

            ENDDO
          ENDDO
   
        CALL ppm_util_fft_step_fd( data_in, ndata(:,idom), FFT_x, info)

         iopt = ppm_param_dealloc
         CALL  ppm_alloc(data_in, ndata(:,idom), iopt, info)

         DO i=1, ndata(1,idom)
           DO j=1, ndata(2,idom)
#if  __DIM == __SFIELD
              DATA_fv_com(i,j,k) = FFT_x(i,j)
#elif __DIM == __VFIELD
              DATA_fv_com(n,i,j,k) = FFT_x(i,j)
#endif


<<<<<<< ppm_util_fft_backward_2d.f
           ENDDO
         ENDDO
       ENDDO
#if __DIM == __VFIELD
      ENDDO
#endif     
=======
#if   __KIND == __SINGLE_PRECISION
      CALL  csfft(isign_fft, Nx_out, scale_fft, data_in(1,j), &
     &                                 data_out(1,j), table, work, isys)
#elif __KIND == __DOUBLE_PRECISION
      CALL   zdfft(isign_fft, Nx_out, scale_fft, data_in(1,j), &
     &                                 data_out(1,j), table, work, isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      incx = 1
      incy = 1
      CALL  cfft(isign_fft, Nx_out, scale_fft, data_in(1,j), incx,        &
     &              data_out(1,j), incy, table, lda_table(1), work, &
     &              lda_work(1), isys)

#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      incx = 1
      incy = 1
      CALL  zfft(isign_fft, Nx_out, scale_fft, data_in(1,j), incx,        &
     &              data_out(1,j), incy, table, lda_table(1), work, &
     &              lda_work(1), isys)
>>>>>>> 1.13

      !-------------------------------------------------------------------------
      !  Transpose x-direction and y-direction
      !-------------------------------------------------------------------------
      mesh_ids_tmp(1) = mesh_id_xpen_complex
      mesh_ids_tmp(2) = mesh_id_ypen
#if __DIM == __SFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,topo_ids_tmp, mesh_ids_tmp, info)
#elif __DIM == __VFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,lda_fv, topo_ids_tmp, mesh_ids_tmp, info)
#endif
      DO k=1,ppm_nsublist(ppm_field_topoid)
          idom = ppm_isublist(k,ppm_field_topoid)     
          lda(1)=2
          lda(2)=idom
          iopt = ppm_param_alloc_fit
          CALL ppm_alloc(ndata_trans,lda,iopt, info)
          
          ndata_trans(1,idom)=ndata(2,idom)
          ndata_trans(2,idom)=ndata(1,idom)
          
          CALL ppm_alloc(data_com,ndata_trans(:,idom),iopt, info)
#if __DIM == __VFIELD
         DO n=1,lda_fv
#endif
         DO i=1, ndata(1,idom)
              DO j=1, ndata(2,idom)
#if __DIM == __SFIELD
                  data_com(j,i)= DATA_fv_com(i,j,k)
#elif __DIM == __VFIELD
                  data_com(j,i)= DATA_fv_com(n,i,j,k)
#endif
              ENDDO
           ENDDO
           
      !-------------------------------------------------------------------------
      !  FFT - Transformation in y-direction
      !-------------------------------------------------------------------------
      CALL ppm_util_fft_step_fd(data_com,ndata_trans(:,idom),FFT_xy,info) 
      iopt = ppm_param_dealloc
      CALL ppm_alloc(data_com,ndata_trans(:,idom),iopt, info)
      
     














<<<<<<< ppm_util_fft_backward_2d.f
=======
#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft_c2r(Plan, MBRank,Nx_out, MBHowmany, &
     &  data_in(1,1), MBiEmbed,MBIstride,MBiDist, &
     &  data_out(1,1),MBoEmbed,MBIstride,MBoDist,FFTW_ESTIMATE)
      CALL sfftw_execute(Plan)
      CALL sfftw_destroy_plan(Plan)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft_c2r(Plan, MBRank,Nx_out, MBHowmany, &
     &  data_in(1,1), MBiEmbed,MBIstride,MBiDist, &
     &  data_out(1,1),MBoEmbed,MBIstride,MBoDist,FFTW_ESTIMATE)
      CALL dfftw_execute(Plan)
      CALL dfftw_destroy_plan(Plan)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      CALL sfftw_plan_many_dft(Plan, MBRank,Nx_out, MBHowmany, &
     &  data_in(1,1), MBiEmbed,MBIstride,MBiDist, &
     &  data_out(1,1),MBoEmbed,MBIstride,MBoDist, &
     &  FFTW_BACKWARD, FFTW_ESTIMATE)
      CALL sfftw_execute(Plan)
      CALL sfftw_destroy_plan(Plan)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      CALL dfftw_plan_many_dft(Plan, MBRank,Nx_out, MBHowmany, &
     &  data_in(1,1), MBiEmbed,MBIstride,MBiDist, &
     &  data_out(1,1),MBoEmbed,MBIstride,MBoDist, &
     &  FFTW_BACKWARD, FFTW_ESTIMATE)
      CALL dfftw_execute(Plan)
      CALL dfftw_destroy_plan(Plan)
#endif
>>>>>>> 1.13


<<<<<<< ppm_util_fft_backward_2d.f
      !-------------------------------------------------------------------------
      !  FFT - Backward Transformation in y-direction
      !-------------------------------------------------------------------------
      CALL ppm_util_fft_step_bd(FFT_xy, ndata_trans(:,idom), data_com,info) 
=======
#endif
#endif

>>>>>>> 1.13

      !-------------------------------------------------------------------------
      !  Transpose y-direction and x-direction
      !-------------------------------------------------------------------------
<<<<<<< ppm_util_fft_backward_2d.f
          DO i=1, ndata(1,idom)
              DO j=1, ndata(2,idom)
#if __DIM == __SFIELD
                  DATA_fv_com(i,j,k) = data_com(j,i)
#elif __DIM == __VFIELD
                  DATA_fv_com(n,i,j,k) = data_com(j,i)
#endif
              ENDDO
          ENDDO
#if __DIM == __VFIELD
        ENDDO
#endif
          iopt = ppm_param_dealloc
          CALL ppm_alloc(data_com, lda, iopt,info)

      ENDDO ! end of do loop over k=1,nsublist

      topo_id    = topo_ids_tmp(1)
      topo_ids_tmp(1)= topo_ids_tmp(2)
      topo_ids_tmp(2)= topo_id
      mesh_id    = mesh_ids_tmp(1)
      mesh_ids_tmp(1)= mesh_ids_tmp(2)
      mesh_ids_tmp(2)= mesh_id

#if __DIM == __SFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,topo_ids_tmp, mesh_ids_tmp, info)
#elif __DIM == __VFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,lda_fv,topo_ids_tmp, mesh_ids_tmp, info)
#endif

      DO k=1,ppm_nsublist(ppm_field_topoid)
          idom = ppm_isublist(k,ppm_field_topoid)
#if __DIM == __VFIELD
       DO n=1, lda_fv
#endif
           DO i=1, ndata_xpen_complex(1,idom)
              DO j=1, ndata_xpen_complex(2,idom)
#if __DIM == __SFIELD
                  FFT_x(i,j)= DATA_fv_com(i,j,k)
#elif __DIM == __VFIELD
                  FFT_x(i,j)= DATA_fv_com(n,i,j,k)
#endif
              ENDDO
          ENDDO

      !-------------------------------------------------------------------------
      !  FFT - Backward Transformation in y-direction
      !-------------------------------------------------------------------------
      lda(1) = ndata_xpen_complex(1,idom)
      lda(2) = ndata_xpen_complex(2,idom)
      CALL ppm_util_fft_step_bd(FFT_x,lda,Result,info) 

      !-------------------------------------------------------------------------
      ! Correct Inverse by problem size factor 1/(Nx*Ny)     
      ! Subtract 1 to fit ppm convention
      !-------------------------------------------------------------------------
        rN = 1/dble((Nm(1)-1)*(Nm(2)-1))
         DO i=1,ndata_xpen_complex(1,idom)
           DO j=1,ndata_xpen_complex(2,idom)
#if __DIM == __SFIELD
              DATA_fv(i,j,k)= Result(i,j)*rN
#elif __DIM == __VFIELD
              DATA_fv(n,i,j,k)= Result(i,j)*rN
#endif
           ENDDO
         ENDDO
#if __DIM == __VFIELD
       ENDDO
#endif
      ENDDO ! end of do loop k=1,nsublist

      !-------------------------------------------------------------------------
      ! Map to original topology if not x-pencil topology
      !-------------------------------------------------------------------------
      IF(.NOT.Its_xpencil_topo) THEN
         topo_ids_tmp(1) = topo_ids_tmp(2)
         topo_ids_tmp(2) = topo_id_user
         mesh_ids_tmp(1) = mesh_id_xpen
         mesh_ids_tmp(2) = mesh_id_user

#if __DIM == __SFIELD
         CALL ppm_fdsolver_map(DATA_fv, topo_ids_tmp, mesh_ids_tmp, info)
#elif __DIM == __VFIELD
         CALL ppm_fdsolver_map(DATA_fv, lda_fv, topo_ids_tmp, mesh_ids_tmp, info)
#endif
      ENDIF 

      !-------------------------------------------------------------------------
      ! Deallocate memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(DATA_fv_com, lda_DATA_fv_com, iopt,info)
      CALL ppm_alloc(ndata,lda,iopt,info)
      CALL ppm_alloc(ndata_trans,lda,iopt, info)
      CALL ppm_alloc(data_in,lda,iopt, info)
      CALL ppm_alloc(FFT_x,lda,iopt, info)
      CALL ppm_alloc(FFT_xy,lda,iopt, info)
      CALL ppm_alloc(min_sub,lda,iopt, info)
      CALL ppm_alloc(max_sub,lda,iopt, info)
      CALL ppm_alloc(cost,lda,iopt, info)
      CALL ppm_alloc(istart,lda,iopt, info)
      CALL ppm_alloc(istart_xpen_complex,lda,iopt, info)
      CALL ppm_alloc(istart_ypen,lda,iopt, info)
      CALL ppm_alloc(istart_trans,lda,iopt, info)
      CALL ppm_alloc(ndata,lda,iopt, info)
      CALL ppm_alloc(ndata_xpen_complex,lda,iopt, info)
      CALL ppm_alloc(ndata_ypen,lda,iopt, info)
      CALL ppm_alloc(ndata_trans,lda,iopt, info)
      CALL ppm_alloc(sub2proc,lda,iopt, info)
      CALL ppm_alloc(isublist,lda,iopt, info)
#endif
=======
      DO j=1,Ny_out
          data_out(lda(1),j) = data_out(1,j)
      ENDDO     
>>>>>>> 1.13

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_fft_backward_2d',t0,info)
      RETURN

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_backward_2d_ss
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_backward_2d_sd
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_backward_2d_vs
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_backward_2d_vd
#endif
#endif

