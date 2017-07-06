      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_util_fft_step_bd_3d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs a single FFT step (i.e. pencil
      !                 or slab) using the precomputed plans from
      !                 ppm_util_fft_init.
      !
      !  Input        : data_in(:,:,:)  (F) data to be transformed 
      !                                
      !  Input/output : lda(:)          (I) size of data               
      !
      !  Output       : data_out(:,:,:) (F) transformed data
      !                 info            (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_fft_step_bd_3d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:02  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if   __CASE == __SLAB

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_slab_bd_3ds(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_slab_bd_3dd(data_in,lda,data_out,info)
#endif

#else

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_step_bd_3ds(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_step_bd_3dd(data_in,lda,data_out,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_util_fft_step_bd_3dc(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_util_fft_step_bd_3dcc(data_in,lda,data_out,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z
      SUBROUTINE ppm_util_fft_z_bd_3dc(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX_Z
      SUBROUTINE ppm_util_fft_z_bd_3dcc(data_in,lda,data_out,info)
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
      USE ppm_module_data_util_fft
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc

      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND ==__SINGLE_PRECISION_COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif   __KIND ==__SINGLE_PRECISION_COMPLEX_Z
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND ==__DOUBLE_PRECISION_COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_double
#elif __KIND ==__DOUBLE_PRECISION_COMPLEX_Z 
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
#ifdef  __FFTW
      INCLUDE "fftw3.f"
#endif

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      COMPLEX(MK), DIMENSION(:,:,:) , INTENT(IN   ) :: data_in
      INTEGER, DIMENSION(:)         , INTENT(INOUT) :: lda
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION 
      REAL(MK), DIMENSION(:,:,:)    , POINTER       :: data_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:) , POINTER       :: data_out 
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z | __KIND == __DOUBLE_PRECISION_COMPLEX_Z
      COMPLEX(MK), DIMENSION(:,:,:) , POINTER       :: data_out 
#endif
      INTEGER                       , INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                         :: t0
      INTEGER                          :: i,j,k,iopt
      INTEGER                          :: Nx_in,Ny_in,Nz_in
      INTEGER                          :: Nx_out,Ny_out,Nz_out
#ifdef __FFTW
      INTEGER*8                        :: Plan
      INTEGER                          :: mbistride,mbrank,mbidist,mbiembed 
      INTEGER                          :: mboembed,mbhowmany,mbodist
#endif
#ifdef __MATHKEISAN
      INTEGER                          :: isign_fft,isys
#if __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      INTEGER                          :: incx,incy
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z | __KIND == __DOUBLE_PRECISION_COMPLEX_Z
      INTEGER                          :: incx,incy
#endif
      REAL(MK)                         :: scale_fft
      REAL(MK), DIMENSION(:),POINTER   :: table,work
      INTEGER, DIMENSION(1)            :: lda_table,lda_work
#endif

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_fft_step_bd_3d',t0,info)

#if  !(defined(__FFTW) | defined(__MATHKEISAN))
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
#ifndef __FFTW
      CALL ppm_error(ppm_err_nofftw,'ppm_util_fft_step_bd_3d',  &
     &    'PPM was compiled without fftw support',__LINE__,info)
#endif
#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_util_fft_step_bd_3d',  &
     &    'PPM was compiled without MATHKEISAN support',__LINE__,info)
#endif
      GOTO 9999      
#else

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (SIZE(lda,1) .LT. 3) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_step_bd_3d',  &
     &            'lda must be at least of size 3',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(1) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_step_bd_3d',  &
     &            'mesh size: lda(1) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(2) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_step_bd_3d',  &
     &            'mesh size: lda(2) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(3) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_step_bd_3d',  &
     &            'mesh size: lda(3) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate result array
      !-------------------------------------------------------------------------
      Nx_in = lda(1)
      Ny_in = lda(2)
      Nz_in = lda(3)
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      Nx_out = (Nx_in-1)*2
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      Nx_out = Nx_in-1       
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z | __KIND == __DOUBLE_PRECISION_COMPLEX_Z
      Nx_out = Nx_in-1       
#endif

      Ny_out=Ny_in
      Nz_out=Nz_in

      lda(1)=Nx_out +1 ! to fit ppm-convention
      lda(2)=Ny_out
      lda(3)=Nz_out
      CALL ppm_alloc(data_out,lda,ppm_param_alloc_fit,info )
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_step_bd_3d',     &
     &        'fft result DATA_OUT',__LINE__,info)
          GOTO 9999
      ENDIF

#ifdef __MATHKEISAN
      !-------------------------------------------------------------------------
      !  NEC version - Use MathKeisan Library
      !-------------------------------------------------------------------------

#if __CASE == __SLAB
      !-------------------------------------------------------------------------
      !  not implemented yet
      !-------------------------------------------------------------------------
      info = ppm_error_fatal
      CALL ppm_error(ppm_err_alloc,'ppm_util_fft_step_bd_3d',     &
     &    'version not implemented',__LINE__,info)
      GOTO 9999
#else

      !-------------------------------------------------------------------------
      !  Allocate working storage
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      lda_work    = 4*Nx_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND== __DOUBLE_PRECISION_COMPLEX
      lda_work(1) = 6*Nx_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z | __KIND== __DOUBLE_PRECISION_COMPLEX_Z
      lda_work(1) = 6*Nx_out
#endif
      CALL ppm_alloc(work,lda_work,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_step_bd_3d',     &
     &        'work not allocated',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Forward FFT 
      !-------------------------------------------------------------------------
      scale_fft = 1
      isign_fft = 1

      DO k=1,Nz_in
        DO j=1,Ny_in
#if   __KIND == __SINGLE_PRECISION
           CALL csfft(isign_fft,Nx_out,scale_fft,data_in(1,j,k), &
     &                                data_out(1,j,k),table_bd_s,work,isys)
#elif __KIND == __DOUBLE_PRECISION
           CALL zdfft(isign_fft,Nx_out,scale_fft,data_in(1,j,k), &
     &                                data_out(1,j,k),table_bd_d,work,isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
           CALL cfft(isign_fft,Nx_out,scale_fft,data_in(1,j,k),incx, &
     &               data_out(1,j,k),incy,table_fd_c_y,lda_table_y, & 
     &               work,lda_work(1),isys)

#elif __KIND == __DOUBLE_PRECISION_COMPLEX
           CALL zfft(isign_fft,Nx_out,scale_fft,data_in(1,j,k),incx, &
     &               data_out(1,j,k),incy,table_fd_cc_y,lda_table_y, & 
     &               work,lda_work(1),isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z
           CALL cfft(isign_fft,Nx_out,scale_fft,data_in(1,j,k),incx, &
     &               data_out(1,j,k),incy,table_fd_c_z,lda_table_z, & 
     &               work,lda_work(1),isys)

#elif __KIND == __DOUBLE_PRECISION_COMPLEX_Z
           CALL zfft(isign_fft,Nx_out,scale_fft,data_in(1,j,k),incx, &
     &               data_out(1,j,k),incy,table_fd_cc_z,lda_table_z, & 
     &               work,lda_work(1),isys)
#endif
         ENDDO      
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate Memory
      !-------------------------------------------------------------------------
      CALL ppm_alloc(work,lda_work,ppm_param_dealloc,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_step_bd_3d',   &
     &        'Mathkeisan WORK memory',__LINE__,info)
      ENDIF
#endif

#else

      !-------------------------------------------------------------------------
      !  FFTW version
      !-------------------------------------------------------------------------
#if __CASE == __SLAB

#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_execute_dft_c2r(Plan_slab_bd_s,data_in(1,1,1),data_out(1,1,1))
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_execute_dft_c2r(Plan_slab_bd_d,data_in(1,1,1),data_out(1,1,1))
#endif

#else

#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_execute_dft_c2r(Plan_bd_s,data_in(1,1,1),data_out(1,1,1))
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_execute_dft_c2r(Plan_bd_d,data_in(1,1,1),data_out(1,1,1))
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      CALL sfftw_execute_dft(Plan_bd_c_y,data_in(1,1,1),data_out(1,1,1))
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      CALL dfftw_execute_dft(Plan_bd_cc_y,data_in(1,1,1),data_out(1,1,1))
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z
      CALL sfftw_execute_dft(Plan_bd_c_z,data_in(1,1,1),data_out(1,1,1))
#elif __KIND == __DOUBLE_PRECISION_COMPLEX_Z
      CALL dfftw_execute_dft(Plan_bd_cc_z,data_in(1,1,1),data_out(1,1,1))
#endif

#endif
#endif
#endif

      !-------------------------------------------------------------------------
      !  Copy margin to conform with PPM convention
      !-------------------------------------------------------------------------
      DO j=1,Ny_out
         DO k=1,Nz_out
            data_out(lda(1),j,k) = data_out(1,j,k)
         ENDDO
      ENDDO     

#if __CASE == __SLAB
      DO i=1,lda(1)
         DO k=1,lda(3)
            data_out(i,lda(2),k) = data_out(i,1,k)
         ENDDO
      ENDDO     
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_fft_step_bd_3d',t0,info)
      RETURN

#if   __CASE == __SLAB

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_slab_bd_3ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_slab_bd_3dd
#endif

#else

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_step_bd_3ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_step_bd_3dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_util_fft_step_bd_3dc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_util_fft_step_bd_3dcc
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z
      END SUBROUTINE ppm_util_fft_z_bd_3dc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX_Z
      END SUBROUTINE ppm_util_fft_z_bd_3dcc
#endif
#endif

