      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_util_fft_step_bd_2d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs a single FFT step (i.e. pencil
      !                 or slab) using the precomputed plans from
      !                 ppm_util_fft_init.
      !
      !  Input        : data_in(:,:)   (F) data to be transformed 
      !
      !  Input/output : lda(:)         (I) size of data
      !
      !  Output       : data_out(:,:)  (F) transformed data
      !                 info           (I) return status. =0 if no error.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_fft_step_bd_2d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:02  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_step_bd_2ds(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_step_bd_2dd(data_in,lda,data_out,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_util_fft_step_bd_2dc(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_util_fft_step_bd_2dcc(data_in,lda,data_out,info)
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
#elif __KIND == __DOUBLE_PRECISION | __KIND ==__DOUBLE_PRECISION_COMPLEX 
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
      COMPLEX(MK), DIMENSION(:,:)   , INTENT(IN   ) :: data_in
      INTEGER, DIMENSION(:)         , INTENT(INOUT) :: lda
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:)      , POINTER       :: data_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:)   , POINTER       :: data_out 
#endif
      INTEGER                       , INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0
      INTEGER                                 :: i,j,iopt
      INTEGER                                 :: Nx_in, Ny_in
      INTEGER                                 :: Nx_out, Ny_out
#ifdef __FFTW
      INTEGER*8                               :: Plan
#endif
#ifdef __MATHKEISAN
      INTEGER                                 :: isign_fft,isys
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      INTEGER                                 :: incx, incy
#endif
      REAL(MK)                                :: scale_fft
      REAL(MK), DIMENSION(:),POINTER          :: table, work
      INTEGER, DIMENSION(1)                   :: lda_table, lda_work
#endif

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_fft_step_bd_2d',t0,info)

#if  !(defined(__FFTW) | defined(__MATHKEISAN))
      !-------------------------------------------------------------------------
      !  Error if library support is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
#ifndef __FFTW
      CALL ppm_error(ppm_err_nofftw,'ppm_util_fft_step_bd_2d',  &
     &    'PPM was compiled without fftw support',__LINE__,info)
#endif
#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_util_fft_step_bd_2d',  &
     &    'PPM was compiled without MATHKEISAN support',__LINE__,info)
#endif
      GOTO 9999      
#else

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (SIZE(lda,1) .LT. 2) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_step_bd_2d',  &
     &            'lda must be at least of size 2',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(1) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_step_bd_2d',  &
     &            'mesh size: lda(1) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(2) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_step_bd_2d',  &
     &            'mesh size: lda(2) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate result array
      !-------------------------------------------------------------------------
      Nx_in = lda(1)
      Ny_in = lda(2)
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      Nx_out = (Nx_in-1)*2
#elif __KIND == __SINGLE_PRECISION_COMPLEX |__KIND == __DOUBLE_PRECISION_COMPLEX
      Nx_out = Nx_in -1
#endif
      Ny_out = Ny_in
      lda(1) = Nx_out+1
      lda(2) = Ny_out
      CALL ppm_alloc(data_out,lda,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_step_bd_2d',     &
     &        'fft result DATA_OUT',__LINE__,info)
          GOTO 9999
      ENDIF
     
#ifdef __MATHKEISAN
      !-------------------------------------------------------------------------
      !  NEC version - Use MathKeisan Library
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      lda_work(1) = 4*Nx_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND== __DOUBLE_PRECISION_COMPLEX
      lda_work(1) = 6*Nx_out
#endif
      CALL ppm_alloc(work,lda_work,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_step_bd_2d',     &
     &        'work storage for Mathkeisan',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Forward FFT 
      !-------------------------------------------------------------------------
      scale_fft = 1
      isign_fft = 1
      DO j=1,Ny_in
#if   __KIND == __SINGLE_PRECISION
          CALL csfft(isign_fft,Nx_out,scale_fft,data_in(1,j), &
     &                                 data_out(1,j),table_bd_s,work,isys)
#elif __KIND == __DOUBLE_PRECISION

          CALL zdfft(isign_fft,Nx_out,scale_fft,data_in(1,j), &
     &                                 data_out(1,j),table_bd_d,work,isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
          CALL cfft(isign_fft,Nx_out,scale_fft,data_in(1,j),incx,     &
     &              data_out(1,j),incy,table_fd_c_y,lda_table_y,work, &
     &              lda_work(1),isys)

#elif __KIND == __DOUBLE_PRECISION_COMPLEX
          CALL zfft(isign_fft,Nx_out,scale_fft,data_in(1,j),incx,     &
     &              data_out(1,j),incy,table_fd_c_y,lda_table_y,work, &
     &              lda_work(1),isys)
#endif
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate Memory
      !-------------------------------------------------------------------------
      CALL ppm_alloc(work,lda_work,ppm_param_dealloc,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_bd_2d',    &
     &         'Mathkeisan WORK memory',__LINE__,info)
      ENDIF

#else

      !-------------------------------------------------------------------------
      !  FFTW version
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_execute_dft_c2r(Plan_bd_s,data_in(1,1),data_out(1,1))
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_execute_dft_c2r(Plan_bd_d,data_in(1,1),data_out(1,1))
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      CALL sfftw_execute_dft(Plan_bd_c_y,data_in(1,1),data_out(1,1))
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      CALL dfftw_execute_dft(Plan_bd_cc_y,data_in(1,1),data_out(1,1))
#endif

#endif
#endif

      !-------------------------------------------------------------------------
      !  Copy margin to conform with PPM convention
      !-------------------------------------------------------------------------
      DO j=1,Ny_out
          data_out(lda(1),j) = data_out(1,j)
      ENDDO     

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_fft_step_bd_2d',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_step_bd_2ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_step_bd_2dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_util_fft_step_bd_2dc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_util_fft_step_bd_2dcc
#endif

