      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_util_fft_pencil_fd_2d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs Fast Fourier Transform using
      !                 FFTW in the first (x) dimension
      !
      !  Input        : data_in(:,:)   (F) data to be transformed
      !
      !  Input/output : lda(:)         (I) size of data
      !
      !
      !  Output       : data_out(:,:)  (F) transformed data
      !                 info           (I) return status. =0 if no error.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_fft_pencil_fd_2d.f,v $
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
      SUBROUTINE ppm_util_fft_pencil_fd_2ds(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_pencil_fd_2dd(data_in,lda,data_out,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_util_fft_pencil_fd_2dc(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_util_fft_pencil_fd_2dcc(data_in,lda,data_out,info)
#endif

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
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
      ! input data
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:)      , INTENT(IN   ) :: data_in
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:)       , INTENT(IN   ) :: data_in
#endif
      ! size of the array
      INTEGER, DIMENSION(:)         , INTENT(INOUT) :: lda
      ! output data, fast fourier transformed
      COMPLEX(MK), DIMENSION(:,:)   , POINTER       :: data_out
      INTEGER                       , INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
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
      INTEGER*8                               :: Plan
      INTEGER                            :: mbistride, mbrank, mbidist, mbiembed
      INTEGER                            :: mboembed, mbhowmany, mbodist
#endif

#ifdef __MATHKEISAN
      ! MATHKEISAN variables for MathKeisan FFTs
      INTEGER                                 :: isign_fft, isys
#if __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      INTEGER                                 :: incx, incy
#endif
      ! scale of the transformation
      REAL(MK)                                :: scale_fft
      ! working storage
      REAL(MK), DIMENSION(:),POINTER          :: table, work
      ! the size of the working storage
      INTEGER, DIMENSION(1)                   :: lda_table, lda_work
#endif


      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_fft_pencil_fd_2d',t0,info)

#ifdef __MATHKEISAN
      ! MATHKEISAN variables for MathKeisan FFTs
      NULLIFY(table,work)
#endif

#if  !(defined(__FFTW) | defined(__MATHKEISAN))

      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error

#ifndef __FFTW
      CALL ppm_error(ppm_err_nofftw,'ppm_util_fft_pencil_fd_2d',  &
     &    'PPM was compiled without fftw support',__LINE__,info)
#endif

#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_util_fft_pencil_fd_2d',  &
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
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_pencil_fd_2d',  &
     &            'lda must be at least of size 2',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(1) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_pencil_fd_2d',  &
     &            'mesh size: lda(1) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(2) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_pencil_fd_2d',  &
     &            'mesh size: lda(2) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      ! subtract 1 to fit ppm-convention
      Nx_in = lda(1)-1
      Ny_in = lda(2)

      !-------------------------------------------------------------------------
      !  Allocate result array
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      Nx_out = Nx_in/2 + 1
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      Nx_out = Nx_in
#endif

      Ny_out=Ny_in

#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      lda(1) = Nx_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      lda(1) = Nx_out+1      ! to fit ppm-convention
#endif
      lda(2)=Ny_out
      CALL ppm_alloc(data_out,lda,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_pencil_fd_2d',     &
     &        'fft result DATA_OUT',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  FFT - Transform in x-direction
      !-------------------------------------------------------------------------


      !-------------------------------------------------------------------------
      !  NEC version - Use MathKeisan Library
      !-------------------------------------------------------------------------
#ifdef __MATHKEISAN

      !-------------------------------------------------------------------------
      !NOTE: This fft routines does not currently work  for  arbitrary n.
      ! Number of data points, of the	form n = 2**p  *  3**q	*
      !	  5**r,	with p,	q, r >=	0
      !-------------------------------------------------------------------------


      !-------------------------------------------------------------------------
      !  Allocate working storage
      !-------------------------------------------------------------------------

      lda_table(1) = 2*Nx_in + 64
      CALL ppm_alloc(table,lda_table,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_pencil_fd_2d',     &
     &        'table not allocated',__LINE__,info)
          GOTO 9999
      ENDIF

      lda_work(1) = 4*Nx_in
      CALL ppm_alloc(work,lda_work,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_pencil_fd_2d',     &
     &        'work not allocated',__LINE__,info)
          GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      scale_fft = 1
      isign_fft = 0
      j     = 1

#if   __KIND == __SINGLE_PRECISION
      CALL  scfft(isign_fft, Nx_in, scale_fft, data_in(1,j), &
     &                                 data_out(1,j), table, work, isys)
#elif __KIND == __DOUBLE_PRECISION
      CALL  dzfft(isign_fft, Nx_in, scale_fft, data_in(1,j), &
     &                                 data_out(1,j), table, work, isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      incx = 1
      incy = 1
      CALL  cfft(isign_fft, Nx_in, scale_fft, data_in(1,j), incx,        &
     &              data_out(1,j), incy, table, lda_table(1), work, &
     &              lda_work(1), isys)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      incx = 1
      incy = 1
      CALL  zfft(isign_fft, Nx_in, scale_fft, data_in(1,j), incx,        &
     &              data_out(1,j), incy, table, lda_table(1), work, &
     &              lda_work(1), isys)
#endif

      !-------------------------------------------------------------------------
      !  Forward FFT
      !-------------------------------------------------------------------------

      isign_fft = -1

      DO j=1,Ny_in

#if   __KIND == __SINGLE_PRECISION
          CALL  scfft(isign_fft, Nx_in, scale_fft, data_in(1,j), &
     &                                 data_out(1,j), table, work, isys)
#elif __KIND == __DOUBLE_PRECISION
          CALL  dzfft(isign_fft, Nx_in, scale_fft, data_in(1,j), &
     &                                 data_out(1,j), table, work, isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
          CALL  cfft(isign_fft, Nx_in, scale_fft, data_in(1,j), incx,        &
     &              data_out(1,j), incy, table, lda_table(1), work, &
     &              lda_work(1), isys)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
          CALL  zfft(isign_fft, Nx_in, scale_fft, data_in(1,j), incx,        &
     &              data_out(1,j), incy, table, lda_table(1), work, &
     &              lda_work(1), isys)
#endif

      ENDDO


      !-------------------------------------------------------------------------
      ! Deallocate working storage
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(table, lda, iopt,info)
      CALL ppm_alloc(work,lda,iopt,info)





#else

      !-------------------------------------------------------------------------
      !  FFTW version for LINUX,...
      !-------------------------------------------------------------------------


      MBRank    = 1
      MBHowmany = Ny_in
      MBiEmbed  = -1
      MBoEmbed  = -1
      MBIstride = 1
      MBiDist    = UBOUND(data_in, 1)
      MBoDist    = UBOUND(data_out,1)


#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft_r2c(Plan, MBRank,Nx_in, MBHowmany, &
           & data_in(1,1), MBiEmbed,MBIstride,MBiDist, &
           & data_out(1,1),MBoEmbed,MBIstride,MBoDist,FFTW_ESTIMATE)
          CALL sfftw_execute(Plan)
          CALL sfftw_destroy_plan(Plan)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft_r2c(Plan, MBRank,Nx_in, MBHowmany, &
           & data_in(1,1), MBiEmbed,MBIstride,MBiDist, &
           & data_out(1,1),MBoEmbed,MBIstride,MBoDist,FFTW_ESTIMATE)
          CALL dfftw_execute(Plan)
          CALL dfftw_destroy_plan(Plan)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      CALL sfftw_plan_many_dft(Plan, MBRank,Nx_in, MBHowmany, &
           & data_in(1,1), MBiEmbed,MBIstride,MBiDist, &
           & data_out(1,1),MBoEmbed,MBIstride,MBoDist, &
           & FFTW_FORWARD, FFTW_ESTIMATE)
          CALL sfftw_execute(Plan)
          CALL sfftw_destroy_plan(Plan)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      CALL dfftw_plan_many_dft(Plan, MBRank,Nx_in, MBHowmany, &
           & data_in(1,1), MBiEmbed,MBIstride,MBiDist, &
           & data_out(1,1),MBoEmbed,MBIstride,MBoDist, &
           & FFTW_FORWARD, FFTW_ESTIMATE)
          CALL dfftw_execute(Plan)
          CALL dfftw_destroy_plan(Plan)
#endif




#endif
#endif

#if __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      !-------------------------------------------------------------------------
      !  Copy margin to conform with PPM convention
      !-------------------------------------------------------------------------
      DO j=1,Ny_out
            data_out(lda(1),j) = data_out(1,j)
      ENDDO
#endif


      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_fft_pencil_fd_2d',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_pencil_fd_2ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_pencil_fd_2dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_util_fft_pencil_fd_2dc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_util_fft_pencil_fd_2dcc
#endif
