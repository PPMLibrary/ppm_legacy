      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_util_fft_forward_3d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs Fast Fourier Transform using
      !                 FFTW in the first (x) dimension
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
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_forward_3ds(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_forward_3dd(data_in,lda,data_out,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_util_fft_forward_3dc(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_util_fft_forward_3dcc(data_in,lda,data_out,info)
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
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:,:)    , INTENT(IN   ) :: data_in
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:) , INTENT(IN   ) :: data_in
#endif
      ! size of array
      INTEGER, DIMENSION(3)           , INTENT(INOUT) :: lda
      ! output data, fast fourier transformed
      COMPLEX(MK), DIMENSION(:,:,:)   , POINTER       :: data_out
      INTEGER                         , INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                :: t0
      ! counters
      INTEGER                                 :: i,j,k,iopt
      ! size of the data_in
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out
      INTEGER                                 :: Nx_out, Ny_out, Nz_out

#ifdef __FFTW
      ! FFTW Plan
      INTEGER*8                          :: Plan
      INTEGER                            :: mbistride, mbrank, mbidist, mbiembed
      INTEGER                            :: mboembed, mbhowmany, mbodist
#endif

#ifdef __MATHKEISAN
      ! MATHKEISAN variables for MathKeisan FFTs
      INTEGER                                 :: isign_fft,isys
#if __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      INTEGER                                 :: incx, incy
#endif
      !scale of the transformation
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
      CALL substart('ppm_util_fft_forward_3d',t0,info)

#ifdef __MATHKEISAN
      NULLIFY(table,work)
#endif

#if  !(defined(__FFTW) | defined(__MATHKEISAN))

      !-------------------------------------------------------------------------
      !  Error if FFT library support is not available
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

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (SIZE(lda,1) .LT. 3) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_forward_3d',  &
     &            'lda must be at least of size 3',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(1) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_forward_3d',  &
     &            'mesh size: lda(1) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(2) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_forward_3d',  &
     &            'mesh size: lda(2) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(3) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_forward_3d',  &
     &            'mesh size: lda(3) must be >0' ,__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF


      ! to fit ppm-convention
      Nx_in=lda(1)-1

      Ny_in=lda(2)
      Nz_in=lda(3)

      !-------------------------------------------------------------------------
      !  Allocate result array
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      Nx_out = Nx_in/2 + 1
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      Nx_out = Nx_in
#endif
      Ny_out = Ny_in
      Nz_out = Nz_in

#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      lda(1) = Nx_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      lda(1) = Nx_out+1      ! to fit ppm-convention
#endif


      lda(2) = Ny_out
      lda(3) = Nz_out
      CALL ppm_alloc(data_out,lda,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_forward_3d',     &
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
      !  Allocate working storage
      !-------------------------------------------------------------------------

      lda_table = 2*Nx_in + 64
      CALL ppm_alloc(table,lda_table,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_forward_2d',     &
     &        'table not allocated',__LINE__,info)
          GOTO 9999
      ENDIF

      lda_work = 4*Nx_in
      CALL ppm_alloc(work,lda_work,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_forward_2d',     &
     &        'work not allocated',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------

      scale_fft = 1
      isign_fft = 0
      j         = 1
      k         = 1


#if   __KIND == __SINGLE_PRECISION
      CALL  scfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), &
     &                                 data_out(1,j,k), table, work, isys)
#elif __KIND == __DOUBLE_PRECISION
      CALL  dzfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), &
     &                                 data_out(1,j,k), table, work, isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      incx = 1
      incy = 1
      CALL  cfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), incx, &
     &               data_out(1,j,k), incy, table, lda_table(1), &
     &              work, lda_work(1),isys)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      incx = 1
      incy = 1
      CALL  zfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), incx, &
     &               data_out(1,j,k), incy, table, lda_table(1), &
     &              work, lda_work(1),isys)
#endif


      !-------------------------------------------------------------------------
      !  Forward FFT
      !-------------------------------------------------------------------------

      isign_fft = -1

      DO k=1,Nz_in
         DO j=1,Ny_in


#if   __KIND == __SINGLE_PRECISION
            CALL  scfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), &
     &                                 data_out(1,j,k), table, work, isys)
#elif __KIND == __DOUBLE_PRECISION
            CALL  dzfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), &
     &                                 data_out(1,j,k), table, work, isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
            CALL  cfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), incx, &
     &               data_out(1,j,k), incy, table, lda_table(1), &
     &              work, lda_work(1),isys)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
            CALL  zfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), incx, &
     &               data_out(1,j,k), incy, table, lda_table(1), &
     &               work, lda_work(1),isys)
#endif

         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate Memory
      !-------------------------------------------------------------------------
      CALL ppm_alloc(table,lda_table,ppm_param_dealloc,info)
      CALL ppm_alloc(work,lda_work,ppm_param_dealloc,info)


#else

      !-------------------------------------------------------------------------
      !  FFTW version for LINUX,...
      !-------------------------------------------------------------------------


      MBRank    = 1
      MBHowmany = Ny_in*Nz_in
      MBiEmbed  = -1
      MBoEmbed  = -1
      MBIstride = 1
      MBiDist    = UBOUND(data_in, 1)
      MBoDist    = UBOUND(data_out,1)
#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft_r2c(Plan, MBRank,Nx_in, MBHowmany, &
           & data_in(1,1,1), MBiEmbed,MBIstride,MBiDist, &
           & data_out(1,1,1),MBoEmbed,MBIstride,MBoDist,FFTW_ESTIMATE)
      CALL sfftw_execute(Plan)
      CALL sfftw_destroy_plan(Plan)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft_r2c(Plan, MBRank,Nx_in, MBHowmany, &
           & data_in(1,1,1), MBiEmbed,MBIstride,MBiDist, &
           & data_out(1,1,1),MBoEmbed,MBIstride,MBoDist,FFTW_ESTIMATE)
      CALL dfftw_execute(Plan)
      CALL dfftw_destroy_plan(Plan)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      CALL sfftw_plan_many_dft(Plan, MBRank,Nx_in, MBHowmany, &
           & data_in(1,1,1), MBiEmbed,MBIstride,MBiDist, &
           & data_out(1,1,1),MBoEmbed,MBIstride,MBoDist, &
           & FFTW_FORWARD, FFTW_ESTIMATE)
      CALL sfftw_execute(Plan)
      CALL sfftw_destroy_plan(Plan)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      CALL dfftw_plan_many_dft(Plan, MBRank,Nx_in, MBHowmany, &
           & data_in(1,1,1), MBiEmbed,MBIstride,MBiDist, &
           & data_out(1,1,1),MBoEmbed,MBIstride,MBoDist, &
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
         DO k=1,Nz_out
            data_out(lda(1),j,k) = data_out(1,j,k)
         ENDDO
      ENDDO
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_fft_forward_3d',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_forward_3ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_forward_3dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_util_fft_forward_3dc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_util_fft_forward_3dcc
#endif
