      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_fft_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Finalizes the FFT module and deallocates all memory.
      !
      !  Input        :  
      !
      !  Input/output : 
      !
      !  Output       : info    (I) Return status. 0 on success.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_fft_finalize.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_util_fft_finalize(info)

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data_util_fft
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_alloc

      IMPLICIT NONE
#ifdef  __FFTW
      INCLUDE "fftw3.f"
#endif

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                    , INTENT(  OUT)   :: info

      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)                        :: t0
      INTEGER, DIMENSION(3)                        :: lda
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_fft_finalize',t0,info)
     
      !-------------------------------------------------------------------------
      !  Deallocate Memory of the FFT solver
      !-------------------------------------------------------------------------
#ifdef __FFTW
      CALL sfftw_destroy_plan(Plan_fd_s)
      CALL sfftw_destroy_plan(Plan_fd_c_y)
      CALL sfftw_destroy_plan(Plan_fd_c_z)
      CALL sfftw_destroy_plan(Plan_bd_s)
      CALL sfftw_destroy_plan(Plan_bd_c_y)
      CALL sfftw_destroy_plan(Plan_bd_c_z)

      CALL sfftw_destroy_plan(Plan_slab_fd_s)
      CALL sfftw_destroy_plan(Plan_slab_bd_s)

      CALL dfftw_destroy_plan(Plan_fd_d)
      CALL dfftw_destroy_plan(Plan_fd_cc_y)
      CALL dfftw_destroy_plan(Plan_fd_cc_z)
      CALL dfftw_destroy_plan(Plan_bd_d)
      CALL dfftw_destroy_plan(Plan_bd_cc_y)
      CALL dfftw_destroy_plan(Plan_bd_cc_z)


      CALL dfftw_destroy_plan(Plan_slab_fd_d)
      CALL dfftw_destroy_plan(Plan_slab_bd_d)
#endif

#ifdef __MATHKEISAN
      iopt = ppm_param_dealloc

      CALL ppm_alloc(table_fd_s,lda,iopt,info)
      CALL ppm_alloc(table_fd_d,lda,iopt,info)
      CALL ppm_alloc(table_fd_c_y,lda,iopt,info)
      CALL ppm_alloc(table_fd_cc_y,lda,iopt,info)
      CALL ppm_alloc(table_fd_c_z,lda,iopt,info)
      CALL ppm_alloc(table_fd_cc_z,lda,iopt,info)

      CALL ppm_alloc(table_bd_s,lda,iopt,info)
      CALL ppm_alloc(table_bd_d,lda,iopt,info)
      CALL ppm_alloc(table_bd_c_y,lda,iopt,info)
      CALL ppm_alloc(table_bd_cc_y,lda,iopt,info)
      CALL ppm_alloc(table_bd_c_z,lda,iopt,info)
      CALL ppm_alloc(table_bd_cc_z,lda,iopt,info)


      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_finalize',   &
     &    'Failed to deallocate memory.',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Later: Destroy the temporary topologies here. As soon as we have
      !         destrictible topologies.
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_fft_finalize',t0,info)
      RETURN

      END SUBROUTINE ppm_util_fft_finalize

