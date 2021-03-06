      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_finalize
      !-------------------------------------------------------------------------
      !
<<<<<<< ppm_fdsolver_finalize.f
      !  Purpose      : finalizes fieldsolver by destroying all FFT-plans 
=======
      !  Purpose      : 
      !                 finalizes fieldsolver by destroying all FFT-plans     
>>>>>>> 1.5
      !
      !  Input        :  
      !
      !  Input/output :
      !
      !  Output       : info       (I) Return status. 0 on success.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_finalize.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.4  2005/08/15 08:31:59  hiebers
      !  Exchanged declaration REAL(8):: t0 by REAL(ppm_kind_double)
      !
      !  Revision 1.3  2005/02/18 08:02:12  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.2  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.1  2005/02/16 10:23:10  hiebers
      !  initial implementation
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_fdsolver_finalize(info)

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
<<<<<<< ppm_fdsolver_finalize.f
      USE ppm_module_data
      USE ppm_module_data_fieldsolver
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_alloc
=======

      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_data_fieldsolver

>>>>>>> 1.5

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
      CALL substart('ppm_fdsolver_finalize',t0,info)

      !-------------------------------------------------------------------------
      !  Finalize the FFT module
      !-------------------------------------------------------------------------
<<<<<<< ppm_fdsolver_finalize.f
      CALL ppm_util_fft_finalize(info)
      IF (info .NE. ppm_param_success) THEN
=======
  


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
          WRITE(mesg,'(A)') 'could not deallocate memory'
>>>>>>> 1.5
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_fdsolver_finalize',   &
     &        'Failed to finalize the FFT module.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fdsolver_finalize',t0,info)
      RETURN

      END SUBROUTINE ppm_fdsolver_finalize

