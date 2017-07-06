      !-------------------------------------------------------------------------
      !  Module       :           ppm_module_util_fft_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Finalizes the FFT module and frees all memory.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_fft_finalize.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_util_fft_finalize

         !----------------------------------------------------------------------
         !  Define interface to the routine
         !----------------------------------------------------------------------
         INTERFACE ppm_util_fft_finalize
            MODULE PROCEDURE ppm_util_fft_finalize
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "ppm_util_fft_finalize.f"

      END MODULE ppm_module_util_fft_finalize
