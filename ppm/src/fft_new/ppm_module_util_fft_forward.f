      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_util_fft_forward
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 2D and 3D parallel FFTs. This does the complete
      !                 FFT from time domain to frequency domain.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_fft_forward.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION           1
#define __DOUBLE_PRECISION           2
#define __SINGLE_PRECISION_COMPLEX   5
#define __DOUBLE_PRECISION_COMPLEX   6
#define __SFIELD                     7
#define __VFIELD                     8

      MODULE ppm_module_util_fft_forward

         !----------------------------------------------------------------------
         !  Define interface to ppm_util_fft_forward
         !----------------------------------------------------------------------
         INTERFACE ppm_util_fft_forward
            MODULE PROCEDURE ppm_util_fft_forward_2d_ss
            MODULE PROCEDURE ppm_util_fft_forward_2d_sd
            MODULE PROCEDURE ppm_util_fft_forward_2d_vs
            MODULE PROCEDURE ppm_util_fft_forward_2d_vd

            MODULE PROCEDURE ppm_util_fft_forward_3d_ss
            MODULE PROCEDURE ppm_util_fft_forward_3d_sd
            MODULE PROCEDURE ppm_util_fft_forward_3d_vs
            MODULE PROCEDURE ppm_util_fft_forward_3d_vd
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#define __DIM __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_util_fft_forward_2d.f"
#include "ppm_util_fft_forward_3d.f"
#undef   __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_util_fft_forward_2d.f"
#include "ppm_util_fft_forward_3d.f"
#undef __KIND
#undef  __DIM

#define __DIM __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_util_fft_forward_2d.f"
#include "ppm_util_fft_forward_3d.f"
#undef   __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_util_fft_forward_2d.f"
#include "ppm_util_fft_forward_3d.f"
#undef __KIND
#undef  __DIM

      END MODULE ppm_module_util_fft_forward

