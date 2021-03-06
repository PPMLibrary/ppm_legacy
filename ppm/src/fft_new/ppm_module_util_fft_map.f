      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_util_fft_map
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains the routines for mapping the 
      !                 data during an FFT.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_fft_map.f,v $
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
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __COMPLEX          3
#define __DOUBLE_COMPLEX   4

#define __SFIELD            9
#define __VFIELD           10


   MODULE ppm_module_util_fft_map

         !----------------------------------------------------------------------
         !  Define interface to ppm_util_fft_map
         !----------------------------------------------------------------------
         INTERFACE ppm_util_fft_map
            MODULE PROCEDURE ppm_util_fft_map_2d_sca_s
            MODULE PROCEDURE ppm_util_fft_map_2d_sca_d
            MODULE PROCEDURE ppm_util_fft_map_2d_sca_c
            MODULE PROCEDURE ppm_util_fft_map_2d_sca_cc
            MODULE PROCEDURE ppm_util_fft_map_2d_vec_s
            MODULE PROCEDURE ppm_util_fft_map_2d_vec_d
            MODULE PROCEDURE ppm_util_fft_map_2d_vec_c
            MODULE PROCEDURE ppm_util_fft_map_2d_vec_cc
            MODULE PROCEDURE ppm_util_fft_map_3d_sca_s
            MODULE PROCEDURE ppm_util_fft_map_3d_sca_d
            MODULE PROCEDURE ppm_util_fft_map_3d_sca_c
            MODULE PROCEDURE ppm_util_fft_map_3d_sca_cc
            MODULE PROCEDURE ppm_util_fft_map_3d_vec_s
            MODULE PROCEDURE ppm_util_fft_map_3d_vec_d
            MODULE PROCEDURE ppm_util_fft_map_3d_vec_c
            MODULE PROCEDURE ppm_util_fft_map_3d_vec_cc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_util_fft_map_2d.f"
#include "ppm_util_fft_map_3d.f"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_fft_map_2d.f"
#include "ppm_util_fft_map_3d.f"
#undef  __KIND

#define __KIND __COMPLEX
#include "ppm_util_fft_map_2d.f"
#include "ppm_util_fft_map_3d.f"
#undef  __KIND

#define __KIND __DOUBLE_COMPLEX
#include "ppm_util_fft_map_2d.f"
#include "ppm_util_fft_map_3d.f"
#undef  __KIND
#undef  __DIM


#define __DIM __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_util_fft_map_2d.f"
#include "ppm_util_fft_map_3d.f"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_fft_map_2d.f"
#include "ppm_util_fft_map_3d.f"
#undef  __KIND

#define __KIND __COMPLEX
#include "ppm_util_fft_map_2d.f"
#include "ppm_util_fft_map_3d.f"
#undef  __KIND

#define __KIND __DOUBLE_COMPLEX
#include "ppm_util_fft_map_2d.f"
#include "ppm_util_fft_map_3d.f"
#undef  __KIND
#undef  __DIM

      END MODULE ppm_module_util_fft_map
