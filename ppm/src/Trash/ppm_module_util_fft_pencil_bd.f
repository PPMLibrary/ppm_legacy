      !-------------------------------------------------------------------------
      !  Module       :          ppm_module_util_fft_pencil_bd
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the fft 
      !                  routines that only do one space direction (pencil),
      !                  but the full, stand-alone, thing:
      !                  init/compute/finalize. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_fft_pencil_bd.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:02  ivos
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
#define __SINGLE_PRECISION_COMPLEX_Z 7
#define __DOUBLE_PRECISION_COMPLEX_Z 8

#define __SLAB                      10

      MODULE ppm_module_util_fft_pencil_bd

         !----------------------------------------------------------------------
         !  Define interface 
         !----------------------------------------------------------------------
         INTERFACE ppm_util_fft_bd_full
            MODULE PROCEDURE ppm_util_fft_pencil_bd_2ds
            MODULE PROCEDURE ppm_util_fft_pencil_bd_2dd
            MODULE PROCEDURE ppm_util_fft_pencil_bd_2dc
            MODULE PROCEDURE ppm_util_fft_pencil_bd_2dcc

            MODULE PROCEDURE ppm_util_fft_pencil_bd_3ds
            MODULE PROCEDURE ppm_util_fft_pencil_bd_3dd
            MODULE PROCEDURE ppm_util_fft_pencil_bd_3dc
            MODULE PROCEDURE ppm_util_fft_pencil_bd_3dcc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_util_fft_pencil_bd_2d.f"
#include "ppm_util_fft_pencil_bd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_fft_pencil_bd_2d.f"
#include "ppm_util_fft_pencil_bd_3d.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_util_fft_pencil_bd_2d.f"
#include "ppm_util_fft_pencil_bd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_util_fft_pencil_bd_2d.f"
#include "ppm_util_fft_pencil_bd_3d.f"
#undef __KIND

      END MODULE ppm_module_util_fft_pencil_bd
