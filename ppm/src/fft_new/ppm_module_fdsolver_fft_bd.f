#ifdef __XLF
@PROCESS NOHOT
#endif
      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_fdsolver_fft_bd
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the backward
      !                 stand-alone fft routines of the fieldsolver.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_fdsolver_fft_bd.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2005/02/16 12:07:35  hiebers
      !  initial implementation
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

      MODULE ppm_module_fdsolver_fft_bd

         !----------------------------------------------------------------------
         !  Define interface to ppm_fdsolver_fft_bd
         !----------------------------------------------------------------------
         INTERFACE ppm_fdsolver_fft_bd
            MODULE PROCEDURE ppm_fdsolver_fft_bd_2ds
            MODULE PROCEDURE ppm_fdsolver_fft_bd_2dd
            MODULE PROCEDURE ppm_fdsolver_fft_bd_2dc
            MODULE PROCEDURE ppm_fdsolver_fft_bd_2dcc

            MODULE PROCEDURE ppm_fdsolver_fft_bd_3ds
            MODULE PROCEDURE ppm_fdsolver_fft_bd_3dd
            MODULE PROCEDURE ppm_fdsolver_fft_bd_3dc
            MODULE PROCEDURE ppm_fdsolver_fft_bd_3dcc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_fft_bd_2d.f"
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_fft_bd_2d.f"
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_fdsolver_fft_bd_2d.f"
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_fdsolver_fft_bd_2d.f"
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

      END MODULE ppm_module_fdsolver_fft_bd
