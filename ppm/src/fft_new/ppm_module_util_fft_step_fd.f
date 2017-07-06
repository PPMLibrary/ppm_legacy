      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_util_fft_step_fd
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the fft 
      !                 routines that only do one direction (step) using
      !                 the init/finalize.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_fft_step_fd.f,v $
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
#define __SINGLE_PRECISION_COMPLEX_Z 7
#define __DOUBLE_PRECISION_COMPLEX_Z 8
#define __SLAB                      10

      MODULE ppm_module_util_fft_step_fd

         !----------------------------------------------------------------------
         !  Define interface 
         !----------------------------------------------------------------------
!        INTERFACE ppm_util_fft_slab_fd
!           MODULE PROCEDURE ppm_util_fft_slab_fd_3ds
!           MODULE PROCEDURE ppm_util_fft_slab_fd_3dd
!        END INTERFACE

!        INTERFACE ppm_util_fft_z_fd
!           MODULE PROCEDURE ppm_util_fft_z_fd_3dc
!           MODULE PROCEDURE ppm_util_fft_z_fd_3dcc
!        END INTERFACE

!        INTERFACE ppm_util_fft_step_fd
!           MODULE PROCEDURE ppm_util_fft_step_fd_2ds
!           MODULE PROCEDURE ppm_util_fft_step_fd_2dd
!           MODULE PROCEDURE ppm_util_fft_step_fd_2dc
!           MODULE PROCEDURE ppm_util_fft_step_fd_2dcc

!           MODULE PROCEDURE ppm_util_fft_step_fd_3ds
!           MODULE PROCEDURE ppm_util_fft_step_fd_3dd
!           MODULE PROCEDURE ppm_util_fft_step_fd_3dc
!           MODULE PROCEDURE ppm_util_fft_step_fd_3dcc
!        END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
!        CONTAINS

#define __KIND __SINGLE_PRECISION
!include "ppm_util_fft_step_fd_2d.f"
!include "ppm_util_fft_step_fd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
!include "ppm_util_fft_step_fd_2d.f"
!include "ppm_util_fft_step_fd_3d.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
!include "ppm_util_fft_step_fd_2d.f"
!include "ppm_util_fft_step_fd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
!include "ppm_util_fft_step_fd_2d.f"
!include "ppm_util_fft_step_fd_3d.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX_Z
!include "ppm_util_fft_step_fd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX_Z
!include "ppm_util_fft_step_fd_3d.f"
#undef __KIND

#define __CASE __SLAB
#define __KIND __SINGLE_PRECISION
!include "ppm_util_fft_step_fd_3d.f"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
!include "ppm_util_fft_step_fd_3d.f"
#undef  __KIND
#undef  __CASE

      END MODULE ppm_module_util_fft_step_fd
