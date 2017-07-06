      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_util_fft_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Initialization for the 2D and 3D FFTs.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_fft_init.f,v $
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
#define __2D                         7
#define __3D                         8
#define __SFIELD                     9
#define __VFIELD                     10
#define __SLAB                       11
#define __PENCIL                     12

      MODULE ppm_module_util_fft_init

         !----------------------------------------------------------------------
         !  Define interface to ppm_util_fft_init
         !----------------------------------------------------------------------
         INTERFACE ppm_util_fft_init
            MODULE PROCEDURE ppm_util_fft_init_2d_sca_s
            MODULE PROCEDURE ppm_util_fft_init_2d_sca_d
            MODULE PROCEDURE ppm_util_fft_init_2d_vec_s
            MODULE PROCEDURE ppm_util_fft_init_2d_vec_d

            MODULE PROCEDURE ppm_util_fft_init_3d_sca_s
            MODULE PROCEDURE ppm_util_fft_init_3d_sca_d
            MODULE PROCEDURE ppm_util_fft_init_3d_vec_s
            MODULE PROCEDURE ppm_util_fft_init_3d_vec_d
         END INTERFACE
         INTERFACE ppm_util_fft_init_pen
            MODULE PROCEDURE ppm_util_fft_init_pen_2d_sca_s
            MODULE PROCEDURE ppm_util_fft_init_pen_2d_sca_d
            MODULE PROCEDURE ppm_util_fft_init_pen_2d_vec_s
            MODULE PROCEDURE ppm_util_fft_init_pen_2d_vec_d

            MODULE PROCEDURE ppm_util_fft_init_pen_3d_sca_s
            MODULE PROCEDURE ppm_util_fft_init_pen_3d_sca_d
            MODULE PROCEDURE ppm_util_fft_init_pen_3d_vec_s
            MODULE PROCEDURE ppm_util_fft_init_pen_3d_vec_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __CASE __SLAB

#define __KIND __SINGLE_PRECISION
#define __MESH_DIM __2D
#define __DIM __SFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#define __DIM __VFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#undef  __MESH_DIM

#define __MESH_DIM __3D
#define __DIM __SFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#define __DIM __VFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#undef __KIND
#undef  __MESH_DIM

#define __KIND __DOUBLE_PRECISION
#define __MESH_DIM __2D
#define __DIM __SFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#define __DIM __VFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#undef  __MESH_DIM

#define __MESH_DIM __3D
#define __DIM __SFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#define __DIM __VFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#undef __KIND
#undef  __MESH_DIM

#undef __CASE
#define __CASE __PENCIL

#define __KIND __SINGLE_PRECISION
#define __MESH_DIM __2D
#define __DIM __SFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#define __DIM __VFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#undef  __MESH_DIM

#define __MESH_DIM __3D
#define __DIM __SFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#define __DIM __VFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#undef __KIND
#undef  __MESH_DIM

#define __KIND __DOUBLE_PRECISION
#define __MESH_DIM __2D
#define __DIM __SFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#define __DIM __VFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#undef  __MESH_DIM

#define __MESH_DIM __3D
#define __DIM __SFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#define __DIM __VFIELD
#include "ppm_util_fft_init.f"
#undef  __DIM
#undef __KIND

#undef __CASE

      END MODULE ppm_module_util_fft_init
