      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_find_neigh
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the
      !                 decomposition routines. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_find_neigh_old.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:02  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 08:56:11  ivos
      !  Renamed.
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

      MODULE ppm_module_find_neigh_old

         !----------------------------------------------------------------------
         !  Define interface to neighbor finding routine
         !----------------------------------------------------------------------
         INTERFACE ppm_find_neigh_old
            MODULE PROCEDURE ppm_find_neigh_old_s
            MODULE PROCEDURE ppm_find_neigh_old_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_find_neigh_old.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_find_neigh_old.f"
#undef __KIND

      END MODULE ppm_module_find_neigh_old
