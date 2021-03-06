      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_tree_done
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the
      !                 routine that determines when a tree is done.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_tree_done.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/09/22 10:32:07  ivos
      !  Initial implementation.
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

      MODULE ppm_module_tree_done

         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_tree_done
            MODULE PROCEDURE ppm_tree_done_s
            MODULE PROCEDURE ppm_tree_done_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_tree_done.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_tree_done.f"
#undef __KIND

      END MODULE ppm_module_tree_done
