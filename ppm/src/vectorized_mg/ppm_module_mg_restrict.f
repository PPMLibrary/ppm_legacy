      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_mg_restrict
      !-------------------------------------------------------------------------
      !
      !
      ! Purpose       : multigrid module
      !               
      !
      ! Remarks       :
      !
      ! References    : 
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_mg_restrict.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:02  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/09/22 18:36:53  kotsalie
      !  MG new version
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------


#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __INTEGER          3
#define __LOGICAL          4
#define __2D               7
#define __3D               8
#define __SFIELD           9
#define __VFIELD          10

MODULE ppm_module_mg_restrict   
  !--------------------------------------------------------------------------
  !Modules
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! INCLUDE THE SOURCES
  !-----------------------------------------------------------------------------

CONTAINS

#define __DIM __SFIELD
#define __MESH_DIM __2D
#define __KIND __SINGLE_PRECISION
#include "vectorized_mg/ppm_mg_restrict.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "vectorized_mg/ppm_mg_restrict.f"
#undef __KIND
#undef __MESH_DIM

#define __MESH_DIM __3D
#define __KIND __SINGLE_PRECISION
#include "vectorized_mg/ppm_mg_restrict.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "vectorized_mg/ppm_mg_restrict.f"
#undef __KIND
#undef __MESH_DIM
#undef __DIM


#define __DIM __VFIELD
#define __MESH_DIM __2D
#define __KIND __SINGLE_PRECISION
#include "vectorized_mg/ppm_mg_restrict.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "vectorized_mg/ppm_mg_restrict.f"
#undef __KIND
#undef __MESH_DIM

#define __MESH_DIM __3D
#define __KIND __SINGLE_PRECISION
#include "vectorized_mg/ppm_mg_restrict.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "vectorized_mg/ppm_mg_restrict.f"
#undef __KIND
#undef __MESH_DIM
#undef __DIM



END MODULE ppm_module_mg_restrict


