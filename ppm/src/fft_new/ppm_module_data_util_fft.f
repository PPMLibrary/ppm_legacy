      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_data_util_fft
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all variables that are
      !                 global to the util_fft routines.
      !
      !  Remarks      : The data structures are initialized in
      !                 ppm_util_fft_init, used in ppm_util_fft_* and
      !                 destroyed in ppm_util_fft_finalize.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_util_fft.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_util_fft

          !---------------------------------------------------------------------
          !  Includes
          !---------------------------------------------------------------------
#include "ppm_define.h"

          !---------------------------------------------------------------------
          !  Modules
          !---------------------------------------------------------------------
          USE ppm_module_data,ONLY:ppm_kind_single,ppm_kind_double
          PRIVATE :: ppm_kind_single,ppm_kind_double

#ifdef __FFTW
          !---------------------------------------------------------------------
          !  Plans for FFTW
          !---------------------------------------------------------------------
          INTEGER*8 Plan_fd_s,     Plan_fd_d
          INTEGER*8 Plan_slab_fd_s,Plan_slab_fd_d
          INTEGER*8 Plan_fd_c_y,   Plan_fd_cc_y
          INTEGER*8 Plan_fd_c_z,   Plan_fd_cc_z
          INTEGER*8 Plan_bd_s,     Plan_bd_d
          INTEGER*8 Plan_slab_bd_s,Plan_slab_bd_d
          INTEGER*8 Plan_bd_c_y,   Plan_bd_cc_y
          INTEGER*8 Plan_bd_c_z,   Plan_bd_cc_z
#endif

#ifdef __MATHKEISAN
          !---------------------------------------------------------------------
          !  Work memory for Mathkeisan FFT
          !---------------------------------------------------------------------
          REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_fd_s    => NULL()
          REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_bd_s    => NULL()
          REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_fd_d    => NULL()
          REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_bd_d    => NULL()
          REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_fd_c_y  => NULL()
          REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_bd_c_y  => NULL()
          REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_fd_cc_y => NULL()
          REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_bd_cc_y => NULL()
          REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_fd_c_z  => NULL()
          REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_bd_c_z  => NULL()
          REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_fd_cc_z => NULL()
          REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_bd_cc_z => NULL()
#endif
          !---------------------------------------------------------------------
          !  Table sizes
          !---------------------------------------------------------------------
          INTEGER, DIMENSION(1)         :: lda_table,lda_table_y,lda_table_z

          !---------------------------------------------------------------------
          !  IDs of temporary topologies and meshes
          !---------------------------------------------------------------------
          INTEGER, DIMENSION(3)         :: topo_ids_tmp
          INTEGER, DIMENSION(4)         :: mesh_ids_tmp

          !---------------------------------------------------------------------
          !  Flags if the current topology already is a pencil or slab type
          !---------------------------------------------------------------------
          LOGICAL                       :: Its_xpencil_topo,Its_xyslab_topo
      END MODULE ppm_module_data_util_fft
