      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_fft_map_3d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : maps the data from topology topo_ids(1)) / mesh 
      !                 (mesh_ids(1)) to topology (topo_ids(2))/mesh 
      !                 (mesh_ids(2)) during an FFT.
      !
      !  Input        : topo_ids(2)        (I) first: current topology   
      !                                        second: destination topology   
      !                 mesh_ids(2)        (I) first: current mesh of the data
      !                                        second: destination mesh   
      !
      !  Input/output : data_f(:,:,:,:,:)  (F)  data to be mapped
      !
      !  Output       : 
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_fft_map_3d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:02  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_map_3d_sca_s(data_fv, topo_ids, mesh_ids,&
     &  ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_map_3d_sca_d(data_fv, topo_ids, mesh_ids,&
     &  ghostsize, info)
#elif __KIND == __COMPLEX
      SUBROUTINE ppm_util_fft_map_3d_sca_c(data_fv, topo_ids, mesh_ids,&
     &  ghostsize, info)
#elif __KIND == __DOUBLE_COMPLEX
      SUBROUTINE ppm_util_fft_map_3d_sca_cc(data_fv, topo_ids, mesh_ids,&
     &  ghostsize, info)
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_map_3d_vec_s(data_fv, lda, topo_ids, mesh_ids, &
     &  ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_map_3d_vec_d(data_fv, lda, topo_ids, mesh_ids, &
     &  ghostsize, info)
#elif __KIND == __COMPLEX
      SUBROUTINE ppm_util_fft_map_3d_vec_c(data_fv, lda, topo_ids, mesh_ids, &
     &  ghostsize, info)
#elif __KIND == __DOUBLE_COMPLEX
      SUBROUTINE ppm_util_fft_map_3d_vec_cc(data_fv, lda, topo_ids, mesh_ids, &
     & ghostsize, info)
#endif
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_map_field
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid

      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:,:,:),          POINTER   :: data_fv
#elif __KIND == __COMPLEX | __KIND == __DOUBLE_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:),       POINTER   :: data_fv
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:,:,:,:),        POINTER   :: data_fv
#elif __KIND == __COMPLEX | __KIND == __DOUBLE_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:,:),     POINTER   :: data_fv
#endif
      INTEGER,                          INTENT(IN)     :: lda
#endif   
      INTEGER, DIMENSION(2),            INTENT(IN)     :: topo_ids
      INTEGER, DIMENSION(2),            INTENT(IN)     :: mesh_ids
      INTEGER, DIMENSION(3),            INTENT(IN)     :: ghostsize
      INTEGER              ,            INTENT(  OUT)  :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0
      INTEGER                                 :: k, i, j
      INTEGER                                 :: from_topo, to_topo
      INTEGER                                 :: from_mesh, to_mesh
#if   __DIM == __SFIELD
      INTEGER, PARAMETER                      :: lda = 1
#endif   
      INTEGER                                 :: maptype
      LOGICAL                                 :: valid
      CHARACTER(LEN=ppm_char)                 :: mesg

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_fft_map_3d',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (ppm_internal_topoid(topo_ids(1)) .NE. ppm_field_topoid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_map_3d',  &
     &            'Passed topology is not the current topology',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (topo_ids(2).GE.0) THEN
              CALL ppm_check_topoid(ppm_param_id_user,topo_ids(2),valid,info)
              IF (.NOT. valid) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_util_fft_map_3d',  &
     &                'Topology ID (from_topo) is invalid!',__LINE__,info)
                  GOTO 9999
              ENDIF           
          ENDIF
          IF (mesh_ids(1) .GT. 0) THEN
                CALL ppm_check_meshid(ppm_param_id_user,mesh_ids(1),     &
     &              ppm_internal_topoid(topo_ids(1)),valid,info)
                IF (.NOT. valid) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_util_fft_map_3d',  &
     &                  'Mesh ID (from_mesh) invalid!',__LINE__,info)
                    GOTO 9999
                ENDIF
          ENDIF
          IF (mesh_ids(2) .GT. 0) THEN
                CALL ppm_check_meshid(ppm_param_id_user,mesh_ids(2),     &
     &             ppm_internal_topoid(topo_ids(2)),valid,info)
                IF (.NOT. valid) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_util_fft_map_3d',  &
     &                  'Mesh ID (to_mesh) invalid!',__LINE__,info)
                    GOTO 9999
                ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Define source and destination 
      !-------------------------------------------------------------------------
      from_topo = topo_ids(1)
      to_topo   = topo_ids(2)
      from_mesh = mesh_ids(1)
      to_mesh   = mesh_ids(2)

      IF (ppm_debug .GT. 1) THEN
         WRITE(mesg,'(A,I5,A,I5)' )'Mapping from topo ',from_topo,   &
     &       ', mesh ',from_mesh 
         CALL ppm_write(ppm_rank,'ppm_util_fft_map',mesg,j)      

         WRITE(mesg,'(A,I5,A,I5)' )'          to topo ',to_topo,   &
     &       ', mesh ',to_mesh 
         CALL ppm_write(ppm_rank,'ppm_util_fft_map',mesg,j)      
      ENDIF

      !-------------------------------------------------------------------------
      !  Map fields
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
      maptype = ppm_param_map_global
      CALL ppm_map_field(DATA_fv,to_topo,from_mesh,to_mesh,ghostsize,   &
     &     maptype,info)
      maptype = ppm_param_map_push
      CALL ppm_map_field(DATA_fv,to_topo,from_mesh,to_mesh,ghostsize,   &
     &     maptype,info)
      maptype = ppm_param_map_send
      CALL ppm_map_field(DATA_fv,to_topo,from_mesh,to_mesh,ghostsize,   &
     &     maptype,info)
      maptype = ppm_param_map_pop
      CALL ppm_map_field(DATA_fv,to_topo,from_mesh,to_mesh,ghostsize,   &
     &     maptype,info)

#elif __DIM == __VFIELD
      maptype = ppm_param_map_global
      CALL ppm_map_field(DATA_fv,lda,to_topo,from_mesh,to_mesh,   &
     &     ghostsize,maptype,info)
      maptype = ppm_param_map_push
      CALL ppm_map_field(DATA_fv,lda,to_topo,from_mesh,to_mesh,   &
     &     ghostsize,maptype,info)
      maptype = ppm_param_map_send
      CALL ppm_map_field(DATA_fv,lda,to_topo,from_mesh,to_mesh,   &
     &     ghostsize,maptype,info)
      maptype = ppm_param_map_pop
      CALL ppm_map_field(DATA_fv,lda,to_topo,from_mesh,to_mesh,   &
     &     ghostsize,maptype,info)
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_fft_map_3d',t0,info)
      RETURN

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_map_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_map_3d_sca_d
#elif __KIND == __COMPLEX
      END SUBROUTINE ppm_util_fft_map_3d_sca_c
#elif __KIND == __DOUBLE_COMPLEX
      END SUBROUTINE ppm_util_fft_map_3d_sca_cc
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_map_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_map_3d_vec_d
#elif __KIND == __COMPLEX
      END SUBROUTINE ppm_util_fft_map_3d_vec_c
#elif __KIND == __DOUBLE_COMPLEX
      END SUBROUTINE ppm_util_fft_map_3d_vec_cc
#endif
#endif

