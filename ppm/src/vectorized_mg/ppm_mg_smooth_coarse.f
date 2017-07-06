!-----------------------------------------------------------------------
!  Subroutine   :            ppm_mg_smooth_coarse
!-----------------------------------------------------------------------
!  Purpose      : In this routine we compute the corrections for
!                 the function based on the Gauss-Seidel iteration
!
!
!  Input        : nsweep      (I) number of iterations(sweeps)
!  Input/output :
!
!  Output       : info        (I) return status. 0 upon success
!
!  Remarks      :
!
!  References   :
!
!  Revisions    :
!-------------------------------------------------------------------------
!  $Log: ppm_mg_smooth_coarse.f,v $
!  Revision 1.1.1.1  2007/07/13 10:19:02  ivos
!  CBL version of the PPM library
!
!  Revision 1.7  2005/01/04 09:45:29  kotsalie
!  ghostsize=2
!
!  Revision 1.6  2004/11/05 18:09:49  kotsalie
!  FINAL FEATURE BEFORE TEST.I DO NOT USE MASKS
!
!  Revision 1.4  2004/10/29 15:59:31  kotsalie
!  RED BLACK SOR
!
!  Revision 1.3  2004/09/28 14:05:31  kotsalie
!  Changes concernig 4th order finite differences
!
!  Revision 1.2  2004/09/23 12:16:49  kotsalie
!  Added USE statement
!
!  Revision 1.1  2004/09/22 18:42:39  kotsalie
!  MG new version
!
!
!------------------------------------------------------------------------
!  Parallel Particle Mesh Library (PPM)
!  Institute of Computational Science
!  ETH Zentrum, Hirschengraben 84
!  CH-8092 Zurich, Switzerland
!-------------------------------------------------------------------------

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_2D_sca_s(nsweep,mlev,c1,c2,c3,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_2D_sca_d(nsweep,mlev,c1,c2,c3,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_3D_sca_s(nsweep,mlev,c1,c2,c3,c4,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_3D_sca_d(nsweep,mlev,c1,c2,c3,c4,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_2D_vec_s(nsweep,mlev,c1,c2,c3,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_2D_vec_d(nsweep,mlev,c1,c2,c3,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_3D_vec_s(nsweep,mlev,c1,c2,c3,c4,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_3D_vec_d(nsweep,mlev,c1,c2,c3,c4,info)
#endif
#endif
#endif

        !----------------------------------------------------------------------
        !  Includes
        !----------------------------------------------------------------------
#include "ppm_define.h"

        !-------------------------------------------------------------------
        !  Modules
        !--------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_data_mg
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_map_field_ghost
        USE ppm_module_data_mesh
        USE ppm_module_write



        IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#else
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
        !-------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------
        INTEGER,                   INTENT(IN)      ::  nsweep
        INTEGER,                   INTENT(IN)      ::  mlev
#if  __MESH_DIM == __2D
        REAL(MK),                  INTENT(IN)      ::  c1,c2,c3
#elif __MESH_DIM == __3D
        REAL(MK),                  INTENT(IN)      ::  c1,c2,c3,c4
#endif
        INTEGER,                   INTENT(INOUT)   ::  info
        !---------------------------------------------------------------------
        !  Local variables
        !---------------------------------------------------------------------
        CHARACTER(LEN=256) :: cbuf
        INTEGER                                    ::  i,j,isub,color
        REAL(MK)                                   ::  c11,c22,c33,c44
        INTEGER                                    ::  ilda,isweep,count
        INTEGER                                    ::  k,idom
        REAL(MK)                                   ::  x,y
        REAL(MK)                                   ::  omega
        INTEGER,DIMENSION(1)                       ::  ldu1
#if __MESH_DIM == __2D
        INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
        INTEGER,DIMENSION(3)                       ::  ldl3,ldu3
#endif
#if __MESH_DIM == __3D
        INTEGER,DIMENSION(5)                       ::  ldl5,ldu5
        INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
#endif
        INTEGER                                    ::  iopt,iface,topoid
        REAL(MK)                                   ::  t0
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_2d_sca_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_2d_sca_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_3d_sca_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_3d_sca_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_2d_vec_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_2d_vec_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_3d_vec_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_3d_vec_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#endif
#endif

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
        REAL(MK),DIMENSION(:,:,:),POINTER :: uc_dummy
#elif __MESH_DIM == __3D
        REAL(MK),DIMENSION(:,:,:,:),POINTER :: uc_dummy
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
        REAL(MK),DIMENSION(:,:,:,:),POINTER :: uc_dummy
#elif __MESH_DIM == __3D
        REAL(MK),DIMENSION(:,:,:,:,:),POINTER :: uc_dummy
#endif
#endif


#if __DIM == __SFIELD
#if __MESH_DIM == __2D
     REAL(MK) :: moldu
#elif __MESH_DIM == __3D
     REAL(MK) :: moldu
#endif
#elif  __DIM == __VFIELD
#if __MESH_DIM == __2D
     REAL(MK),DIMENSION(:),POINTER :: moldu
#elif __MESH_DIM == __3D
     REAL(MK),DIMENSION(:),POINTER :: moldu
#endif
#endif


#if __DIM == __SFIELD
#if __MESH_DIM == __2D
        REAL(MK),DIMENSION(:,:),POINTER :: tuc
#elif __MESH_DIM == __3D
       REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
      REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
#elif __MESH_DIM == __3D
      REAL(MK),DIMENSION(:,:,:,:),POINTER :: tuc
#endif
#endif

#if __MESH_DIM == __2D
        LOGICAL,DIMENSION(:,:),POINTER :: mask_red
        LOGICAL,DIMENSION(:,:),POINTER :: mask_black
#elif __MESH_DIM == __3D
       LOGICAL,DIMENSION(:,:,:),POINTER :: mask_red
       LOGICAL,DIMENSION(:,:,:),POINTER :: mask_black
#endif


#if __KIND == __SINGLE_PRECISION
      omega=omega_s
#elif __KIND == __DOUBLE_PRECISION
      omega=omega_d
#endif

        !-----------------------------------------------------------------------
        !Externals
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !Initialize
        !-----------------------------------------------------------------------

        CALL substart('ppm_mg_smooth_coarse',t0,info)
        IF (l_print) THEN
         WRITE (cbuf,*) 'SMOOTHER entering ','mlev:',mlev
         CALL PPM_WRITE(ppm_rank,'mg_smooth',cbuf,info)
        ENDIF

        !-----------------------------------------------------------------------
        !  Check arguments
        !-----------------------------------------------------------------------
        IF (ppm_debug .GT. 0) THEN
          IF (nsweep.LT.1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
     &            'nsweep must be >=1',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (mlev.LE.1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
     &            'level must be >1',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c1.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
     &            'Factor c1 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c2.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
     &            'Factor c2 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c3.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
     &            'Factor c3 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#if __MESH_DIM == __3D
          IF (c4.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
     &            'Factor c4 must be >0',__LINE__,info)
              GOTO 9999
         ENDIF
#endif
        ENDIF

#if __DIM ==__VFIELD
        NULLIFY(moldu)
#endif
        NULLIFY(uc_dummy)

        !-----------------------------------------------------------------------
        !Definition of necessary variables and allocation of arrays
        !-----------------------------------------------------------------------
        topoid=ppm_field_topoid


#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        mgfield=>mgfield_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
        mgfield=>mgfield_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        mgfield=>mgfield_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
        mgfield=>mgfield_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        mgfield=>mgfield_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
        mgfield=>mgfield_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        mgfield=>mgfield_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
        mgfield=>mgfield_3d_vec_d
#endif
#endif
#endif


#if  __DIM == __SFIELD
#if  __MESH_DIM == __2D

        !-----------------------------------------------------------------------
        !Implementation
        !-----------------------------------------------------------------------

            iopt = ppm_param_alloc_fit
            ldl3(1) = 1-ghostsize(1)
            ldl3(2) = 1-ghostsize(2)
            ldl3(3) = 1
            ldu3(1) = max_node(1,mlev)+ghostsize(1)
            ldu3(2) = max_node(2,mlev)+ghostsize(2)
            ldu3(3) = nsubs
            CALL ppm_alloc(uc_dummy,ldl3,ldu3,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'uc_dummy',__LINE__,info)
            GOTO 9999
            ENDIF


        count = 0
            iopt = ppm_param_alloc_fit
            ldl3(1) = 1-ghostsize(1)
            ldl3(2) = 1-ghostsize(2)
            ldl3(3) = 1
            ldu3(1) = max_node(1,mlev)+ghostsize(1)
            ldu3(2) = max_node(2,mlev)+ghostsize(2)
            ldu3(3) = nsubs
            CALL ppm_alloc(mask_dummy_2d,ldl3,ldu3,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'mask_dummy_2d',__LINE__,info)
            GOTO 9999
            ENDIF
        DO isweep=1,nsweep
           DO color=0,1

              DO isub=1,nsubs

                 IF (color.EQ.0) THEN
                    mask_red=>mgfield(isub,mlev)%mask_red
                    mask_dummy_2d(:,:,&
     &                            isub)=mask_red(:,:)
                 ELSE
                    mask_black=>mgfield(isub,mlev)%mask_black
                    mask_dummy_2d(:,:,&
     &                             isub)=mask_black(:,:)
                 ENDIF
                 tuc=>mgfield(isub,mlev)%uc
                 uc_dummy(:,:,isub)=tuc(:,:)


              ENDDO!DO isub

              !-----------------------------------------------------------------
              !Communicate red(even) if color==0 or communicate black(odd)
              !if color==1
              !-----------------------------------------------------------------


              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                    ghostsize,ppm_param_map_ghost_get,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_push,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_send,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                          ghostsize,ppm_param_map_pop,info,mask_dummy_2d)



              DO isub=1,nsubs
                 tuc=>mgfield(isub,mlev)%uc
                          tuc(:,:)=uc_dummy(&
     &                         :,:,isub)
                 DO j=start(2,isub,mlev),stop(2,isub,mlev)
                    DO i=start(1,isub,mlev)+mod(j+color,2),stop(1,isub,mlev),2
                          mgfield(isub,mlev)%uc(i,j) = c1*(&
     &                                   (mgfield(isub,mlev)%uc(i-1,j)+ &
     &                                mgfield(isub,mlev)%uc(i+1,j))*c2 + &
     &                                 (mgfield(isub,mlev)%uc(i,j-1)+&
     &                                  mgfield(isub,mlev)%uc(i,j+1))*c3-&
     &                                         mgfield(isub,mlev)%fc(i,j))

                    ENDDO
                 ENDDO
              ENDDO
              IF (isweep.EQ.nsweep) THEN
               IF (color.EQ.1) THEN

                 DO isub=1,nsubs
                    mask_red=>mgfield(isub,mlev)%mask_red
                    mask_dummy_2d(:,:,&
     &                            isub)=mask_red(:,:)

                 tuc=>mgfield(isub,mlev)%uc
                 uc_dummy(:,:,isub)=tuc(:,:)
                ENDDO
               ENDIF
              ENDIF

             ENDDO!DO color

             IF (isweep.EQ.nsweep) THEN

              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                    ghostsize,ppm_param_map_ghost_get,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_push,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_send,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                          ghostsize,ppm_param_map_pop,info,mask_dummy_2d)


              DO isub=1,nsubs
                 tuc=>mgfield(isub,mlev)%uc
                          tuc(:,:)=uc_dummy(&
     &                         :,:,isub)
              ENDDO
            ENDIF


           ENDDO!DO nsweep



            iopt = ppm_param_dealloc
            ldl3(1) = 1-ghostsize(1)
            ldl3(2) = 1-ghostsize(2)
            ldl3(3) = 1
            ldu3(1) = max_node(1,mlev)+ghostsize(1)
            ldu3(2) = max_node(2,mlev)+ghostsize(2)
            ldu3(3) = nsubs
            CALL ppm_alloc(uc_dummy,ldl3,ldu3,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'uc_dummy',__LINE__,info)
            GOTO 9999
            ENDIF

#elif __MESH_DIM == __3D

        !-----------------------------------------------------------------------
        !Implementation
        !-----------------------------------------------------------------------

            iopt = ppm_param_alloc_fit
            ldl4(1) = 1-ghostsize(1)
            ldl4(2) = 1-ghostsize(2)
            ldl4(3) = 1-ghostsize(3)
            ldl4(4) = 1
            ldu4(1) = max_node(1,mlev)+ghostsize(1)
            ldu4(2) = max_node(2,mlev)+ghostsize(2)
            ldu4(3) = max_node(3,mlev)+ghostsize(3)
            ldu4(4) = nsubs
            CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'uc_dummy',__LINE__,info)
            GOTO 9999
            ENDIF





        DO isweep=1,nsweep
           DO color=0,1


              DO isub=1,nsubs

                 tuc=>mgfield(isub,mlev)%uc
                   DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                    DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                     DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                       uc_dummy(i,j,k,isub)=tuc(i,j,k)
                     ENDDO
                    ENDDO
                   ENDDO

              ENDDO!DO isub


              !-----------------------------------------------------------------
              !Communicate red(even) if color==0 or communicate black(odd)
              !if color==1
              !-----------------------------------------------------------------


              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                    ghostsize,ppm_param_map_ghost_get,info)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_push,info)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_send,info)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                          ghostsize,ppm_param_map_pop,info)



              DO isub=1,nsubs
                 tuc=>mgfield(isub,mlev)%uc
                   DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                    DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                     DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                         tuc(i,j,k)=uc_dummy(i,j,k,isub)
                     ENDDO
                    ENDDO
                  ENDDO


                 DO k=start(3,isub,mlev),stop(3,isub,mlev)
                    DO j=start(2,isub,mlev),stop(2,isub,mlev)
                       DO i=start(1,isub,mlev)+mod(j+k+color,2),stop(1,isub,mlev),2

                            moldu=tuc(i,j,k)

                             mgfield(isub,mlev)%uc(i,j,k) = moldu+&
     &                             omega*(&
     &                             c1*((mgfield(isub,mlev)%uc(i-1,j,k)+ &
     &                            mgfield(isub,mlev)%uc(i+1,j,k))*c2 + &
     &                                 (mgfield(isub,mlev)%uc(i,j-1,k)+&
     &                            mgfield(isub,mlev)%uc(i,j+1,k))*c3 + &
     &                           (mgfield(isub,mlev)%uc(i,j,k-1)+&
     &                            mgfield(isub,mlev)%uc(i,j,k+1))*c4 - &
     &                                    mgfield(isub,mlev)%fc(i,j,k))&
     &                            -moldu)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO!subs

                  IF (isweep.EQ.nsweep) THEN

                    IF (color.EQ.1) THEN
                     DO isub=1,nsubs

                      tuc=>mgfield(isub,mlev)%uc

                   DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                    DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                     DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                       uc_dummy(i,j,k,isub)=tuc(i,j,k)
                     ENDDO
                    ENDDO
                   ENDDO

                    ENDDO
                    ENDIF
                  ENDIF

          ENDDO!DO color

              IF (isweep.EQ.nsweep) THEN


              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                    ghostsize,ppm_param_map_ghost_get,info)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_push,info)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_send,info)
              CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
     &                          ghostsize,ppm_param_map_pop,info)



              ENDIF

              DO isub=1,nsubs
                 tuc=>mgfield(isub,mlev)%uc

                   DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                    DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                     DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                         tuc(i,j,k)=uc_dummy(i,j,k,isub)
                     ENDDO
                    ENDDO
                  ENDDO

              ENDDO
        ENDDO!Do isweep

            iopt = ppm_param_dealloc
            ldl4(1) = 1-ghostsize(1)
            ldl4(2) = 1-ghostsize(2)
            ldl4(3) = 1-ghostsize(3)
            ldl4(4) = 1
            ldu4(1) = max_node(1,mlev)+ghostsize(1)
            ldu4(2) = max_node(2,mlev)+ghostsize(2)
            ldu4(3) = max_node(3,mlev)+ghostsize(3)
            ldu4(4) = nsubs
            CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'uc_dummy',__LINE__,info)
            GOTO 9999
            ENDIF
#endif
#elif __DIM == __VFIELD
#if  __MESH_DIM == __2D

        !-----------------------------------------------------------------------
        !Implementation
        !-----------------------------------------------------------------------

            iopt = ppm_param_alloc_fit
            ldl4(1) = 1
            ldl4(2) = 1-ghostsize(1)
            ldl4(3) = 1-ghostsize(2)
            ldl4(4) = 1
            ldu4(1) = vecdim
            ldu4(2) = max_node(1,mlev)+ghostsize(1)
            ldu4(3) = max_node(2,mlev)+ghostsize(2)
            ldu4(4) = nsubs
            CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'uc_dummy',__LINE__,info)
            GOTO 9999
            ENDIF


            count = 0

            iopt = ppm_param_alloc_fit
            ldl3(1) = 1-ghostsize(1)
            ldl3(2) = 1-ghostsize(2)
            ldl3(3) = 1
            ldu3(1) = max_node(1,mlev)+ghostsize(1)
            ldu3(2) = max_node(2,mlev)+ghostsize(2)
            ldu3(3) = nsubs
            CALL ppm_alloc(mask_dummy_2d,ldl3,ldu3,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'mask_dummy_2d',__LINE__,info)
            GOTO 9999
            ENDIF
        DO isweep=1,nsweep
           DO color=0,1

              DO isub=1,nsubs

                 IF (color.EQ.0) THEN
                    mask_red=>mask_red
                    mask_dummy_2d(:,:,&
     &                            isub)=mgfield(isub,mlev)%mask_red(:,:)
                 ELSE
                    mask_black=>mgfield(isub,mlev)%mask_black
                    mask_dummy_2d(:,:,&
     &                             isub)=mask_black(:,:)
                 ENDIF
                 tuc=>mgfield(isub,mlev)%uc
                 uc_dummy(:,:,:,isub)=tuc(:,:,:)

              ENDDO!DO isub

              !-----------------------------------------------------------------
              !Communicate red(even) if color==0 or communicate black(odd)
              !if color==1
              !-----------------------------------------------------------------

              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                    ghostsize,ppm_param_map_ghost_get,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_push,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_send,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                          ghostsize,ppm_param_map_pop,info,mask_dummy_2d)



              DO isub=1,nsubs
                 tuc=>mgfield(isub,mlev)%uc
                 tuc(:,:,:)=uc_dummy(&
     &                         :,:,:,isub)
                 DO j=start(2,isub,mlev),stop(2,isub,mlev)
                    DO i=start(1,isub,mlev)+mod(j+color,2),stop(1,isub,mlev),2
                     DO ilda=1,vecdim
                          mgfield(isub,mlev)%uc(ilda,i,j) = c1*(&
     &                                   (mgfield(isub,mlev)%uc(ilda,i-1,j)+ &
     &                                mgfield(isub,mlev)%uc(ilda,i+1,j))*c2 + &
     &                                 (mgfield(isub,mlev)%uc(ilda,i,j-1)+&
     &                                  mgfield(isub,mlev)%uc(ilda,i,j+1))*c3-&
     &                                         mgfield(isub,mlev)%fc(ilda,i,j))

                     ENDDO
                    ENDDO
                 ENDDO
              ENDDO
                   IF (isweep.EQ.nsweep) THEN
                    IF (color.EQ.1) THEN

                     DO isub=1,nsubs
                      mask_red=>mask_red
                      mask_dummy_2d(:,:,&
     &                            isub)=mgfield(isub,mlev)%mask_red(:,:)

                      tuc=>mgfield(isub,mlev)%uc
                      uc_dummy(:,:,:,isub)=tuc(:,:,:)
                     ENDDO
                    ENDIF
                   ENDIF




           ENDDO!DO color

             IF (isweep.EQ.nsweep) THEN
              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                    ghostsize,ppm_param_map_ghost_get,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_push,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_send,info,mask_dummy_2d)
              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                          ghostsize,ppm_param_map_pop,info,mask_dummy_2d)


              DO isub=1,nsubs
                 tuc=>mgfield(isub,mlev)%uc
                 tuc(:,:,:)=uc_dummy(&
     &                         :,:,:,isub)
              ENDDO
             ENDIF

        ENDDO!DO nsweep



            iopt = ppm_param_dealloc
            ldl4(1) = 1
            ldl4(2) = 1-ghostsize(1)
            ldl4(3) = 1-ghostsize(2)
            ldl4(4) = 1
            ldu4(1) = vecdim
            ldu4(2) = max_node(1,mlev)+ghostsize(1)
            ldu4(3) = max_node(2,mlev)+ghostsize(2)
            ldu4(4) = nsubs
            CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'uc_dummy',__LINE__,info)
            GOTO 9999
            ENDIF
#elif __MESH_DIM == __3D

        !-----------------------------------------------------------------------
        !Implementation
        !-----------------------------------------------------------------------

            iopt = ppm_param_alloc_fit
            ldl5(1) = 1
            ldl5(2) = 1-ghostsize(1)
            ldl5(3) = 1-ghostsize(2)
            ldl5(4) = 1-ghostsize(3)
            ldl5(5) = 1
            ldu5(1) = vecdim
            ldu5(2) = max_node(1,mlev)+ghostsize(1)
            ldu5(3) = max_node(2,mlev)+ghostsize(2)
            ldu5(4) = max_node(3,mlev)+ghostsize(3)
            ldu5(5) = nsubs
            CALL ppm_alloc(uc_dummy,ldl5,ldu5,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'uc_dummy',__LINE__,info)
            GOTO 9999
            ENDIF





           iopt = ppm_param_alloc_fit
           ldu1(1)=vecdim
           CALL ppm_alloc(moldu,ldu1,iopt,info)
           IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'moldu',__LINE__,info)
            GOTO 9999
           ENDIF



        DO isweep=1,nsweep

           DO color=0,1

            DO isub=1,nsubs
                 !--------------------------------------------------------------
                 !Impose boundaries on even if color=0 or odd if color=1
                 !--------------------------------------------------------------

                 tuc=>mgfield(isub,mlev)%uc

                   DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                    DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                     DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                       uc_dummy(1,i,j,k,isub)=tuc(1,i,j,k)
                       uc_dummy(2,i,j,k,isub)=tuc(2,i,j,k)
                       uc_dummy(3,i,j,k,isub)=tuc(3,i,j,k)
                     ENDDO
                    ENDDO
                   ENDDO
              ENDDO!DO isub




              !-----------------------------------------------------------------
              !Communicate red(even) if color==0 or communicate black(odd)
              !if color==1
              !-----------------------------------------------------------------


              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                    ghostsize,ppm_param_map_ghost_get,info)
              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_push,info)
              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                         ghostsize,ppm_param_map_send,info)
              CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
     &                          ghostsize,ppm_param_map_pop,info)


              DO isub=1,nsubs
                 tuc=>mgfield(isub,mlev)%uc

                   DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                    DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                     DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                         tuc(1,i,j,k)=uc_dummy(1,i,j,k,isub)
                         tuc(2,i,j,k)=uc_dummy(2,i,j,k,isub)
                         tuc(3,i,j,k)=uc_dummy(3,i,j,k,isub)
                     ENDDO
                    ENDDO
                  ENDDO



                 DO k=start(3,isub,mlev),stop(3,isub,mlev)
                    DO j=start(2,isub,mlev),stop(2,isub,mlev)
                       DO i=start(1,isub,mlev)+mod(j+k+color,2),stop(1,isub,mlev),2

                        moldu(1) = tuc(1,i,j,k)
                        moldu(2) = tuc(2,i,j,k)
                        moldu(3) = tuc(3,i,j,k)



                             mgfield(isub,mlev)%uc(1,i,j,k) = moldu(1)+&
     &                             omega*(&
     &                             c1*((mgfield(isub,mlev)%uc(1,i-1,j,k)+ &
     &                            mgfield(isub,mlev)%uc(1,i+1,j,k))*c2 + &
     &                                 (mgfield(isub,mlev)%uc(1,i,j-1,k)+&
     &                            mgfield(isub,mlev)%uc(1,i,j+1,k))*c3 + &
     &                           (mgfield(isub,mlev)%uc(1,i,j,k-1)+&
     &                            mgfield(isub,mlev)%uc(1,i,j,k+1))*c4 - &
     &                            mgfield(isub,mlev)%fc(1,i,j,k))&
     &                            -moldu(1))



                             mgfield(isub,mlev)%uc(2,i,j,k) = moldu(2)+&
     &                             omega*(&
     &                             c1*((mgfield(isub,mlev)%uc(2,i-1,j,k)+ &
     &                            mgfield(isub,mlev)%uc(2,i+1,j,k))*c2 + &
     &                                 (mgfield(isub,mlev)%uc(2,i,j-1,k)+&
     &                            mgfield(isub,mlev)%uc(2,i,j+1,k))*c3 + &
     &                           (mgfield(isub,mlev)%uc(2,i,j,k-1)+&
     &                            mgfield(isub,mlev)%uc(2,i,j,k+1))*c4 - &
     &                            mgfield(isub,mlev)%fc(2,i,j,k))&
     &                            -moldu(2))

                             mgfield(isub,mlev)%uc(3,i,j,k) = moldu(3)+&
     &                             omega*(&
     &                             c1*((mgfield(isub,mlev)%uc(3,i-1,j,k)+ &
     &                            mgfield(isub,mlev)%uc(3,i+1,j,k))*c2 + &
     &                                 (mgfield(isub,mlev)%uc(3,i,j-1,k)+&
     &                            mgfield(isub,mlev)%uc(3,i,j+1,k))*c3 + &
     &                           (mgfield(isub,mlev)%uc(3,i,j,k-1)+&
     &                            mgfield(isub,mlev)%uc(3,i,j,k+1))*c4 - &
     &                            mgfield(isub,mlev)%fc(3,i,j,k))&
     &                            -moldu(3))

                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO!subs


                IF (isweep.EQ.nsweep) THEN
                   IF (color.EQ.1) THEN
                    DO isub=1,nsubs


                      tuc=>mgfield(isub,mlev)%uc

                   DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                    DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                     DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                       uc_dummy(1,i,j,k,isub)=tuc(1,i,j,k)
                       uc_dummy(2,i,j,k,isub)=tuc(2,i,j,k)
                       uc_dummy(3,i,j,k,isub)=tuc(3,i,j,k)
                     ENDDO
                    ENDDO
                   ENDDO
                    ENDDO!subs

                   ENDIF
                  ENDIF


          ENDDO!DO color


         IF (isweep.EQ.nsweep) THEN

          CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
       &             ghostsize,ppm_param_map_ghost_get,info)
          CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
       &                   ghostsize,ppm_param_map_push,info)
          CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
       &                   ghostsize,ppm_param_map_send,info)
          CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
       &                ghostsize,ppm_param_map_pop,info)



                  DO isub=1,nsubs
                   tuc=>mgfield(isub,mlev)%uc

                   DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                    DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                     DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                         tuc(1,i,j,k)=uc_dummy(1,i,j,k,isub)
                         tuc(2,i,j,k)=uc_dummy(2,i,j,k,isub)
                         tuc(3,i,j,k)=uc_dummy(3,i,j,k,isub)
                     ENDDO
                    ENDDO
                  ENDDO
                  ENDDO
         ENDIF

        ENDDO!Do isweep

            iopt = ppm_param_dealloc
            ldl5(1) = 1
            ldl5(2) = 1-ghostsize(1)
            ldl5(3) = 1-ghostsize(2)
            ldl5(4) = 1-ghostsize(3)
            ldl5(5) = 1
            ldu5(1) = vecdim
            ldu5(2) = max_node(1,mlev)+ghostsize(1)
            ldu5(4) = max_node(2,mlev)+ghostsize(2)
            ldu5(4) = max_node(3,mlev)+ghostsize(3)
            ldu5(5) = nsubs
            CALL ppm_alloc(uc_dummy,ldl5,ldu5,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'uc_dummy',__LINE__,info)
            GOTO 9999
            ENDIF

           iopt = ppm_param_dealloc
           ldu1(1)=vecdim
           CALL ppm_alloc(moldu,ldu1,iopt,info)
           IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'GSsolv',    &
      &                       'moldu',__LINE__,info)
            GOTO 9999
           ENDIF

#endif
#endif


        !----------------------------------------------------------------------
        !  Return
        !----------------------------------------------------------------------
9999    CONTINUE
        CALL substop('ppm_mg_smooth_coarse',t0,info)
        RETURN
#if __DIM == __SFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_coarse_2D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_coarse_2D_sca_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_coarse_3D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_coarse_3D_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_coarse_2D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_coarse_2D_vec_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_coarse_3D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_coarse_3D_vec_d
#endif
#endif
#endif




