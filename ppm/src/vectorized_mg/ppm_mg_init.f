      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_mg_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine initializes the solver for
      !                 2D and 3D problems
      !
      !  Input        :  equation   (I)  :  KIND OF EQUATION TO BE SOLVED
      !                                     FOR THE MOMENT ONLY POISSON
      !                  order      (I)  :  ORDER OF FINITE DIFFERENCES
      !                                     NOW SECOND. THE GHOSTSIZE IS
      !                                     AUTOMATICALLY ADJUSTED
      !                  smoother   (I)  :  NOW GAUSS-SEIDEL
      !
      !                  [lda]     (I)   : LEADING DIMENSION, ONLY TO BE
      !                                    GIVEN FOR VECTOR CASES
      !
      !                  ibcdef     (I)  : ARRAY OF BOUNDARY CONDITION
      !
      !
      !                  bcvalue   (F)   : ARRAY WHERE THE VALUES OF THE BC
      !                                    ARE STORED.IN CASE OF PERIODIC
      !                                    JUST GIVE ANY KIND OF VALUE
      !
      !                  EPSU      (F)   : STOPPING CRITERIUM. DETAIL:SHOULD
      !                                    BE SCALED WITH THE MAXIMUM VALUE           !                                    OF THE RHS.
      !
      !                  limlev    (I)    :Number of levels that the user
      !                                    wants to coarse.
      !
      !                  wcycle    (L)    : TRUE if the user wants W-cycle.
      !                                    OTHERWISE FALSE
      !                  lprint    (L)    : TRUE IF YOU WANT TO DUMP OUT
      !                                     INFORMATION
      !
      !
      !
      !  Input/output :
      !
      !  Output       : info       (I) return status. 0 upon success.
      !
      !  Remarks      :  PLEASE PAY ATTENTION THAT IN ORDER TO DIVIDE
      !                  FURTHER A MESH IT SHOULD BE DIVISIBLE WITH 2.
      !                  IF YOU WANT TO SOLVE DIFFERENT EQUATIONS
      !                  THE WHOLE MACHINERY SHOULD BE CALLED TWICE.
      !                  ALSO THE SOLVER IS NOW PROGRAMMED FOR THE POISSON
      !                  PROBLEM. A FUTURE IMPROVEMENT WOULD BE
      !                  TO USE A GENERAL STENCIL.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_mg_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:02  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2005/01/04 09:47:45  kotsalie
      !  ghostsize=2 for scalar case
      !
      !  Revision 1.4  2004/10/29 15:59:14  kotsalie
      !  RED BLACK SOR FOR 3d vec case. 2d will soon follow.
      !
      !  Revision 1.3  2004/09/28 14:04:49  kotsalie
      !  Changes concerning 4th order finite differences
      !
      !  Revision 1.2  2004/09/23 09:38:30  kotsalie
      !  added details in the header
      !
      !  Revision 1.1  2004/09/22 18:27:09  kotsalie
      !  MG new version
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if __DIM == __SFIELD
#if __MESH_DIM  == __2D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_init_2d_sca_s(equation,iorder,smoother,ibcdef,&
     &          bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_init_2d_sca_d(equation,iorder,smoother,ibcdef,&
     &                          bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#endif
#elif  __MESH_DIM  == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_init_3d_sca_s(equation,iorder,smoother,ibcdef,&
     &                         bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_init_3d_sca_d(equation,iorder,smoother,ibcdef,&
     &                         bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM  == __2D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_init_2d_vec_s(equation,iorder,smoother,lda,ibcdef,&
     &    bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_init_2d_vec_d(equation,iorder,smoother,lda,ibcdef,&
     &   bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#endif
#elif  __MESH_DIM  == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_init_3d_vec_s(equation,iorder,smoother,lda,ibcdef,&
     &              bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_init_3d_vec_d(equation,iorder,smoother,lda,ibcdef,&
     &                    bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#endif
#endif
#endif

        !----------------------------------------------------------------------
        !  Includes
        !----------------------------------------------------------------------
#include "ppm_define.h"

        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_data_mesh
        USE ppm_module_data_mg
        USE ppm_module_alloc
        USE ppm_module_mg_alloc
        USE ppm_module_error
        USE ppm_module_mesh_derive
        USE ppm_module_substart
        USE ppm_module_substop

        IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#else
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER, INTENT(IN)                                :: equation
        INTEGER, INTENT(IN)                                :: iorder
        INTEGER, INTENT(IN)                                :: smoother
#if __DIM == __VFIELD
        INTEGER,              INTENT(IN)                   ::  lda
#endif
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
        INTEGER,DIMENSION(:)                               ::  ibcdef
        REAL(MK),DIMENSION(:,:)                            ::  bcvalue
#elif __MESH_DIM == __3D
        INTEGER,DIMENSION(:)                               ::  ibcdef
        REAL(MK),DIMENSION(:,:,:)                          ::  bcvalue
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
        INTEGER,DIMENSION(:,:)                               ::  ibcdef
        REAL(MK),DIMENSION(:,:,:)                            ::  bcvalue
#elif __MESH_DIM == __3D
        INTEGER,DIMENSION(:,:)                               ::  ibcdef
        REAL(MK),DIMENSION(:,:,:,:)                          ::  bcvalue
#endif
#endif

        INTEGER,  INTENT(IN)                               :: mesh_id
        REAL(MK),INTENT(IN)                                :: EPSU
        INTEGER,INTENT(IN)                                 :: limlev
        LOGICAL,INTENT(IN)                                 :: wcycle
        LOGICAL,INTENT(IN)                                 :: lprint
        REAL(MK),INTENT(IN)                                :: omega
        INTEGER, INTENT(OUT)                               :: info
        !--------------------------------------------------------------------
        !  Local variables
        !-----------------------------------------------------------------------
        REAL(MK)                             :: t0
        INTEGER                              :: meshid,mlev
        INTEGER                              :: idom
        INTEGER                              ::  count,ilda,iface
        INTEGER                              :: i,j,k
        INTEGER                              :: kk
#if __MESH_DIM == __2D
        INTEGER                              :: dir
#endif
        INTEGER                              :: iter1,iter2,ix,iy
        INTEGER                              :: newmeshid,lmesh_id
        INTEGER , DIMENSION(1)               :: ldu1
        INTEGER , DIMENSION(2)               :: ldu2,ldl2
        INTEGER , DIMENSION(3)               :: ldu3,ldl3
#if __MESH_DIM == __3D
        INTEGER                              :: dir1,dir2,jj,iz
        INTEGER , DIMENSION(4)               :: ldu4,ldl4
#endif
        INTEGER , DIMENSION(ppm_dim)         :: Nml
        REAL(MK), DIMENSION(ppm_dim)         :: min_phys,max_phys
        INTEGER                              :: iopt,topoid


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
REAL(MK),DIMENSION(:,:),POINTER :: tuc
REAL(MK),DIMENSION(:,:),POINTER :: terr
#elif __MESH_DIM == __3D
REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
REAL(MK),DIMENSION(:,:,:),POINTER :: terr
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
REAL(MK),DIMENSION(:,:,:),POINTER :: terr
#elif __MESH_DIM == __3D
REAL(MK),DIMENSION(:,:,:,:),POINTER :: tuc
REAL(MK),DIMENSION(:,:,:,:),POINTER :: terr
#endif
#endif



        !-----------------------------------------------------------------------
        !  Externals
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !  Initialize
        !-----------------------------------------------------------------------

        CALL substart('ppm_mg_init',t0,info)


        !-----------------------------------------------------------------------
        !  Check arguments
        !-----------------------------------------------------------------------
          IF (ppm_debug.GT.0) THEN
#if __DIM == __VFIELD
          IF (lda.LE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_poiss_mg_init',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
          IF (EPSU.LE.0.0_MK) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_poiss_mg_init',  &
    &            'EPSU must be >0',__LINE__,info)
             GOTO 9999
          ENDIF
        ENDIF

        !---------------------------------------------------------------------
        ! Definition of necessary variables and allocation of arrays
        !---------------------------------------------------------------------
#if __DIM == __SFIELD
        vecdim = 1
#elif __DIM == __VFIELD
        vecdim = lda
#endif
        w_cycle=wcycle
        l_print=lprint

        topoid = ppm_field_topoid
        nsubs  = ppm_nsublist(topoid)
        PRINT *,'nsub:',nsubs
        meshid = ppm_meshid(topoid)%internal(mesh_id)
        lmesh_id = mesh_id


#if    __KIND == __SINGLE_PRECISION
        min_phys(:)=ppm_min_physs(:,topoid)
        max_phys(:)=ppm_max_physs(:,topoid)
        EPSU_s = EPSU
        omega_s=omega
#elif  __KIND == __DOUBLE_PRECISION
        min_phys(:)=ppm_min_physd(:,topoid)
        max_phys(:)=ppm_max_physd(:,topoid)
        EPSU_d = EPSU
        omega_d=omega
#endif
#if __MESH_DIM == __2D
        Nml(1) = ppm_cart_mesh(meshid,topoid)%Nm(1)
        Nml(2) = ppm_cart_mesh(meshid,topoid)%Nm(2)
        maxlev = INT(log10(Nml(1)*Nml(2)*REAL(ppm_nproc,MK))/log10(2.0_MK))
        IF (maxlev.GT.limlev) THEN
         maxlev=limlev
        ENDIF
#if __KIND == __SINGLE_PRECISION
        dx_s = (max_phys(1)-min_phys(1))/(Nml(1)-1)
        dy_s = (max_phys(2)-min_phys(2))/(Nml(2)-1)
        rdx2_s  = 1/(dx_s*dx_s)
        rdy2_s  = 1/(dy_s*dy_s)
#elif __KIND == __DOUBLE_PRECISION
        dx_d = (max_phys(1)-min_phys(1))/(Nml(1)-1)
        dy_d = (max_phys(2)-min_phys(2))/(Nml(2)-1)

        rdx2_d  = 1/(dx_d*dx_d)
        rdy2_d  = 1/(dy_d*dy_d)

#endif
#elif __MESH_DIM == __3D
        Nml(1) = ppm_cart_mesh(meshid,topoid)%Nm(1)
        Nml(2) = ppm_cart_mesh(meshid,topoid)%Nm(2)
        Nml(3) = ppm_cart_mesh(meshid,topoid)%Nm(3)
        maxlev = INT(log10(Nml(1)*Nml(2)*Nml(3)* &
     &           REAL(ppm_nproc,MK))/log10(2.0_MK))

        IF (maxlev.GT.limlev) THEN
         maxlev=limlev
        ENDIF
#if __KIND == __SINGLE_PRECISION
        dx_s = (max_phys(1)-min_phys(1))/(Nml(1)-1)
        dy_s = (max_phys(2)-min_phys(2))/(Nml(2)-1)
        dz_s = (max_phys(3)-min_phys(3))/(Nml(3)-1)
        rdx2_s = 1/(dx_s*dx_s)
        rdy2_s = 1/(dy_s*dy_s)
        rdz2_s = 1/(dz_s*dz_s)
#elif __KIND == __DOUBLE_PRECISION
        dx_d = (max_phys(1)-min_phys(1))/(Nml(1)-1)
        dy_d = (max_phys(2)-min_phys(2))/(Nml(2)-1)
        dz_d = (max_phys(3)-min_phys(3))/(Nml(3)-1)
        rdx2_d = 1/(dx_d*dx_d)
        rdy2_d = 1/(dy_d*dy_d)
        rdz2_d = 1/(dz_d*dz_d)
#endif
#endif


#if __DIM == __SFIELD
        iopt = ppm_param_alloc_fit
        ldu1(1) = 2*ppm_dim
        CALL ppm_alloc(bcdef_sca,ldu1,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
     &                   'Boundary condiotions',__LINE__,info)
           GOTO 9999
        ENDIF

        bcdef_sca(:)=ibcdef(:)
#elif __DIM == __VFIELD
        iopt = ppm_param_alloc_fit
        ldu2(1) = vecdim
        ldu2(2) = 2*ppm_dim
        CALL ppm_alloc(bcdef_vec,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
     &                   'Boundary condiotions',__LINE__,info)
           GOTO 9999
        ENDIF

        bcdef_vec(:,:)=ibcdef(:,:)
#endif


        iopt = ppm_param_alloc_fit
        ldu1(1) = ppm_dim
        CALL ppm_alloc(ghostsize,ldu1,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
     &                   'ghostsize',__LINE__,info)
           GOTO 9999
        ENDIF

        IF (iorder.EQ.ppm_param_order_2) THEN
         ghostsize(:)=1
         order=iorder
        ELSEIF (iorder.EQ.ppm_param_order_4) THEN
         ghostsize(:)=2
         order=iorder
        ENDIF

        iopt = ppm_param_alloc_fit
        ldu1(1) = ppm_dim
        CALL ppm_alloc(factor,ldu1,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
     &                   'factor',__LINE__,info)
           GOTO 9999
        ENDIF

        iopt = ppm_param_alloc_fit
        ldu1(1) = maxlev
        CALL ppm_alloc(meshid_g,ldu1,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
      &                  'meshid_g',__LINE__,info)
           GOTO 9999
        ENDIF

        iopt = ppm_param_alloc_fit
        ldu1(1) = maxlev
        CALL ppm_alloc(mesh_id_g,ldu1,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
      &                  'mesh_id_g',__LINE__,info)
           GOTO 9999
        ENDIF

        iopt = ppm_param_alloc_fit
        ldu3(1) = ppm_dim
        ldu3(2) = nsubs
        ldu3(3) = maxlev
        CALL ppm_alloc(start,ldu3,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
     &             'starting indices when updating the field',__LINE__,info)
           GOTO 9999
        ENDIF


        iopt = ppm_param_alloc_fit
        ldu3(1) = ppm_dim
        ldu3(2) = nsubs
        ldu3(3) = maxlev
        CALL ppm_alloc(stop,ldu3,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'stopping indices when updating the field',__LINE__,info)
           GOTO 9999
        ENDIF

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_2d_sca_s,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_2d_sca_s

#elif __KIND == __DOUBLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_2d_sca_d,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF
        mgfield => mgfield_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_3d_sca_s,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
                &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_3d_sca_s

#elif __KIND == __DOUBLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_3d_sca_d,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_3d_sca_d

#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_2d_vec_s,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_2d_vec_s

#elif __KIND == __DOUBLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_2d_vec_d,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF
        mgfield => mgfield_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_3d_vec_s,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
                &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_3d_vec_s

#elif __KIND == __DOUBLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_3d_vec_d,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_3d_vec_d

#endif
#endif
#endif




        iopt = ppm_param_alloc_fit
        ldu2(1) = 2*ppm_dim
        ldu2(2) = nsubs
        CALL ppm_alloc(lboundary,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the boundary alloc.',__LINE__,info)
           GOTO 9999
        ENDIF

        iopt = ppm_param_alloc_fit
        ldu2(1) = ppm_dim
        ldu2(2) = maxlev
        CALL ppm_alloc(max_node,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with a maximum number alloc.',__LINE__,info)
           GOTO 9999
        ENDIF

        max_node(:,:)=0

        lboundary(:,:)=.FALSE.
        start(:,:,:)=1


        !-----------------------------------------------------------------------
        ! Derive coarser meshes
        !-----------------------------------------------------------------------

        DO mlev=1,maxlev


#if __MESH_DIM == __2D

           !--------------------------------------------------------------------
           ! Go through the subs, define the stopping indices on each mesh,
           ! check and store if it is on the boundary, allocate the
           ! multigrid fields, pass the boundary values.
           !--------------------------------------------------------------------
           DO i=1,nsubs
              idom=ppm_isublist(i,topoid)

               stop(:,i,mlev)= ppm_cart_mesh(meshid,topoid)%nnodes(:,idom)

              DO j=1,ppm_dim
                 IF (max_node(j,mlev).LT.stop(j,i,mlev)) THEN
                    max_node(j,mlev)=stop(j,i,mlev)
                 ENDIF
              ENDDO



              !-----------------------------------------------------------------
              ! Allocate the function correction, the restricted errors,
              ! the residuals and the values on the boundary on each level.
              !----------------------------------------------------------------
#if __DIM == __SFIELD
              iopt = ppm_param_alloc_fit
              ldl2(1) = 1-ghostsize(1)
              ldl2(2) = 1-ghostsize(2)
              ldu2(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu2(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl2,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF


              tuc=>mgfield(i,mlev)%uc
              tuc(:,:)=0.0_MK

              iopt = ppm_param_alloc_fit
              ldu2(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)
              ldu2(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

               mgfield(i,mlev)%fc(:,:)=0.0_MK

              iopt = ppm_param_alloc_fit
              ldl2(1) = 1-ghostsize(1)
              ldl2(2) = 1-ghostsize(2)
              ldu2(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu2(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)

              CALL ppm_alloc(mgfield(i,mlev)%err,ldl2,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              terr=>mgfield(i,mlev)%err
              terr(:,:)=0.0_MK
#elif __DIM == __VFIELD
              iopt = ppm_param_alloc_fit
              ldl3(1) = 1
              ldl3(2) = 1-ghostsize(1)
              ldl3(3) = 1-ghostsize(2)
              ldu3(1) = vecdim
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF


              tuc=>mgfield(i,mlev)%uc
              tuc(:,:,:)=0.0_MK

              iopt = ppm_param_alloc_fit
              ldu3(1) = vecdim
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

               mgfield(i,mlev)%fc(:,:,:)=0.0_MK

              iopt = ppm_param_alloc_fit
              ldl3(1) = 1
              ldl3(2) = 1-ghostsize(1)
              ldl3(3) = 1-ghostsize(2)
              ldu3(1) = vecdim
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)

              CALL ppm_alloc(mgfield(i,mlev)%err,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              terr=>mgfield(i,mlev)%err
              terr(:,:,:)=0.0_MK
#endif
              iopt = ppm_param_alloc_fit
              ldl2(1) = 1-ghostsize(1)
              ldl2(2) = 1-ghostsize(2)
              ldu2(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu2(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%mask_red,ldl2,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the mask  alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              iopt = ppm_param_alloc_fit
              ldl2(1) = 1-ghostsize(1)
              ldl2(2) = 1-ghostsize(2)
              ldu2(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu2(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%mask_black,ldl2,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with mask alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF


              !----------------------------------------------------------------
              !Filling the mask for communication (red black)
              !----------------------------------------------------------------
              DO iy=1-ghostsize(2),ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
                 DO ix=1-ghostsize(1),ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)

                    IF (MOD(ix+iy,2).EQ.0) THEN

                       mgfield(i,mlev)%mask_red(ix,iy)=.TRUE.
                       mgfield(i,mlev)%mask_black(ix,iy)=.FALSE.

                    ELSE

                       mgfield(i,mlev)%mask_red(ix,iy)   = .FALSE.
                       mgfield(i,mlev)%mask_black(ix,iy) = .TRUE.

                    ENDIF
                 ENDDO
              ENDDO

           ENDDO!DO 1,nsubs



#elif __MESH_DIM == __3D


           DO i=1,nsubs

              idom=ppm_isublist(i,topoid)
              stop(:,i,mlev) = ppm_cart_mesh(meshid,topoid)%nnodes(:,idom)

              DO j=1,ppm_dim
                 IF (max_node(j,mlev).LT.stop(j,i,mlev)) THEN
                    max_node(j,mlev)=stop(j,i,mlev)
                 ENDIF
              ENDDO

              IF (ppm_subs_bc(1,i,topoid).EQ.1) THEN

                 lboundary(1,i)=.TRUE.


              ELSEIF (ppm_subs_bc(3,i,topoid).EQ.1) THEN

                 lboundary(3,i)=.TRUE.

              ELSEIF (ppm_subs_bc(2,i,topoid).EQ.1) THEN

                 lboundary(2,i)=.TRUE.


              ELSEIF (ppm_subs_bc(4,i,topoid).EQ.1) THEN

                 lboundary(4,i)=.TRUE.

              ELSEIF (ppm_subs_bc(5,i,topoid).EQ.1) THEN

                 lboundary(5,i)=.TRUE.


              ELSEIF (ppm_subs_bc(6,i,topoid).EQ.1) THEN

                 lboundary(6,i)=.TRUE.


              ENDIF


              !----------------------------------------------------------------
              ! Allocate the function correction, the restricted errors and the
              !residuals on each level.
              !----------------------------------------------------------------

#if __DIM == __SFIELD
              iopt = ppm_param_alloc_fit
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1-ghostsize(3)
              ldu3(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              tuc=>mgfield(i,mlev)%uc
              tuc(:,:,:)=0.0_MK

              iopt = ppm_param_alloc_fit
              ldu3(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              mgfield(i,mlev)%fc(:,:,:)=0.0_MK



              iopt = ppm_param_alloc_fit
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1-ghostsize(3)
              ldu3(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%err,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              terr=>mgfield(i,mlev)%err
              terr(:,:,:)=0.0_MK

#elif __DIM == __VFIELD
              iopt = ppm_param_alloc_fit
              ldl4(1) = 1
              ldl4(2) = 1-ghostsize(1)
              ldl4(3) = 1-ghostsize(2)
              ldl4(4) = 1-ghostsize(3)
              ldu4(1) = vecdim
              ldu4(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu4(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu4(4) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl4,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              tuc=>mgfield(i,mlev)%uc
              tuc(:,:,:,:)=0.0_MK


              iopt = ppm_param_alloc_fit
              ldu4(1) = vecdim
              ldu4(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)
              ldu4(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)
              ldu4(4) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              mgfield(i,mlev)%fc(:,:,:,:)=0.0_MK


              iopt = ppm_param_alloc_fit
              ldl4(1) = 1
              ldl4(2) = 1-ghostsize(1)
              ldl4(3) = 1-ghostsize(2)
              ldl4(4) = 1-ghostsize(3)
              ldu4(1) = vecdim
              ldu4(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu4(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu4(4) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%err,ldl4,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              terr=>mgfield(i,mlev)%err
              terr(:,:,:,:)=0.0_MK

#endif


              iopt = ppm_param_alloc_fit
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1-ghostsize(3)
              ldu3(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%mask_red,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the mask  alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              iopt = ppm_param_alloc_fit
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1-ghostsize(3)
              ldu3(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%mask_black,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
                      &        'Problem with mask alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF


              !----------------------------------------------------------------
              !Filling the mask for communication (red black)
              !-----------------------------------------------------------------
              DO iz=1-ghostsize(3),&
               ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)

                 DO iy=1-ghostsize(2),&
                   & ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
                    DO ix=1-ghostsize(1),&
                   &  ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)

                       IF (MOD(ix+iy+iz,2).EQ.0) THEN

                          mgfield(i,mlev)%mask_red(ix,iy,iz)=.TRUE.
                          mgfield(i,mlev)%mask_black(ix,iy,iz)=.FALSE.
                          !mgfield(i,mlev)%mask_black(ix,iy,iz)=.TRUE.

                       ELSE

                          mgfield(i,mlev)%mask_red(ix,iy,iz)   = .FALSE.
                          mgfield(i,mlev)%mask_black(ix,iy,iz) = .TRUE.

                      ENDIF
                    ENDDO
                 ENDDO
              ENDDO



           ENDDO!DO i=1,nsubs

#endif


           factor(:)=2
           mesh_id_g(mlev)=lmesh_id
           meshid_g(mlev)=meshid
           newmeshid=-1


           IF (mlev.LT.maxlev) THEN
            CALL ppm_mesh_derive(topoid,meshid,ppm_param_mesh_coarsen,factor,&
     &                          newmeshid,info)


            lmesh_id = newmeshid
            meshid = ppm_meshid(topoid)%internal(lmesh_id)

           ENDIF
        ENDDO!DO mlev=1,maxlev


        !----------------------------------------------------------------------
        !  Return
        !----------------------------------------------------------------------
9999    CONTINUE
        CALL substop('ppm_mg_init',t0,info)
        RETURN
#if    __DIM       == __SFIELD
#if    __MESH_DIM  == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_init_2d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_init_2d_sca_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_init_3d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_init_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if    __MESH_DIM  == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_init_2d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_init_2d_vec_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_init_3d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_init_3d_vec_d
#endif
#endif
#endif
