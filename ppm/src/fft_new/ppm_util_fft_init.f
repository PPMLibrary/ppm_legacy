      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_fft_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Initializes the FFT module by setting up all the 
      !                 temporary topologies that are needed, creating the 
      !                 FFT plans, and allocating all memory. This must be 
      !                 called before ppm_util_fft_forward/backward.
      !
      !  Input        :  DATA_fv(:,:,:,:) (F) field data
      !                  lda_fv      (I) size of leading dimension in the 
      !                                  vector case
      !                  mesh_id_user(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology(xpencils)
      !                                  topo_ids(2) second  topology(ypencils)
      !                                  topo_ids(3) third   topology(zpencils)
      !                                         (3D only!!)
      !                                 
      !                  mesh_ids(3) (I) mesh IDs on which the FFTs are to be
      !                                  performed:
      !                                  mesh_ids(1) first mesh 
      !                                         (xpencils,real)
      !                                  mesh_ids(2) second mesh 
      !                                         (xpencils,complex)
      !                                  mesh_ids(3) third mesh 
      !                                         (ypencils,complex)
      !                                  mesh_ids(4) forth mesh 
      !                                         (zpencils,complex) (3D only!!)
      !                  ghostsize(3) (I)ghostsize
      !
      !  Input/output : 
      !
      !  Output       : info          (I) Return status. 0 on success.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_fft_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if   __CASE == __SLAB

#if   __DIM == __SFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_init_2d_sca_s(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_init_2d_sca_d(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_init_3d_sca_s(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_init_3d_sca_d(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#endif
#endif
#endif

#if   __DIM == __VFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_init_2d_vec_s(DATA_fv,lda_fv,mesh_id_user, &
     &   topo_ids, mesh_ids,ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_init_2d_vec_d(DATA_fv,lda_fv,mesh_id_user, &
     &  topo_ids,mesh_ids,ghostsize,info)
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_init_3d_vec_s(DATA_fv,lda_fv,mesh_id_user, &
     &  topo_ids,mesh_ids,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_init_3d_vec_d(DATA_fv,lda_fv,mesh_id_user, &
     & topo_ids,mesh_ids,ghostsize, info)
#endif
#endif
#endif

#else

#if   __DIM == __SFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_init_pen_2d_sca_s(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_init_pen_2d_sca_d(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_init_pen_3d_sca_s(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_init_pen_3d_sca_d(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#endif
#endif
#endif

#if   __DIM == __VFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_init_pen_2d_vec_s(DATA_fv,lda_fv,mesh_id_user, &
     &   topo_ids, mesh_ids,ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_init_pen_2d_vec_d(DATA_fv,lda_fv,mesh_id_user, &
     &  topo_ids,mesh_ids,ghostsize,info)
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_fft_init_pen_3d_vec_s(DATA_fv,lda_fv,mesh_id_user, &
     &  topo_ids,mesh_ids,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_fft_init_pen_3d_vec_d(DATA_fv,lda_fv,mesh_id_user, &
     & topo_ids,mesh_ids,ghostsize, info)
#endif
#endif
#endif
#endif

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_data_util_fft
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_mktopo
      USE ppm_module_error

      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
#ifdef  __FFTW
      INCLUDE "fftw3.f"
#endif

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
#if   __MESH_DIM == __2D
      REAL(MK), DIMENSION(:,:,:),          POINTER   :: data_fv
#elif __MESH_DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:),        POINTER   :: data_fv
#endif
#endif

#if   __DIM == __VFIELD
#if   __MESH_DIM == __2D
      REAL(MK), DIMENSION(:,:,:,:),        POINTER   :: data_fv
      INTEGER,                             INTENT(IN):: lda_fv
#elif __MESH_DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:,:),      POINTER   :: data_fv
      INTEGER,                             INTENT(IN):: lda_fv
#endif
#endif
      INTEGER                      , INTENT(IN)      :: mesh_id_user

#if   __MESH_DIM == __2D
      INTEGER, DIMENSION(2)        , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(3)        , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(2)        , INTENT(IN   )   :: ghostsize
#elif __MESH_DIM == __3D
      INTEGER, DIMENSION(3)        , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(4)        , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(3)        , INTENT(IN   )   :: ghostsize
#endif
      INTEGER                      , INTENT(  OUT)   :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)                   :: t0
      INTEGER                                 :: k,i,j
      CHARACTER(LEN=ppm_char)                 :: mesg
#if   __MESH_DIM == __2D
      REAL(MK), DIMENSION(:,:),     POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:),  POINTER   :: data_comp,data_compl
      INTEGER, DIMENSION(2)                   :: lda
#elif __MESH_DIM == __3D
      REAL(MK), DIMENSION(:,:,:),   POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:,:),POINTER   :: data_comp,data_compl
      INTEGER, DIMENSION(3)                   :: lda
#endif
      INTEGER,DIMENSION(1)                    :: MB_in
      INTEGER                                 :: Nx_in,Ny_in,Nz_in
      INTEGER                                 :: Nx_out,Ny_out,Nz_out
#ifdef __FFTW
      INTEGER                                 :: mbistride,mbrank,mbidist
      INTEGER, DIMENSION(1)                   :: mbiembed,mboembed
      INTEGER                                 :: mbhowmany,mbodist
      INTEGER, DIMENSION(2)                   :: iembed_slab,oembed_slab
#endif
#ifdef __MATHKEISAN
      INTEGER                                 :: isign_fft,scale_fft,isys
      INTEGER                                 :: incx,incy
      INTEGER                                 :: isys
      REAL(MK), DIMENSION(:),POINTER          :: work
#endif
      REAL(MK), DIMENSION(:,:), POINTER       :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp,asign
      REAL(MK), DIMENSION(3  )                :: min_phys,max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id,mesh_id
      INTEGER                                 :: yhmax,zhmax
      INTEGER                                 :: mesh_id_internal
      INTEGER                                 :: mesh_id_xpen,mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      INTEGER                                 :: mesh_id_slab
      INTEGER                                 :: mesh_id_slab_complex
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart,istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen,istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata,ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen,ndata_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen,ndata_slab
      INTEGER , DIMENSION(:  ), POINTER       :: sub2proc
      INTEGER , DIMENSION(:  ), POINTER       :: isublist
      INTEGER                                 :: nsublist,idom
      INTEGER                                 :: dim,n
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_user,topo_id_internal
      INTEGER                                 :: topo_id_xpen,topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER                                 :: topo_id_slab
      INTEGER, DIMENSION(3)                   :: Nm,Nm_com,Nm_poisson
      INTEGER, DIMENSION(2)                   :: Nm_slab

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_fft_init',t0,info)

#if  !(defined(__FFTW) | defined(__MATHKEISAN))
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
#ifndef __FFTW
      CALL ppm_error(ppm_err_nofftw,'ppm_util_fft_init',  &
     &    'PPM was compiled without fftw support',__LINE__,info)
#endif
#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_util_fft_init',  &
     &    'PPM was compiled without MATHKEISAN support',__LINE__,info)
#endif
      GOTO 9999   
#else

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
           IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_util_fft_init',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ppm_field_topoid .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_fft_init',  &
     &            'No field topology has been defined so far',__LINE__,info)
              GOTO 9999
          ENDIF         
      ENDIF

      !----------------------------------------------------------------------
      !  Allocate the isublist array
      !----------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      lda(1) = ppm_nsublist(ppm_field_topoid)
      CALL ppm_alloc(isublist,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          WRITE (mesg,'(A,I10,A)') 'allocating ',      &
     &                ppm_nsublist(ppm_field_topoid) ,' isublist failed'
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',      &
     &        mesg,__LINE__,info)
          GOTO 9999
      ENDIF
      lda(1) = ppm_dim
      lda(2) = ppm_nsubs(ppm_field_topoid)
      CALL ppm_alloc(ndata,lda,iopt,info)
      CALL ppm_alloc(ndata_slab,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          WRITE (mesg,'(A,I10,A)') 'allocating ',      &
     &                ppm_nsublist(ppm_field_topoid) ,' ndata failed'
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',      &
     &        mesg,__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize variables
      !-------------------------------------------------------------------------
      asign           = ppm_param_assign_internal
      topo_id_internal= ppm_field_topoid
      topo_id_user    = ppm_user_topoid(topo_id_internal)
      mesh_id_internal= ppm_meshid(topo_id_internal)%internal(mesh_id_user)
      nsublist        = ppm_nsublist(topo_id_internal)
      isublist        = ppm_isublist(:,topo_id_internal)
      Nm(1)           = ppm_cart_mesh(mesh_id_internal, topo_id_internal)%Nm(1)
      Nm(2)           = ppm_cart_mesh(mesh_id_internal, topo_id_internal)%Nm(2)
      Nm(3)           = ppm_cart_mesh(mesh_id_internal, topo_id_internal)%Nm(3)
      Nm_com(1)       = Nm(1)/2+1
      Nm_com(2)       = Nm(2)
      Nm_com(3)       = Nm(3)
      ndata           = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes
      bcdef(1)        = ppm_param_bcdef_periodic 
      bcdef(2)        = ppm_param_bcdef_periodic 
      bcdef(3)        = ppm_param_bcdef_periodic 
      bcdef(4)        = ppm_param_bcdef_periodic 
      bcdef(5)        = ppm_param_bcdef_periodic 
      bcdef(6)        = ppm_param_bcdef_periodic 
      Npart           = 0
      ! size of dummy arrays
      lda             = 1

      !-------------------------------------------------------------------------
      !  Allocate memory
      !-------------------------------------------------------------------------
#ifdef __MATHKEISAN
      CALL ppm_alloc(work,1,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',     &
     &        'Mathkeisan work memory WORK',__LINE__,info)
          GOTO 9999
      ENDIF
#endif
      CALL ppm_alloc(data_real,lda,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',     &
     &        'data_real',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(data_comp,lda,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',     &
     &        'data_comp',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(data_compl,lda,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',     &
     &        'data_compl',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the min_phys and max_phys
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      min_phys(1) = ppm_min_physs(1,topo_id_internal) 
      min_phys(2) = ppm_min_physs(2,topo_id_internal) 
      min_phys(3) = ppm_min_physs(3,topo_id_internal) 
      max_phys(1) = ppm_max_physs(1,topo_id_internal)
      max_phys(2) = ppm_max_physs(2,topo_id_internal)
      max_phys(3) = ppm_max_physs(3,topo_id_internal)
#elif __KIND == __DOUBLE_PRECISION
      min_phys(1) = ppm_min_physd(1,topo_id_internal) 
      min_phys(2) = ppm_min_physd(2,topo_id_internal) 
      min_phys(3) = ppm_min_physd(3,topo_id_internal) 
      max_phys(1) = ppm_max_physd(1,topo_id_internal)
      max_phys(2) = ppm_max_physd(2,topo_id_internal)
      max_phys(3) = ppm_max_physd(3,topo_id_internal)
#endif
      length_phys(1) = max_phys(1) - min_phys(1)
      length_phys(2) = max_phys(2) - min_phys(2)
      length_phys(3) = max_phys(3) - min_phys(3)

      IF (ppm_debug .GT. 1) THEN
         WRITE(mesg,'(A,3F15.3)' ) 'minimal extent', min_phys(1), min_phys(2), &
     &        min_phys(3)
         CALL ppm_write(ppm_rank,'ppm_util_fft_init',mesg,j)
         WRITE(mesg,'(A,3F15.3)' ) 'maximal extent', max_phys(1), max_phys(2), &
     &        max_phys(3)
         CALL ppm_write(ppm_rank,'ppm_util_fft_init',mesg,j)
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if the input topology already is an x-pencil or xy-slab
      !-------------------------------------------------------------------------
      Its_xpencil_topo = .TRUE.
      Its_xyslab_topo  = .TRUE.

      DO k=1,ppm_nsublist(topo_id_internal)
          idom = ppm_isublist(k,topo_id_internal)

#if   __KIND == __SINGLE_PRECISION
         length(1)=  ppm_max_subs(1,idom,topo_id_internal) - &
     &                                 ppm_min_subs(1,idom,topo_id_internal) 
         length(2)=  ppm_max_subs(2,idom,topo_id_internal) - &
     &                                 ppm_min_subs(2,idom,topo_id_internal) 
         IF( abs(length(1) - length_phys(1)).GT.(ppm_myepss) ) THEN 
            Its_xpencil_topo=.FALSE.
            Its_xyslab_topo =.FALSE.
         ENDIF
         IF( abs(length(2) - length_phys(2)).GT.(ppm_myepss) ) THEN 
            Its_xyslab_topo =.FALSE.
         ENDIF
#elif __KIND == __DOUBLE_PRECISION
         length(1)=ppm_max_subd(1,idom,topo_id_internal) - &
     &                                 ppm_min_subd(1,idom,topo_id_internal) 
         length(2)=ppm_max_subd(2,idom,topo_id_internal) - &
     &                                 ppm_min_subd(2,idom,topo_id_internal) 
         IF( abs(length(1) - length_phys(1)).GT.(ppm_myepsd) ) THEN 
            Its_xpencil_topo=.FALSE.
            Its_xyslab_topo =.FALSE.
         ENDIF
         IF( abs(length(2) - length_phys(2)).GT.(ppm_myepsd) ) THEN 
            Its_xyslab_topo =.FALSE.
         ENDIF
#endif
      ENDDO

      IF (ppm_debug .GT. 0) THEN
         IF(Its_xyslab_topo) THEN
            WRITE(mesg,'(A)' ) 'XY slab topology' 
            CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         
         ELSE
            WRITE(mesg,'(A)' ) 'Not XY slab topology' 
            CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         ENDIF
         IF(Its_xpencil_topo) THEN
            WRITE(mesg,'(A)' ) 'X pencil topology' 
            CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         
         ELSE
            WRITE(mesg,'(A)' ) 'Not X pencil topology' 
            CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Define the temporary topologies for the transforms
      !-------------------------------------------------------------------------
#if __CASE == __SLAB
      !-------------------------------------------------------------------------
      !  XY-SLABS
      !-------------------------------------------------------------------------
      IF(Its_xyslab_topo) THEN
         topo_id_xpen         = topo_id_user
         mesh_id_xpen         = mesh_id_user
         topo_id_zpen         = topo_ids(2)
         mesh_id_xpen_complex = mesh_ids(2)
         mesh_id_zpen         = mesh_ids(3)
      ELSE
         topo_id_xpen         = topo_ids(1)
         topo_id_zpen         = topo_ids(2)
         mesh_id_xpen         = mesh_ids(1)
         mesh_id_xpen_complex = mesh_ids(2)
         mesh_id_zpen         = mesh_ids(3)
      ENDIF

      IF (ppm_debug .GT. 1) THEN
         WRITE(mesg,'(A)' ) '  ID             topo  mesh' 
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A)' ) '-----------------------------------'
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'Original        ',topo_id_user, mesh_id_user  
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'XZ Slab        ', topo_id_xpen, mesh_id_xpen
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'XZ Slab  Complex', topo_id_xpen, &
     &       mesh_id_xpen_complex
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'Z Pencil Complex', topo_id_zpen, mesh_id_zpen
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
      ENDIF

      IF (.NOT.Its_xyslab_topo) THEN
          decomp = ppm_param_decomp_xy_slab
          CALL ppm_mktopo(xp,Npart,Nm,decomp,asign, min_phys,max_phys,      &
     &                    bcdef,ghostsize,topo_id_xpen,mesh_id_xpen,min_sub,&
     &                    max_sub,cost,sub2proc,nsubs,isublist,             &
     &                    nsublist,istart,ndata,info)

         topo_ids_tmp(1) = topo_id_user
         topo_ids_tmp(2) = topo_id_xpen

         mesh_ids_tmp(1) = mesh_id_user
         mesh_ids_tmp(2) = mesh_id_xpen
      ENDIF
#else
      !-------------------------------------------------------------------------
      !  X-PENCIL
      !-------------------------------------------------------------------------
      IF(Its_xpencil_topo) THEN
         topo_id_xpen         = topo_id_user
         mesh_id_xpen         = mesh_id_user
         topo_id_ypen         = topo_ids(2)
         mesh_id_xpen_complex = mesh_ids(2)
         mesh_id_ypen         = mesh_ids(3)
#if   __MESH_DIM == __3D
         topo_id_zpen         = topo_ids(3)
         mesh_id_zpen         = mesh_ids(4)
#endif
      ELSE
         topo_id_xpen         = topo_ids(1)
         topo_id_ypen         = topo_ids(2)
         mesh_id_xpen         = mesh_ids(1)
         mesh_id_xpen_complex = mesh_ids(2)
         mesh_id_ypen         = mesh_ids(3)
#if   __MESH_DIM == __3D
         topo_id_zpen         = topo_ids(3)
         mesh_id_zpen         = mesh_ids(4)
#endif
      ENDIF

      IF (ppm_debug .GT. 1) THEN
         WRITE(mesg,'(A)' ) '  ID             topo  mesh' 
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A)' ) '-----------------------------------'
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'Original        ',topo_id_user, mesh_id_user  
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'X Pencil        ', topo_id_xpen, mesh_id_xpen
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'X Pencil Complex', topo_id_xpen, &
     &       mesh_id_xpen_complex
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'Y Pencil Complex', topo_id_ypen, mesh_id_ypen
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'Z Pencil Complex', topo_id_zpen, mesh_id_zpen
         CALL ppm_write(ppm_rank,'ppm_util_fft_forward_3d',mesg,j)
      ENDIF

      IF (.NOT.Its_xpencil_topo) THEN
          decomp = ppm_param_decomp_xpencil
          CALL ppm_mktopo(xp,Npart,Nm,decomp,asign,min_phys,max_phys,       &
     &                    bcdef,ghostsize,topo_id_xpen,mesh_id_xpen,min_sub,&
     &                    max_sub,cost,sub2proc,nsubs,isublist,             &
     &                    nsublist,istart,ndata,info)

         topo_ids_tmp(1) = topo_id_user
         topo_ids_tmp(2) = topo_id_xpen

         mesh_ids_tmp(1) = mesh_id_user
         mesh_ids_tmp(2) = mesh_id_xpen
      ENDIF

      !-------------------------------------------------------------------------
      !  Y-PENCIL
      !-------------------------------------------------------------------------
      decomp = ppm_param_decomp_ypencil
      CALL ppm_mktopo(xp,Npart,Nm_com,decomp,asign,min_phys,max_phys,   &
     &                bcdef,ghostsize,topo_id_ypen,mesh_id_ypen,min_sub,&
     &                max_sub,cost,sub2proc,nsubs,isublist,             &
     &                nsublist,istart_ypen,ndata_ypen,info)
#endif

#if __MESH_DIM == __3D
      !-------------------------------------------------------------------------
      !  Z-PENCIL
      !-------------------------------------------------------------------------
      decomp = ppm_param_decomp_zpencil
      CALL ppm_mktopo(xp,Npart,Nm_com,decomp,asign,min_phys,max_phys,   &
     &                bcdef,ghostsize,topo_id_zpen,mesh_id_zpen,min_sub,&
     &                max_sub,cost,sub2proc,nsubs,isublist,             &
     &                nsublist,istart_zpen,ndata_zpen,info)
#endif

#if   __CASE == __SLAB
      !-------------------------------------------------------------------------
      !  Create Plan for xy slab topology
      !-------------------------------------------------------------------------
      idom = isublist(nsublist)
      ! substract 1 to fit ppm-convention
      Nx_in=ndata_slab(1,idom)-1
      Ny_in=ndata_slab(2,idom)-1
#if __MESH_DIM == __3D
      Nz_in=ndata_slab(3,idom)
#endif
      Nx_out = Nx_in/2 + 1
      Ny_out = Ny_in 
      Nz_out = Nz_in

      !-------------------------------------------------------------------------
      !  FFTW 
      !-------------------------------------------------------------------------
#ifdef __FFTW
#if __MESH_DIM == __3D
      MBRank    = 2
      MBIstride = 1
      MBHowmany = Nz_in
      lda(1)    = Nx_in+1
      lda(2)    = Ny_in+1
      lda(3)    = Nz_in
 
      iopt      = ppm_param_alloc_fit
      CALL ppm_alloc(data_real,lda,iopt,info)
      data_real = 1.0_MK

      lda(1) = Nx_out
      lda(2) = Ny_out + 1
      lda(3) = Nz_out
      CALL ppm_alloc(data_comp,lda,iopt,info)
      data_comp = 1.0_MK

      Nm_slab(1) = Nx_in
      Nm_slab(2) = Ny_in

      iEmbed_slab(1)  = Nx_in+1
      iEmbed_slab(2)  = Ny_in+1
      oEmbed_slab(1)  = Nx_out 
      oEmbed_slab(2)  = Ny_out + 1
      MBiDist         = (Nx_in+1) *  (Ny_in+1)
      MBoDist         =  Nx_out * (Ny_out+1)
#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft_r2c(Plan_slab_fd_s,MBRank,Nm_slab(1),MBHowmany,&
     &     data_real(1,1,1),iEmbed_slab(1),MBIstride,MBiDist, &
     &     data_comp(1,1,1),oEmbed_slab(1),MBIstride,MBoDist,FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft_r2c(Plan_slab_fd_d,MBRank,Nm_slab,MBHowmany, &
     &     data_real(1,1,1),iEmbed_slab,MBIstride,MBiDist, &
     &     data_comp(1,1,1),oEmbed_slab,MBIstride,MBoDist,FFTW_MEASURE)
#endif

      oEmbed_slab(1)  = Nx_in+1
      oEmbed_slab(2)  = Ny_in+1
      iEmbed_slab(1)  = Nx_out 
      iEmbed_slab(2)  = Ny_out + 1
      MBoDist         = (Nx_in+1) *  (Ny_in+1)
      MBiDist         =  Nx_out * (Ny_out+1)

#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft_c2r(Plan_slab_bd_s,MBRank,Nm_slab,MBHowmany, &
     &     data_comp(1,1,1),iEmbed_slab,MBIstride,MBiDist, &
     &     data_real(1,1,1),oEmbed_slab,MBIstride,MBoDist,FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft_c2r(Plan_slab_bd_d,MBRank,Nm_slab,MBHowmany, &
     &     data_comp(1,1,1), iEmbed_slab,MBIstride,MBiDist, &
     &     data_real(1,1,1), oEmbed_slab,MBIstride,MBoDist,FFTW_MEASURE)
#endif

       CALL ppm_alloc(data_real,lda,ppm_param_dealloc,info)
#endif         
#endif

#else
      !-------------------------------------------------------------------------
      !  Create Plan for xpencil topology
      !-------------------------------------------------------------------------
      idom = isublist(nsublist)
      ! substract 1 to fit ppm-convention
      Nx_in=ndata(1,idom)-1
      Ny_in=ndata(2,idom)
#if __MESH_DIM == __3D
      Nz_in=ndata(3,idom)
#endif
      Nx_out = Nx_in/2 + 1
      Ny_out = Ny_in
      Nz_out = Nz_in

#ifdef __FFTW
      !-------------------------------------------------------------------------
      !  FFTW 
      !-------------------------------------------------------------------------
#if __MESH_DIM == __2D
      MBRank      = 1
      MBiEmbed(1) = -1
      MBoEmbed(1) = -1
      MBIstride   = 1
      MBHowmany   = Ny_in

      iopt   = ppm_param_alloc_fit
      lda(1) = Nx_in+1
      lda(2) = Ny_in
      CALL ppm_alloc(data_real,lda,iopt,info)

      lda(1) = Nx_out
      lda(2) = Ny_out
      CALL ppm_alloc(data_comp,lda,iopt,info)

      MBiDist    = Nx_in+1
      MBoDist    = Nx_out
      MB_in (1)  = Nx_in
#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft_r2c(Plan_fd_s,MBRank,MB_in,MBHowmany, &
     &     data_real(1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_comp(1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
      MBiDist    = Nx_out
      MBoDist    = Nx_in+1
      CALL sfftw_plan_many_dft_c2r(Plan_bd_s,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_real(1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft_r2c(Plan_fd_d,MBRank,MB_in,MBHowmany, &
     &     data_real(1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_comp(1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
      MBiDist    = Nx_out
      MBoDist    = Nx_in+1
      CALL dfftw_plan_many_dft_c2r(Plan_bd_d,MBRank,MB_in,MBHowmany, &
     &     data_real(1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_comp(1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
#endif
#endif

#if __MESH_DIM == __3D
      MBRank      = 1
      MBiEmbed(1) = -1
      MBoEmbed(1) = -1
      MBIstride   = 1
      MBHowmany   = Ny_in*Nz_in

      iopt   = ppm_param_alloc_fit
      lda(1) = Nx_in+1
      lda(2) = Ny_in
      lda(3) = Nz_in
      CALL ppm_alloc(data_real,lda,iopt,info)

      lda(1) = Nx_out
      lda(2) = Ny_out
      lda(3) = Nz_out
      CALL ppm_alloc(data_comp,lda,iopt,info)

      MBiDist    = Nx_in+1
      MBoDist    = Nx_out
      MB_in (1)  = Nx_in
#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft_r2c(Plan_fd_s,MBRank,MB_in,MBHowmany, &
     &     data_real(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_comp(1,1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
      MBiDist    = Nx_out
      MBoDist    = Nx_in+1
      CALL sfftw_plan_many_dft_c2r(Plan_bd_s,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_real(1,1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft_r2c(Plan_fd_d,MBRank,MB_in,MBHowmany, &
     &     data_real(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_comp(1,1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
      MBiDist    = Nx_out
      MBoDist    = Nx_in+1
      CALL dfftw_plan_many_dft_c2r(Plan_bd_d,MBRank,MB_in,MBHowmany, &
     &     data_real(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_comp(1,1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
#endif
#endif
#endif

#ifdef __MATHKEISAN
      !-------------------------------------------------------------------------
      !  MATHKEISAN 
      !-------------------------------------------------------------------------
      lda_table = 2*Nx_in + 64
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(table_fd_s,lda_table,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_s,lda_table,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',     &
     &        'table_fd_s not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#elif __KIND == __DOUBLE_PRECISION
      CALL ppm_alloc(table_fd_d,lda_table,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_d,lda_table,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',     &
     &        'table_fd_d not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      scale_fft = 1
      isign_fft = 0
#if __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      CALL  scfft(isign_fft,Nx_in,scale_fft,data_real(1,1), &
     &                                 data_comp(1,1),table_fd_s,work,isys)
      CALL  csfft(isign_fft,Nx_in,scale_fft,data_comp(1,1), &
     &                                 data_real(1,1),table_bd_s,work,isys)

#elif __KIND == __DOUBLE_PRECISION
      CALL  dzfft(isign_fft,Nx_in,scale_fft,data_real(1,1), &
     &                                 data_comp(1,1),table_fd_d,work,isys)
      CALL  zdfft(isign_fft,Nx_in,scale_fft,data_comp(1,1), &
     &                                 data_real(1,1),table_bd_d,work,isys)
#endif
#endif

#if __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      CALL  scfft(isign_fft,Nx_in,scale_fft,data_real(1,1,1), &
     &                                 data_comp(1,1,1),table_fd_s,work,isys)
      CALL  csfft(isign_fft,Nx_in,scale_fft,data_comp(1,1,1), &
     &                                 data_real(1,1,1),table_bd_s,work,isys)
#elif __KIND == __DOUBLE_PRECISION
      CALL  dzfft(isign_fft,Nx_in,scale_fft,data_real(1,1,1), &
     &                                 data_comp(1,1,1),table_fd_d,work,isys)
      CALL  zdfft(isign_fft,Nx_in,scale_fft,data_comp(1,1,1), &
     &                                 data_real(1,1,1),table_bd_d,work,isys)
#endif
#endif
#endif

      !-------------------------------------------------------------------------
      !  Create Plan/ Table for ypencil topology
      !-------------------------------------------------------------------------
      idom = isublist(nsublist)
      Nx_in=ndata_ypen(2,idom)-1
      Ny_in=ndata_ypen(1,idom)
#if __MESH_DIM == __3D
      Nz_in=ndata_ypen(3,idom)
#endif
      Nx_out = Nx_in   
      Ny_out = Ny_in
      Nz_out = Nz_in

#ifdef __FFTW
      !-------------------------------------------------------------------------
      !  FFTW 
      !-------------------------------------------------------------------------
#if __MESH_DIM == __2D
      MBRank      = 1
      MBiEmbed(1) = -1
      MBoEmbed(1) = -1
      MBIstride   = 1
      MBiDist     = Nx_in+1
      MBoDist     = Nx_out+1
      MB_in(1)    = Nx_in
      MBHowmany   = Ny_in

      lda(1) = Nx_in+1
      lda(2) = Ny_in
      CALL ppm_alloc(data_comp,lda,iopt,info)
      CALL ppm_alloc(data_compl,lda,iopt,info)

#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft(Plan_fd_c_y,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_FORWARD, FFTW_MEASURE)
      CALL sfftw_plan_many_dft(Plan_bd_c_y,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_BACKWARD, FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft(Plan_fd_cc_y,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_FORWARD, FFTW_MEASURE)
      CALL dfftw_plan_many_dft(Plan_bd_cc_y,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_BACKWARD, FFTW_MEASURE)
#endif
#endif

#if __MESH_DIM == __3D
      MBRank      = 1
      MBiEmbed(1) = -1
      MBoEmbed(1) = -1
      MBIstride   = 1
      MBiDist     = Nx_in+1
      MBoDist     = Nx_out+1
      MB_in(1)    = Nx_in
      MBHowmany   = Ny_in*Nz_in

      lda(1) = Nx_in+1
      lda(2) = Ny_in
      lda(3) = Nz_in
      CALL ppm_alloc(data_comp,lda,iopt,info)
      CALL ppm_alloc(data_compl,lda,iopt,info)

#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft(Plan_fd_c_y,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_FORWARD, FFTW_MEASURE)
      CALL sfftw_plan_many_dft(Plan_bd_c_y,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_BACKWARD, FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft(Plan_fd_cc_y,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_FORWARD, FFTW_MEASURE)
      CALL dfftw_plan_many_dft(Plan_bd_cc_y,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_BACKWARD, FFTW_MEASURE)
#endif
#endif
#endif

#ifdef __MATHKEISAN
      !-------------------------------------------------------------------------
      !  MATHKEISAN 
      !-------------------------------------------------------------------------
      lda_table_y = 2*Nx_in + 64
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(table_fd_c_y,lda_table_y,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_c_y,lda_table_y,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',     &
     &        'table_fd_s not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#elif __KIND == __DOUBLE_PRECISION
      CALL ppm_alloc(table_fd_cc_y,lda_table_y,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_cc_y,lda_table_y,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',     &
     &        'table_fd_d not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      scale_fft = 1
      isign_fft = 0
      incx      = 1
      incy      = 1
#if __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      CALL  cfft(isign_fft,Nx_in,scale_fft,data_comp(1,1),incx, &
     &            data_compl(1,1),incy,table_fd_c_y,lda_table_y,work,1,isys)
#elif __KIND == __DOUBLE_PRECISION
      CALL  zfft(isign_fft,Nx_in,scale_fft,data_comp(1,1),incx, &
     &            data_compl(1,1),incy,table_fd_cc_y,lda_table_y,work,1,isys)
#endif
#endif

#if __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      CALL  cfft(isign_fft,Nx_in,scale_fft,data_comp(1,1,1),incx, &
     &            data_compl(1,1,1),incy,table_fd_c_y,lda_table_y,work,1,isys)
#elif __KIND == __DOUBLE_PRECISION
      CALL  zfft(isign_fft,Nx_in,scale_fft,data_comp(1,1,1),incx, &
     &            data_compl(1,1,1),incy,table_fd_cc_y,lda_table_y,work,1,isys)
#endif
#endif
#endif
#endif 

#if __MESH_DIM == __3D
      !-------------------------------------------------------------------------
      !  Create Plan for zpencil topology
      !-------------------------------------------------------------------------
      idom   = isublist(nsublist)
      Nx_in  = ndata_zpen(3,idom)-1
      Ny_in  = ndata_zpen(2,idom)
      Nz_in  = ndata_zpen(1,idom)
      Nx_out = Nx_in     
      Ny_out = Ny_in
      Nz_out = Nz_in

#ifdef __FFTW
      !-------------------------------------------------------------------------
      !  FFTW 
      !-------------------------------------------------------------------------
      MBRank      = 1
      MBiEmbed(1) = -1
      MBoEmbed(1) = -1
      MBIstride   = 1
      MBiDist     = Nx_in+1
      MBoDist     = Nx_out+1
      MB_in(1)    = Nx_in
      MBHowmany   = Ny_in*Nz_in

      lda(1) = Nx_in+1
      lda(2) = Ny_in
      lda(3) = Nz_in
      CALL ppm_alloc(data_comp,lda,iopt,info)
      CALL ppm_alloc(data_compl,lda,iopt,info)

#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft(Plan_fd_c_z,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_FORWARD, FFTW_MEASURE)
      CALL sfftw_plan_many_dft(Plan_bd_c_z,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_BACKWARD, FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft(Plan_fd_cc_z,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_FORWARD, FFTW_MEASURE)
      CALL dfftw_plan_many_dft(Plan_bd_cc_z,MBRank,MB_in,MBHowmany, &
     &     data_comp(1,1,1),MBiEmbed(1),MBIstride,MBiDist, &
     &     data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
     &     FFTW_BACKWARD, FFTW_MEASURE)
#endif
#endif

#ifdef __MATHKEISAN
      !-------------------------------------------------------------------------
      !  MATHKEISAN 
      !-------------------------------------------------------------------------
      lda_table_z = 2*Nx_in + 64
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(table_fd_c_z,lda_table_z,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_c_z,lda_table_z,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',     &
     &        'table_fd_s not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#elif __KIND == __DOUBLE_PRECISION
      CALL ppm_alloc(table_fd_cc_z,lda_table_z,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_cc_z,lda_table_z,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_fft_init',     &
     &        'table_fd_d not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      scale_fft = 1
      isign_fft = 0
      incx      = 1
      incy      = 1
#if   __KIND == __SINGLE_PRECISION
      CALL  cfft(isign_fft,Nx_in,scale_fft,data_comp(1,1,1),incx, &
     &            data_compl(1,1,1),incy,table_fd_c_z,lda_table_z,work,1,isys)
#elif __KIND == __DOUBLE_PRECISION
      CALL  zfft(isign_fft,Nx_in,scale_fft,data_comp(1,1),incx, &
     &            data_compl(1,1,1),incy,table_fd_cc_z,lda_table_z,work,1,isys)
#endif
#endif
#endif

      !-------------------------------------------------------------------------
      !  Deallocate memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(isublist,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_init',   &
     &        'temporary isublist',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ndata,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_init',   &
     &        'temporary ndata',__LINE__,info)
      ENDIF
#ifdef __MATHKEISAN
      CALL ppm_alloc(work,1,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_init',   &
     &        'temporary work array for Mathkeisan',__LINE__,info)
      ENDIF
#endif
      CALL ppm_alloc(data_real,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_init',   &
     &        'temporary real_data',__LINE__,info)
      ENDIF
      CALL ppm_alloc(data_comp,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_init',   &
     &        'temporary data_comp',__LINE__,info)
      ENDIF
      CALL ppm_alloc(data_compl,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_fft_init',   &
     &        'temporary data_compl',__LINE__,info)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_fft_init',t0,info)
      RETURN

#if   __CASE == __SLAB

#if   __DIM == __SFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_2d_sca_s 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_2d_sca_d 
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_3d_sca_s 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_3d_sca_d 
#endif
#endif
#endif

#if   __DIM == __VFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_2d_vec_s 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_3d_vec_d
#endif
#endif
#endif

#else

#if   __DIM == __SFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_pen_2d_sca_s 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_pen_2d_sca_d 
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_pen_3d_sca_s 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_pen_3d_sca_d 
#endif
#endif
#endif

#if   __DIM == __VFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_pen_2d_vec_s 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_pen_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_pen_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_fft_init_pen_3d_vec_d
#endif
#endif
#endif

#endif

