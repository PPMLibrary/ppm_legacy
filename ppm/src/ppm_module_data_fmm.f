      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_data_fmm
      !-------------------------------------------------------------------------
      !
      ! Purpose       : fast mulipole method, data
      !
      !
      ! Remarks       :
      !
      ! References    :
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_fmm.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:52  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.12  2006/06/29 10:28:36  pchatela
      !  Added vector strengths support
      !
      !  Revision 1.11  2006/06/16 07:50:28  hiebers
      !  added list of topo IDs
      !
      !  Revision 1.10  2005/09/19 13:03:30  polasekb
      !  code cosmetics
      !
      !  Revision 1.9  2005/09/12 13:30:19  polasekb
      !  added ppm_subid
      !
      !  Revision 1.8  2005/08/11 15:13:33  polasekb
      !  added maxboxcost
      !
      !  Revision 1.7  2005/08/08 13:34:44  polasekb
      !  removed fmm_prec
      !
      !  Revision 1.6  2005/08/04 16:02:46  polasekb
      !  addes some new data
      !
      !  Revision 1.5  2005/07/29 12:36:51  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.4  2005/07/27 21:11:54  polasekb
      !  added totalmass (again)
      !
      !  Revision 1.3  2005/07/25 14:28:32  polasekb
      !  added some constants for the spherical harmonics
      !
      !  Revision 1.2  2005/06/02 13:55:16  polasekb
      !  removed totalmass
      !
      !  Revision 1.1  2005/05/27 08:04:09  polasekb
      !  initial implementation
      !
      !
      !	 Revision 0 2004/11/11 4:04:15 polasekb
      !  Start.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

MODULE ppm_module_data_fmm
      !-------------------------------------------------------------------------
      !Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data,ONLY:ppm_kind_single,ppm_kind_double
      PRIVATE :: ppm_kind_single,ppm_kind_double

      !-------------------------------------------------------------------------
      ! Define Initialization-FLAG
      !-------------------------------------------------------------------------
      LOGICAL                                       :: ppm_fmm_initialized = .FALSE.

      !-------------------------------------------------------------------------
      ! Define order of expansion
      !-------------------------------------------------------------------------
      INTEGER                                       :: order

      !-------------------------------------------------------------------------
      ! Define radius of tree boxes
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:),POINTER    :: radius_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:),POINTER    :: radius_d => NULL()

      !-------------------------------------------------------------------------
      ! Define expansions of all tree boxes
      ! 1st index: boxid
      ! 2nd/3rd index: expansion
      !-------------------------------------------------------------------------
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:)  ,POINTER :: expansion_s_sf => NULL()
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:)  ,POINTER :: expansion_d_sf => NULL()
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:),POINTER :: expansion_s_vf => NULL()
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:),POINTER :: expansion_d_vf => NULL()

      !-------------------------------------------------------------------------
      ! Define center of mass of tree boxes
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:,:),POINTER   :: centerofbox_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:,:),POINTER   :: centerofbox_d => NULL()

      !-------------------------------------------------------------------------
      ! Define totalmass of tree boxes
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:),POINTER     :: totalmass_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:),POINTER     :: totalmass_d => NULL()

      !-------------------------------------------------------------------------
      ! Store tree output in data file
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! Define min_box, minimum extent of tree boxes
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:,:),POINTER   :: min_box_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:,:),POINTER   :: min_box_d => NULL()

      !-------------------------------------------------------------------------
      ! Define max_box, maximum extent of tree boxes
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:,:),POINTER   :: max_box_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:,:),POINTER   :: max_box_d => NULL()

      !-------------------------------------------------------------------------
      ! Define nbox, total number of boxes
      !-------------------------------------------------------------------------
      INTEGER                                        :: nbox

      !-------------------------------------------------------------------------
      ! Define nchld, number of children of the box
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:),POINTER                  :: nchld

      !-------------------------------------------------------------------------
      ! Define lhbx, pointer to first and last point in box (in lpdx)
      ! 1st index: 1 or 2, first and last
      ! 2nd index: box id
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:,:),POINTER               :: lhbx

      !-------------------------------------------------------------------------
      ! Define lpdx, permutation of xp, particles ordered according to tree
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:),POINTER                :: lpdx

      !-------------------------------------------------------------------------
      ! Define boxcost
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:),POINTER   :: boxcost_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:),POINTER   :: boxcost_d => NULL()

      !-------------------------------------------------------------------------
      ! Define parent, the partent of the box
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:),POINTER                 :: parent => NULL()

      !-------------------------------------------------------------------------
      ! Define   child, the child ids of the box
      ! 1st index: child number (1-8 in octtree)
      ! 2nd index: box id
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:,:),POINTER               :: child => NULL()

      !-------------------------------------------------------------------------
      ! Define blevel, level of box
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:),POINTER                 :: blevel => NULL()

      !-------------------------------------------------------------------------
      ! Define nbpl, number of boxes per level
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:),POINTER                 :: nbpl => NULL()

      !-------------------------------------------------------------------------
      ! Define nlevel, total number of levels
      !-------------------------------------------------------------------------
      INTEGER                                       :: nlevel

      !-------------------------------------------------------------------------
      ! Define list of topo ids
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:),POINTER                 :: topoidlist => NULL()

      !-------------------------------------------------------------------------
      ! End of Tree data definition
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! Define boxid, box id of sub id
      ! 1st index: sub id
      ! 2nd index: user topology id
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:,:),POINTER               :: ppm_boxid => NULL()

      !-------------------------------------------------------------------------
      ! Define subid. sub id of box id
      ! 1st index: box id
      ! 2nd index: user topology id
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:,:),POINTER               :: ppm_subid => NULL()

      !-------------------------------------------------------------------------
      ! Define boxpart, which particle is in which box
      !-------------------------------------------------------------------------
      INTEGER,DIMENSION(:),POINTER               :: boxpart => NULL()

      !-------------------------------------------------------------------------
      ! Define maxboxcost, the maximum nr of particles per box
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single)                      :: maxboxcost_s
      REAL(ppm_kind_double)                      :: maxboxcost_d

      !-------------------------------------------------------------------------
      ! Data for spherical harmonics
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Define Anm
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:,:),POINTER :: Anm_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: Anm_d => NULL()

      !-------------------------------------------------------------------------
      ! Define sqrtfac
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:,:),POINTER :: sqrtfac_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: sqrtfac_d => NULL()

      !-------------------------------------------------------------------------
      ! Define Cnm
      !-------------------------------------------------------------------------
      COMPLEX(ppm_kind_single),DIMENSION(:,:),POINTER   :: Cnm_s_sf => NULL()
      COMPLEX(ppm_kind_double),DIMENSION(:,:),POINTER   :: Cnm_d_sf => NULL()
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:),POINTER :: Cnm_s_vf => NULL()
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:),POINTER :: Cnm_d_vf => NULL()

      !-------------------------------------------------------------------------
      ! Define Inner
      !-------------------------------------------------------------------------
      COMPLEX(ppm_kind_single),DIMENSION(:,:),POINTER :: Inner_s => NULL()
      COMPLEX(ppm_kind_double),DIMENSION(:,:),POINTER :: Inner_d => NULL()

      !-------------------------------------------------------------------------
      ! Define Outer
      !-------------------------------------------------------------------------
      COMPLEX(ppm_kind_single),DIMENSION(:,:),POINTER :: Outer_s => NULL()
      COMPLEX(ppm_kind_double),DIMENSION(:,:),POINTER :: Outer_d => NULL()

      !-------------------------------------------------------------------------
      ! Define Ynm
      !-------------------------------------------------------------------------
      COMPLEX(ppm_kind_single),DIMENSION(:,:),POINTER :: Ynm_s => NULL()
      COMPLEX(ppm_kind_double),DIMENSION(:,:),POINTER :: Ynm_d => NULL()

      !-------------------------------------------------------------------------
      ! Define Pnm
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:,:),POINTER :: Pnm_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: Pnm_d => NULL()

      !-------------------------------------------------------------------------
      ! Define fracfac
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:  ),POINTER :: fracfac_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:  ),POINTER :: fracfac_d => NULL()

      !-------------------------------------------------------------------------
      ! Define rho
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:  ),POINTER :: rho_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:  ),POINTER :: rho_d => NULL()

      !-------------------------------------------------------------------------
      ! Define theta
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:  ),POINTER :: theta_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:  ),POINTER :: theta_d => NULL()

      !-------------------------------------------------------------------------
      ! Define phi
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:  ),POINTER :: phi_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:  ),POINTER :: phi_d => NULL()

      !-------------------------------------------------------------------------
      ! Define fac
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single),DIMENSION(:  ),POINTER :: fac_s => NULL()
      REAL(ppm_kind_double),DIMENSION(:  ),POINTER :: fac_d => NULL()

END MODULE ppm_module_data_fmm
