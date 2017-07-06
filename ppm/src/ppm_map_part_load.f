      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_map_part_load
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine loads the internally stored particle 
      !                 mapping.
      !
      !  Input        : 
      !                                    
      !  Output       : info         (I) : return status, 0 on success
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_load.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2006/10/10 21:10:59  walther
      !  *** empty log message ***
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_part_load(info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_state
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  This routine is ALL integer, except the timing variable t0
      !  We pick double as the precision - no real reason for that 
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER    :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(OUT)  :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(1) :: ldu
      INTEGER               :: i,j,k
      INTEGER               :: iopt
      REAL(MK)              :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_load',t0,info)

      !-------------------------------------------------------------------------
      !  load the map type (is used in _send() line 229)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_map_type_state 

      !-------------------------------------------------------------------------
      !  load the ppm_nsendlist: the number of processers to send to
      !-------------------------------------------------------------------------
      ppm_nsendlist = ppm_nsendlist_state

      !-------------------------------------------------------------------------
      !  load the ppm_nsendbuffer: the size of the send buffer - here zero
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ppm_nsendbuffer_state 

      !-------------------------------------------------------------------------
      !  load the ppm_buffer_set: the current number of sets stored - here zero
      !-------------------------------------------------------------------------
      ppm_buffer_set = ppm_buffer_set_state 

      !-------------------------------------------------------------------------
      !  load the ppm_psendbuffer:
      !  first allocate memory for it
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist + 1
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)

      !-------------------------------------------------------------------------
      !  then load it: pointer to the first particle in the send buffer
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsendlist+1
         ppm_psendbuffer(k) = ppm_psendbuffer_state(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  load the ppm_buffer2part
      !  first allocate memory for it
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = SIZE(ppm_buffer2part_state)
      CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
 
      !-------------------------------------------------------------------------
      !  then load it: pointer to the particle id of the j-th entry in the 
      !  send buffer
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsendlist
         DO j=ppm_psendbuffer(k),ppm_psendbuffer(k+1)-1
            ppm_buffer2part(j) = ppm_buffer2part_state(j)
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Store the ppm_nrecvlist and ppm_sendlist 
      !-------------------------------------------------------------------------
      ppm_nrecvlist = ppm_nrecvlist_state 

      ppm_nsendlist = ppm_nsendlist_state 

      !-------------------------------------------------------------------------
      !  load the receive list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      DO k=1,ppm_nrecvlist
         ppm_irecvlist(k) = ppm_irecvlist_state(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  load the send list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      DO k=1,ppm_nsendlist
         ppm_isendlist(k) = ppm_isendlist_state(k) 
      ENDDO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_load',t0,info)
      RETURN
      END SUBROUTINE ppm_map_part_load
