!>
!!  Subroutine to initialize the CBL case for HDCP2
!!
!!
!! @par Revision History
!! - first version by Anurag Dipankar , MPIM, (2012-12-12)
!! @par Literature
!! -
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_nh_torus_exp
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_io_units,            ONLY: find_next_free_unit
  USE mo_physical_constants,  ONLY: rd, cpd, p0ref, cvd_o_rd, rd_o_cpd, &
     &                              grav, alv, vtmpc1
  !SBr, CHa: add th_prof
  USE mo_nh_testcases_nml,    ONLY: u_cbl, v_cbl, th_cbl, th_prof, psfc_cbl, &
                                    bubctr_x, bubctr_y, nh_test_name
  USE mo_nh_wk_exp,           ONLY: bub_amp, bub_ver_width, bub_hor_width, bubctr_z
  USE mo_model_domain,        ONLY: t_patch
  USE mo_math_constants,      ONLY: pi, rad2deg, pi_2
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_nh_ref
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_math_utilities,      ONLY: plane_torus_distance
  USE mo_sync,                ONLY: sync_patch_array, SYNC_C
  USE mo_nh_init_utils,       ONLY: init_w
  USE mo_run_config,          ONLY: iqv, iqc
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_e
  USE mo_satad,               ONLY: spec_humi, sat_pres_water
  USE mo_les_config,          ONLY: les_config
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_les_utilities,       ONLY: vert_intp_linear_1d

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: init_nh_state_cbl, cbl_stevens_fluxes, init_nh_state_rico, &
            sfcflx_uniform, init_torus_with_sounding, init_warm_bubble

  !add parameters for Patton setup, SBr, CHa
  !DEFINED PARAMETERS (Patton et al. 2005 JAS) for init_nh_state_cbl:
  REAL(wp), PARAMETER :: zh_inv1     = 775._wp      !bottom height of capping inversion 
  REAL(wp), PARAMETER :: zh_inv2     = 875._wp      !top height of capping inversion

!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
!
  !>
  !! Initialization of prognostic state vector for the nh CBL test case 
  !!  without moisture
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_state_cbl( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
    &                           ptr_int, ptr_metrics)

    TYPE(t_patch),TARGET,  INTENT(IN)   ::  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   ::  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT)::  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT)::  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   ::  ptr_metrics      
    TYPE(t_nh_ref),        INTENT(INOUT)::  ptr_nh_ref

    REAL(wp) :: z_exner_h(1:nproma,ptr_patch%nlev+1), z_help(1:nproma) 
    REAL(wp) :: zvn1, zvn2, zu, zv, zt00, zh00, ex_sfc
    INTEGER  :: jc,jk,jb,i_startblk,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e,npromz_e
    INTEGER  :: nlev, nlevp1                  !< number of full and half levels
    INTEGER  :: nlen, jcn, jbn, jg, ntropo

    REAL(wp), PARAMETER :: zh0     = 0._wp      !< height (m) above which temperature increases
    REAL(wp), PARAMETER :: lambda  = 1500._wp   !moist height from Stevens(2007)
    REAL(wp), PARAMETER :: dtdz_st = 0.04_wp    !< theta lapse rate in stratosphere (T>0!)
    REAL(wp), PARAMETER :: z_tropo = 11000._wp  !height tropopause
    REAL(wp), PARAMETER :: rh_sfc  = 0.8_wp    !RH at surface [1], default 0.8

  !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e
    npromz_e = ptr_patch%npromz_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    !patch id
    jg = ptr_patch%id

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = psfc_cbl
    ex_sfc   = (psfc_cbl/p0ref)**rd_o_cpd
    les_config(jg)%psfc = psfc_cbl

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp
    !SBr,CHa: set qv to its minimal value
    ptr_nh_prog%tracer(:,:,:,iqv) = 0.e-12_wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      !Tracers
      IF(.NOT.les_config(jg)%is_dry_cbl)THEN
        DO jk = 1, nlev
          ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = rh_sfc * spec_humi(sat_pres_water(th_cbl(1)),psfc_cbl) * &
                    EXP(-ptr_metrics%z_mc(1:nlen,jk,jb)/lambda)

        END DO
      END IF

      ntropo = 0
      DO jk = nlev, 1, -1
         ! init potential temperature
         !z_help(1:nlen) = th_cbl(1) + max(0._wp, (ptr_metrics%z_mc(1:nlen,jk,jb)-zh0)*th_cbl(2))

         !SBr, Cha: add th_prof
         if (th_prof == 1) then
         z_help(1:nlen) = th_cbl(1) + max(0._wp, (ptr_metrics%z_mc(1:nlen,jk,jb)-zh0)*th_cbl(2))
         
         elseif (th_prof == 2) then
            if (ptr_metrics%z_mc(1,jk,jb) > zh_inv2 ) then
               z_help(1:nlen) = th_cbl(1) + th_cbl(2)*(zh_inv2 - zh_inv1) + max(0._wp, (ptr_metrics%z_mc(1:nlen,jk,jb)-zh_inv2)*th_cbl(3))
            elseif ((ptr_metrics%z_mc(1,jk,jb) > zh_inv1) .and. (ptr_metrics%z_mc(1,jk,jb) <= zh_inv2)) then 
               z_help(1:nlen) = th_cbl(1) + max(0._wp, (ptr_metrics%z_mc(1:nlen,jk,jb)-zh_inv1)*th_cbl(2))
            elseif ((ptr_metrics%z_mc(1,jk,jb) <= zh_inv1)) then
               z_help(1:nlen) = th_cbl(1)
            endif !z_ifc

         else
            z_help(1:nlen) = th_cbl(1)
         endif !th_prof

         ! constant temperature above tropopause
         if ((ptr_metrics%z_mc(1,jk,jb) > z_tropo) .and. (ntropo == 0)) then
            ntropo = 1
            zt00   = z_help(1)
            zh00   = ptr_metrics%z_mc(1,jk,jb)
         endif
         if (ptr_metrics%z_mc(1,jk,jb) > z_tropo) then
            z_help(1:nlen) = zt00 + (ptr_metrics%z_mc(1:nlen,jk,jb)-zh00) * dtdz_st
         endif

         ! virtual potential temperature
         ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO

      !Get hydrostatic exner at the surface using surface pressure 
      z_exner_h(1:nlen,nlevp1) = ex_sfc
 
      !Get exner at full levels starting from exner at surface
      DO jk = nlev, 1, -1
         !exner at next half level after surface
         z_exner_h(1:nlen,jk) = z_exner_h(1:nlen,jk+1) - grav/cpd *     &
                                ptr_metrics%ddqz_z_full(1:nlen,jk,jb)/ &
                                ptr_nh_prog%theta_v(1:nlen,jk,jb)
        
         !exner at main levels
         ptr_nh_prog%exner(1:nlen,jk,jb) = 0.5_wp * &
                                     (z_exner_h(1:nlen,jk)+z_exner_h(1:nlen,jk+1))
      END DO

      DO jk = 1 , nlev
         ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO !jk

    ENDDO !jb

!--------------------------------------------------------------------------------
    !Mean wind 
!--------------------------------------------------------------------------------
    i_startblk = ptr_patch%edges%start_blk(2,1)
    DO jb = i_startblk , nblks_e
     CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)
     DO jk = 1 , nlev 
      DO jc = i_startidx, i_endidx

        !Torus geometry is flat so zu is only function of height which is same for all cells
        !But it is kept varyign with jc,jb to introduce topography lateron
        jcn  =   ptr_patch%edges%cell_idx(jc,jb,1)
        jbn  =   ptr_patch%edges%cell_blk(jc,jb,1)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)

        zvn1 =  zu * ptr_patch%edges%primal_normal_cell(jc,jb,1)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(jc,jb,1)%v2      
 
        jcn  =   ptr_patch%edges%cell_idx(jc,jb,2)
        jbn  =   ptr_patch%edges%cell_blk(jc,jb,2)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
      
        zvn2 =  zu * ptr_patch%edges%primal_normal_cell(jc,jb,2)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(jc,jb,2)%v2      

        ptr_nh_prog%vn(jc,jk,jb) = ptr_int%c_lin_e(jc,1,jb)*zvn1 + &
                                   ptr_int%c_lin_e(jc,2,jb)*zvn2

        ptr_nh_ref%vn_ref(jc,jk,jb) = ptr_nh_prog%vn(jc,jk,jb)
      END DO
     END DO
    END DO     
    
    !W wind and reference
    CALL init_w(ptr_patch, ptr_int, ptr_nh_prog%vn, ptr_metrics%z_ifc, ptr_nh_prog%w)
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)
    ptr_nh_ref%w_ref = ptr_nh_prog%w


  END SUBROUTINE init_nh_state_cbl

  !>
  !! Initialization of prognostic state vector for the nh GATE test case 
  !!  with moisture
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_torus_with_sounding ( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
                                        ptr_int, ptr_metrics)

    TYPE(t_patch),TARGET,  INTENT(IN)   ::  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   ::  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT)::  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT)::  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   ::  ptr_metrics      
    TYPE(t_nh_ref),        INTENT(INOUT)::  ptr_nh_ref

    INTEGER  :: je,jk,jb,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e,npromz_e, jg
    INTEGER  :: nlev, nlevp1                        !< number of full and half levels
    INTEGER  :: nlen, i_rcstartlev, jcn, jbn, ist

    REAL(wp), DIMENSION(ptr_patch%nlev) :: theta_in, qv_in, u_in, v_in
    REAL(wp) :: zvn1, zvn2, zu, zv, psfc_in, ex_sfc
    REAL(wp) :: z_exner_h(1:nproma,ptr_patch%nlev+1), z_help(1:nproma) 

    CHARACTER(len=max_char_length), PARAMETER :: &
       &  routine = 'mo_nh_torus_exp:init_torus_with_sounding'
  !-------------------------------------------------------------------------
    
    !Read the sounding file
    CALL  read_ext_profile (ptr_metrics%z_mc(2,:,2), theta_in, qv_in, u_in, v_in, &
                            psfc_in)

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e
    npromz_e = ptr_patch%npromz_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    !patch id
    jg = ptr_patch%id

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = psfc_in
    ex_sfc   = (psfc_in/p0ref)**rd_o_cpd
    les_config(jg)%psfc = psfc_in

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = qv_in(jk)
      END DO

      DO jk = nlev, 1, -1
         z_help(1:nlen) = theta_in(jk)
         ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO

      !Get hydrostatic exner at the surface using surface pressure 
      z_exner_h(1:nlen,nlevp1) = ex_sfc
 
      !Get exner at full levels starting from exner at surface
      DO jk = nlev, 1, -1
         !exner at next half level after surface
         z_exner_h(1:nlen,jk) = z_exner_h(1:nlen,jk+1) - grav/cpd *     &
                                ptr_metrics%ddqz_z_full(1:nlen,jk,jb)/ &
                                ptr_nh_prog%theta_v(1:nlen,jk,jb)
        
         !exner at main levels
         ptr_nh_prog%exner(1:nlen,jk,jb) = 0.5_wp * &
                                     (z_exner_h(1:nlen,jk)+z_exner_h(1:nlen,jk+1))
      END DO

      DO jk = 1 , nlev
        ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO !jk

    ENDDO !jb

    !CALL diagnose_pres_temp (ptr_metrics, ptr_nh_prog,ptr_nh_prog, ptr_nh_diag,     &
    !                         ptr_patch, opt_calc_pres=.TRUE., opt_calc_temp=.TRUE.)

    !Mean wind 
    DO jb = 1 , nblks_e
     CALL get_indices_e( ptr_patch, jb, 1, nblks_e, i_startidx, i_endidx, grf_bdywidth_e+1)
     DO jk = 1 , nlev 
      DO je = i_startidx, i_endidx

        jcn  =   ptr_patch%edges%cell_idx(je,jb,1)
        jbn  =   ptr_patch%edges%cell_blk(je,jb,1)
        zu   =   u_in(jk)
        zv   =   v_in(jk)

        zvn1 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,1)%v2      
 
        jcn  =   ptr_patch%edges%cell_idx(je,jb,2)
        jbn  =   ptr_patch%edges%cell_blk(je,jb,2)
        zu   =   u_in(jk)
        zv   =   v_in(jk)
      
        zvn2 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,2)%v2      

        ptr_nh_prog%vn(je,jk,jb) = ptr_int%c_lin_e(je,1,jb)*zvn1 + &
                                   ptr_int%c_lin_e(je,2,jb)*zvn2

        ptr_nh_ref%vn_ref(je,jk,jb) = ptr_nh_prog%vn(je,jk,jb)
      END DO
     END DO
    END DO     
    
    !W wind and reference
    CALL init_w(ptr_patch, ptr_int, ptr_nh_prog%vn, ptr_metrics%z_ifc, ptr_nh_prog%w)
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)
    ptr_nh_ref%w_ref  = ptr_nh_prog%w


  END SUBROUTINE init_torus_with_sounding

!-------------------------------------------------------------------------
  
  !>
  !! Initialization of prognostic state vector for the nh RICO test case 
  !!  with moisture
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_state_rico( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
    &                           ptr_int, ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(IN)   :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   :: &
      &  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   :: &
      &  ptr_metrics                          !< NH metrics state
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    REAL(wp) :: rho_sfc, z_help(1:nproma), zvn1, zvn2, zu, zv
    INTEGER  :: je,jk,jb,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e,npromz_e
    INTEGER  :: nlev, nlevp1                        !< number of full and half levels
    INTEGER  :: nlen, i_rcstartlev, jcn, jbn, jg

    !DEFINED PARAMETERS (RICO case):
    REAL(wp), PARAMETER :: zpsfc   = 101540._wp   !< surface pressure
    REAL(wp), PARAMETER :: zh1     = 740._wp      !< height (m) above which temperature increases
    REAL(wp), PARAMETER :: ztsfc   = 297.9_wp
    REAL(wp), PARAMETER :: zh2     = 3260._wp     !< moist height for RICO
    REAL(wp), PARAMETER :: zh3     = 15000._wp    !< height for extrapolated RICO profiles
    REAL(wp), PARAMETER :: zh4     = 30000._wp    !< height for extrapolated RICO profiles
    REAL(wp), PARAMETER :: zh5     = 60000._wp    !< height for extrapolated RICO profiles

  !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e
    npromz_e = ptr_patch%npromz_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    !patch id
    jg = ptr_patch%id

    !Set some reference density    
    rho_sfc = zpsfc / (rd * les_config(jg)%sst)

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zpsfc

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      !Tracers
      DO jk = 1, nlev
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = min((0.016_wp - 0.0022_wp * ptr_metrics%z_mc(1:nlen,jk,jb)/zh1), &
                                              (0.0138_wp - 0.0114_wp * (ptr_metrics%z_mc(1:nlen,jk,jb)-zh1)/(zh2 - zh1)))
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = max(ptr_nh_prog%tracer(1:nlen,jk,jb,iqv),             &
                                              (0.0024_wp - 0.0006_wp * (ptr_metrics%z_mc(1:nlen,jk,jb) - zh2)/(4000._wp - zh2)))
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = max(ptr_nh_prog%tracer(1:nlen,jk,jb,iqv), 3e-6_wp)                            
      END DO

      DO jk = 1, nlev
       ! init potential temperature
       z_help(1:nlen) = ztsfc + max(0._wp, (ptr_metrics%z_mc(1:nlen,jk,jb)-zh1)*19.1_wp/(4000._wp-zh1))
       z_help(1:nlen) = max(z_help(1:nlen), ztsfc +(ptr_metrics%z_mc(1:nlen,jk,jb)-zh3)*502._wp/(zh4-zh3))
       z_help(1:nlen) = max(z_help(1:nlen), ztsfc +(ptr_metrics%z_mc(1:nlen,jk,jb)-zh4)*2600._wp/(zh5-zh4))

       ! virtual potential temperature
       ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO

      !Get hydrostatic pressure and exner at lowest level
      ptr_nh_diag%pres(1:nlen,nlev,jb) = zpsfc - rho_sfc * ptr_metrics%geopot(1:nlen,nlev,jb)
      ptr_nh_prog%exner(1:nlen,nlev,jb) = (ptr_nh_diag%pres(1:nlen,nlev,jb)/p0ref)**rd_o_cpd 

      !Get exner at other levels
      DO jk = nlev-1, 1, -1
         z_help(1:nlen) = 0.5_wp * ( ptr_nh_prog%theta_v(1:nlen,jk,jb) +  &
                                     ptr_nh_prog%theta_v(1:nlen,jk+1,jb) )
   
         ptr_nh_prog%exner(1:nlen,jk,jb) = ptr_nh_prog%exner(1:nlen,jk+1,jb) &
            &  -grav/cpd*ptr_metrics%ddqz_z_half(1:nlen,jk+1,jb)/z_help(1:nlen)
      END DO

      DO jk = 1 , nlev
        ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO !jk

    ENDDO !jb

    !Mean wind 
    DO jb = 1 , nblks_e
     CALL get_indices_e( ptr_patch, jb, 1, nblks_e, i_startidx, i_endidx, grf_bdywidth_e+1)
     DO jk = 1 , nlev 
      DO je = i_startidx, i_endidx

        !Torus geometry is flat so zu is only function of height which is same for all cells
        !But it is kept varyign with jc,jb to introduce topography lateron
        jcn  =   ptr_patch%edges%cell_idx(je,jb,1)
        jbn  =   ptr_patch%edges%cell_blk(je,jb,1)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zu   =   min(zu, -2._wp)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)

        zvn1 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,1)%v2      
 
        jcn  =   ptr_patch%edges%cell_idx(je,jb,2)
        jbn  =   ptr_patch%edges%cell_blk(je,jb,2)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zu   =   min(zu, -2._wp)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
      
        zvn2 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,2)%v2      

        ptr_nh_prog%vn(je,jk,jb) = ptr_int%c_lin_e(je,1,jb)*zvn1 + &
                                   ptr_int%c_lin_e(je,2,jb)*zvn2

        ptr_nh_ref%vn_ref(je,jk,jb) = ptr_nh_prog%vn(je,jk,jb)
      END DO
     END DO
    END DO     
    
    !W wind and reference
    CALL init_w(ptr_patch, ptr_int, ptr_nh_prog%vn, ptr_metrics%z_ifc, ptr_nh_prog%w)
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)
    ptr_nh_ref%w_ref = ptr_nh_prog%w
    
  END SUBROUTINE init_nh_state_rico
  
!-------------------------------------------------------------------------


  SUBROUTINE cbl_stevens_fluxes( t_l, qv_l, p_l, rho_l, tsk, shfx, lhfx )

  !-------------------------------------------------------------------------
  ! Calculate sensible and latent heat fluxes from buoyancy flux for Stevens 
  ! (2007) case. (code from Wayne Angevine 2013) 
  !
  ! Variable explanations
  ! tsk  = skin temperature
  ! zp0  = surface pressure (Pa)
  ! qv1d = water vapor mixing ratio in lowest level
  ! th1d = potential temperature in lowest level
  ! kqfx = exchange coefficient for water vapor (used in qfx calculation below)
  !
  ! Constants for saturation calculation:
  ! svp1  = 0.6112
  ! svp2  = 17.67
  ! svp3  = 29.65
  !-------------------------------------------------------------------------

    ! INPUT PARAMETERS:
    REAL(wp), INTENT(IN) ::    t_l       , & ! temperature at lowest level [K]
     &                         qv_l      , & ! moisture at lowest level    [kg/kg]
     &                         p_l       , & ! pressure at lowest level    [Pa]
     &                         rho_l         ! rho at lowest level         [kg/m3]

    ! INPUT/OUTPUT PARAMETERS:
    REAL(wp), INTENT(INOUT) :: tsk           ! skin temperature            [K]

    ! OUTPUT PARAMETERS:
    REAL(wp), INTENT(OUT) ::   shfx      , & ! sensible heat flux (+ down) [W/m2]
      &                        lhfx          ! latent heat flux   (+ down) [W/m2]

    ! LOCAL VARIABLES:
    REAL(wp) :: Beta0, Vs, mav, qsfc, qsfc_air, kqfx, khfx, th_l

    REAL(wp), PARAMETER :: zp0     = 101540._wp !p0ref !< surface pressure

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
    Beta0 = 7e-4_wp       ! Fixed specified buoyancy flux (25 W/m^2)
   !Beta0 = 4.2e-4_wp     ! Fixed specified buoyancy flux (15 W/m^2)
   !Beta0 = 11.2e-4_wp    ! Fixed specified buoyancy flux (40 W/m^2)
    Vs    = 0.01_wp       ! Specified velocity scale
    mav   = 0.90_wp       ! Surface moisture availability

  ! Calculate saturation mixing ratio at SST (from previous timestep)
  ! The rest of the buoyancy flux goes to sensible heat
  ! Calculate surface saturated q and q in air at surface 
   !e1=svp1*exp(svp2*(tsk-tmelt)/(tsk-svp3))                       
   !qsfc=rd/rv*e1/((zp0/1000.)-e1)
    qsfc     = spec_humi(sat_pres_water(tsk),zp0)
    qsfc_air = qsfc * mav

    th_l =  t_l * (p0ref/p_l)**rd_o_cpd 

    ! Calculate hfx,qfx, and SST to keep buoyancy flux constant
    ! Could calculate moisture flux first, but should be OK either way
    kqfx = Vs * (qsfc_air - qv_l)   ! exchange coefficient for water vapor
    khfx = Beta0 * th_l/grav - 0.608_wp * th_l * kqfx
    tsk  = khfx / Vs + th_l

    ! Convert units
    shfx = - rho_l * cpd * khfx
    lhfx = - rho_l * alv * kqfx

  END SUBROUTINE cbl_stevens_fluxes

!-------------------------------------------------------------------------
!
! This subroutine creates a simple two valued field for the sensible heat flux
! and the water vapor flux.  The domain is simply divided in two with the 
! division determined by the longitude value given.  on each side of the 
! division the sensible and latent heat fluxes have different values.
!
  SUBROUTINE sfcflx_uniform(ptr_patch, shflux_sfc, sh_high, sh_low, qvflux_sfc,   &
    & qv_high, qv_low, wallLonDeg)
    TYPE(t_patch),TARGET,  INTENT(IN) :: ptr_patch  !< patch on which computation is performed
    REAL(wp), INTENT(in) :: wallLonDeg
    REAL(wp), INTENT(out):: shflux_sfc(:,:) ! sensible heat flux [W/m2]
    REAL(wp), INTENT(out):: qvflux_sfc(:,:) ! Water vapor flux at sfc [kg/kg]
    Real(wp), INTENT(in) :: sh_high   ! upper value of sensible heat flux
    Real(wp), INTENT(in) :: sh_low    ! lower value of sensible heat flux
    Real(wp), INTENT(in) :: qv_high   ! upper value of vapor flux
    Real(wp), INTENT(in) :: qv_low    ! lower value of vapor flux

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg

!-------------------------------------------------------------------------

    all_cells => ptr_patch%cells%ALL

    !Add horizontal variation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        lat_deg = ptr_patch%cells%center(jc,jb)%lat * rad2deg
        lon_deg = ptr_patch%cells%center(jc,jb)%lon * rad2deg

        IF((lon_deg) >= wallLonDeg) THEN
          shflux_sfc(jc,jb) = sh_high
          qvflux_sfc(jc,jb) = qv_high
        ELSE
          shflux_sfc(jc,jb) = sh_low
          qvflux_sfc(jc,jb) = qv_low
        ENDIF

      END DO
    END DO

  END SUBROUTINE sfcflx_uniform
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  ! read sounding from external file and then interpolate
  ! to model levels 
  ! Sounding file is assumed to be in this format:
  ! ps,no_vert_levels
  ! z(m) theta(k) qv(kg/kg) u(m/s) v(m/s)
  ! where the top row is close to surface/or surface and bottom row
  ! is near the top
  !!
  SUBROUTINE  read_ext_profile(z_in, theta_in, qv_in, u_in, v_in, psfc_in)
  
    REAL(wp),  INTENT(IN)  :: z_in(:)
    REAL(wp),  INTENT(OUT) :: theta_in(:)
    REAL(wp),  INTENT(OUT) :: qv_in(:)
    REAL(wp),  INTENT(OUT) :: u_in(:)
    REAL(wp),  INTENT(OUT) :: v_in(:)
    REAL(wp),  INTENT(OUT) :: psfc_in
  
    REAL(wp), ALLOCATABLE, DIMENSION(:):: zs, ths, qvs, us, vs
    CHARACTER(len=max_char_length),PARAMETER :: routine  = &
         &   'mo_nh_torus_exp:read_ext_profile'
  
    INTEGER :: ist, iunit
    INTEGER :: jk, klev 
  
    !-------------------------------------------------------------------------
  
    CALL message(TRIM(routine), 'READING FROM SOUNDING!')
    
    !open file again to read data this time
    iunit = find_next_free_unit(10,100)
    OPEN (unit=iunit,file='sound_in', action='READ', status='OLD', IOSTAT=ist) 
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), 'open verticaling sound file failed')
    ENDIF
  
    !Read the header : ps,klev
    READ (iunit,*,IOSTAT=ist)psfc_in,klev

    !now allocate
    ALLOCATE(zs(klev),ths(klev),us(klev),vs(klev),qvs(klev))
    zs = 0.0_wp; ths = 0._wp; us = 0._wp; vs = 0._wp; qvs = 0._wp

    DO jk = klev,1,-1 
      READ (iunit,*,IOSTAT=ist) zs(jk),ths(jk),qvs(jk),us(jk),vs(jk)
      IF(ist/=success)THEN
        CALL finish (TRIM(routine), 'reading souding file failed')
      ENDIF
    END DO

    CLOSE(iunit)

    !Check if the file is written in descending order
    IF(zs(1) < zs(klev)) &
         CALL finish (TRIM(routine), 'Writing souding data in descending order!')

    !Now perform interpolation to grid levels assuming:
    !a) linear interpolation
    !b) Beyond the last Z level the values are linearly extrapolated 
    !c) Assuming model grid is flat-NOT on sphere
  
    CALL vert_intp_linear_1d(zs,ths,z_in,theta_in)
    CALL vert_intp_linear_1d(zs,qvs,z_in,qv_in)
    CALL vert_intp_linear_1d(zs,us,z_in,u_in)
    CALL vert_intp_linear_1d(zs,vs,z_in,v_in)
  
    DEALLOCATE(zs, ths, qvs, us, vs)
  
  
  END SUBROUTINE  read_ext_profile

  !>
  !! Initialization of prognostic state vector for the warm bubble experiment
  !! on a torus. It is simplified form of Wesimann Klemp testcase from
  !! G. H. Bryan and J. M. Fritsch, "A benchmark simulation for moist
  !! nonhydrostatic numerical models", MWR, 2002
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_warm_bubble( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
    &                          ptr_int, ptr_metrics)

    TYPE(t_patch),TARGET,  INTENT(IN)   ::  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   ::  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT)::  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT)::  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   ::  ptr_metrics      
    TYPE(t_nh_ref),        INTENT(INOUT)::  ptr_nh_ref

    REAL(wp) :: z_exner_h(1:nproma,ptr_patch%nlev+1), z_help(1:nproma) 
    REAL(wp) :: ex_sfc, x_loc(3), x_c(3), psfc_in, dis, inv_th0, pres_new
    REAL(wp) :: th_new, qv_new, qc_new, th_old, th_ptb, temp_new
    REAL(wp), DIMENSION(ptr_patch%nlev) :: theta_in, qv_in, qc_in, tmp
    INTEGER  :: jc,jk,jb,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c
    INTEGER  :: nlev, nlevp1                  
    INTEGER  :: nlen, jg, itr

    REAL(wp), DIMENSION(3) :: x_bubble 
    CHARACTER(len=max_char_length),PARAMETER :: routine  = &
         &   'mo_nh_torus_exp:init_warm_bubble'
  !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !Note that this souding is created from a matlab code that 
    !iterates through the Eqs. 25,26,and 34 of Bryan and Fritsch's paper
    !given theta_e, qt, and surface pressure. The source code is 
    !in icon-aes-and/scripts/preprocessing/ named init.f90 and findzero.m
    !which creates a sound_** file that is then used to read in below
    !
    !One can also use the fortran source code from Bryan's website to generate
    !the initial condition
    !-------------------------------------------------------------------------
    !Read the sounding file
    CALL  read_ext_profile (ptr_metrics%z_mc(2,:,2), theta_in, qv_in, qc_in, tmp, &
                            psfc_in)
  
    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    !patch id
    jg = ptr_patch%id

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = psfc_in
    ex_sfc   = (psfc_in/p0ref)**rd_o_cpd
    les_config(jg)%psfc = psfc_in

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = qv_in(jk) 
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) = qc_in(jk)
      END DO

      DO jk = nlev, 1, -1
         z_help(1:nlen) = theta_in(jk)
         ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO

      !Get hydrostatic exner at the surface using surface pressure 
      z_exner_h(1:nlen,nlevp1) = ex_sfc
 
      !Get exner at full levels starting from exner at surface
      DO jk = nlev, 1, -1
         !exner at next half level after surface
         z_exner_h(1:nlen,jk) = z_exner_h(1:nlen,jk+1) - grav/cpd *     &
                                ptr_metrics%ddqz_z_full(1:nlen,jk,jb)/ &
                                ptr_nh_prog%theta_v(1:nlen,jk,jb)
        
         !exner at main levels
         ptr_nh_prog%exner(1:nlen,jk,jb) = 0.5_wp * &
                                     (z_exner_h(1:nlen,jk)+z_exner_h(1:nlen,jk+1))
      END DO

    ENDDO !jb

    !Zero wind
    ptr_nh_prog%vn = 0._wp
    ptr_nh_prog%w  = 0._wp

    !--------------------------------------------------------------------------
    !Add perturbation to theta_v and iterate to get balanced thermodynamic state
    !--------------------------------------------------------------------------

    !Bubble center, note that torus domain has center in the middle
    !Uses bubctr_lon and bubctr_lat to represent x,y in torus
    x_bubble = (/bubctr_x,bubctr_y,bubctr_z/)

    !First non-dimensionalize bubble ceter
    x_c(1) = x_bubble(1) / bub_hor_width
    x_c(2) = x_bubble(2) / bub_hor_width
    x_c(3) = x_bubble(3) / bub_ver_width


    !the th0 ratio
    inv_th0 = 1._wp / 300._wp

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jc = 1 , nlen
        DO jk = 1 , nlev
          x_loc(1) = ptr_patch%cells%cartesian_center(jc,jb)%x(1)/bub_hor_width
          x_loc(2) = ptr_patch%cells%cartesian_center(jc,jb)%x(2)/bub_hor_width
          x_loc(3) = ptr_metrics%z_mc(jc,jk,jb)/bub_ver_width
            
          IF(nh_test_name.eq.'2D_BUBBLE')THEN
           x_c(2)   = x_loc(2)
          END IF

          dis = plane_torus_distance(x_loc,x_c,ptr_patch%geometry_info)

          IF(dis < 1._wp)THEN
            th_ptb =  bub_amp * cos(pi_2*dis)**2 * inv_th0
            qv_new = qv_in(jk)
            qc_new = qc_in(jk)
            pres_new = p0ref*ptr_nh_prog%exner(jc,jk,jb)**(cpd/rd)
  
            DO itr = 1 , 20 
             th_new = ( th_ptb + 1._wp ) * ptr_nh_prog%theta_v(jc,jk,jb) / &
                       (1._wp + vtmpc1 * qv_new - qc_new)
             temp_new = th_new * ptr_nh_prog%exner(jc,jk,jb)
             qv_new   = spec_humi(sat_pres_water(temp_new),pres_new)
             qc_new   = qv_in(jk) + qc_in(jk) - qv_new

             IF(qc_new<0._wp)CALL finish(TRIM(routine), 'qc < 0')
            END DO

            !assign values to proper prog vars
            ptr_nh_prog%tracer(jc,jk,jb,iqv) = qv_new 
            ptr_nh_prog%tracer(jc,jk,jb,iqc) = qc_new
            ptr_nh_prog%theta_v(jc,jk,jb) = th_new * ( 1._wp + vtmpc1*qv_new - qc_new ) 

          END IF
            
        END DO
      END DO
    END DO 


    !calculate some of the variables
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jk = 1 , nlev
        ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO !jk
    ENDDO !jb


  END SUBROUTINE init_warm_bubble

END MODULE mo_nh_torus_exp
