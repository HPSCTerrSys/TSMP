!>
!! Defines the artificial testcases for the nonhydrostatic atmospheric model.
!! 
!! 
!! @par Revision History
!! Initial release by Almut Gassmann (2008-03-18)
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
MODULE mo_nh_testcases  
!-------------------------------------------------------------------------  
!  
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2006                       
!  
!-------------------------------------------------------------------------  
!  
!  
!  
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, finish, message_text
  USE mo_nh_testcases_nml
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, inwp, icosmo, iedmf
  USE mo_grid_config,          ONLY: lplane, n_dom, l_limited_area
  USE mo_model_domain,         ONLY: t_patch
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_math_constants,       ONLY: pi
  USE mo_math_utilities,       ONLY: gc2cc, t_cartesian_coordinates, &
                                   & t_geographical_coordinates,     &
                                   & arc_length
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: ltransport, iforcing, iqv
  USE mo_extpar_config,        ONLY: itopo
    
  USE mo_dynamics_config,      ONLY: nnow, nnew, lcoriolis
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_physical_constants,   ONLY: grav, cpd, rd, cvd_o_rd, p0ref

  USE mo_nonhydro_types,       ONLY: t_nh_state
  USE mo_nonhydro_state,       ONLY: duplicate_prog_state
  
  USE mo_intp_data_strc,       ONLY: t_int_state
  USE mo_nwp_lnd_types,        ONLY: t_lnd_state
  USE mo_lnd_nwp_config,       ONLY: isub_seaice
  USE mo_nh_pa_test,           ONLY: init_nh_state_prog_patest
  USE mo_nh_df_test,           ONLY: init_nh_state_prog_dftest
  USE mo_nh_hs_test,           ONLY: init_nh_state_prog_held_suarez
  USE mo_nh_jabw_exp,          ONLY: init_nh_topo_jabw, init_nh_state_prog_jabw,  & 
                                   & init_passive_tracers_nh_jabw, init_nh_inwp_tracers
  USE mo_nh_mrw_exp,           ONLY: init_nh_topo_mrw, init_nh_state_prog_mrw,    &
                                   & init_nh_prog_mwbr_const, mount_half_width                     
  USE mo_nh_wk_exp,            ONLY: init_nh_topo_wk, init_nh_env_wk,             &
                                   & init_nh_buble_wk                       
  USE mo_nh_bb13_exp,          ONLY: init_nh_env_bb13, init_nh_bubble_bb13                       
  USE mo_nh_dcmip_gw,          ONLY: init_nh_dcmip_gw, init_nh_gw_analyt
  USE mo_nh_dcmip_hadley,      ONLY: init_nh_dcmip_hadley         
  USE mo_nh_dcmip_schaer,      ONLY: init_nh_prog_dcmip_schaer,                   &
                                   & init_nh_topo_dcmip_schaer
  USE mo_nh_dcmip_rest_atm,   ONLY : init_nh_topo_dcmip_rest_atm,                 &
                                   & init_nh_prog_dcmip_rest_atm  
  USE mo_nh_dcmip_tc,          ONLY: init_nh_dcmip_tc
  USE mo_nh_dcmip_bw,          ONLY: init_nh_dcmip_bw
  USE mo_nh_dcmip_terminator,  ONLY: init_nh_dcmip_terminator
  USE mo_nh_lim_area_testcases,ONLY: init_nh_atmo_ana_nconstlayers,               &
                                   & init_nh_anaprof_uv, init_nh_topo_ana,        &
                                   & itype_atmo_ana, init_nh_atmo_ana_poly
  USE mo_nh_prog_util,         ONLY: nh_prog_add_random
  USE mo_random_util,          ONLY: add_random_noise_global
  USE mo_grid_geometry_info,   ONLY: planar_torus_geometry
  USE mo_nh_rce_exp,           ONLY: init_nh_state_rce_glb
  USE mo_nh_torus_exp,         ONLY: init_nh_state_cbl, init_nh_state_rico, &
                                     init_torus_with_sounding, init_warm_bubble
  USE mo_nh_tpe_exp,           ONLY: init_nh_state_prog_TPE
  
  IMPLICIT NONE  
  
  PRIVATE
  
  PUBLIC :: init_nh_testtopo,init_nh_testcase

  ! !DEFINED PARAMETERS for jablonowski williamson: 
  !  The rest of the needed parameters are define in mo_nh_jabw_exp
  REAL(wp), PARAMETER :: ps0      = 1.e5_wp     !< surface pressure (Pa)
  
  ! !DEFINED PARAMETERS for mountain induced Rossby wave train:
  REAL(wp), PARAMETER :: pres_sp  = 93000.0_wp  !< pressure surface at the south pole
  REAL(wp), PARAMETER :: temp_mrw = 288._wp     !< temperature of isothermal atmosphere

  ! !DEFINED PARAMETERS for APE
  REAL(wp), PARAMETER :: zp_ape   = 101325._wp  !< surface pressure
  REAL(wp), PARAMETER :: ztmc_ape = 25.006_wp   !< total moisture content 

  CONTAINS  
  
!-------------------------------------------------------------------------
!
!
  !>
  !! Initialize topography for nonhydrostatic artificial testcases.
  !! 
  !! Initialize topography for nonhydrostatic artificial testcases
  !! 
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-03-19)
  !! Modification by Daniel Reinert, DWD (2010-07-15)
  !! - moved initialization of topography into new subroutine 
  !!   init_nh_testtopo, which is called after the domain-decomposition.
  !!   (because of possible conflicts with the external-data type)
  !! 
  SUBROUTINE init_nh_testtopo (p_patch, ext_data)
!
! !INPUT VARIABLES:
  TYPE(t_patch),TARGET,  INTENT(INOUT) :: p_patch(n_dom)
  TYPE(t_external_data), INTENT(INOUT) :: ext_data(n_dom)

  INTEGER        :: jg, jc, jb, nlen
  INTEGER        :: nblks_c, npromz_c
  REAL(wp)       :: z_lon, z_lat, z_dist
  TYPE(t_geographical_coordinates) :: z_x2_geo
  TYPE(t_cartesian_coordinates)    :: z_x1_cart, z_x2_cart
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = '(mo_nh_testcases) init_nh_testtopo:'
 
  LOGICAL       :: l_modified

!-----------------------------------------------------------------------

  ! Initialize topography to zero if idealized topo is used
  IF ( itopo == 0 ) THEN
    DO jg = 1, n_dom
      ext_data(jg)%atm%topography_c(1:nproma,1:p_patch(jg)%nblks_c) = 0.0_wp
    ENDDO
  ENDIF

  SELECT CASE (nh_test_name)

  CASE ('zero', 'HS_jw')
    
    IF(nh_test_name=='HS_jw') CALL message(TRIM(routine),'running the Held-Suarez test')

  CASE ('schaer')
 
    !IF(.NOT.lplane) CALL finish(TRIM(routine),'Schaer test case only for lplane=True')

    ! At present the mountain is at position lat=0,lon=0 (given in meters)
    z_x2_geo%lon = 0.0_wp
    z_x2_geo%lat = 0.0_wp

    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_c
        IF (jb /=  p_patch(jg)%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
          IF(p_patch(jg)%geometry_info%geometry_type==planar_torus_geometry)THEN
            z_lon = p_patch(jg)%cells%cartesian_center(jc,jb)%x(1)
          ELSE
            z_lon = p_patch(jg)%cells%center(jc,jb)%lon*torus_domain_length/pi*0.5_wp
          END IF
          z_dist = z_lon-z_x2_geo%lon
          ext_data(jg)%atm%topography_c(jc,jb) = 250.0_wp &
          & * EXP(-(z_dist/5000.0_wp)**2)*((COS(pi*z_dist/4000.0_wp))**2)
        ENDDO
      ENDDO 
    ENDDO 

  CASE ('bell')

    ! At present the mountain is at position lat=0,lon=0 (given in meters)
    z_x2_geo%lon = 0.0_wp
    z_x2_geo%lat = 0.0_wp

    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_c
        IF (jb /=  p_patch(jg)%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
          IF (lplane) THEN
            z_lat  = 0.0_wp!p_patch(jg)%cells%center(jc,jb)%lat*torus_domain_length/pi*0.5_wp
            z_lon  = p_patch(jg)%cells%center(jc,jb)%lon*torus_domain_length/pi*0.5_wp
            z_dist = SQRT((z_lat-z_x2_geo%lat)**2+(z_lon-z_x2_geo%lon)**2)
          ELSE
            z_x1_cart    = gc2cc(p_patch(jg)%cells%center(jc,jb))
            z_x2_cart    = gc2cc(z_x2_geo)
            z_dist       = arc_length(z_x1_cart,z_x2_cart)
          ENDIF
          ext_data(jg)%atm%topography_c(jc,jb) = mount_height/ &
                 (1.0_wp+ (z_dist/mount_half_width)**2)**1.5_wp
        ENDDO
      ENDDO 
    ENDDO 

  CASE ('jabw', 'jabw_s')

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_c
     npromz_c  = p_patch(jg)%npromz_c

     CALL init_nh_topo_jabw ( p_patch(jg),ext_data(jg)%atm%topography_c, nblks_c, npromz_c)
    END DO

  CASE ('jabw_m')  

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_c
     npromz_c  = p_patch(jg)%npromz_c

     CALL init_nh_topo_jabw ( p_patch(jg),ext_data(jg)%atm%topography_c, nblks_c, npromz_c, &
                          & opt_m_height = mount_height, opt_m_half_width = mount_half_width )
    END DO

  CASE ('dcmip_bw_11')
    ! itopo == 0 --> The topography is initialized to 0 at the beginning of this subroutine
    CALL message(TRIM(routine),'running DCMIP2016 baroclinic wave test case 1-1')

  CASE ('mrw_nh', 'mrw2_nh' , 'mwbr_const')

   CALL message(TRIM(routine),'running mrw, setting topography')

   l_modified = .FALSE.
   IF (nh_test_name=='mrw2_nh'  ) THEN
     l_modified = .TRUE.
   ENDIF

   DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_c
     npromz_c  = p_patch(jg)%npromz_c

     CALL init_nh_topo_mrw ( p_patch(jg),ext_data(jg)%atm%topography_c, nblks_c, npromz_c, l_modified) 
   ENDDO

   CALL message(TRIM(routine),'topography is initialised ')

  CASE ('wk82')  

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_c
     npromz_c  = p_patch(jg)%npromz_c

     CALL init_nh_topo_wk ( p_patch(jg),ext_data(jg)%atm%topography_c, nblks_c, npromz_c)
    END DO

  CASE ('bb13')

    ! Testcase Baldauf, Brdar (2013) QJRMS (linear gravity/sound wave expansion in a channel)

    CALL message(TRIM(routine), "no orography for testcase bb13")

    ! no orography:
    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_c

        IF (jb /=  p_patch(jg)%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF

        DO jc = 1, nlen
          ext_data(jg)%atm%topography_c(jc,jb) = 0.0_wp
        ENDDO
      ENDDO 
    ENDDO 

  CASE ('PA')
   ! The topography ist initialized in "init_nh_state_prog_patest"
    CALL message(TRIM(routine),'running the Pure 3D-Advection test.')

  CASE ('DF1')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 1.')

  CASE ('DF2')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 2.')

  CASE ('DF3')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 3.')

  CASE ('DF4')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 4.')

  CASE ('HS_nh')
    ! The topography has been initialized to 0 
    CALL message(TRIM(routine),'running the Held-Suarez test')
    
  CASE ('APE_nwp')
    ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Aqua-Planet Experiment with non-hydrostatic atm. dynamics and NWP physics')
  
  CASE ('APE_echam')
    ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Aqua-Planet Experiment with non-hydrostatic atm. dynamics and ECHAM physics')
  
  CASE ('APE_nh')
    ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Aqua-Planet Experiment with non-hydrostatic atm. dynamics')

  CASE ('APEc_nh')
    ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running coupled Aqua-Planet Experiment with non-hydrostatic atm. dynamics')

  CASE ('TPEo')

    ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Terra-Planet Experiment with ECHAM physics')
    IF ( itopo == 0 ) THEN
      CALL message(TRIM(routine), 'using zero topography for TPEc experiment')
    END IF

  CASE ('TPEc')

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Terra-Planet Experiment with ECHAM physics')
    IF ( itopo == 0 ) THEN
      CALL message(TRIM(routine), 'using zero topography for TPEc experiment')
    END IF
  
  CASE ('g_lim_area')

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_c
     npromz_c  = p_patch(jg)%npromz_c

     CALL init_nh_topo_ana ( p_patch(jg), lplane, ext_data(jg)%atm%topography_c, nblks_c, npromz_c)

    END DO
    CALL message(TRIM(routine),'running g_lim_area')

  CASE ('dcmip_pa_12')

    ! The topography has been initialized to 0 
    CALL message(TRIM(routine),'running the dcmip_pa_12 (PA with Hadley-like circulation) test')

  CASE ('dcmip_gw_31')

    ! The topography has been initialized to 0 
    CALL message(TRIM(routine),'running the dcmip_gw_31 (small planet gravity wave) test')

  CASE ('dcmip_gw_32')

    ! The topography has been initialized to 0 
    CALL message(TRIM(routine),'running the dcmip_gw_32 (analyt. small planet gravity wave) test')

  CASE ('dcmip_rest_200')

    DO jg = 1, n_dom 

     CALL init_nh_topo_dcmip_rest_atm ( p_patch(jg), ext_data(jg)%atm%topography_c, ext_data(jg)%atm%fis  )

    END DO

    CALL message(TRIM(routine),'running the dcmip_rest_200 (steady state at rest dcmip) test')

!!$  CASE ('dcmip_mw_2x')
!DR topography_v no longer read in available. If needed, it should be 
!DR interpolated from topography_c.
!!$
!!$    DO jg = 1, n_dom 
!!$
!!$     CALL init_nh_topo_dcmip_schaer ( p_patch(jg),  ext_data(jg)%atm%topography_c,  &
!!$                          & ext_data(jg)%atm%topography_v, ext_data(jg)%atm%fis  )
!!$    END DO
!!$    CALL message(TRIM(routine),'running the dcmip_mw_2x (schaer-type dcmip) test')

  CASE ('dcmip_tc_51')
    ! itopo == 0 --> The topography is initialized to 0 at the begining of this subroutine
    CALL message(TRIM(routine),'running DCMIP tropical cyclone testcase 51')

  CASE ('dcmip_tc_52')
    ! itopo == 0 --> The topography is initialized to 0 at the begining of this subroutine
    CALL message(TRIM(routine),'running DCMIP tropical cyclone testcase 52')

  CASE ('CBL')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'CBL case is only for plane torus!')

   ! SBr,CHa
    DO jg = 1, n_dom
        ext_data(jg)%atm%alb_dif(:,:) = 0.2_wp !0.15_wp !0.2_wp
        ext_data(jg)%atm%albuv_dif(:,:) = 0.2_wp !0.15_wp !0.2_wp
        ext_data(jg)%atm%albni_dif(:,:) = 0.2_wp !0.15_wp !0.2_wp
        ext_data(jg)%atm%fr_land(:,:) = 1._wp
        ext_data(jg)%atm%fr_land_smt(:,:) = 1._wp
    END DO

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Convective Boundary Layer Experiment')

  !SBr,CHa: catchment-scale circulation case
  CASE ('CBL_ccs')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'CBL case is only for plane torus!')

   ! SBr,CHa
    DO jg = 1, n_dom
        ext_data(jg)%atm%alb_dif(:,:) = 0.2_wp !0.15_wp !0.2_wp
        ext_data(jg)%atm%albuv_dif(:,:) = 0.2_wp !0.15_wp !0.2_wp
        ext_data(jg)%atm%albni_dif(:,:) = 0.2_wp !0.15_wp !0.2_wp
        ext_data(jg)%atm%fr_land(:,:) = 1._wp
        ext_data(jg)%atm%fr_land_smt(:,:) = 1._wp
    END DO

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Convective Boundary Layer Experiment')

  CASE ('2D_BUBBLE', '3D_BUBBLE')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'2D BUBBLE case is only for plane torus!')

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running 2D Warm Bubble test case')

  CASE ('RCE_glb')

   ! Running Radiative Convective Equilibrium testcase
   CALL message(TRIM(routine),'running ICON in RCE on a global domain')

  CASE ('RICO')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'RICO case is only for plane torus!')

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Rain in the Culumus Over the Ocean LES Experiment')

  CASE ('RCE','GATE')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'To initialize with sounding is only for torus!')

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running LES with sounding')

  CASE DEFAULT

    CALL finish(routine,'wrong input for nh_test_name')

  END SELECT


  END SUBROUTINE init_nh_testtopo



!-------------------------------------------------------------------------
!
!
  !>
  !! Defines nonhydrostatic artificial initial conditions.
  !! 
  !! Initializes meteorological fields
  !! 
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-04-14)
  !! 
  SUBROUTINE init_nh_testcase (p_patch, p_nh_state, p_int, p_lnd_state, ext_data, ntl)
!
! !INPUT VARIABLES:
  TYPE(t_patch),TARGET,  INTENT(INOUT) :: p_patch(n_dom)
  TYPE(t_int_state),     INTENT(   IN) :: p_int(n_dom)
  TYPE(t_lnd_state),     INTENT(INOUT) :: p_lnd_state(n_dom)
  TYPE(t_external_data), INTENT(INOUT) :: ext_data(n_dom)
  INTEGER :: ntl
  TYPE(t_nh_state), TARGET, INTENT(INOUT):: p_nh_state(n_dom)

  INTEGER        :: jg, je, jc, jb, jk, jt,   &
                    nlen, nblks_e, npromz_e,  nblks_c, npromz_c
  INTEGER        :: nlev, nlevp1        !< number of full and half levels

  TYPE(t_nh_state), POINTER       :: p_nhdom
                            
  REAL(wp)              :: p_sfc_jabw  ! surface pressure for the jabw test case, 
                                       ! standard values is 100000 Pa   
  REAL(wp)              :: global_moist
  REAL(wp) :: z_help

  LOGICAL  :: l_hydro_adjust, l_moist
  
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine =  &
                                   '(mo_nh_testcases) init_nh_testcase:' 

!-----------------------------------------------------------------------

  SELECT CASE (nh_test_name)

  CASE ('jabw', 'jabw_s', 'jabw_m', 'HS_jw') ! Jablonowski test

  CALL message(TRIM(routine),'Jablonowski test')
  IF ( iforcing == inwp ) THEN
    CALL message(TRIM(routine),' iforcing == inwp')
  ELSE
    CALL message(TRIM(routine),'Attention: iforcing /= inwp')
  ENDIF

  IF (nh_test_name == "jabw_s" .OR. nh_test_name == "jabw_m") THEN
   jw_up = 0.0_wp
  END IF  
   p_sfc_jabw = ps0

  DO jg = 1, n_dom

    CALL   init_nh_state_prog_jabw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                   & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                   & p_int(jg),                                   &
                                   & p_sfc_jabw,jw_up )


    IF ( ltransport .AND. iforcing /= inwp ) THEN   ! passive tracers

       CALL init_passive_tracers_nh_jabw (p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                         & rotate_axis_deg, tracer_inidist_list, p_sfc_jabw) 

    END IF

    IF ( ltransport .AND. iforcing == inwp ) THEN 
     IF ( atm_phy_nwp_config(jg)%inwp_gscp /= 0 .OR.&
                   &                 atm_phy_nwp_config(jg)%inwp_convection /= 0  ) THEN   !


      CALL init_nh_inwp_tracers ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                & rh_at_1000hpa, qv_max, l_rediag=.TRUE. )




     ELSE

        p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,:) = 0.0_wp   

     END IF

    END IF

    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

  ENDDO !jg

  CALL message(TRIM(routine),'End setup Jablonowski test')


  CASE ('dcmip_bw_11')

    DO jg = 1, n_dom
      CALL init_nh_dcmip_bw (p_patch(jg),                   &
        &                    p_nh_state(jg)%prog(nnow(jg)), &
        &                    p_nh_state(jg)%diag,           &
        &                    p_int(jg),                     &
        &                    p_nh_state(jg)%metrics         )
    END DO  ! jg
    CALL message(TRIM(routine),'End setup baroclinic wave (dcmip_bw_11) test')


  CASE ('mrw_nh', 'mrw2_nh')

   CALL message(TRIM(routine),'MRW test')
  
   l_hydro_adjust = .TRUE.
   l_moist = .FALSE.

   DO jg = 1, n_dom

     IF ( iforcing == inwp ) THEN
       CALL message(TRIM(routine),' iforcing == inwp')     
       IF ( atm_phy_nwp_config(jg)%inwp_gscp /= 0 .OR.&
                      &                 atm_phy_nwp_config(jg)%inwp_convection /= 0  ) THEN 
         l_moist = .TRUE.
       END IF
     ENDIF


     IF (.NOT. l_moist) THEN

       CALL   init_nh_state_prog_mrw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                        &
                                     & ext_data(jg)%atm%topography_c,              &
                                     & p_nh_state(jg)%metrics,                     &
                                     & p_int(jg), l_hydro_adjust, iforcing, l_moist)

    ELSE

       CALL   init_nh_state_prog_mrw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                        &
                                     & ext_data(jg)%atm%topography_c,              &
                                     & p_nh_state(jg)%metrics,                     &
                                     & p_int(jg), l_hydro_adjust, iforcing, l_moist , &
                                     & opt_rh_at_1000hpa= rh_at_1000hpa,           &
                                     & opt_qv_max=qv_max                           )
    END IF
   
    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))  
 
   ENDDO !jg

   CALL message(TRIM(routine),'End setup MRW test')


  CASE ('mwbr_const')

   CALL message(TRIM(routine),'mwbr_const test case')
  
   l_hydro_adjust = .TRUE.
   l_moist = .FALSE.

   DO jg = 1, n_dom

     IF ( iforcing == inwp ) THEN
       CALL message(TRIM(routine),' iforcing == inwp')     
       IF ( atm_phy_nwp_config(jg)%inwp_gscp /= 0 .OR.&
                      &                 atm_phy_nwp_config(jg)%inwp_convection /= 0  ) THEN 
         l_moist = .TRUE.
       END IF
     ENDIF 

     IF (.NOT. l_moist) THEN

       CALL   init_nh_prog_mwbr_const ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                      & p_nh_state(jg)%diag,                        &
                                      & ext_data(jg)%atm%topography_c,              &
                                      & p_nh_state(jg)%metrics,                     &
                                      & p_int(jg), l_hydro_adjust, iforcing, l_moist)

     ELSE
  
       CALL   init_nh_prog_mwbr_const ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                      & p_nh_state(jg)%diag,                        &
                                      & ext_data(jg)%atm%topography_c,              &
                                      & p_nh_state(jg)%metrics,                     &
                                      & p_int(jg), l_hydro_adjust, iforcing, l_moist,&
                                      & opt_rh_at_1000hpa= rh_at_1000hpa,           &
                                      & opt_qv_max=qv_max                           )

    END IF

    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg))) 
 
   ENDDO !jg

   CALL message(TRIM(routine),'End setup mwbr_const test')


  CASE ('zero','bell','schaer')

  ! For the moment we think of a given Brunt Vaisala frequency and a given      
  ! zonal wind. The lplane and the lcorio=F options are assumed

  DO jg = 1, n_dom
    p_nhdom   => p_nh_state(jg)
    nblks_e   = p_patch(jg)%nblks_e
    npromz_e  = p_patch(jg)%npromz_e

    ! number of vertical levels
    nlev   = p_patch(jg)%nlev
    nlevp1 = p_patch(jg)%nlevp1

    DO jt = 1, ntl 
      ! normal wind
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
           nlen = nproma
        ELSE
           nlen = npromz_e
        ENDIF
        DO jk = 1, nlev
          DO je = 1, nlen
            p_nh_state(jg)%prog(jt)%vn(je,jk,jb) = nh_u0 &
            !(p_nhdom%metrics%geopot(1,jk,1)/grav/1000.0_wp+5.0_wp)& !shear
            *p_patch(jg)%edges%primal_normal(je,jb)%v1

            ! copy vn to reference state vector (needed by Rayleigh damping mechanism)
            p_nh_state(jg)%ref%vn_ref(je,jk,jb)  &
                      = p_nh_state(jg)%prog(jt)%vn(je,jk,jb)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    nblks_c   = p_patch(jg)%nblks_c
    npromz_c  = p_patch(jg)%npromz_c
    ! scalars (all is dry!)
    DO jt = 1, ntl 
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF
        DO jk = 1, nlev
          DO jc = 1, nlen
            z_help=(nh_brunt_vais/grav)**2*p_nhdom%metrics%geopot(jc,jk,jb)
            ! profile of theta is explicitly given
            p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb) = nh_t0*EXP(z_help)
          ENDDO
        ENDDO
        DO jk = nlev, 1, -1
          DO jc = 1, nlen
            IF (jk == nlev) THEN
              IF (nh_brunt_vais /= 0.0_wp) THEN
                ! for exner(nlev) (lowermost level): analytical exner
                z_help = (nh_brunt_vais/grav)**2*p_nhdom%metrics%geopot(jc,jk,jb)
                p_nh_state(jg)%prog(jt)%exner(jc,jk,jb) =    &
                  (grav/nh_brunt_vais)**2/nh_t0/cpd*(EXP(-z_help)-1.0_wp)+1.0_wp
              ELSE
                p_nh_state(jg)%prog(jt)%exner(jc,jk,jb) = &
                  1.0_wp-p_nhdom%metrics%geopot(jc,jk,jb)/cpd/nh_t0
              ENDIF
            ELSE ! other levels are hydrostatically balanced with respect to model numerics
              z_help=0.5_wp*(p_nh_state(jg)%prog(jt)%theta_v(jc,jk  ,jb) &
                           + p_nh_state(jg)%prog(jt)%theta_v(jc,jk+1,jb))
              p_nh_state(jg)%prog(jt)%exner(jc,jk,jb) = &
              & p_nh_state(jg)%prog(jt)%exner(jc,jk+1,jb) &
              & -grav/cpd*p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb)/z_help
            ENDIF
          ENDDO
        ENDDO
        DO jk = 1, nlev
          DO jc = 1, nlen

       !     ! perturbation in theta_v for gravity test case
       !     p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)=&
       !              p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)+ 0.01_wp*sin(&
       !              pi*p_nhdom%metrics%geopot(jc,jk,jb)/grav/10000.0_wp)&
       !            /(1.0_wp+(p_patch(jg)%cells%center(jc,jb)%lon*30.0/pi)**2)
       !            !Das ist fuer dx=500 mit 600 Punkten (auf 2pi verteilt)

       !     ! perturbation in theta_v for Straka test case
       !     z_help = SQRT( ( p_patch(jg)%cells%center(jc,jb)%lon*6.4_wp/pi)**2 &
       !     &      +((p_nhdom%metrics%z_mc(jc,jk,jb)-3000.0_wp)/2000.0_wp)**2)
       !     IF (z_help<=1.0_wp) THEN
       !       p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)=&
       !       &   p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb) &
       !       &   -15.0_wp*(COS(pi*z_help)+1.0_wp)*0.5_wp&
       !       &   /p_nh_state(jg)%prog(jt)%exner(jc,jk,jb)
       !     ENDIF

            ! exner and theta_v are given, so rho is deduced...
            p_nh_state(jg)%prog(jt)%rho(jc,jk,jb) = &
            &        (p_nh_state(jg)%prog(jt)%exner(jc,jk,jb)**cvd_o_rd)*p0ref/rd &
            &       /p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)

          ENDDO
        ENDDO
        DO jk = 1, nlevp1
          p_nh_state(jg)%prog(jt)%w(1:nlen,jk,jb) = 0.0_wp

          ! copy w to reference state vector (needed by Rayleigh damping mechanism)
          p_nh_state(jg)%ref%w_ref(je,jk,jb) = &
              p_nh_state(jg)%prog(jt)%w(je,jk,jb)
        ENDDO
      ENDDO
    ENDDO
  ENDDO


  CASE ('PA')  ! pure advection test case, no mountain

    DO jg = 1, n_dom
      CALL init_nh_state_prog_patest(p_patch(jg),p_int(jg),             &
        &                       p_nh_state(jg)%prog(nnow(jg)),          &
        &                       p_nh_state(jg)%diag,ext_data(jg),       &
        &                       p_nh_state(jg)%metrics,rotate_axis_deg, &
        &                       linit_tracer_fv, tracer_inidist_list )

      CALL init_nh_state_prog_patest(p_patch(jg),p_int(jg),             &
        &                       p_nh_state(jg)%prog(nnew(jg)),          &
        &                       p_nh_state(jg)%diag,ext_data(jg),       &
        &                       p_nh_state(jg)%metrics,rotate_axis_deg, &
        &                       linit_tracer_fv, tracer_inidist_list )

    ENDDO !jg


  CASE ('DF1', 'DF2', 'DF3', 'DF4')  ! 2D deformational flow test case, no mountain

    DO jg = 1, n_dom

      CALL init_nh_state_prog_dftest(p_patch(jg),                         &
           &                         p_nh_state(jg)%prog(nnow(jg)),       &
           &                         p_nh_state(jg)%diag,                 &
           &                         p_int(jg), ext_data(jg),             &
           &                         rotate_axis_deg, nh_test_name,       &
           &                         linit_tracer_fv, tracer_inidist_list )

      CALL init_nh_state_prog_dftest(p_patch(jg),                         &
           &                         p_nh_state(jg)%prog(nnew(jg)),       &
           &                         p_nh_state(jg)%diag,                 &
           &                         p_int(jg), ext_data(jg),             &
           &                         rotate_axis_deg, nh_test_name,       &
           &                         linit_tracer_fv, tracer_inidist_list )

    ENDDO !jg


  CASE ('HS_nh')  ! Held-Suarez test case, no mountain, isothermal atmosphere

    DO jg = 1, n_dom

      ! number of vertical levels
      nlev   = p_patch(jg)%nlev

      CALL init_nh_state_prog_held_suarez(p_patch(jg),p_nh_state(jg)%prog(nnow(jg)),  &
        &                            p_nh_state(jg)%diag,ext_data(jg),           &
        &                            p_nh_state(jg)%metrics)

      IF (lhs_nh_vn_ptb) THEN
          CALL nh_prog_add_random( p_patch(jg),           & ! input
               & p_nh_state(jg)%prog(nnow(jg))%vn(:,:,:), & ! in and out
               & "edge", hs_nh_vn_ptb_scale, 1, nlev ) ! input

          CALL message(TRIM(routine),'Initial state used in the &
               & Held-Suarez test: random noised added to the normal wind')
      END IF

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

    ENDDO !jg


  CASE ('APE_nwp', 'APE_echam', 'APE_nh', 'APEc_nh')  ! Aqua-Planet Experiment, no mountain

    p_sfc_jabw   = zp_ape          ! Pa
    global_moist = ztmc_ape        ! kg/m**2 total moisture content
    jw_up = 1._wp

    DO jg = 1, n_dom
    
      CALL   init_nh_state_prog_jabw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                     & p_int(jg),                                   &
                                     & p_sfc_jabw,jw_up )
    
      IF ( ltransport ) THEN   !
    
        CALL init_nh_inwp_tracers ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                  & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                  & rh_at_1000hpa, qv_max, l_rediag=.TRUE.,  &
                                  & opt_global_moist=global_moist)
    
      END IF
    
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    
    ENDDO !jg

    CALL message(TRIM(routine),'End setup non-hydrostatic APE test (APE_nwp, APE_echam, APE_nh, APEc_nh)')

  CASE ('TPEc', 'TPEo')  ! Terra-Planet Experiment

    jw_up = 1._wp

    DO jg = 1, n_dom
    
      CALL init_nh_state_prog_TPE(p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%diag, &
                                  ext_data(jg), p_nh_state(jg)%metrics,                            &
                                  rh_at_1000hpa, qv_max, tpe_moist, tpe_psfc, tpe_temp)

      ! why do we call this?   
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    
    ENDDO !jg

    CALL message(TRIM(routine),'End setup TPEc test')

  CASE ('wk82')

   CALL message(TRIM(routine),'wk82 test')
  
   l_hydro_adjust = .TRUE.

   DO jg = 1, n_dom

         ! initialize environment atmosphere
  
    CALL   init_nh_env_wk ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                 &
                                     & p_nh_state(jg)%metrics,              &
                                     & p_int(jg), l_hydro_adjust )
         ! add perturbation to theta and recalculate theta_v and rho 
    CALL init_nh_buble_wk ( p_patch(jg),p_nh_state(jg)%metrics,            &
                                     & p_nh_state(jg)%prog(nnow(jg)),      &
                                     & p_nh_state(jg)%diag )
         

    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

   ENDDO !jg

   CALL message(TRIM(routine),'End setup wk82 test')

  
  CASE ('bb13')

    CALL message(TRIM(routine), 'Baldauf, Brdar (2013) QJRMS test (linear gravity/sound waves in a channel)')
  
    l_hydro_adjust = .TRUE.

    DO jg = 1, n_dom

         ! initialize environment atmosphere
      CALL init_nh_env_bb13   ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                 &
                                     & p_nh_state(jg)%metrics,              &
                                     & p_int(jg), l_hydro_adjust  )
         ! add perturbation to theta and recalculate theta_v and rho 
      CALL init_nh_bubble_bb13 ( p_patch(jg), p_nh_state(jg)%metrics,       &
                                     & p_nh_state(jg)%prog(nnow(jg)),       &
                                     & p_nh_state(jg)%diag )
         

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

    ENDDO !jg

    CALL message(TRIM(routine),'End setup Baldauf, Brdar (2013) test')

  
  CASE ('g_lim_area')

   CALL message(TRIM(routine),'g_lim_area test')
  
   l_hydro_adjust = .TRUE.

   IF (.NOT. l_limited_area  .OR. lcoriolis) THEN

     WRITE(message_text,'(a)') &
             & 'For g_lim_area test case l_limited_area must &
             & be .TRUE. and lcoriolis must be .FALSE.'
            CALL finish  (routine, TRIM(message_text))
   END IF

   DO jg = 1, n_dom

         ! initialize environment atmosphere
    SELECT CASE (itype_atmo_ana)

    CASE(1)
    
      CALL  init_nh_atmo_ana_nconstlayers( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                 &
                                     & p_nh_state(jg)%metrics,              &
                                     & p_int(jg),l_hydro_adjust  )
    CASE(2) 
      CALL  init_nh_atmo_ana_poly( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                 &
                                     & p_nh_state(jg)%metrics,              &
                                     & p_int(jg),l_hydro_adjust  )

    CASE default
      WRITE(message_text,'(a)') &
             & 'You should define a valid option for    &
             & itype_atmo_ana for the &
             & g_lim_area test case'
            CALL finish  (routine, TRIM(message_text))

    END SELECT
 
         ! initialize wind
    CALL   init_nh_anaprof_uv  ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg))%vn, &
                               & p_nh_state(jg)%prog(nnow(jg))%w,               &
                               & p_nh_state(jg)%metrics, p_int(jg) ) 
         
    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

   ENDDO !jg

   CALL message(TRIM(routine),'End setup g_lim_area test')


  CASE ('dcmip_pa_12')

    CALL message(TRIM(routine),'setup dcmip_pa_12 (PA with Hadley-like circulation) test')

    DO jg = 1, n_dom
      CALL init_nh_dcmip_hadley( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)),  &
        &                        p_nh_state(jg)%diag, p_int(jg),              &
        &                        p_nh_state(jg)%metrics )

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    ENDDO

    CALL message(TRIM(routine),'End setup dcmip_pa_12 test')


  CASE ('dcmip_gw_31')

    CALL message(TRIM(routine),'setup dcmip_gw_31 (gravity waves on small planet) test')

    DO jg = 1, n_dom
      CALL init_nh_dcmip_gw( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
        &                    p_nh_state(jg)%diag, p_nh_state(jg)%metrics )
    ENDDO

    CALL message(TRIM(routine),'End setup dcmip_gw_31 test')


  CASE ('dcmip_gw_32')

    CALL message(TRIM(routine),'setup dcmip_gw_32 (analyt. gravity waves on small planet) test')

    DO jg = 1, n_dom
      CALL init_nh_gw_analyt( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)),           &
        &                    p_nh_state(jg)%diag, p_nh_state(jg)%metrics, p_int(jg) )
    ENDDO

    CALL message(TRIM(routine),'End setup dcmip_gw_32 test')


  CASE ('dcmip_rest_200')

    CALL message(TRIM(routine),'setup dcmip_rest_200 (steady state at rest) test')
  
    l_hydro_adjust = .TRUE.

    IF ( lcoriolis) THEN

      WRITE(message_text,'(a)') &
             & 'For dcmip_rest_200 test case  lcoriolis must be .FALSE.'
            CALL finish  (routine, TRIM(message_text))
    END IF

    DO jg = 1, n_dom
      CALL init_nh_prog_dcmip_rest_atm( p_patch(jg),                              &
        &                    p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%diag,  &
        &                    p_nh_state(jg)%metrics, p_int(jg),l_hydro_adjust )
    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    ENDDO

    CALL message(TRIM(routine),'End setup dcmip_rest_200 test')


  CASE ('dcmip_mw_2x')

    CALL message(TRIM(routine),'setup dcmip_mw_2x (schaer-type on small planet) test')
  
    l_hydro_adjust = .FALSE.

    IF ( lcoriolis) THEN

      WRITE(message_text,'(a)') &
             & 'For dcmip_mw_2x test case  lcoriolis must be .FALSE.'
            CALL finish  (routine, TRIM(message_text))
    END IF

    DO jg = 1, n_dom
      CALL init_nh_prog_dcmip_schaer( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
        &                    p_nh_state(jg)%diag, p_nh_state(jg)%ref,             &
        &                    p_nh_state(jg)%metrics, p_int(jg),l_hydro_adjust )
    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    ENDDO

    CALL message(TRIM(routine),'End setup dcmip_mw_2x test')


  CASE ('dcmip_tc_51','dcmip_tc_52')

    ! 'dcmip_tc_51' and 'dcmip_tc_52' have the same initial state.

    DO jg = 1, n_dom

      CALL init_nh_dcmip_tc ( p_patch(jg),                   &
        &                     p_nh_state(jg)%prog(nnow(jg)), &
        &                     p_nh_state(jg)%diag,           &
        &                     p_nh_state(jg)%metrics,        &
        &                     p_int(jg)                     )

      CALL duplicate_prog_state( p_nh_state(jg)%prog(nnow(jg)), &
        &                        p_nh_state(jg)%prog(nnew(jg)) )

    END DO !jg

    CALL message(TRIM(routine),'End setup dcmip_tc_51/52')


  CASE ('CBL')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'CBL case is only for plane torus!')

    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev
      CALL init_nh_state_cbl ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%ref,  &
                      & p_nh_state(jg)%diag, p_int(jg), p_nh_state(jg)%metrics )
 
      call add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%w(:,:,:),         &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=w_perturb )   

      call add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=th_perturb )   

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    END DO !jg

    CALL message(TRIM(routine),'End setup CBL test')

  !SBr,CHa: catchment-scale circulation case
  CASE ('CBL_ccs')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'CBL case is only for plane torus!')

    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev
      CALL init_nh_state_cbl ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%ref,  &
                      & p_nh_state(jg)%diag, p_int(jg), p_nh_state(jg)%metrics )
 
      call add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%w(:,:,:),         &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=w_perturb )   

      call add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=th_perturb )   

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    END DO !jg

    CALL message(TRIM(routine),'End setup CBL test')

  CASE ('RCE_glb')

    ! u,v,w are initialized to zero.  exner and rho are similar/identical to CBL
    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev
      CALL init_nh_state_rce_glb ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%ref,  &
                      & p_nh_state(jg)%diag, p_int(jg), p_nh_state(jg)%metrics )

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=th_perturb )   

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
!
      CALL message(TRIM(routine),'End setup global RCE test')
    END DO !jg

  CASE ('RICO')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'RICO case is only for plane torus!')

    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev
      CALL init_nh_state_rico ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%ref,  &
                      & p_nh_state(jg)%diag, p_int(jg), p_nh_state(jg)%metrics )

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%w(:,:,:),         &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=w_perturb )   

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=th_perturb )   
 
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    END DO !jg

    CALL message(TRIM(routine),'End setup RICO test')

  CASE ('RCE','GATE') !to initialize with sounding

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'To initizialize with sounding is only for torus!')

    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev

      CALL init_torus_with_sounding ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                 p_nh_state(jg)%ref, p_nh_state(jg)%diag, p_int(jg), p_nh_state(jg)%metrics )
 
      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%w(:,:,:),         &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=w_perturb )   

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=th_perturb )   
                    
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    END DO !jg

    CALL message(TRIM(routine),'End init with sounding')

  CASE ('2D_BUBBLE', '3D_BUBBLE') !to initialize with sounding

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'2D warm bubble case is only for torus!')

    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev

      CALL init_warm_bubble ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                 p_nh_state(jg)%ref, p_nh_state(jg)%diag, p_int(jg), p_nh_state(jg)%metrics )

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    END DO !jg

    CALL message(TRIM(routine),'End initilization of 2D warm bubble')
  END SELECT

 END SUBROUTINE init_nh_testcase


 
END MODULE mo_nh_testcases
