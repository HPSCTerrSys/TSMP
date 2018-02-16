!!>
!! mo_surface_les
!!
!! Surface calculations for les physics using Businger Dyer relationship
!!
!! @author Anurag Dipankar, MPI-M
!!
!! @par Revision History
!! Initial release by Anurag Dipankar, MPI-M (2013-03-07)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_surface_les

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text
  USE mo_nonhydro_types,      ONLY: t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: msg_level
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_impl_constants    ,  ONLY: min_rlcell_int
  USE mo_sync,                ONLY: SYNC_C, sync_patch_array, global_sum_array
  USE mo_physical_constants,  ONLY: cpd, p0ref, grav, alv, rd, rgrav, rd_o_cpd, vtmpc1
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag 
  USE mo_satad,               ONLY: spec_humi, sat_pres_water
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_les_config,          ONLY: les_config
  USE mo_math_constants,      ONLY: pi_2, ln2
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_data_turbdiff,       ONLY: akt, alpha0
  USE mo_turbdiff_config,     ONLY: turbdiff_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: surface_conditions

  !Parameters for surface layer parameterizations: From Zeng_etal 1997 J. Clim
  REAL(wp), PARAMETER :: bsm = 5.0_wp  !Businger Stable Momentum
  REAL(wp), PARAMETER :: bum = 16._wp  !Businger Untable Momentum
  REAL(wp), PARAMETER :: bsh = 5.0_wp  !Businger Stable Heat
  REAL(wp), PARAMETER :: buh = 16._wp  !Businger Untable Heat

  !Parameters for surface parameterizations from COSMO docs
  REAL(wp), PARAMETER :: beta_10 = 0.042_wp
  REAL(wp), PARAMETER :: h_10    = 10._wp
  REAL(wp), PARAMETER :: zh_max  = 0.1_wp

  !Parameters for RICO case
  REAL(wp), PARAMETER :: c_m = 0.001229_wp  
  REAL(wp), PARAMETER :: c_h = 0.001094_wp
  REAL(wp), PARAMETER :: c_q = 0.001133_wp
  REAL(wp), PARAMETER :: th0_rico = 298.5_wp
  REAL(wp), PARAMETER :: psfc = 101540._wp

  CHARACTER(len=12)  :: str_module = 'surface_les'  ! Output of module for 1 line debug
  INTEGER            :: idt_src    = 4           ! Determines level of detail for 1 line debug

  CONTAINS


  !>
  !! surface_conditions
  !!------------------------------------------------------------------------
  !! Calculate surface temperature and moisture given the fluxes using Businger
  !! Dyer relationships .OR. vice versa. All calculations are done at cell center
  !!  
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-06)
  SUBROUTINE  surface_conditions(p_nh_metrics, p_patch, p_nh_diag, p_int, &
                                 p_prog_lnd_now, p_prog_lnd_new, p_diag_lnd, &
                                 prm_diag, theta, qv)

    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch    !< single patch
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_int_state), INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_lnd_prog),  INTENT(in)        :: p_prog_lnd_now!<land prog state 
    TYPE(t_lnd_prog),  INTENT(inout)     :: p_prog_lnd_new!<land prog state 
    TYPE(t_lnd_diag),  INTENT(inout)     :: p_diag_lnd    !<land diag state 
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars
    REAL(wp),          INTENT(in)        :: theta(:,:,:)  !pot temp  
    REAL(wp),          INTENT(in)        :: qv(:,:,:)     !spec humidity

    REAL(wp) :: rhos, obukhov_length, z_mc, ustar, mwind
    REAL(wp) :: zrough, exner, var(nproma,p_patch%nblks_c), theta_nlev, qv_nlev
    REAL(wp) :: theta_sfc, shfl, lhfl, umfl, vmfl, bflx1, bflx2, theta_sfc1, diff
    REAL(wp) :: RIB, zh, tcn_mom, tcn_heat, t_sfc, ex_sfc, inv_bus_mom
    REAL(wp) :: ustar_mean
    REAL(wp) :: pres_sfc(nproma,p_patch%nblks_c)
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc, isidx, isblk, rl
    INTEGER :: nlev, jg, itr, jkp1
    
    CHARACTER(len=*), PARAMETER :: routine = 'mo_surface_les:surface_conditions'

    IF (msg_level >= 15) &
         CALL message(TRIM(routine), '')

    jg = p_patch%id

    IF(les_config(jg)%psfc < 0._wp)THEN
      !use pressure from dynamics
      pres_sfc = p_nh_diag%pres_sfc
    ELSE
      !use imposed pressure 
      pres_sfc = les_config(jg)%psfc
    END IF 

    ! number of vertical levels
    nlev = p_patch%nlev
    jk   = nlev
    jkp1 = jk+1

    i_nchdom   = MAX(1,p_patch%n_childdom)
 
    !loop indices: exclude halo 
    rl_start   = grf_bdywidth_c+1 
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)         


    SELECT CASE(les_config(jg)%isrfc_type)

    !Default case with TERRA
    !Most of what follows has been taken from mo_nwp_turbtrans_interface 
    !for COSMO turbulence

    CASE(1)

     !For now the LES scheme uses the surface fluxes from land directly. Exchange coefficients
     !are calculated by calling turbtrans using GME approach. See comment on mo_interface_les
     !where it is called. Status as on 11.09.2013 (AD)

     !change sign of surface fluxes
     prm_diag%shfl_s  = - prm_diag%shfl_s
     prm_diag%lhfl_s  = - prm_diag%lhfl_s
     prm_diag%umfl_s  = - prm_diag%umfl_s
     prm_diag%vmfl_s  = - prm_diag%vmfl_s
 
    !Prescribed latent/sensible heat fluxes: get ustar and surface temperature / moisture
    !Ideally, one should do an iteration to get ustar corresponding to given fluxes

    !Fixed surface pressure to be comparable to incompressible LES models. Pseudo density
    !(constant in time) that is used here is fixed to the initial value so that the flux 
    !is fixed in density units. 
    CASE(2)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,exner,zrough,mwind,z_mc,RIB,rhos, &
!$OMP            ustar,obukhov_length,theta_sfc),ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx

            !Surface exner
            exner = EXP( rd_o_cpd*LOG(pres_sfc(jc,jb)/p0ref) )

            !Roughness length
            zrough = prm_diag%gz0(jc,jb) * rgrav

            !Get reference surface pot. temperature
            !First time step t_g takes value assigned in nwp_phy_init          
            theta_sfc = p_prog_lnd_now%t_g(jc,jb) / exner

            !Mean wind at nlev
            mwind  = MAX( les_config(jg)%min_sfc_wind,SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) )
           
            !Z height
            z_mc   = p_nh_metrics%z_mc(jc,jk,jb) - p_nh_metrics%z_ifc(jc,jkp1,jb)

            !Now diagnose friction velocity (ustar)
            IF(les_config(jg)%ufric<0._wp)THEN
              !Bulk Richardson no at first model level
              RIB = grav * (theta(jc,jk,jb)-theta_sfc)*(z_mc-zrough)/(theta_sfc*mwind**2)
              !first guess
              ustar = SQRT( diag_ustar_sq(z_mc,zrough,RIB,mwind) )
              DO itr = 1 , 5
                !"-" sign in the begining because ustar*thstar = -shflx
                obukhov_length = - theta_sfc*ustar**3/(akt*grav*les_config(jg)%shflx)
                ustar = mwind / businger_mom(zrough,z_mc,obukhov_length)  
              END DO

            ELSE
              ustar = les_config(jg)%ufric
              obukhov_length = -theta_sfc*ustar**3/(akt*grav*les_config(jg)%shflx)
            END IF

            !New temperature
            theta_sfc   = theta(jc,jk,jb) + les_config(jg)%shflx / ustar * &
                          businger_heat(zrough,z_mc,obukhov_length) 

            p_prog_lnd_new%t_g(jc,jb) = theta_sfc * exner

            !Get surface qv
            p_diag_lnd%qv_s(jc,jb) = qv(jc,jk,jb) + les_config(jg)%lhflx / ustar * &
                                     businger_heat(zrough,z_mc,obukhov_length) 

            !rho at surface: no qc at suface
            rhos   =  pres_sfc(jc,jb)/( rd * &
                     p_prog_lnd_new%t_g(jc,jb)*(1._wp+vtmpc1*p_diag_lnd%qv_s(jc,jb)) )  

            !Get surface fluxes
            prm_diag%shfl_s(jc,jb)  = les_config(jg)%shflx * rhos * cpd
            prm_diag%lhfl_s(jc,jb)  = les_config(jg)%lhflx * rhos * alv
            prm_diag%umfl_s(jc,jb)  = ustar**2 * rhos * p_nh_diag%u(jc,jk,jb) / mwind
            prm_diag%vmfl_s(jc,jb)  = ustar**2 * rhos * p_nh_diag%v(jc,jk,jb) / mwind

         END DO  
      END DO
!$OMP END DO 
!$OMP END PARALLEL

    !Prescribed buoyancy flux and transfer coefficient at surface to get a uniform SST (Stevens 2007 JAS)
    !It uses fixed transfer coefficient and assumes that q_s is saturated
    !Fixed surface pressure to be comparable to incompressible LES models. Pseudo density
    !(constant in time) that is used here is fixed to the initial value so that the flux 
    !is fixed in density units. 
    CASE(3)

      !Get mean theta and qv at first model level and surface pressure

      var(:,:) = theta(:,jk,:)
      WHERE(.NOT.p_patch%cells%decomp_info%owner_mask(:,:)) var(:,:) = 0._wp
      theta_nlev =  global_sum_array(var)/REAL(p_patch%n_patch_cells_g,wp)

      var(:,:) = qv(:,jk,:)
      WHERE(.NOT.p_patch%cells%decomp_info%owner_mask(:,:)) var(:,:) = 0._wp
      qv_nlev =  global_sum_array(var)/REAL(p_patch%n_patch_cells_g,wp)

      !Iterate to get surface temperature given buoyancy flux:following UCLA-LES

      !Note that t_g(:,:) is uniform so the first index is used below
      rl    = grf_bdywidth_c+1
      isblk = p_patch%cells%start_blk(rl,1)
      isidx = p_patch%cells%start_idx(rl,1)
       
      ex_sfc  = 1._wp !EXP( rd_o_cpd*LOG(pres_sfc(isidx,isblk)/p0ref) )
      t_sfc     = p_prog_lnd_now%t_g(isidx,isblk) 
      theta_sfc = t_sfc / ex_sfc       
      
      diff = 1._wp
      itr = 0 
      DO WHILE (diff > 1.e-6_wp .AND. itr < 10)
         bflx1 = les_config(jg)%tran_coeff*( (theta_sfc-theta_nlev)+vtmpc1* &
                 theta_nlev*(spec_humi(sat_pres_water(t_sfc),pres_sfc(isidx,isblk))- &
                 qv_nlev) )*grav/theta_nlev

         theta_sfc1 = theta_sfc + 0.1_wp 
         t_sfc      = theta_sfc1 * ex_sfc

         bflx2 = les_config(jg)%tran_coeff*( (theta_sfc1-theta_nlev)+vtmpc1* &
                 theta_nlev*(spec_humi(sat_pres_water(t_sfc),pres_sfc(isidx,isblk))- &
                 qv_nlev) )*grav/theta_nlev

         theta_sfc = theta_sfc1 + 0.1_wp*(les_config(jg)%bflux-bflx1)/(bflx2-bflx1) 
         t_sfc     = theta_sfc * ex_sfc
         
         diff = ABS(1._wp - theta_sfc/theta_sfc1)
         itr = itr + 1         
      END DO               
       

      !WRITE(message_text,'(i4,f14.6,f14.7)')itr,diff,bflx2
      !CALL message('FINAL ITR, RESID, AND BFLX:',message_text )

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,zrough,z_mc,mwind,RIB,ustar, &
!$OMP            obukhov_length,rhos),ICON_OMP_RUNTIME_SCHEDULE 
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx

            !Mean wind at nlev
            mwind  = MAX( les_config(jg)%min_sfc_wind, SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) )
           
            z_mc   = p_nh_metrics%z_mc(jc,jk,jb) - p_nh_metrics%z_ifc(jc,jkp1,jb)

            !Now diagnose friction velocity (ustar)
            IF(les_config(jg)%ufric<0._wp)THEN
              !Roughness length
              zrough = prm_diag%gz0(jc,jb) * rgrav

              !Bulk Richardson no at first model level
              RIB = grav * (theta(jc,jk,jb)-theta_sfc) * (z_mc-zrough)/(theta_sfc*mwind**2)
              ustar = SQRT( diag_ustar_sq(z_mc,zrough,RIB,mwind) )
              DO itr = 1 , 5
                obukhov_length = - theta_sfc*ustar**3/(akt*grav* les_config(jg)%tran_coeff*&
                                   (theta_sfc-theta(jc,jk,jb)))
                ustar = mwind / businger_mom(zrough,z_mc,obukhov_length)  
              END DO

            ELSE
              ustar = les_config(jg)%ufric
            END IF

            !Surface temperature
            p_prog_lnd_new%t_g(jc,jb) = t_sfc

            !Get surface qv 
            p_diag_lnd%qv_s(jc,jb) = spec_humi(sat_pres_water(p_prog_lnd_new%t_g(jc,jb)),pres_sfc(jc,jb))

            !rho at surface: no qc at suface
            rhos   =  pres_sfc(jc,jb)/( rd * &
                     p_prog_lnd_new%t_g(jc,jb)*(1._wp+vtmpc1*p_diag_lnd%qv_s(jc,jb)) )  

            prm_diag%shfl_s(jc,jb) = rhos*cpd*les_config(jg)%tran_coeff* &
                                     (theta_sfc-theta(jc,jk,jb))
            prm_diag%lhfl_s(jc,jb) = rhos*alv*les_config(jg)%tran_coeff* &
                                     (p_diag_lnd%qv_s(jc,jb)-qv(jc,jk,jb))
            prm_diag%umfl_s(jc,jb)  = ustar**2 * rhos * p_nh_diag%u(jc,jk,jb) / mwind
            prm_diag%vmfl_s(jc,jb)  = ustar**2 * rhos * p_nh_diag%v(jc,jk,jb) / mwind

         END DO  
      END DO
!$OMP END DO 
!$OMP END PARALLEL

   !Rico case
    CASE(4)

      !RICO case - bulk aerodynamic formulation with fixed exchange coef.
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
                              
            !Get surface qv and temperature
            p_diag_lnd%qv_s(jc,jb) = spec_humi(sat_pres_water(les_config(jg)%sst),psfc)
            p_prog_lnd_new%t_g(jc,jb) = les_config(jg)%sst
            
            !Mean wind at nlev
            mwind  = SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) 
 
            shfl  =   c_h * mwind * (th0_rico-theta(jc,jk,jb))
            lhfl  =   c_q * mwind * (p_diag_lnd%qv_s(jc,jb)-qv(jc,jk,jb))
            umfl  =   c_m * mwind * p_nh_diag%u(jc,jk,jb)
            vmfl  =   c_m * mwind * p_nh_diag%v(jc,jk,jb)
                                 
            !Surface density 
            rhos   =  psfc/( rd * &
                      p_prog_lnd_new%t_g(jc,jb)*(1._wp+vtmpc1*p_diag_lnd%qv_s(jc,jb)) )  
                      
            !Get surface fluxes                       
            prm_diag%shfl_s(jc,jb)  = shfl * rhos * cpd
            prm_diag%lhfl_s(jc,jb)  = lhfl * rhos * alv
            prm_diag%umfl_s(jc,jb)  = umfl * rhos
            prm_diag%vmfl_s(jc,jb)  = vmfl * rhos

        END DO  
      END DO
 
    !Fix SST case
    CASE(5)

!   Get roughness length * grav           
    IF(turbdiff_config(jg)%lconst_z0 .AND. turbdiff_config(jg)%const_z0 <= 0._wp)THEN
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx           
           mwind = MAX( les_config(jg)%min_sfc_wind, SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) )
           var(jc,jb) = SQRT( MAX(0._wp,prm_diag%tcm(jc,jb)) ) * mwind
        END DO         
      END DO

      WHERE(.NOT.p_patch%cells%decomp_info%owner_mask(:,:)) var(:,:) = 0._wp
      ustar_mean =  global_sum_array(var)/REAL(p_patch%n_patch_cells_g,wp)
      
      prm_diag%gz0(:,:) = MAX(0.001_wp,0.016_wp*ustar_mean**2)   
    END IF


!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,zrough,theta_sfc,mwind,z_mc, &
!$OMP            RIB,tcn_mom,tcn_heat,rhos,itr,shfl,lhfl,bflx1,ustar,   &
!$OMP            obukhov_length,inv_bus_mom),ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx

           !Roughness length
           zrough = prm_diag%gz0(jc,jb) * rgrav

           !Get surface pot. temperature and humidity
           p_prog_lnd_new%t_g(jc,jb) = les_config(jg)%sst
           theta_sfc = p_prog_lnd_new%t_g(jc,jb) / EXP( rd_o_cpd*LOG(pres_sfc(jc,jb)/p0ref) )

           p_diag_lnd%qv_s(jc,jb) = spec_humi(sat_pres_water(les_config(jg)%sst),pres_sfc(jc,jb))

           !rho at surface: no qc at suface
           rhos   =  pres_sfc(jc,jb)/( rd * &
                     p_prog_lnd_new%t_g(jc,jb)*(1._wp+vtmpc1*p_diag_lnd%qv_s(jc,jb)) )  

           mwind = MAX( les_config(jg)%min_sfc_wind, SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) ) 

           !Z height to be used as a reference height in surface layer
           z_mc   = p_nh_metrics%z_mc(jc,jk,jb) - p_nh_metrics%z_ifc(jc,jkp1,jb)

           !First guess for tch and tcm using bulk approach
           RIB = grav * (theta(jc,jk,jb)-theta_sfc) * (z_mc-zrough) / (theta_sfc * mwind**2)

           tcn_mom             = (akt/LOG(z_mc/zrough))**2
           prm_diag%tcm(jc,jb) = tcn_mom * stability_function_mom(RIB,z_mc/zrough,tcn_mom)

           tcn_heat            = akt**2/(LOG(z_mc/zrough)*LOG(z_mc/zrough))
           prm_diag%tch(jc,jb) = tcn_heat * stability_function_heat(RIB,z_mc/zrough,tcn_heat)

           !now iterate
           DO itr = 1 , 5
              shfl = prm_diag%tch(jc,jb)*mwind*(theta_sfc-theta(jc,jk,jb))
              lhfl = prm_diag%tch(jc,jb)*mwind*(p_diag_lnd%qv_s(jc,jb)-qv(jc,jk,jb))
              bflx1= shfl + vtmpc1 * theta_sfc * lhfl
              ustar= SQRT(prm_diag%tcm(jc,jb))*mwind
             
              obukhov_length = -ustar**3 * theta_sfc * rgrav / (akt * bflx1)

              inv_bus_mom = 1._wp / businger_mom(zrough,z_mc,obukhov_length)
              prm_diag%tch(jc,jb) = inv_bus_mom / businger_heat(zrough,z_mc,obukhov_length)

              prm_diag%tcm(jc,jb) = inv_bus_mom * inv_bus_mom
           END DO

           !Get surface fluxes
           prm_diag%shfl_s(jc,jb) = rhos*cpd*prm_diag%tch(jc,jb)*mwind*(theta_sfc-theta(jc,jk,jb))
           prm_diag%lhfl_s(jc,jb) = rhos*alv*prm_diag%tch(jc,jb)*mwind*(p_diag_lnd%qv_s(jc,jb)-qv(jc,jk,jb))
           prm_diag%umfl_s(jc,jb) = rhos*prm_diag%tcm(jc,jb)*mwind*p_nh_diag%u(jc,jk,jb) 
           prm_diag%vmfl_s(jc,jb) = rhos*prm_diag%tcm(jc,jb)*mwind*p_nh_diag%v(jc,jk,jb) 
           
         END DO
      END DO   
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !No fluxes
    CASE(0)

      prm_diag%shfl_s  = 0._wp
      prm_diag%lhfl_s  = 0._wp 
      prm_diag%umfl_s  = 0._wp 
      prm_diag%vmfl_s  = 0._wp 

  END SELECT 


  !Sync is required for mom fluxes 
  CALL sync_patch_array(SYNC_C, p_patch, prm_diag%umfl_s)
  CALL sync_patch_array(SYNC_C, p_patch, prm_diag%vmfl_s)


  END SUBROUTINE surface_conditions

  !>
  !! factor_heat
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile:
  !! Louis (1979) A Parametirc model of vertical eddy fluxes in the atmosphere 
  !! and R. B. Stull's book
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-06)
  FUNCTION businger_heat(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
     REAL(wp) :: factor, zeta, lamda, psi
     REAL(wp) :: zeta0, lamda0, psi0

     IF(L > 0._wp)THEN !Stable
       zeta   = z1/L 
       zeta0  = z0/L 
       IF(zeta > 1._wp)THEN !Zeng etal 1997 J. Clim
         psi    = -bsh*LOG(zeta) - zeta + 1
         psi0   = -bsh*LOG(zeta0) - zeta0 + 1
         factor = (LOG(L/z0) + bsh - psi + psi0  ) / akt
       ELSE
         psi    = -bsh*zeta
         psi0   = -bsh*zeta0
         factor = (LOG(z1/z0) - psi + psi0) / akt
       END IF
     ELSEIF(L < 0._wp)THEN !unstable
       zeta   = z1/L 
       zeta0  = z0/L 
       lamda  = SQRT(1._wp - buh*zeta)  
       lamda0 = SQRT(1._wp - buh*zeta0)  
       psi    = 2._wp * ( LOG(1._wp+lamda) - ln2 )
       psi0   = 2._wp * ( LOG(1._wp+lamda0) - ln2 )
       factor = (LOG(z1/z0) - psi + psi0) / akt
     ELSE !Neutral
       factor = LOG(z1/z0) / akt
     END IF 

  END FUNCTION businger_heat 
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  FUNCTION phi_heat(z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z1, L
     REAL(wp) :: factor, zeta, lamda

     IF(L > 0._wp)THEN !Stable
       zeta   = z1/L 
       IF(zeta > 1._wp)THEN
         factor = bsh + zeta
       ELSE     
         lamda = bsh*zeta
         factor = 1._wp + lamda
       END IF
     ELSEIF(L < 0._wp)THEN !unstable
       zeta   = z1/L 
       lamda  = SQRT(1._wp - buh*zeta)  
       factor = 1._wp / lamda
     ELSE !neutral
       factor = 1._wp
     END IF 

  END FUNCTION phi_heat 

  !>
  !! factor_mom
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile: 
  !! Louis (1979) A Parametirc model of vertical eddy fluxes in the atmosphere 
  !! and R. B. Stull's book
  !!------------------------------------------------------------------------
  FUNCTION businger_mom(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
     REAL(wp) :: factor, zeta, psi, lamda
     REAL(wp) :: zeta0, psi0, lamda0

     IF(L > 0._wp)THEN !Stable
       zeta  = z1/L 
       zeta0 = z0/L 
       IF(zeta > 1._wp)THEN !Zeng etal 1997 J. Clim
         psi    = -bsm*LOG(zeta) - zeta + 1
         psi0   = -bsm*LOG(zeta0) - zeta0 + 1
         factor = (LOG(L/z0) + bsh - psi + psi0  ) / akt
       ELSE
         psi  = -bsm*zeta
         psi0 = -bsm*zeta0
         factor = ( LOG(z1/z0) - psi + psi0 ) / akt
       END IF
     ELSEIF(L < 0._wp)THEN !unstable
       zeta   = z1/L 
       zeta0  = z0/L 
       lamda  = SQRT(SQRT(1._wp - bum*zeta))  
       lamda0 = SQRT(SQRT(1._wp - bum*zeta0))  

       psi    = 2._wp * LOG(1._wp+lamda) + LOG(1._wp+lamda*lamda) - &
                2._wp * ATAN(lamda) + pi_2 - 3._wp*ln2

       psi0   = 2._wp * LOG(1._wp+lamda0) + LOG(1._wp+lamda0*lamda0) - &
                2._wp * ATAN(lamda0) + pi_2 - 3._wp*ln2

       factor = ( LOG(z1/z0) - psi + psi0 ) / akt
     ELSE !neutral
       factor = LOG(z1/z0) / akt
     END IF

  END FUNCTION businger_mom 
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  FUNCTION phi_mom(z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z1, L
     REAL(wp) :: factor, zeta, lamda

     IF(L > 0._wp)THEN !Stable
       zeta   = z1/L 
       IF(zeta > 1._wp)THEN
         factor = bsm + zeta
       ELSE
         factor = 1._wp + bsm * zeta
       END IF
     ELSEIF(L < 0._wp)THEN !unstable 
       zeta   = z1/L 
       lamda  = SQRT(SQRT(1._wp - bum*zeta))  
       factor = 1._wp / lamda
     ELSE !neutral
       factor  = 1._wp
     END IF

  END FUNCTION phi_mom 

  !>
  !! diagnose ustar
  !!------------------------------------------------------------------------
  FUNCTION diag_ustar_sq(h, z0, RIB, wind) RESULT(ustar_sq)
     REAL(wp), INTENT(IN) :: h, z0, RIB, wind

     REAL(wp) :: ustar_sq, tcn_mom

     tcn_mom  = (akt/LOG(h/z0))**2
     ustar_sq = tcn_mom*wind**2*stability_function_mom(RIB,h/z0,tcn_mom)
        
  END FUNCTION diag_ustar_sq
 
  !>
  !! diagnose wstar
  !!------------------------------------------------------------------------
  FUNCTION diag_wstar_sq(h, z0, RIB, wind) RESULT(wstar_sq)
     REAL(wp), INTENT(IN) :: h, z0, RIB, wind

     REAL(wp) :: wstar_sq, tcn_heat, zh
  
     zh = MIN(z0, zh_max)

     IF(RIB < 0._wp)THEN
       tcn_heat = akt**2/(LOG(h/z0)*LOG(h/zh))
       wstar_sq = ( tcn_heat * wind**3 * ABS(RIB) * &
                    stability_function_heat(RIB,h/zh,tcn_heat) )**(2._wp/3._wp)
     ELSE
       wstar_sq = 0._wp
     END IF      
        
  END FUNCTION diag_wstar_sq
 
  !>
  !! init_zrough
  !!------------------------------------------------------------------------
  FUNCTION init_gz(h, wind) RESULT(gz)
     REAL(wp), INTENT(IN) :: h, wind

     REAL(wp) :: z01, z02, gz
 

     z01 = alpha0 * wind**2 / ( (1._wp/beta_10) + LOG(h/h_10)/akt )**2
     !z02 = ( alpha0 * wind**2 * ABS(RIB) )**1.5_wp / (5 * SQRT(grav*h))

     !gz = MAX(z01, z02) 
     gz = z01
        
  END FUNCTION init_gz

  !>
  !! stability_function_mom
  !! Taken from COSMO docs and Holstag & Boville 1992 
  !!------------------------------------------------------------------------
  FUNCTION stability_function_mom(RIB, hz0, tc) RESULT(stab_fun)
     REAL(wp), INTENT(IN) :: RIB, hz0, tc

     REAL(wp) :: stab_fun, hz0_fac
 
     IF(RIB.GE.0._wp)THEN
       !Cosmo
       !stab_fun = 1._wp / ( 1._wp + 10._wp*RIB/SQRT(1._wp+5*RIB) ) 

       !H&B
       stab_fun = 1._wp / ( 1._wp + 10._wp*RIB*(1._wp+8._wp*RIB) ) 
     ELSE
       hz0_fac = ( hz0**(1._wp/3._wp) - 1._wp )**1.5_wp
       !for water surface (z0/h)**(1/3)<<1 giving hz0_fac=SQRT(h/z0)
       !Generally it is explicitly written for water surface but i don't
       !see any reason to do that.
       stab_fun = 1._wp + 10._wp*ABS(RIB)/(1._wp + 75._wp*tc*hz0_fac*SQRT(ABS(RIB)))
     END IF 
        
  END FUNCTION stability_function_mom
  !>
  !! stability_function_heat
  !!------------------------------------------------------------------------
  FUNCTION stability_function_heat(RIB, hzh, tc) RESULT(stab_fun)
     REAL(wp), INTENT(IN) :: RIB, hzh, tc

     REAL(wp) :: stab_fun, hzh_fac
 
     IF(RIB.GE.0._wp)THEN
       !Cosmo
       !stab_fun = 1._wp / ( 1._wp + 15._wp*RIB*SQRT(1._wp+5*RIB) ) 
      
       !H&B
       stab_fun = 1._wp / ( 1._wp + 10._wp*RIB*(1._wp+8._wp*RIB) ) 
     ELSE
       hzh_fac = ( hzh**(1._wp/3._wp) - 1._wp )**1.5_wp
       !for water surface (zh/h)**(1/3)<<1 giving hzh_fac=SQRT(h/zh)
       !Generally it is explicitly written for water surface but i don't
       !see any reason to do that.
       stab_fun = 1._wp + 15._wp*ABS(RIB)/(1._wp + 75._wp*tc*hzh_fac*SQRT(ABS(RIB)))
     END IF 
  END FUNCTION stability_function_heat

!-------------------------------------------------------------------------------

    
END MODULE mo_surface_les




