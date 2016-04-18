SUBROUTINE send_fld_2cos

!---------------------------------------------------------------------
! Description:
!  This routine sends the fluxes from CLM3.5 to COSMO
!
! References:
!  CEREFACS/ETH: E. Maisonnave, Edoward Davin
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in CLM3.5, Initial release
! 2.1        2012/09/18 Markus Uebel 
!   Inclusion of CO2 coupling (photosynthesis rate)
! 3.1        2013/02/01 Prabhakar Shrestha
!   nee used for CO2, albd and albi allocated for direct/diffuse albedos
!   Included aerodynamic resistance and surface temperature/moisture
! 3.1        2013/07/23 Fabian Gasper
!   Bug fix in albt_gcell for sending with multiple threads
!
!   This gives 2 options for coupling COSMO and CLM i.e either flux or transfer coefficients
! @VERSION@    @DATE@     <Your name>
!  <Modification comments>         
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:
 
USE nanMod
USE clm_atmlnd   ,         ONLY : clm_l2a, atm_l2a, clm_mapl2a
USE clmtype      ,         ONLY : clm3, nameg
USE subgridAveMod,         ONLY : p2g, c2g
USE clm_varpar  ,          ONLY : numrad
USE decompMod   ,          ONLY : get_proc_global, get_proc_bounds, adecomp
USE spmdGathScatMod ,      ONLY : gather_data_to_master
USE shr_kind_mod ,         ONLY : r8 => shr_kind_r8
USE spmdMod      ,         ONLY : masterproc
USE clm_time_manager ,     ONLY :                            &
                                  get_nstep,                 &! return timestep number
                                  dtime                       ! timestep in second
USE oas_clm_vardef
!
USE netcdf                 !CPS 

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local Variables:
INTEGER            :: g              ! indices
! processor bounds indices
INTEGER            :: numg           ! total number of gridcells across all processors
INTEGER            :: numl           ! total number of landunits across all processors
INTEGER            :: numc           ! total number of columns across all processors
INTEGER            :: nump           ! total number of pfts across all processors
INTEGER            :: begg,endg      ! local beg/end gridcells gdc
INTEGER            :: begl,endl      ! local beg/end landunits
INTEGER            :: begc,endc      ! local beg/end columns
INTEGER            :: begp,endp      ! local beg/end pfts

INTEGER, PARAMETER ::   jps_taux   =  1    !  zonal wind stress
INTEGER, PARAMETER ::   jps_tauy   =  2    !  meridional wind stress
INTEGER, PARAMETER ::   jps_lat    =  3    !  total latent heat flux (W/m**2)
INTEGER, PARAMETER ::   jps_sens   =  4    !  total sensible heat flux (W/m**2)
INTEGER, PARAMETER ::   jps_ir     =  5    !  emitted infrared (longwave) radiation (W/m**2)
INTEGER, PARAMETER ::   jps_albd   =  6    !  direct albedo
INTEGER, PARAMETER ::   jps_albi   =  7    !  diffuse albedo
!MU (18.09.2012)
INTEGER, PARAMETER ::   jps_co2fl  =  8    !  net CO2 flux (now only photosynthesis rate) (umol CO2 m-2s-1)
!MU (18.09.2012)
INTEGER, PARAMETER ::   jps_ram1   =  9    !  Aerodynamic Resistance (s/m). !CPS
INTEGER, PARAMETER ::   jps_rah1   =  10    !  Aerodynamic Resistance (s/m). !CPS
INTEGER, PARAMETER ::   jps_raw1   =  11    !  Aerodynamic Resistance (s/m). !CPS
INTEGER, PARAMETER ::   jps_tsf1   =  12    !  Surface Temperature (K)  !CPS
INTEGER, PARAMETER ::   jps_qsf1   =  13    !  Surface Humidity (kg/kg) !CPS
!MU (12.04.13)
INTEGER, PARAMETER ::   jps_fpsn   =  14   !  photosynthesis rate (umol CO2 m-2s-1)
INTEGER, PARAMETER ::   jps_fplres =  15   !  plant respiration (umol CO2 m-2s-1)
!MU (12.04.13)

INTEGER            ::   isec, info, jn, ier, jj, ji, g1      ! temporary integer

REAL(KIND=r8), POINTER  :: taux_pft(:)
REAL(KIND=r8), POINTER  :: taux_gcell(:)
REAL(KIND=r8), POINTER  :: tauy_pft(:)
REAL(KIND=r8), POINTER  :: tauy_gcell(:)
REAL(KIND=r8), POINTER  :: eflx_lh_tot_pft(:)
REAL(KIND=r8), POINTER  :: eflx_lh_tot_gcell(:)
REAL(KIND=r8), POINTER  :: eflx_sh_tot_pft(:)
REAL(KIND=r8), POINTER  :: eflx_sh_tot_gcell(:)
REAL(KIND=r8), POINTER  :: eflx_lwrad_out_pft(:)
REAL(KIND=r8), POINTER  :: eflx_lwrad_out_gcell(:)
REAL(KIND=r8), POINTER  :: albd_pft(:,:)
REAL(KIND=r8), POINTER  :: albd_gcell(:,:)
REAL(KIND=r8), POINTER  :: albi_pft(:,:)
REAL(KIND=r8), POINTER  :: albi_gcell(:,:)
!MU (18.09.2012)
REAL(KIND=r8), POINTER  :: co2fl_pft(:)
REAL(KIND=r8), POINTER  :: co2fl_gcell(:)
!MU (18.09.2012)
REAL(KIND=r8), POINTER  :: ram1_pft(:)                !CPS
REAL(KIND=r8), POINTER  :: ram1_gcell(:)              !CPS
REAL(KIND=r8), POINTER  :: rah1_pft(:)                !CPS
REAL(KIND=r8), POINTER  :: rah1_gcell(:)              !CPS
REAL(KIND=r8), POINTER  :: raw1_pft(:)                !CPS
REAL(KIND=r8), POINTER  :: raw1_gcell(:)              !CPS
REAL(KIND=r8), POINTER  :: tsf1_pft(:)                !CPS
REAL(KIND=r8), POINTER  :: tsf1_gcell(:)              !CPS
REAL(KIND=r8), POINTER  :: qsf1_pft(:)                !CPS
REAL(KIND=r8), POINTER  :: qsf1_gcell(:)              !CPS
!MU (12.04.13)
REAL(KIND=r8), POINTER  :: fpsn_pft(:)
REAL(KIND=r8), POINTER  :: fpsn_gcell(:)
REAL(KIND=r8), POINTER  :: fplres_pft(:)
REAL(KIND=r8), POINTER  :: fplres_gcell(:)
!MU (12.04.13)

REAL(KIND=r8), ALLOCATABLE      :: snd_field(:)


INTEGER                         :: status, ncid, ncvarid, dimids(3) !CPS Debug Outputs
    
!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine send_fld2cos 
!------------------------------------------------------------------------------

 ! Set pointers into derived type (pft level)

 taux_pft             => clm3%g%l%c%p%pmf%taux
 tauy_pft             => clm3%g%l%c%p%pmf%tauy
 eflx_lh_tot_pft      => clm3%g%l%c%p%pef%eflx_lh_tot
 eflx_sh_tot_pft      => clm3%g%l%c%p%pef%eflx_sh_tot
 eflx_lwrad_out_pft   => clm3%g%l%c%p%pef%eflx_lwrad_out
 albd_pft             => clm3%g%l%c%p%pps%albd
 albi_pft             => clm3%g%l%c%p%pps%albi
#ifdef COUP_OAS_COS
 co2fl_pft            => clm3%g%l%c%p%pcf%fco2    !CPS net flux(+=atm)
 ram1_pft             => clm3%g%l%c%p%pps%ram1    !CPS
 rah1_pft             => clm3%g%l%c%p%pps%rah1    !CPS
 raw1_pft             => clm3%g%l%c%p%pps%raw1    !CPS
 tsf1_pft             => clm3%g%l%c%p%pes%t_sf    !CPS
 qsf1_pft             => clm3%g%l%c%p%pes%q_sf    !CPS
!MU (12.04.13)
 fpsn_pft             => clm3%g%l%c%p%pcf%fpsn
 fplres_pft           => clm3%g%l%c%p%pcf%fplres
!MU (12.04.13)
#endif 
 ! Set pointers into derived type (gridcell level)

 taux_gcell           => clm_l2a%taux
 tauy_gcell           => clm_l2a%tauy
 eflx_lh_tot_gcell    => clm_l2a%eflx_lh_tot
 eflx_sh_tot_gcell    => clm_l2a%eflx_sh_tot
 eflx_lwrad_out_gcell => clm_l2a%eflx_lwrad_out
 albd_gcell           => clm_l2a%albd
 albi_gcell           => clm_l2a%albi
#ifdef COUP_OAS_COS
 co2fl_gcell          => clm_l2a%nee    !CPS
 ram1_gcell           => clm_l2a%ram1   !CPS
 rah1_gcell           => clm_l2a%rah1   !CPS
 raw1_gcell           => clm_l2a%raw1   !CPS
 tsf1_gcell           => clm_l2a%t_sf   !CPS
 qsf1_gcell           => clm_l2a%q_sf   !CPS
!MU (12.04.13)
 fpsn_gcell           => clm_l2a%fpsn
 fplres_gcell         => clm_l2a%fplres
!MU (12.04.13)
#endif
 ! Get total global number of grid cells, landunits, columns and pfts

 CALL get_proc_global(numg,numl,numc,nump)
 CALL get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)


     ALLOCATE(snd_field(begg:endg), stat=ier)
     IF (ier /= 0) THEN
        CALL prism_abort_proto( ncomp_id, 'send_fld_2cos', 'Failure in allocating snd_field' )
        RETURN
     ENDIF


isec = dtime * ( get_nstep() - 1 )

 ! Do flux averaging from pft to gridcell

 ! wind stress in x direction
 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, taux_pft, taux_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

 ! wind stress in y direction
 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, tauy_pft, tauy_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

 ! Latent Heat Flux
 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, eflx_lh_tot_pft, eflx_lh_tot_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

 ! Sensible Heat Flux
 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, eflx_sh_tot_pft, eflx_sh_tot_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

 ! Longwave Upward Radiation
 CALL  p2g(begp, endp, begc, endc, begl, endl, begg, endg, eflx_lwrad_out_pft, eflx_lwrad_out_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

 ! albedo (direct beam)
 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, albd_pft, albd_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

 ! albedo (diffuse)
 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, albi_pft, albi_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

#ifdef COUP_OAS_COS
 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, co2fl_pft, co2fl_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, ram1_pft, ram1_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
  
 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, rah1_pft, rah1_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, raw1_pft, raw1_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
 
 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, tsf1_pft, tsf1_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, qsf1_pft, qsf1_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

!MU (12.04.13)
 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, fpsn_pft, fpsn_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

 CALL p2g(begp, endp, begc, endc, begl, endl, begg, endg, fplres_pft, fplres_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
!MU (12.04.13)
#endif 
!CPS WRITE(6,*) "oasclm: send_fld_2cos: p2g done"   
 

 ! Create remapping matrix between land and atmosphere model
 CALL clm_mapl2a(clm_l2a, atm_l2a)




snd_field(:)=atm_l2a%taux
IF( ssnd(jps_taux)%laction )  CALL oas_clm_snd( jps_taux, isec,snd_field(:),begg,endg, info )

snd_field(:)=atm_l2a%tauy
IF( ssnd(jps_tauy)%laction )  CALL oas_clm_snd( jps_tauy, isec,snd_field(:),begg,endg, info )

snd_field(:)=atm_l2a%eflx_lh_tot
IF( ssnd(jps_lat)%laction )  CALL oas_clm_snd( jps_lat, isec,snd_field(:),begg,endg, info )

snd_field(:)=atm_l2a%eflx_sh_tot
IF( ssnd(jps_sens)%laction )  CALL oas_clm_snd( jps_sens, isec,snd_field(:),begg,endg, info )

snd_field(:)=atm_l2a%eflx_lwrad_out
IF( ssnd(jps_ir)%laction )  CALL oas_clm_snd( jps_ir, isec,snd_field(:),begg,endg, info )

snd_field(:)= 0.5_r8*(atm_l2a%albd(:,1)+atm_l2a%albd(:,2))
IF( ssnd(jps_albd)%laction )  CALL oas_clm_snd( jps_albd, isec,snd_field(:),begg,endg, info )

snd_field(:)= 0.5_r8*(atm_l2a%albi(:,2)+atm_l2a%albi(:,2))
IF( ssnd(jps_albi)%laction )  CALL oas_clm_snd( jps_albi, isec,snd_field(:),begg,endg, info )

#ifdef COUP_OAS_COS
snd_field(:)=atm_l2a%nee
IF( ssnd(jps_co2fl)%laction )  CALL oas_clm_snd( jps_co2fl, isec,snd_field(:),begg,endg, info )

snd_field(:)=atm_l2a%ram1
IF( ssnd(jps_ram1)%laction )  CALL oas_clm_snd( jps_ram1, isec,snd_field(:),begg,endg, info )

snd_field(:)=atm_l2a%rah1
IF( ssnd(jps_rah1)%laction )  CALL oas_clm_snd( jps_rah1, isec,snd_field(:),begg,endg, info )

snd_field(:)=atm_l2a%raw1
IF( ssnd(jps_raw1)%laction )  CALL oas_clm_snd( jps_raw1, isec,snd_field(:),begg,endg, info )

snd_field(:)=atm_l2a%t_sf
IF( ssnd(jps_tsf1)%laction )  CALL oas_clm_snd( jps_tsf1, isec,snd_field(:),begg,endg, info )

snd_field(:)=atm_l2a%q_sf
IF( ssnd(jps_qsf1)%laction )  CALL oas_clm_snd( jps_qsf1, isec,snd_field(:),begg,endg, info )

!MU (12.04.13)
snd_field(:)=atm_l2a%fpsn
IF( ssnd(jps_fpsn)%laction )  CALL oas_clm_snd( jps_fpsn, isec,snd_field(:),begg,endg, info )

snd_field(:)=atm_l2a%fplres
IF( ssnd(jps_fplres)%laction )  CALL oas_clm_snd( jps_fplres, isec,snd_field(:),begg,endg, info )
!MU (12.04.13)

#endif


  DEALLOCATE(snd_field)


!CPS WRITE(6,*) 'oasclm: send complete'

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------


END SUBROUTINE send_fld_2cos
