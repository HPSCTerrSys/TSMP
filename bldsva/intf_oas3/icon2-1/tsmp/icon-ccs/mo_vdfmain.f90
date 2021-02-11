!>
!! Top for single VDF call for EDMF DUALM
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2011-9-30)
!!   (IFS cycle CY36R1_DUALM_M8b)
!!
!! Modifications by Dmitrii Mironov, DWD (2016-08-08)
!! - Changes related to the use of prognostic the sea-ice albedo.
!!
!!-----------------------------------------------------------------------------
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!-----------------------------------------------------------------------------

MODULE mo_vdfmain
 
  PUBLIC :: vdfmain

CONTAINS

!! !OPTIONS XOPT(HSFUN)
!! #ifdef RS6K
!! @PROCESS NOSTRICT
!! #endif
SUBROUTINE VDFMAIN ( CDCONF , &
 & KIDIA  , KFDIA  , KLON   , KLEV   , KLEVS  , KSTEP  , KTILES , &
 & KTRAC  , KLEVSN , KLEVI  , KDHVTLS, KDHFTLS, KDHVTSS, KDHFTSS, &
 & KDHVTTS, KDHFTTS, KDHVTIS, KDHFTIS, &
 & PTSPHY , KTVL   , KTVH   , KCNT   , PCVL   , PCVH   , PSIGFLT, &
 & PUM1   , PVM1   , PTM1   , PQM1   , PLM1   , PIM1   , PAM1   , PCM1   , &
 & PAPHM1 , PAPM1  , PGEOM1 , PGEOH  , PTSKM1M, PTSAM1M, PWSAM1M, &
 & PSSRFL , PSLRFL , PEMIS  , PHRLW  , PHRSW  , &
 & PTSNOW , PTICE  , &
 & PHLICE , PTLICE , PTLWML , &  
 & PSST   , KSOTY  , PFRTI  , PALBTI , PWLMX  , &
 & PCHAR  , PUCURR , PVCURR , PTSKRAD, PCFLX  , &
 & PSOTEU , PSOTEV , PSOBETA, PVERVEL, &
 ! OUTPUT
 & PZ0M   , PZ0H   , &
 & PVDIS  , PVDISG , PDISGW3D,PAHFLEV, PAHFLSB, PFWSB  , PBIR   , PVAR   , &
 & PU10M  , PV10M  , PT2M   , PD2M   , PQ2M   , PZINV  , PBLH   , KHPBLN , KVARTOP, &
 & PSSRFLTI,PEVAPSNW,PGUST  , PWUAVG , LDNODECP,KPBLTYPE,PLDIFF , &
 & PFPLVL , PFPLVN , PFHPVL , PFHPVN , &
 ! DIAGNOSTIC OUTPUT
 & PEXTR2 , KFLDX2 , PEXTRA , KLEVX  , KFLDX  , LLDIAG , &
 ! OUTPUT TENDENCIES
 & PTE    , PQE    , PLE    , PIE    , PAE    , PVOM   , PVOL   , &
 & PTENC  , PTSKE1 , &
 ! UPDATED FIELDS FOR TILES
 & PUSTRTI, PVSTRTI, PAHFSTI, PEVAPTI, PTSKTI , &
 ! OUTPUT FLUXES
 & PDIFTS , PDIFTQ , PDIFTL , PDIFTI , PSTRTU , PSTRTV , PTOFDU , PTOFDV, &
 & PSTRSOU, PSTRSOV, PKH    , &
!amk
 & LDLAND , &   
!xxx
 ! DDH OUTPUTS
 & PDHTLS , PDHTSS , PDHTTS , PDHTIS &
! TERRA data
 & , ext_data                                                           & !in
 & , jb, jg                                                             & ! -
 & , t_snow_ex, t_snow_mult_ex, t_s_ex, t_g_ex, qv_s_ex                 & !inout
 & , w_snow_ex                                                          & ! -
 & , rho_snow_ex, rho_snow_mult_ex, h_snow_ex, w_i_ex, w_p_ex, w_s_ex   & ! -
 & , t_so_ex, w_so_ex, w_so_ice_ex  &  !, t_2m_ex, u_10m_ex, v_10m_ex   & ! -
 & , freshsnow_ex, snowfrac_lc_ex, snowfrac_ex                          & ! -
 & , wliq_snow_ex, wtot_snow_ex, dzh_snow_ex                            & ! -
 & , prr_con_ex, prs_con_ex, prr_gsp_ex, prs_gsp_ex                     & !in
 & , tch_ex, tcm_ex, tfv_ex                                             & !inout
 & , sobs_ex, thbs_ex, pabs_ex                                          & !in
 & , runoff_s_ex, runoff_g_ex                                           & !inout
 & , t_g, qv_s                                                          & ! -
 & , t_ice, h_ice, t_snow_si, h_snow_si, alb_si                         & ! -
 & , fr_seaice                                                          ) ! -
!***

!**   *VDFMAIN* - DOES THE VERTICAL EXCHANGE OF U,V,SLG,QT BY TURBULENCE.

!     J.F.GELEYN       20/04/82   Original  
!     C.A.BLONDIN      18/12/86
!     A.C.M. BELJAARS  20/10/89   IFS-VERSION (TECHNICAL REVISION OF CY34)
!     A.C.M. BELJAARS  26/03/90   OBUKHOV-L UPDATE 
!     A.C.M. BELJAARS  30/09/98   SURFACE TILES 
!     P. Viterbo       17/05/2000 Surface DDH for TILES
!     D. Salmond       15/10/2001 FULLIMP mods
!     S. Abdalla       27/11/2001 Passing Zi/L to waves
!     A. Beljaars       2/05/2003 New tile coupling     
!     P.Viterbo        24/05/2004 Change surface units
!     M. Ko"hler        3/12/2004 Moist Advection-Diffusion
!     A. Beljaars       4/04/2005 Turb. orogr. drag
!     A. Beljaars      30/09/2005 Include Subgr. Oro. in solver
!     A.Beljaars       31/03/2005 Introduction of ocean current b.c.
!     P. Viterbo       17/06/2005 surf external library
!     M. Ko"hler        6/06/2006 Single Column Model option (LSCMEC) 
!     R. Neggers       24/05/2006 PBL dualM scheme + bimodal cloud scheme
!     G. Balsamo       15/01/2007 Soil type
!     A. Beljaars      27/02/2009 Delete PZIDLWV
!     G. Balsamo       07/04/2008 Lake model (FLAKE)
!     M. Janiskova/G. Balsamo 21/10/2008 Move VDFDIFH5 in VDFDIFH
!     M. Ko"hler/G. Modzinski     Bug correction in optimised qsat computation
!     P.de Rosnay/G.Balsamo   07/03/2009 Offline Jacobians EKF 
 
!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE FOUR
!     PROGNOSTIC VARIABLES U,V,T AND Q DUE TO THE VERTICAL EXCHANGE BY
!     TURBULENT (= NON-MOIST CONVECTIVE) PROCESSES. THESE TENDENCIES ARE
!     OBTAINED AS THE DIFFERENCE BETWEEN THE RESULT OF AN IMPLICIT
!     TIME-STEP STARTING FROM VALUES AT T-1 AND THESE T-1 VALUES. ALL
!     THE DIAGNOSTIC COMPUTATIONS (EXCHANGE COEFFICIENTS, ...) ARE DONE
!      FROM THE T-1 VALUES. AS A BY-PRODUCT THE ROUGHNESS LENGTH OVER SEA
!     IS UPDATED ACCORDINGLY TO THE *CHARNOCK FORMULA. HEAT AND MOISTURE
!     SURFACE FLUXES AND THEIR DERIVATIVES AGAINST TS, WS AND WL
!     (THE LATTER WILL BE LATER WEIGHTED WITH THE SNOW FACTOR IN
!     *VDIFF*), LATER TO BE USED FOR SOIL PROCESSES TREATMENT, ARE ALSO
!     COMPUTED AS WELL AS A STABILITY VALUE TO BE USED AS A DIAGNOSTIC
!     OF THE DEPTH OF THE WELL MIXED LAYER IN CONVECTIVE COMPUTATIONS.

!     INTERFACE.
!     ----------
!          *VDIFF* TAKES THE MODEL VARIABLES AT T-1 AND RETURNS THE VALUES
!     FOR THE PROGNOSTIC TIME T+1 DUE TO VERTICAL DIFFUSION.
!     THE MODEL VARIABLES, THE MODEL DIMENSIONS AND THE DIAGNOSTICS DATA
!     ARE PASSED AS SUBROUTINE ARGUMENTS. CONSTANTS THAT DO NOT CHANGE
!     DURING A MODEL RUN (E.G. PHYSICAL CONSTANTS, SWITCHES ETC.) ARE
!     STORED IN A SINGLE COMMON BLOCK *YOMVDF*, WHICH IS INITIALIZED
!     BY SET-UP ROUTINE *SUVDF*.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLEV*         NUMBER OF LEVELS
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*        NUMBER OF SOIL LAYERS
!    *KSTEP*        CURRENT TIME STEP INDEX
!    *KTILES*       NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT 
!                   OF SURFACE BOUNDARY CONDITION)
!    *KTRAC*        Number of tracers
!    *KLEVSN*       Number of snow layers (diagnostics) 
!    *KLEVI*        Number of sea ice layers (diagnostics)
!    *KDHVTLS*      Number of variables for individual tiles
!    *KDHFTLS*      Number of fluxes for individual tiles
!    *KDHVTSS*      Number of variables for snow energy budget
!    *KDHFTSS*      Number of fluxes for snow energy budget
!    *KDHVTTS*      Number of variables for soil energy budget
!    *KDHFTTS*      Number of fluxes for soil energy budget
!    *KDHVTIS*      Number of variables for sea ice energy budget
!    *KDHFTIS*      Number of fluxes for sea ice energy budget

!    *KTVL*         VEGETATION TYPE FOR LOW VEGETATION FRACTION
!    *KTVH*         VEGETATION TYPE FOR HIGH VEGETATION FRACTION
!    *KSOTY*        SOIL TYPE                                     (1-7)

!    *KCNT*         Index of vdf sub steps.

!     INPUT PARAMETERS (LOGICAL)

!     INPUT PARAMETERS (REAL)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):

!    *PCVL*         LOW VEGETATION COVER                          -  
!    *PCVH*         HIGH VEGETATION COVER                         -  
!    *PSIGFLT*      STANDARD DEVIATION OF FILTERED OROGRAPHY      M
!    *PUM1*         X-VELOCITY COMPONENT                          M/S
!    *PVM1*         Y-VELOCITY COMPONENT                          M/S
!    *PTM1*         TEMPERATURE                                   K
!    *PQM1*         SPECIFIC HUMIDITY                             KG/KG
!    *PLM1*         SPECIFIC CLOUD LIQUID WATER                   KG/KG
!    *PIM1*         SPECIFIC CLOUD ICE                            KG/KG
!    *PAM1*         CLOUD FRACTION                                1
!    *PCM1*         TRACER CONCENTRATION                          KG/KG
!    *PAPM1*        PRESSURE ON FULL LEVELS                       PA
!    *PAPHM1*       PRESSURE ON HALF LEVELS                       PA
!    *PGEOM1*       GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL AT HALF LEVELS                   M2/S2
!    *PTSKM1M*      SKIN TEMPERATURE                              K
!    *PTSAM1M*      SURFACE TEMPERATURE                           K
!    *PWSAM1M*      SOIL MOISTURE ALL LAYERS                      M**3/M**3
!    *PSSRFL*       NET SHORTWAVE RADIATION FLUX AT SURFACE       W/M2
!    *PSLRFL*       NET LONGWAVE RADIATION FLUX AT SURFACE        W/M2
!    *PEMIS*        MODEL SURFACE LONGWAVE EMISSIVITY
!    *PHRLW*        LONGWAVE HEATING RATE                         K/s
!    *PHRSW*        SHORTWAVE HEATING RATE                        K/s
!    *PTSNOW*       SNOW TEMPERATURE                              K
!    *PTICE*        ICE TEMPERATURE (TOP SLAB)                    K
!    *PHLICE*       LAKE ICE THICKNESS                            m
!    *PTLICE*       LAKE ICE TEMPERATURE                          K
!    *PTLWML*       LAKE MEAN WATER TEMPERATURE                   K
!    *PSST*         (OPEN) SEA SURFACE TEMPERATURE                K
!    *PFRTI*        TILE FRACTIONS                                (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!    *PALBTI*       BROADBAND ALBEDO FOR TILE FRACTIONS
!    *PWLMX*        MAXIMUM SKIN RESERVOIR CAPACITY               kg/m**2
!    *PCHAR*        "EQUIVALENT" CHARNOCK PARAMETER
!    *PUCURR*       OCEAN CURRENT X_COMPONENT      
!    *PVCURR*       OCEAN CURRENT Y_COMPONENT      
!    *PTSKRAD*      SKIN TEMPERATURE OF LATEST FULL RADIATION
!                      TIMESTEP                                   K
!    *PCFLX*        TRACER SURFACE FLUX                           kg/(m2 s)
!    *PSOTEU*       Explicit part of U-tendency from subgrid orography scheme    
!    *PSOTEV*       Explicit part of V-tendency from subgrid orography scheme     
!    *PSOBETA*      Implicit part of subgrid orography 

!    *PVERVEL*      VERTICAL VELOCITY

!     INPUT PARAMETERS (LOGICAL):

!     CONTRIBUTIONS TO BUDGETS (OUTPUT,REAL):

!    *PVDIS*        TURBULENT DISSIPATION                         W/M2
!    *PVDISG*       SUBGRID OROGRAPHY DISSIPATION                 W/M2
!    *PAHFLEV*      LATENT HEAT FLUX  (SNOW/ICE FREE PART)        W/M2
!    *PAHFLSB*      LATENT HEAT FLUX  (SNOW/ICE COVERED PART)     W/M2

!     UPDATED PARAMETERS (REAL):

!    *PTE*          TEMPERATURE TENDENCY                          K/S
!    *PQE*          MOISTURE TENDENCY                             KG/(KG S)
!    *PLE*          LIQUID WATER TENDENCY                         KG/(KG S)
!    *PIE*          ICE WATER TENDENCY                            KG/(KG S)
!    *PAE*          CLOUD FRACTION TENDENCY                       1/S)
!    *PVOM*         MERIODINAL VELOCITY TENDENCY (DU/DT)          M/S2
!    *PVOL*         LATITUDE TENDENCY            (DV/DT)          M/S2
!    *PTENC*        TRACER TENDENCY                               KG/(KG S)
!    *PTSKE1*       SKIN TEMPERATURE TENDENCY                     K/S
!    *PZ0M*         AERODYNAMIC ROUGHNESS LENGTH                  M
!    *PZ0H*         ROUGHNESS LENGTH FOR HEAT                     M

!     UPDATED PARAMETERS FOR TILES (REAL): 

!    *PUSTRTI*      SURFACE U-STRESS                              N/M2 
!    *PVSTRTI*      SURFACE V-STRESS                              N/M2
!    *PAHFSTI*      SURFACE SENSIBLE HEAT FLUX                    W/M2
!    *PEVAPTI*      SURFACE MOISTURE FLUX                         KG/M2/S
!    *PTSKTI*       SKIN TEMPERATURE                              K

!     OUTPUT PARAMETERS (REAL):

!    *PFPLVL*       PBL PRECIPITATION FLUX AS RAIN                KG/(M**2*S)
!    *PFPLVN*       PBL PRECIPITATION FLUX AS SNOW                KG/(M**2*S)
!    *PFHPVL*       ENTHALPY FLUX OF PBL PRECIPITATION AS RAIN    J/(M**2*S)
!    *PFHPVN*       ENTHALPY FLUX OF PBL PRECIPITATION AS SNOW    J/(M**2*S)

!    *PLDIFF*       CONTRIB TO PBL CONDENSATE BY PASSIVE CLOUDS   KG/KG

!    *PFWSB*        EVAPORATION OF SNOW                           KG/(M**2*S)
!    *PU10M*        U-COMPONENT WIND AT 10 M                      M/S
!    *PV10M*        V-COMPONENT WIND AT 10 M                      M/S
!    *PT2M*         TEMPERATURE AT 2M                             K
!    *PD2M*         DEW POINT TEMPERATURE AT 2M                   K
!    *PQ2M*         SPECIFIC HUMIDITY AT 2M                       KG/KG
!    *PGUST*        GUST AT 10 M                                  M/S
!    *PBLH*         PBL HEIGHT (dry diagnostic based on Ri#)      M
!    *PZINV*        PBL HEIGHT (moist parcel, not for stable PBL) M
!    *PSSRFLTI*     NET SHORTWAVE RADIATION FLUX AT SURFACE, FOR
!                      EACH TILE                                  W/M2
!    *PEVAPSNW*     EVAPORATION FROM SNOW UNDER FOREST            KG/(M2*S)
!    *PSTRTU*       TURBULENT FLUX OF U-MOMEMTUM            KG*(M/S)/(M2*S)
!    *PSTRTV*       TURBULENT FLUX OF V-MOMEMTUM            KG*(M/S)/(M2*S)
!    *PTOFDU*       TOFD COMP. OF TURBULENT FLUX OF U-MOMEMTUM   KG*(M/S)/(M2*S)
!    *PTOFDV*       TOFD COMP. OF TURBULENT FLUX OF V-MOMEMTUM   KG*(M/S)/(M2*S)
!    *PDIFTS*       TURBULENT FLUX OF HEAT                        W/M2
!    *PDIFTQ*       TURBULENT FLUX OF SPECIFIC HUMIDITY           KG/(M2*S)
!    *PDIFTL*       TURBULENT FLUX OF LIQUID WATER                KG/(M2*S)
!    *PDIFTI*       TURBULENT FLUX OF ICE WATER                   KG/(M2*S)
!    *PSTRSOU*      SUBGRID OROGRAPHY FLUX OF U-MOMEMTUM    KG*(M/S)/(M2*S)
!    *PSTRSOV*      SUBGRID OROGRAPHY FLUX OF V-MOMEMTUM    KG*(M/S)/(M2*S)

!    *PKH*          TURB. DIFF. COEFF. FOR HEAT ABOVE SURF. LAY.  (M2/S)
!                   IN SURFACE LAYER: CH*U                        (M/S)
!    *PDHTLS*       Diagnostic array for tiles (see module yomcdh)
!                      (Wm-2 for energy fluxes, kg/(m2s) for water fluxes)
!    *PDHTSS*       Diagnostic array for snow T (see module yomcdh)
!                      (Wm-2 for fluxes)
!    *PDHTTS*       Diagnostic array for soil T (see module yomcdh)
!                      (Wm-2 for fluxes)
!    *PDHTIS*       Diagnostic array for ice T (see module yomcdh)
!                      (Wm-2 for fluxes)

!     Additional parameters for flux boundary condition (in SCM model):

!    *LSFCFLX*      If .TRUE. flux boundary condtion is used 
!    *REXTSHF*      Specified sensible heat flux [W/m2]
!    *REXTLHF*      Specified latent heat flux [W/m2]

!     METHOD.
!     -------

!          FIRST AN AUXIALIARY VARIABLE CP(Q)T+GZ IS CREATED ON WHICH
!     THE VERTICAL DIFFUSION PROCESS WILL WORK LIKE ON U,V AND Q. THEN
!     ALONG THE VERTICAL AND AT THE SURFACE, EXCHANGE COEFFICIENTS (WITH
!     THE DIMENSION OF A PRESSURE THICKNESS) ARE COMPUTED FOR MOMENTUM
!     AND FOR HEAT (SENSIBLE PLUS LATENT). THE LETTERS M AND H ARE USED
!     TO DISTINGUISH THEM AND THE COMPUTATION IS THE RESULT OF A
!     CONDITIONAL MERGE BETWEEN THE STABLE AND THE UNSTABLE CASE
!     (DEPENDING ON THE SIGN OF THE *RICHARDSON BULK NUMBER).
!          IN THE SECOND PART OF THE ROUTINE THE IMPLICIT LINEAR
!     SYSTEMS FOR U,V FIRST AND T,Q SECOND ARE SOLVED BY A *GAUSSIAN
!     ELIMINATION BACK-SUBSTITUTION METHOD. FOR T AND Q THE LOWER
!     BOUNDARY CONDITION DEPENDS ON THE SURFACE STATE.
!     OVER LAND, TWO DIFFERENT REGIMES OF EVAPORATION PREVAIL:
!     A STOMATAL RESISTANCE DEPENDENT ONE OVER THE VEGETATED PART
!     AND A SOIL RELATIVE HUMIDITY DEPENDENT ONE OVER THE
!     BARE SOIL PART OF THE GRID MESH.
!     POTENTIAL EVAPORATION TAKES PLACE OVER THE SEA, THE SNOW
!     COVERED PART AND THE LIQUID WATER COVERED PART OF THE
!     GRID MESH AS WELL AS IN CASE OF DEW DEPOSITION.
!          FINALLY ONE RETURNS TO THE VARIABLE TEMPERATURE TO COMPUTE
!     ITS TENDENCY AND THE LATER IS MODIFIED BY THE DISSIPATION'S EFFECT
!     (ONE ASSUMES NO STORAGE IN THE TURBULENT KINETIC ENERGY RANGE) AND
!     THE EFFECT OF MOISTURE DIFFUSION ON CP. Z0 IS UPDATED AND THE
!     SURFACE FLUXES OF T AND Q AND THEIR DERIVATIVES ARE PREPARED AND
!     STORED LIKE THE DIFFERENCE BETWEEN THE IMPLICITELY OBTAINED
!     CP(Q)T+GZ AND CP(Q)T AT THE SURFACE.

!     EXTERNALS.
!     ----------

!     *VDFMAIN* CALLS SUCESSIVELY:
!         *SURFEXCDRIVER*
!         *VDFEXCU*
!         *VDFTOFDC*
!         *VDFDIFM*
!         *VDFDIFH*
!         *VDFDIFC*
!         *VDFINCR*
!         *VDFSDRV*
!         *VDFPPCFL*
!         *VDFUPDZ0*

!     REFERENCE.
!     ----------

!          SEE VERTICAL DIFFUSION'S PART OF THE MODEL'S DOCUMENTATION
!     FOR DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
! USE YOMCT0   , ONLY : LSCMEC   ,LSFCFLX  ,REXTSHF  ,REXTLHF
! USE YOEVDF   , ONLY : RVDIFTS
! USE YOMCST   , ONLY : RG       ,RD       ,&
!                     & RCPD     ,RETV     ,RLVTT    ,RLSTT    ,RTT 
! USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
!                     & R4IES    ,R5LES    ,R5IES    ,RVTMP2   ,R5ALVCP  ,&
!                     & R5ALSCP  ,RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,&
!                     & RTICECU  ,RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
! USE YOMJFH   , ONLY : N_VMASS
! USE YOEPHY   , ONLY : LVDFTRAC 
! USE YOMSEKF  , ONLY : N_SEKF_PT, LUSEKF_REF, LUSE_JATM
! amk
! USE YOECLDP   , ONLY : NCLDTOP
! xxx

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
                & RG       ,RD       ,&                               !yomcst
                & RCPD     ,RETV     ,RLVTT    ,RLSTT    ,RTT      ,& ! -
                & R2ES     ,R3LES    ,R3IES    ,R4LES    ,&           !yoethf
                & R4IES    ,R5LES    ,R5IES    ,RVTMP2   ,R5ALVCP  ,& ! -
                & R5ALSCP  ,RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,& ! -
                & RTICECU  ,RTWAT_RTICE_R      ,RTWAT_RTICECU_R    ,& ! -
                & RVDIFTS                                             !yoevdf
USE mo_edmf_param   ,ONLY : &
                & LSCMEC   ,LSFCFLX  ,REXTSHF  ,REXTLHF  ,&           !yomct0
                & LVDFTRAC ,&                                         !yoephy 
                & N_SEKF_PT          ,LUSEKF_REF         ,LUSE_JATM,& !yomsekf
                & N_VMASS  ,&                                         !yomjfh
                & FOEALFA  ,&                                         !fcttre.f
                & ntiles_edmf
USE mo_physical_constants,  ONLY:  p0ref, rd_o_cpd
USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
USE mo_lnd_nwp_config,ONLY: nlev_soil, nlev_snow, ntiles_total, ntiles_water, &
                            isub_water, isub_seaice
USE mo_ext_data_types,ONLY: t_external_data

USE mo_vdfdpbl      ,ONLY : vdfdpbl
USE mo_vdfhghtn     ,ONLY : vdfhghtn
USE mo_vdfexcu      ,ONLY : vdfexcu
USE mo_vdftofdc     ,ONLY : vdftofdc
USE mo_vdfdifm      ,ONLY : vdfdifm
USE mo_vdfdifh      ,ONLY : vdfdifh
USE mo_vdfincr      ,ONLY : vdfincr
USE mo_vdfdifc      ,ONLY : vdfdifc
USE mo_vdffblend    ,ONLY : vdffblend
USE mo_vdfcloud     ,ONLY : vdfcloud
USE mo_surfexcdriver,ONLY : surfexcdriver
USE mo_surfpp       ,ONLY : surfpp

USE mo_run_config,   ONLY : ltestcase
USE mo_nh_testcases_nml, ONLY : nh_test_name
USE mo_nh_torus_exp ,ONLY : cbl_stevens_fluxes

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVSN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTIS 
CHARACTER(LEN=1)  ,INTENT(IN)    :: CDCONF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCNT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIGFLT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCM1(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNOW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTICE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(KLON)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLICE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLWML(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKRAD(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFLX(KLON,KTRAC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOBETA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0M(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0H(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDISG(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISGW3D(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLEV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLSB(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFWSB(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBIR(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVAR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU10M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV10M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZINV(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBLH(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KHPBLN(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSRFLTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEVAPSNW(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLVL(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLVN(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPVL(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPVN(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGUST(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWUAVG(KLON)
LOGICAL           ,INTENT(IN)    :: LDNODECP(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KPBLTYPE(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLDIFF(KLON,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KVARTOP(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PIE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKE1(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTQ(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTL(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTI(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOFDU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOFDV(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRSOU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRSOV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTLS(KLON,KTILES,KDHVTLS+KDHFTLS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTSS(KLON,KLEVSN,KDHVTSS+KDHFTSS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTTS(KLON,KLEVS,KDHVTTS+KDHFTTS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTIS(KLON,KLEVI,KDHVTIS+KDHFTIS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
!amk
LOGICAL                          :: LDLAND(KLON)
!xxx
!          DIAGNOSTIC OUTPUT
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX2, KLEVX, KFLDX
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTR2(KLON,KFLDX2), PEXTRA(KLON,KLEVX,KFLDX)
LOGICAL           ,INTENT(IN)    :: LLDIAG

! TERRA data

INTEGER          ,INTENT(IN)                                               :: &
  jb             ,jg                 
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,0:nlev_snow,ntiles_total) :: &
  t_snow_mult_ex 
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,nlev_snow,ntiles_total)   :: &
  rho_snow_mult_ex  
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  t_g_ex         ,qv_s_ex  
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total)             :: &
  t_snow_ex      ,t_s_ex         ,                                            & 
  w_snow_ex      ,rho_snow_ex    ,h_snow_ex       ,                           &
  w_i_ex         ,w_p_ex         ,w_s_ex
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,0:nlev_soil,ntiles_total) :: &
  t_so_ex             
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,nlev_soil,ntiles_total)   :: &
  w_so_ex        ,w_so_ice_ex          
!REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON)                         :: &
! u_10m_ex       ,v_10m_ex             
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total)             :: &
  freshsnow_ex   ,snowfrac_lc_ex ,snowfrac_ex 
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,nlev_snow,ntiles_total)   :: &
  wliq_snow_ex   ,wtot_snow_ex   ,dzh_snow_ex          
REAL(KIND=JPRB)  ,INTENT(IN)     ,DIMENSION(KLON)                          :: &
  prr_con_ex     ,prs_con_ex     ,prr_gsp_ex     ,prs_gsp_ex           
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  tch_ex         ,tcm_ex         ,tfv_ex               
REAL(KIND=JPRB)  ,INTENT(IN)     ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  sobs_ex        ,thbs_ex        ,pabs_ex              
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total)             :: &
  runoff_s_ex    ,runoff_g_ex        
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON)                          :: &
  t_g            ,qv_s
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON)                          :: &
  t_ice          ,h_ice          ,t_snow_si      ,h_snow_si       ,alb_si 
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON)                          :: &
  fr_seaice
TYPE(t_external_data), INTENT(INOUT)                                       :: &
  ext_data

!*         0.2    LOCAL VARIABLES

REAL(KIND=JPRB) ::    ZDIFTQT(KLON,0:KLEV), ZDIFTSLG(KLON,0:KLEV) 

REAL(KIND=JPRB) ::    ZCPTGZ(KLON,KLEV) , ZCFM(KLON,KLEV)   , ZCFH(KLON,KLEV)   ,&
                    & ZUDIF(KLON,KLEV)  , ZVDIF(KLON,KLEV)  ,&
                    & ZQTDIF(KLON,KLEV) , ZSLGDIF(KLON,KLEV),&
                    & ZSLGM1(KLON,KLEV) , ZQTM1(KLON,KLEV)  , ZQTE(KLON,KLEV)   ,&
                    & ZSLGE(KLON,KLEV)  , ZTOFDC(KLON,KLEV) , ZSOC(KLON,KLEV)
REAL(KIND=JPRB) ::    ZKHFL(KLON)       , ZKQFL(KLON)       , ZKMFL(KLON)  
REAL(KIND=JPRB) ::    ZQEA(KLON,KLEV)   , ZLEA(KLON,KLEV)   , ZIEA(KLON,KLEV)   ,&
                    & ZQTEA(KLON,KLEV)  , ZSLGEA(KLON,KLEV) ,                    &
                    & ZSLGEWODIS(KLON,KLEV) , ZPTEA(KLON,KLEV)
REAL(KIND=JPRB) ::    ZQTEP(KLON,KLEV)  , ZSLGEP(KLON,KLEV), ZDZRHOI

!RN --- variables associated with dualM scheme -------------------------

INTEGER(KIND=JPIM), PARAMETER   ::   IDRAFT = 3    ! nr of updrafts

   !  thermodynamic transport
REAL(KIND=JPRB) ::    ZMFLX(KLON,0:KLEV,IDRAFT)  ,  ZQTUH(KLON,0:KLEV,IDRAFT)   ,&   
                    & ZSLGUH(KLON,0:KLEV,IDRAFT) ,  ZWUH(KLON,0:KLEV,IDRAFT)
   !  momentum transport
REAL(KIND=JPRB) ::    ZMFLXM(KLON,0:KLEV,IDRAFT)              , ZUUH(KLON,0:KLEV,IDRAFT),&
                    & ZVUH(KLON,0:KLEV,IDRAFT) ,&
                    & ZUCURR(KLON)       , ZVCURR(KLON)        , ZTAUX(KLON)    ,&
                    & ZTAUY(KLON)

REAL(KIND=JPRB) ::    ZFRACB(KLON,IDRAFT), ZZPTOP(KLON,IDRAFT) , ZZPLCL(KLON,IDRAFT)
INTEGER(KIND=JPIM) :: IPTOP(KLON,IDRAFT) , IPLCL(KLON,IDRAFT)  , IPLZB(KLON,IDRAFT)  

!RN --------------------------------------------------------------------

REAL(KIND=JPRB) ::    ZZ0MW(KLON)       , ZZ0HW(KLON)       , ZZ0QW(KLON)       ,&
                    & ZBLEND(KLON)      , ZFBLEND(KLON)
REAL(KIND=JPRB) ::    ZZCPTS(KLON)      , ZZQSA(KLON)       , ZZBUOM(KLON)      ,&
                    & ZZZDL(KLON)
REAL(KIND=JPRB) ::    ZTUPD(KLON,KLEV)  , ZQUPD(KLON,KLEV)  , ZLUPD(KLON,KLEV)  ,&
                    & ZIUPD(KLON,KLEV)  , ZQTUPD(KLON,KLEV) , ZLIUPD(KLON,KLEV) ,&
                    & ZSLGUPD(KLON,KLEV), ZAUPD(KLON,KLEV)  
REAL(KIND=JPRB) ::    ZTAUXCG(KLON,KLEV), ZTAUYCG(KLON,KLEV)
REAL(KIND=JPRB) ::    ZSVFLUXCLD(KLON,0:KLEV)               , ZSVFLUXSUB(KLON,0:KLEV),&
                    & ZSVFLUX(KLON,0:KLEV),ZBUOYPOS(KLON)   , ZBUOYNEG(KLON)    ,&
                    & ZDZH(KLON,0:KLEV)

REAL(KIND=JPRB) ::    ZCLDBASE(KLON)    , ZCLDTOP(KLON)     , &
                    & ZRICUI(KLON)      , ZMCU(KLON)        , ZDTHV(KLON)

REAL(KIND=JPRB) ::    ZPFLXUSUM(KLON,0:KLEV)

!RN --- VDF qt variance budget & bimodal cloud scheme variables ---------

REAL(KIND=JPRB) ::    ZVARGEN           , ZVARTRANS         , ZVARDISS          ,&
                    & ZTAU(KLON,KLEV)   , ZTAUNEW           , &
                    & ZDQTDZ(KLON,0:KLEV),ZWQTF             , &
                    & ZWQT2(KLON,0:KLEV), & 
                    & ZQTUP(KLON,KLEV)  , ZSLGUP(KLON,KLEV) , &               
                    & ZQTTEST(KLON,KLEV)   
                    
REAL(KIND=JPRB) ::    ZCLDFRAC(KLON,KLEV),ZQLAV(KLON,KLEV)

!RN ---------------------------------------------------------------------

REAL(KIND=JPRB) ::    ZEXTSHF(KLON)     , ZEXTLHF(KLON)

REAL(KIND=JPRB) ::    ZCPTSTI(KLON,KTILES), ZQSTI(KLON,KTILES)  ,&
                    & ZDQSTI(KLON,KTILES) , ZCSATTI(KLON,KTILES),&
                    & ZCAIRTI(KLON,KTILES), ZCFHTI(KLON,KTILES) ,&
                    & ZCFQTI(KLON,KTILES) , ZAHFLTI(KLON,KTILES),&
                    & ZTSKTIP1(KLON,KTILES)
REAL(KIND=JPRB) ::    ZSTR(KLON,KTILES)   , ZG0(KLON,KTILES)

LOGICAL ::            LLRUNDRY(KLON), LLPBL(KLON,KLEV)

INTEGER(KIND=JPIM) :: ITOP, JD, JK, JL, JT, KCAP

REAL(KIND=JPRB) ::    ZGDPH, ZRHO, ZTMST, ZRG, ZRTMST
LOGICAL ::            LLSFCFLX, LLTERRA
REAL(KIND=JPRB) ::    ZHU1
REAL(KIND=JPRB) ::    ZALFAW(KLON,KLEV)

REAL(KIND=JPRB) ::    ZDETR(KLON,KLEV), ZAWAKE(KLON,KLEV)


!amk ... for moist static energy BL equilibrium M closure
REAL(KIND=JPRB) ::    ZDHPBL(KLON), ZMSCALE(KLON), ZMGEOM, &
                    & ZDT, ZDQ, ZKHFLTOP, ZKQFLTOP, ZSLGENH, ZDMSE, ZMFLXNEW, &
                    & ZMFS(KLON), ZCONS10, ZMFMAX
LOGICAL ::            LMPBLEQU
!xxx
REAL(KIND=JPRB) ::    ZQTENH(KLON,0:KLEV) 

REAL(KIND=JPRB), DIMENSION(KLON,ntiles_total+ntiles_water) ::             &
                    & shfl_soil_t, lhfl_soil_t, shfl_snow_t, lhfl_snow_t, &
                    & shfl_s_t   , lhfl_s_t   , qhfl_s_t                , &
                    & lhfl_bs_t  , rstom_t                   
REAL(KIND=JPRB), DIMENSION(KLON,nlev_soil,ntiles_total+ntiles_water) :: &
  lhfl_pl_t
REAL(KIND=JPRB)       ZTMEAN

!amk testing only!!!
REAL(KIND=JPRB) ::    ZCFNC1  
!xxx

!Parameters for RICO case
REAL(KIND=JPRB), PARAMETER :: c_h = 0.001094_jprb   ! drag coefficient for heat
REAL(KIND=JPRB), PARAMETER :: c_q = 0.001133_jprb   ! drag coefficient for moisture
REAL(KIND=JPRB), PARAMETER :: th0_rico = 298.5_jprb ! SST
REAL(KIND=JPRB), PARAMETER :: psfc = 101540._jprb   ! surface pressure
REAL(KIND=JPRB) mwind, th_l

REAL(KIND=JPRB) ::    ZHOOK_HANDLE

!INTERFACE
! #include "surfexcdriver.h"
! #include "surfpp.h"
!END INTERFACE

! #include "vdfdifh_c.intfb.h"
! #include "vdfdifh.intfb.h"
! #include "vdfdifm_c.intfb.h"
! #include "vdfdifm.intfb.h"
! #include "vdfdifc.intfb.h"
! #include "vdfdpbl.intfb.h"
! #include "vdfexcu.intfb.h"
! #include "vdfhghtn.intfb.h"
! #include "vdfincr.intfb.h"
! #include "vdffblend.intfb.h"
! #include "vdfcloud.intfb.h"
! #include "vdftofdc.intfb.h"

! #include "fcttre.h"  ! replaced by use statement


!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFMAIN',0,ZHOOK_HANDLE)

ZTMST       = PTSPHY
ZRTMST      = 1.0_JPRB/PTSPHY    ! optimization
ZRG         = 1.0_JPRB/RG        !     -"-
LLRUNDRY(:) = .FALSE.  ! option to run dry updrafts with no condensation


!*         1.0  SCM: Fixed fluxes for flux boundary condition ([W/m^2] downward)

IF (LSCMEC) THEN
  LLSFCFLX   = LSFCFLX   ! scm namelist parameters
  ZEXTSHF(:) = REXTSHF
  ZEXTLHF(:) = REXTLHF
ELSE
  LLSFCFLX   = .FALSE.
  ZEXTSHF(:) = 0.0_JPRB
  ZEXTLHF(:) = 0.0_JPRB
ENDIF
!IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
!  LLTERRA = .TRUE.
!ELSE
 LLTERRA = .TRUE.
!ENDIF

!amk  turn on specified surface fluxes everywhere globally
!     (attention: number here and in mo_nwp_conv_interactive.f90)
!LLSFCFLX   = .TRUE.
!ZEXTSHF(:) = -15.0_JPRB
!ZEXTLHF(:) = -60.0_JPRB
!xxx


!*         1.0a First guess of 10m wind (best using a log profile with z0 ???)

DO JL=KIDIA,KFDIA
  PU10M(JL) = PUM1(JL,KLEV)   ! simple guess using lowest level (~10m)
  PV10M(JL) = PVM1(JL,KLEV)   ! attention: u/v10m using in TERRA
ENDDO


!*         1.1  Store initial tendencies for flux calculation
!*              and initialize variable.

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZQEA(JL,JK)=PQE(JL,JK)
    ZLEA(JL,JK)=PLE(JL,JK)
    ZIEA(JL,JK)=PIE(JL,JK)
    ZQTEP(JL,JK)  = 0._JPRB
    ZSLGEP(JL,JK) = 0._JPRB
    PLDIFF(JL,JK) = 0._JPRB
  ENDDO
ENDDO

DO JD=1,IDRAFT
  DO JK=0,KLEV
    DO JL=KIDIA,KFDIA
      ZMFLX(JL,JK,JD)  = 0.0_JPRB
      ZSLGUH(JL,JK,JD) = 0.0_JPRB
      ZQTUH(JL,JK,JD)  = 0.0_JPRB
      ZWUH(JL,JK,JD)   = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO

DO JK=0,KLEV
  DO JL=KIDIA,KFDIA
    PFPLVL(JL,JK) = 0.0_JPRB
    PFPLVN(JL,JK) = 0.0_JPRB
    PFHPVL(JL,JK) = 0.0_JPRB
    PFHPVN(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO



!     ------------------------------------------------------------------

!*         2.     NEW VARIABLES S, SLG, QT
!*                (at initial time level)
!                 ------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

!*         2.1  dry static energy cp(q)*T + gz

    ZCPTGZ(JL,JK)  =PGEOM1(JL,JK) + PTM1(JL,JK)*RCPD * (1.0_JPRB + RVTMP2*PQM1(JL,JK))

!*         2.2  total water and generalized liquid water static energy 
!*              slg = cp*T + gz - Lcond*ql - Ldep*qi

    ZSLGM1(JL,JK) = ZCPTGZ(JL,JK) - RLVTT * PLM1(JL,JK) - RLSTT * PIM1(JL,JK)
    ZSLGE(JL,JK)  = RCPD * ( ( 1.0_JPRB + RVTMP2 * PQM1(JL,JK) ) * PTE(JL,JK)   &  !dcpT/dt
                &                       + RVTMP2 * PTM1(JL,JK)   * PQE(JL,JK) ) &  
                & - RLVTT * PLE(JL,JK) - RLSTT * PIE(JL,JK)                        !dLqli/dt  

    ZQTM1(JL,JK ) = PQM1(JL,JK) + PLM1(JL,JK) + PIM1(JL,JK)
    ZQTE(JL,JK)   = PQE(JL,JK)  + PLE(JL,JK)  + PIE(JL,JK)             !rad+dyn. qt tendency
    
    ZSLGEA(JL,JK) = ZSLGE(JL,JK)
    ZPTEA(JL,JK)  = PTE(JL,JK)
    ZQTEA(JL,JK)  = ZQTE(JL,JK)
    
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         3.  Compute all surface related quantities
!          ------------------------------------------

! initialize skin temperature for surfexcdriver

DO JL=KIDIA,KFDIA
  ztmean = 0.0_jprb
  DO jt=1,ntiles_total
    ztmean = ztmean + t_g_ex(jl,jt) * ext_data%atm%frac_t(jl,jb,jt)
  ENDDO
  IF (SUM(ext_data%atm%frac_t(jl,jb,1:ntiles_total)) > 0.0_JPRB ) THEN
    ztmean = ztmean / SUM(ext_data%atm%frac_t(jl,jb,1:ntiles_total))
  ELSE
    ztmean = t_g_ex(jl,isub_water)       ! set to SST if no land for safety
  ENDIF  

  PTSKTI(JL,1) = PSST(jl)                ! ocean tile  (use SST as input, not skin)
  PTSKTI(JL,2) = t_g_ex(jl,isub_seaice)  ! SEAICE tile
  PTSKTI(JL,3:ntiles_edmf) = ztmean      ! TESSEL/IFS land tiles (take mean land value)
  IF (ext_data%atm%frac_t(jl,jb,isub_water) == 0.0_jprb) THEN
    PTSKTI(JL,1) = t_g_ex(jl,isub_water) ! safety for vupdz0 calculations
  ENDIF
ENDDO

! use sea ice fluxes from previous step for sea ice calculation in the surface interface
DO JL=KIDIA,KFDIA
  shfl_soil_t(jl,isub_seaice) = PAHFSTI(jl,2)
  lhfl_soil_t(jl,isub_seaice) = PEVAPTI(jl,2) * RLSTT
  shfl_snow_t(jl,isub_seaice) = PAHFSTI(jl,2)   !snow over sea ice not used currently
  lhfl_snow_t(jl,isub_seaice) = PEVAPTI(jl,2) * RLSTT
ENDDO

!debug
DO JL=KIDIA,KFDIA
 if ( PTSKM1M(JL) > 400.0 .or. PTSKM1M(JL) < 100.0 ) then
  write(*,*) 'vdfmain1: ', JL, JB, PTSKM1M(JL), PTM1(JL,KLEV), PTM1(JL,KLEV-1)
 endif
ENDDO
!xxxxx

CALL SURFEXCDRIVER( &
  ! TERRA data
   &   ext_data                                                           & !in
   & , jb, jg                                                             & ! -
   & , t_snow_ex, t_snow_mult_ex, t_s_ex, t_g_ex, qv_s_ex                 & !inout
   & , w_snow_ex                                                          & ! -
   & , rho_snow_ex, rho_snow_mult_ex, h_snow_ex, w_i_ex, w_p_ex, w_s_ex   & ! -
   & , t_so_ex, w_so_ex, w_so_ice_ex                                      & ! -
   & , PU10M, PV10M                    &  !, t_2m_ex, u_10m_ex, v_10m_ex  & ! -
   & , freshsnow_ex, snowfrac_lc_ex, snowfrac_ex                          & ! -
   & , wliq_snow_ex, wtot_snow_ex, dzh_snow_ex                            & ! -
   & , prr_con_ex, prs_con_ex, prr_gsp_ex, prs_gsp_ex                     & !in
   & , tch_ex, tcm_ex, tfv_ex                                             & !inout
   & , sobs_ex, thbs_ex, pabs_ex                                          & !in
   & , runoff_s_ex, runoff_g_ex                                           & !inout
   & , t_g, qv_s                                                          & ! -
   & , t_ice, h_ice, t_snow_si, h_snow_si, alb_si                         & ! -
   & , fr_seaice                                                          & !in
   & , shfl_soil_t, lhfl_soil_t, shfl_snow_t, lhfl_snow_t                 & !out
   & , shfl_s_t   , lhfl_s_t   , qhfl_s_t                                 &
   & , lhfl_bs_t  , lhfl_pl_t  , rstom_t    ,                             & !out
  ! standard input
   & CDCONF=CDCONF, &
   & KIDIA=KIDIA, KFDIA=KFDIA, KLON=KLON, KLEVS=KLEVS, KTILES=KTILES, KSTEP=KSTEP, &
   & KLEVSN=KLEVSN, KLEVI=KLEVI, KDHVTLS=KDHVTLS, KDHFTLS=KDHFTLS, &
   & KDHVTSS=KDHVTSS, KDHFTSS=KDHFTSS, KDHVTTS=KDHVTTS, KDHFTTS=KDHFTTS, &
   & KDHVTIS=KDHVTIS, KDHFTIS=KDHFTIS, K_VMASS=N_VMASS, &
   & PTSTEP=PTSPHY, PRVDIFTS=RVDIFTS, &
  ! input data, non-tiled
   & KTVL=KTVL, KTVH=KTVH, PCVL=PCVL, PCVH=PCVH, &
   & PUMLEV=PUM1(:,KLEV), PVMLEV=PVM1(:,KLEV), PTMLEV=PTM1(:,KLEV), &
   & PQMLEV=PQM1(:,KLEV), PAPHMS=PAPHM1(:,KLEV), PAPMS=PAPM1(:,KLEV), PGEOMLEV=PGEOM1(:,KLEV), &
   & PCPTGZLEV=ZCPTGZ(:,KLEV), PSST=PSST, PTSKM1M=PTSKM1M, PCHAR=PCHAR, &
   & PSSRFL=PSSRFL, PSLRFL=PSLRFL, PEMIS=PEMIS, PTICE=PTICE, PTSNOW=PTSNOW, &
   & PHLICE=PHLICE,PTLICE=PTLICE,PTLWML=PTLWML, & 
   & PWLMX=PWLMX, PUCURR=PUCURR, PVCURR=PVCURR, &
  ! input data, soil
   & PTSAM1M=PTSAM1M, PWSAM1M=PWSAM1M, KSOTY=KSOTY, &
  ! input data, tiled
   & PFRTI=PFRTI, PALBTI=PALBTI, &
  ! updated data, tiled
   & PUSTRTI=PUSTRTI, PVSTRTI=PVSTRTI, PAHFSTI=PAHFSTI, PEVAPTI=PEVAPTI, &
   & PTSKTI=PTSKTI, &
  ! updated data, non-tiled
   & PZ0M=PZ0M, PZ0H=PZ0H, &
  ! output data, tiled
   & PSSRFLTI=PSSRFLTI, PQSTI=ZQSTI, PDQSTI=ZDQSTI, PCPTSTI=ZCPTSTI, &
   & PCFHTI=ZCFHTI, PCFQTI=ZCFQTI, PCSATTI=ZCSATTI, PCAIRTI=ZCAIRTI, &
  ! output data, non-tiled
   & PKHLEV=PKH(:,KLEV), PCFMLEV=ZCFM(:,KLEV), PKMFL=ZKMFL, PKHFL=ZKHFL, &
   & PKQFL=ZKQFL, PEVAPSNW=PEVAPSNW, PZ0MW=ZZ0MW, PZ0HW=ZZ0HW, PZ0QW=ZZ0QW, &
   & PBLENDPP=ZBLEND, PCPTSPP=ZZCPTS, PQSAPP=ZZQSA, PBUOMPP=ZZBUOM, &
   & PZDLPP=ZZZDL )
! output data, diagnostics
!dmk   & PDHTLS=PDHTLS, PDHTSS=PDHTSS, PDHTTS=PDHTTS, PDHTIS=PDTHIS &


!amk  overwrite SCM surface fluxes from above calculation
!     (attention: number here and in mo_nwp_conv_interactive.f90)
DO JL=KIDIA,KFDIA

! surface fluxes from TERRA
!      should really be handed over for each tile not mean ???????

  ZEXTSHF(JL) = 0.0_JPRB
  ZEXTLHF(JL) = 0.0_JPRB
  DO JT=1,ntiles_total      ! land points only
    IF (LLTERRA) THEN
      ZEXTSHF(JL) = ZEXTSHF(JL) + ext_data%atm%frac_t(JL,JB,JT) *    &
        ( SHFL_SOIL_T(JL,JT) * (1.0_JPRB - SNOWFRAC_EX(JL,JT)) +  &
          SHFL_SNOW_T(JL,JT) *             SNOWFRAC_EX(JL,JT)  )
      ZEXTLHF(JL) = ZEXTLHF(JL) + ext_data%atm%frac_t(JL,JB,JT) *    &
        ( LHFL_SOIL_T(JL,JT) * (1.0_JPRB - SNOWFRAC_EX(JL,JT)) +  &
          LHFL_SNOW_T(JL,JT) *             SNOWFRAC_EX(JL,JT)  )
    END IF
  ENDDO

! update skin temperature from TERRA and SEAICE

  ztmean = 0.0_jprb
  DO jt=1,ntiles_total
    ztmean = ztmean + t_g_ex(jl,jt) * ext_data%atm%frac_t(jl,jb,jt)
  ENDDO
  IF (SUM(ext_data%atm%frac_t(jl,jb,1:ntiles_total)) > 0.0_JPRB ) THEN
    ztmean = ztmean / SUM(ext_data%atm%frac_t(jl,jb,1:ntiles_total))
  ELSE
    ztmean = t_g_ex(jl,isub_water)       ! set to SST if no land for safety
  ENDIF

  PTSKTI(JL,1) = PSST(jl)                ! ocean tile  (use SST as input, not skin)
  PTSKTI(JL,2) = t_g_ex(jl,isub_seaice)  ! SEAICE tile
  PTSKTI(JL,3:ntiles_edmf) = ztmean      ! TESSEL/IFS land tiles (take mean land value)
  IF (ext_data%atm%frac_t(jl,jb,isub_water) == 0.0_jprb) THEN
    PTSKTI(JL,1) = t_g_ex(jl,isub_water) ! safety for vdfdifh calculations
  ENDIF

! sea ice fluxes (mean flux stored in '_soil' flux)
!  JT = isub_seaice
!  ZEXTSHF(JL) = ZEXTSHF(JL) + ext_data%atm%frac_t(JL,JB,JT) *    &
!    SHFL_SOIL_T(JL,JT)
!  ZEXTLHF(JL) = ZEXTLHF(JL) + ext_data%atm%frac_t(JL,JB,JT) *    &
!    LHFL_SOIL_T(JL,JT)

! normalize to get mean over land points
!  ZEXTSHF(JL) = ZEXTSHF(JL) / SUM(ext_data%atm%frac_t(JL,JB,1:ntiles_total))
!  ZEXTLHF(JL) = ZEXTLHF(JL) / SUM(ext_data%atm%frac_t(JL,JB,1:ntiles_total))

! test: fixed surface fluxes
   !ZEXTSHF(JL) = -15.0_JPRB
   !ZEXTLHF(JL) = -60.0_JPRB
   !ZKMFL(JL)   = 0.0_JPRB  ! - " -      (0 is bad idea!!!)

! optional SCM case definition after Stevens(2007) for dry BL

  IF (ltestcase) THEN
    ZRHO = PAPHM1(JL,KLEV)/( RD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV)) )
    SELECT CASE(nh_test_name)
    CASE ('CBL')
      CALL cbl_stevens_fluxes( PTM1(JL,KLEV), PQM1(JL,KLEV), PAPM1(JL,KLEV), &
                             & ZRHO         , T_G(JL)      , &
                             & ZEXTSHF(JL)  , ZEXTLHF(JL)  )
      LLSFCFLX    = .TRUE.
    !SBr,CHa: catchment-scale circulation case  
    CASE ('CBL_ccs')
      CALL cbl_stevens_fluxes( PTM1(JL,KLEV), PQM1(JL,KLEV), PAPM1(JL,KLEV), &
                             & ZRHO         , T_G(JL)      , &
                             & ZEXTSHF(JL)  , ZEXTLHF(JL)  )
      LLSFCFLX    = .TRUE.
    CASE ('RICO')
      mwind = SQRT(PUM1(JL,KLEV)**2 + PVM1(JL,KLEV)**2) 
      th_l  = PTM1(JL,KLEV) * (p0ref/PAPM1(JL,KLEV))**rd_o_cpd  
      ZEXTSHF(JL) = - ZRHO * RCPD  * c_h * mwind * (th0_rico - th_l)
      ZEXTLHF(JL) = - ZRHO * RLVTT * c_q * mwind * (qv_s(JL) - PQM1(JL,KLEV))
      LLSFCFLX    = .TRUE.
    END SELECT
  ENDIF

ENDDO

!          FLUX BOUNDARY CONDITION

DO JL=KIDIA,KFDIA  
  IF (LDLAND(JL) .OR. LLSFCFLX) THEN
    ZRHO = PAPHM1(JL,KLEV)/( RD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV)) )
    ZKHFL(JL) = ZEXTSHF(JL) / ZRHO / RCPD
    ZKQFL(JL) = ZEXTLHF(JL) / ZRHO / RLVTT
   !ZKMFL(JL) = calculated as mean of TESSEL tiles in surfexcdriver
    IF (LLSFCFLX) THEN      ! ATTENTION: surface momentum flux
      ZKMFL(JL) = 0.0_JPRB  ! needs to be calculated in surfexcdriver!!!
    END IF
  END IF
ENDDO


!     ------------------------------------------------------------------

!*         4.     EXCHANGE COEFFICIENTS
!                 ---------------------

!*         4.4  COMPUTATION OF THE PBL EXTENSION


!          SET PBL HEIGHT-INDEX TO 1

DO JL=KIDIA,KFDIA
  KHPBLN(JL)=1
  KVARTOP(JL)  = 0
ENDDO
ITOP=1  !ITOP is used in some solvers. ITOP=1 means: always integrate over whole atmosphere depth


!*         4.5  BOUNDARY LAYER HEIGHT FOR DIANOSTICS ONLY

CALL VDFDPBL(KIDIA,KFDIA,KLON,KLEV,&
 & PUM1,PVM1,PTM1,PQM1,PGEOM1,&
 & ZKMFL,ZKHFL,ZKQFL,PBLH)  


!*         4.6  ORGANIZED UPDRAFTS

CALL VDFHGHTN (KIDIA   , KFDIA   , KLON    , KLEV    , IDRAFT   , ZTMST   , KSTEP   , &
             & PUM1    , PVM1    , PTM1    , PQM1    , PLM1     , PIM1    , PAM1    , &
             & PAPHM1  , PAPM1   , PGEOM1  , PGEOH   , PVERVEL  , PQE     , PTE     , &
             & ZKMFL   , ZKHFL   , ZKQFL   , ZMFLX   , &
! DIAGNOSTIC OUTPUT
             & PEXTR2  , KFLDX2  , PEXTRA  , KLEVX   , KFLDX    , KCNT    , LLDIAG  , &
!              
             & ZUUH    , ZVUH    , ZSLGUH  , ZQTUH   , ZFRACB   , ZWUH    , &
             & ZZPTOP  , IPTOP   , ZZPLCL  , IPLCL   , IPLZB    , &
             & PWUAVG  , ZRICUI  , ZMCU    , ZDTHV   , &
             & PFPLVL  , PFPLVN  , ZDETR   , &
!amk: for convective preconditioning
             & PVAR    , &
             & LDLAND  , &   
             & PBIR    , LDNODECP, LLRUNDRY, KPBLTYPE, ZWQT2 )


!--- set some PBL heights based on i) PBL type and ii) various updraft properties ---

DO JL=KIDIA,KFDIA
    
  SELECT CASE (KPBLTYPE(JL))
    
      CASE(0)
          !Stable PBL
          PZINV(JL)    = 0.0_JPRB
          KHPBLN(JL)   = KLEV
          ZCLDBASE(JL) = -100._JPRB
          ZCLDTOP(JL)  = -100._JPRB
          KVARTOP(JL)  = 0
     
      CASE(1)
          !Dry convective PBL
          PZINV(JL)    = ZZPTOP(JL,2)
          KHPBLN(JL)   = IPTOP(JL,2)
          ZCLDBASE(JL) = -100._JPRB
          ZCLDTOP(JL)  = -100._JPRB
          KVARTOP(JL)  = IPTOP(JL,2)

      CASE(2)
          !Stratocumulus
          PZINV(JL)    = ZZPTOP(JL,3)
          KHPBLN(JL)   = IPTOP(JL,3)
          ZCLDBASE(JL) = ZZPLCL(JL,3)
          ZCLDTOP(JL)  = ZZPTOP(JL,3)
          KVARTOP(JL)  = IPTOP(JL,3)
          
      CASE(3)
          !Shallow cumulus
          PZINV(JL)    = ZZPLCL(JL,3)
          KHPBLN(JL)   = IPLCL(JL,3)+1   !equivalent to definition of top level: use first half level *BELOW* boundary
          ZCLDBASE(JL) = ZZPLCL(JL,3)
          ZCLDTOP(JL)  = ZZPTOP(JL,3)
          KVARTOP(JL)  = IPTOP(JL,1)
                
      CASE(4)
          !Deep cumulus - only do a dry subcloud ML
          PZINV(JL)    = ZZPTOP(JL,2)
          KHPBLN(JL)   = IPTOP(JL,2)
          ZCLDBASE(JL) = ZZPLCL(JL,1)
          ZCLDTOP(JL)  = ZZPTOP(JL,1)
          KVARTOP(JL)  = IPTOP(JL,2)
    
  END SELECT 
  
ENDDO !JL

IF ( LLDIAG ) THEN
  DO JL=KIDIA,KFDIA
        PEXTRA(JL,11,41) = 0.0_JPRB
        PEXTRA(JL,12,41) = 0.0_JPRB
        PEXTRA(JL,13,41) = 0.0_JPRB
        PEXTRA(JL,14,41) = 0.0_JPRB
        PEXTRA(JL,15,41) = 0.0_JPRB
        PEXTRA(JL,17,41) = 0.0_JPRB
        PEXTRA(JL,18,41) = 0.0_JPRB
        PEXTRA(JL,19,41) = 0.0_JPRB
        PEXTRA(JL,20,41) = 0.0_JPRB
    SELECT CASE (KPBLTYPE(JL))
      CASE (0)
        PEXTRA(JL, 5,41) = PEXTRA(JL, 5,41) + ZTMST
        PEXTRA(JL,15,41) = 1.0_JPRB
      CASE (1)
        PEXTRA(JL, 1,41) = PEXTRA(JL, 1,41) + ZTMST
        PEXTRA(JL,11,41) = 1.0_JPRB
      CASE (2)
        PEXTRA(JL, 2,41) = PEXTRA(JL, 2,41) + ZTMST
        PEXTRA(JL,12,41) = 1.0_JPRB
        PEXTRA(JL, 7,41) = PEXTRA(JL, 7,41) + ZCLDBASE(JL) * ZTMST
        PEXTRA(JL, 8,41) = PEXTRA(JL, 8,41) + ZCLDTOP(JL)  * ZTMST
        PEXTRA(JL,17,41) = ZCLDBASE(JL)
        PEXTRA(JL,18,41) = ZCLDTOP(JL)
      CASE (3)
        PEXTRA(JL, 3,41) = PEXTRA(JL, 3,41) + ZTMST
        PEXTRA(JL,13,41) = 1.0_JPRB
        PEXTRA(JL, 9,41) = PEXTRA(JL, 9,41) + ZCLDBASE(JL) * ZTMST
        PEXTRA(JL,10,41) = PEXTRA(JL,10,41) + ZCLDTOP(JL)  * ZTMST
        PEXTRA(JL,19,41) = ZCLDBASE(JL)
        PEXTRA(JL,20,41) = ZCLDTOP(JL)
      CASE (4)
        PEXTRA(JL, 4,41) = PEXTRA(JL, 4,41) + ZTMST
        PEXTRA(JL,14,41) = 1.0_JPRB
    END SELECT
        PEXTRA(JL, 6,41) = PEXTRA(JL, 6,41) + PZINV(JL)    * ZTMST
        PEXTRA(JL,21,41) = REAL( KPBLTYPE(JL), JPRB )
        PEXTRA(JL,16,41) = PZINV(JL)
   ENDDO
ENDIF


!--- set PBL indicator ---   
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    LLPBL(JL,JK) = JK >= KVARTOP(JL) .AND. KVARTOP(JL)>0
  ENDDO
ENDDO


!*         4.7  EXCHANGE COEFFICIENTS ABOVE THE SURFACE LAYER

CALL VDFEXCU(KIDIA  , KFDIA  , KLON   , KLEV    , IDRAFT  , ZTMST  , PZ0M   , &
           & PHRLW  , PHRSW  , &
           & PUM1   , PVM1   , PTM1   , PQM1    , PLM1    , PIM1   , &
           & PAPHM1 , PAPM1  , PGEOM1 , PGEOH   , ZCPTGZ  , &
! DIAGNOSTIC OUTPUT
            &PEXTR2 , KFLDX2 , PEXTRA , KLEVX   , KFLDX   , &
!              
           & ZKMFL  , ZKHFL  , ZKQFL  , ZCFM    , ZCFH    , ZTAUXCG, ZTAUYCG, &
           & ZRICUI , ZMCU   , ZDTHV  , ZMFLX   , KVARTOP , &
           & PZINV  , KHPBLN , PKH    , ZCLDBASE, ZCLDTOP , KPBLTYPE)  


!amk

!*         4.7.1  Moist static energy MSE mass-flux closure (diagnostic only)
 
  ! turbulent fluxes (surface and cloud base, upward definition)

  DO JL=KIDIA,KFDIA

    ZRHO       = PAPHM1(JL,KLEV-1)/( RD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV)) )
    ZDHPBL(JL) = ZRHO * ( - RCPD * ZKHFL(JL) - RLVTT * ZKQFL(JL) )            !surface MSE flux

    IF ( LLDIAG ) THEN
      PEXTRA(JL,28,41) = ZDHPBL(JL)                                           !surface MSE flux
    ENDIF

!xmkIF ( KPBLTYPE(JL) .NE. 0 ) THEN
    IF ( KHPBLN(JL) .LT. KLEV ) THEN 
      JK         = KHPBLN(JL)
      ZMGEOM     = PGEOM1(JL,JK) - PGEOM1(JL,JK+1)
      ZDT        =(ZCPTGZ(JL,JK) - ZCPTGZ(JL,JK+1)) * (1.0_JPRB/RCPD)  
      ZDQ        = PQM1(JL,JK)   - PQM1(JL,JK+1)
      ZKHFLTOP   = - PKH(JL,JK) * ZDT * RG / ZMGEOM
      ZKQFLTOP   = - PKH(JL,JK) * ZDQ * RG / ZMGEOM 
      ZRHO       = PAPHM1(JL,JK)/( RD*PTM1(JL,JK)*(1.0_JPRB+RETV*PQM1(JL,JK)) )
      ZDHPBL(JL) = ZDHPBL(JL) &
               & + ZRHO * ( - RCPD * ZKHFLTOP  - RLVTT * ZKQFLTOP  )          !LCL MSE flux (entrainment)
    ENDIF

    IF ( LLDIAG ) THEN
      PEXTRA(JL,29,41) = ZDHPBL(JL)                                           !sfc+top
    ENDIF
  ENDDO


  ! integrate dynamical and radiative tendencies (attention: ICON needs that input!!!)

  DO JK=KLEV,1,-1
    DO JL=KIDIA,KFDIA
      IF ( JK >= KHPBLN(JL) ) THEN
        ZDHPBL(JL) = ZDHPBL(JL) + ( RLVTT*PQE(JL,JK) + RCPD*PTE(JL,JK) ) &
                              & * ( PAPHM1(JL,JK) - PAPHM1(JL,JK-1) ) / RG
      ENDIF  
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    IF ( LLDIAG ) THEN
      PEXTRA(JL,30,41) = ZDHPBL(JL)                                           !sfc+top+rad+dyn
    ENDIF
  ENDDO
!xxx


!*         4.8  MASS FLUX MODIFICATIONS (FOR SOLVER STABILITY)

DO JD=2,IDRAFT

  !-- Remove single and/or double massflux layers --
  DO JL=KIDIA,KFDIA
!    IF ( ZMFLX(JL,KLEV-2,JD) < 1.E-40_JPRB ) THEN
    IF ( ZMFLX(JL,KLEV-3,JD) < 1.E-40_JPRB ) THEN
      ZMFLX(JL,KLEV-1,JD) = 0._JPRB
      ZMFLX(JL,KLEV-2,JD) = 0._JPRB
    ENDIF
  ENDDO
  
  !-- Prune massflux-spikes in top PBL layer --
  DO JK=2,KLEV-2
  DO JL=KIDIA,KFDIA
      IF ( ZMFLX(JL,JK-1,JD) < 1.E-40_JPRB .AND. ZMFLX(JL,JK,JD) > 1.E-40_JPRB ) THEN
        ZMFLX(JL,JK,JD) = MIN( ZMFLX(JL,JK,JD) , 1.5_JPRB*ZMFLX(JL,JK+1,JD) )
      ENDIF
  ENDDO
  ENDDO
  
ENDDO


!*         4.9     TURBULENT OROGRAPHIC DRAG COEFFICIENTS

CALL VDFTOFDC(KIDIA,KFDIA,KLON,KLEV,ZTMST,&
 & PUM1,PVM1,PGEOM1,PSIGFLT,&
 & ZTOFDC)  



!     ------------------------------------------------------------------

!*         5.     SOLVE ADVECTION-DIFFUSION EQUATION
!                 ----------------------------------

!*         5.1  MOMENTUM

ZMFLXM(KIDIA:KFDIA,:,:)=ZMFLX(KIDIA:KFDIA,:,:)
ZUCURR(:)    =0._JPRB   ! ocean currents not yet active
ZVCURR(:)    =0._JPRB   ! ...
ZHU1=RVDIFTS*ZTMST
ZSOC(KIDIA:KFDIA,1:KLEV)=PSOBETA(KIDIA:KFDIA,1:KLEV)*ZHU1

!...centered
!CALL VDFDIFM_C (KIDIA,KFDIA, KLON , KLEV  , IDRAFT , ITOP  , &
!              & ZTMST, PUM1 , PVM1 , PAPHM1, ZCFM   , ZMFLXM, ZUUH , ZVUH ,&
!              & ZTOFDC,PSOTEU,PSOTEV,ZSOC  , &
!              & PVOM , PVOL , ZUCURR,ZVCURR, ZUDIF  , ZVDIF , ZTAUX, ZTAUY)

!debug
!DO JL=KIDIA,KFDIA
!  if (zuuh(jl,klev-1,1) > 100.0  .or. zvuh(jl,klev-1,1)   > 100.0 .or. &
!      zcfm(jl,klev)     > 1000.0 .or. zmflxm(jl,klev-1,1) > 1000.0  ) THEN
!    write(*,*) 'vdfmain before vdfdifm', zuuh(jl,klev-1,1), zvuh(jl,klev-1,1), &
!      zcfm(jl,klev), zmflxm(jl,klev-1,1)
!  endif
!ENDDO
!xxxxx

!...fully upstream M (phi_up - phi_bar) term
CALL VDFDIFM (KIDIA, KFDIA, KLON , KLEV  , IDRAFT , ITOP  , &
            & ZTMST, PUM1 , PVM1 , PAPHM1, ZCFM   , ZMFLXM, ZUUH , ZVUH ,&
            & ZTOFDC,PSOTEU,PSOTEV,ZSOC  , &
            & PVOM , PVOL , ZUCURR,ZVCURR, ZUDIF  , ZVDIF , ZTAUX, ZTAUY)


!*         5.2  GENERALIZED LIQUID WATER STATIC ENERGY AND TOTAL WATER

!...centered
!CALL VDFDIFH_C (KIDIA  , KFDIA  , KLON   , KLEV   , IDRAFT , ITOP   , KTILES, &
!            & ZTMST  , ZEXTSHF, ZEXTLHF, LLSFCFLX, &
!            & PFRTI  , PSSRFLTI,PSLRFL , PEMIS  , PEVAPSNW, &
!            & PHLICE , PTLICE , PTLWML , &
!            & ZSLGM1 , PTM1   , PQM1   , ZQTM1  , PAPHM1 , &
!            & ZCFH   , ZCFHTI , ZCFQTI , ZMFLX  , ZSLGUH , ZQTUH  , &
!            & ZSLGDIF, ZQTDIF , ZCPTSTI, ZQSTI  , ZCAIRTI, ZCSATTI, &
!            & ZDQSTI , PTSKTI , PTSKRAD, PTSAM1M(1,1)    , PTSNOW , PTICE  , PSST, &
!            & ZTSKTIP1,ZSLGE  , PTE    , ZQTE, &
!            & PEVAPTI, PAHFSTI, ZAHFLTI, ZSTR   , ZG0)

!...fully upstream M (phi_up - phi_bar) term
CALL VDFDIFH (KIDIA  , KFDIA  , KLON   , KLEV   , IDRAFT , ITOP   , KTILES, &
            & ZTMST  , ZEXTSHF, ZEXTLHF, LLSFCFLX,LLTERRA, LDLAND , &
            & PFRTI  , PSSRFLTI,PSLRFL , PEMIS  , PEVAPSNW, &
            & PHLICE , PTLICE , PTLWML , &
            & ZSLGM1 , PTM1   , PQM1   , ZQTM1  , PAPHM1 , &
            & ZCFH   , ZCFHTI , ZCFQTI , ZMFLX  , ZSLGUH , ZQTUH  , &
            & ZSLGDIF, ZQTDIF , ZCPTSTI, ZQSTI  , ZCAIRTI, ZCSATTI, &
            & ZDQSTI , PTSKTI , PTSKRAD, PTSAM1M(1,1)    , PTSNOW , t_ice , PSST, &
            & ZTSKTIP1,ZSLGE  , PTE    , ZQTE, &
            & PEVAPTI, PAHFSTI, ZAHFLTI, ZSTR   , ZG0)

!DO JL=KIDIA,KFDIA
!  IF ( (SUM(PAHFSTI(JL,:)) == 0.0) .or. (SUM(PEVAPTI(JL,:)) == 0.0) .or. &
!       (PDIFTS(JL,KLEV)    == 0.0) .or. (PDIFTQ(JL,KLEV)    == 0.0) ) THEN
!    write(*,*) 'vdfmain4: ', PDIFTS(JL,KLEV), PDIFTQ(JL,KLEV), PAHFSTI(JL,:), PEVAPTI(JL,:)*RLVTT, &
!      & PFRTI(JL,:)
!  ENDIF
!ENDDO


!*         5.3  INCREMENTATION OF U AND V TENDENCIES, STORAGE OF
!*              THE DISSIPATION, COMPUTATION OF MULTILEVEL FLUXES.

CALL VDFINCR (KIDIA  , KFDIA  , KLON   , KLEV   , ITOP   , ZTMST  , &
            & PUM1   , PVM1   , ZSLGM1 , PTM1   , ZQTM1  , PAPHM1 , PGEOM1 , &
            & ZCFM   , ZTOFDC , PSOTEU , PSOTEV , ZSOC   ,&
            & ZUDIF  , ZVDIF  , PUCURR , PVCURR ,ZSLGDIF, ZQTDIF , &
            & PVOM   , PVOL   , ZSLGE  , ZQTE   , ZSLGEWODIS, &
            & PVDIS  , PVDISG , PSTRTU , PSTRTV , PSTRSOU, PSTRSOV , PTOFDU , PTOFDV, & 
! DIAGNOSTIC OUTPUT
            & PEXTR2 , KFLDX2 , PEXTRA , KLEVX  , KFLDX  , KCNT   , LLDIAG,&
            & PDISGW3D)  


!          5.4  Solve for tracers

IF (LVDFTRAC .AND. KTRAC > 0) THEN 
  CALL VDFDIFC(KIDIA,KFDIA,KLON,KLEV,ITOP,KTRAC,&
             & ZTMST,PCM1,PTENC,PAPHM1,ZCFH,PCFLX)
ENDIF



!     ------------------------------------------------------------------
!
!*        6.     TIME INTEGRATION OF CONSERVED STATE VARIABLES QT AND SLG
!                 --------------------------------------------------------
!
!         Compute the conserved state after rad+dyn *AND* pbl conv+diff.
!         This will be used later to obtain the tendencies of the non-conserved 
!         prognostic variables T and QV.
!
  
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
      
      !calculate updraft precipitation flux divergence (tendencies)
      ZDZRHOI = RG/( PAPHM1(JL,JK)-PAPHM1(JL,JK-1) )
      ZQTEP(JL,JK)  = -( PFPLVL(JL,JK) - PFPLVL(JL,JK-1) ) * ZDZRHOI &
                    & -( PFPLVN(JL,JK) - PFPLVN(JL,JK-1) ) * ZDZRHOI
      ZSLGEP(JL,JK) =   RLVTT * ( PFPLVL(JL,JK) - PFPLVL(JL,JK-1) ) * ZDZRHOI &
                    & + RLSTT * ( PFPLVN(JL,JK) - PFPLVN(JL,JK-1) ) * ZDZRHOI
    ENDDO
  ENDDO


  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      
      !account for updraft precipitation tendencies
      ZQTE(JL,JK)  = ZQTE(JL,JK)  + ZQTEP(JL,JK)
      ZSLGE(JL,JK) = ZSLGE(JL,JK) + ZSLGEP(JL,JK)

      !integrate in time
      ZQTUPD(JL,JK)  = ZQTM1(JL,JK)  + ZQTE(JL,JK)  * ZTMST
      ZSLGUPD(JL,JK) = ZSLGM1(JL,JK) + ZSLGE(JL,JK) * ZTMST
      
      !total specific humidity limiter --- ATTENTION: non-conserving
      ZQTUPD(JL,JK) = MAX( 0._JPRB, ZQTUPD(JL,JK))
      ZQTE(JL,JK)   = (ZQTUPD(JL,JK) - ZQTM1(JL,JK) ) * ZRTMST
      
    ENDDO
  ENDDO



!     ------------------------------------------------------------------

!*         7.     SURFACE FLUXES - TILES
!                 ----------------------
!*         AND    COMPUTE 2M TEMPERATURE AND HUMIDITY, 10M WIND,
!*                  and gustiness

!  Compute wind speed at blending height

CALL VDFFBLEND(KIDIA,KFDIA,KLON,KLEV, &
 & PUM1, PVM1, PGEOM1, PUCURR, PVCURR, ZBLEND, &
 & ZFBLEND)

! Wrap-up computations for the surface and 2T/2D/10U/10V/gustiness computation

!amk  ATTENTION: Surface fluxes need to be aggregated from tiles in SURFPP!!!
!PDIFTS(KIDIA:KFDIA,KLEV) = 0.0_JPRB
!PDIFTQ(KIDIA:KFDIA,KLEV) = 0.0_JPRB
!DO JT=1,KTILES
!  DO JL=KIDIA,KFDIA
!    PDIFTS(JL,KLEV) = PDIFTS(JL,KLEV) + PFRTI(JL,JT) * PAHFSTI(JL,JT)
!    PDIFTQ(JL,KLEV) = PDIFTQ(JL,KLEV) + PFRTI(JL,JT) * PEVAPTI(JL,JT)
!  ENDDO
!ENDDO
!xxx

!debug
DO JL=KIDIA,KFDIA
 if ( PTSKM1M(JL) > 400.0 .or. PTSKM1M(JL) < 100.0  ) then
  write(*,*) 'vdfmain3: ', PTSKM1M(JL), PTM1(JL,KLEV), PTM1(JL,KLEV-1), PFRTI(JL,1)
 endif
ENDDO
!xxxxx

CALL SURFPP( KIDIA=KIDIA,KFDIA=KFDIA,KLON=KLON,KTILES=KTILES, &
 & KDHVTLS=KDHVTLS,KDHFTLS=KDHFTLS, &
 & PTSTEP=PTSPHY, &
! input
 & PFRTI=PFRTI, PAHFLTI=ZAHFLTI, PG0TI=ZG0, &
 & PSTRTULEV=PSTRTU(:,KLEV), PSTRTVLEV=PSTRTV(:,KLEV), PTSKM1M=PTSKM1M, &
 & PUMLEV=PUM1(:,KLEV), PVMLEV=PVM1(:,KLEV), PQMLEV=PQM1(:,KLEV), &
 & PGEOMLEV=PGEOM1(:,KLEV), PCPTSPP=ZZCPTS, PCPTGZLEV=ZCPTGZ(:,KLEV), &
 & PAPHMS=PAPHM1(:,KLEV), PZ0MW=ZZ0MW, PZ0HW=ZZ0HW, PZ0QW=ZZ0QW, &
 & PZDL=ZZZDL, PQSAPP=ZZQSA, PBLEND=ZBLEND, PFBLEND=ZFBLEND, PBUOM=ZZBUOM, &
 & PZ0M=PZ0M, PEVAPSNW=PEVAPSNW,PSSRFLTI=PSSRFLTI, PSLRFL=PSLRFL, PSST=PSST, &
 & PUCURR=PUCURR, PVCURR=PVCURR, &
! updated
 & PAHFSTI=PAHFSTI, PEVAPTI=PEVAPTI, PTSKE1=PTSKE1,PTSKTIP1=ZTSKTIP1, &
! output
 & PDIFTSLEV=PDIFTS(:,KLEV), PDIFTQLEV=PDIFTQ(:,KLEV), PUSTRTI=PUSTRTI, &
 & PVSTRTI=PVSTRTI,  PTSKTI=PTSKTI, PAHFLEV=PAHFLEV, PAHFLSB=PAHFLSB, &
 & PFWSB=PFWSB, PU10M=PU10M, PV10M=PV10M, PT2M=PT2M, PD2M=PD2M, PQ2M=PQ2M, &
 & PGUST=PGUST, &
! output DDH
 & PDHTLS=PDHTLS &
 & )

! store skin temperature in t_g_ex for next step

DO JL=KIDIA,KFDIA
 !ztmean = 0.0_jprb
 !DO jt=3,ntiles_edmf
 !  ztmean = ztmean + PTSKTI(jl,jt) * PFRTI(jl,jt)   ! why not use ZTSKTIP1 (new)?
 !ENDDO
 !IF (SUM(PFRTI(jl,3:ntiles_edmf)) > 0.0_JPRB ) THEN
 !  ztmean = ztmean / SUM(PFRTI(jl,3:ntiles_edmf))
 !ELSE
 !  ztmean = PTSKTI(jl,1)                  ! set to SST if no land for safety
 !ENDIF

! ATTENTION: overwriting t_g_ex highly dangerous!!!!!

  t_g_ex(jl,isub_water)     = PTSKTI(jl,1) !ocean tiles (includes cold skin ..., NOT SST) ???
 !t_g_ex(jl,isub_seaice)    = PTSKTI(jl,2) !sea ice tiles
 !t_g_ex(jl,isub_seaice)    = t_g_ex(jl,isub_seaice) !sea ice tiles
 !t_g_ex(jl,1:ntiles_total) = ztmean       !TERRA tiles (mean)
ENDDO

!DO JL=KIDIA,KFDIA
!  IF ( ABS(PDIFTS(JL,KLEV))          > 1000.0_JPRB  .OR. & 
!       ABS(PDIFTQ(JL,KLEV) * RLVTT ) > 2000.0_JPRB ) THEN
!    write(*,*) 'vdfmain: SHF, LHF ', PDIFTS(JL,KLEV), PDIFTQ(JL,KLEV) * RLVTT 
!  ENDIF
!ENDDO

!amk  diagnostic values have to be set !!!!
! DO JL=KIDIA,KFDIA
!   PU10M(JL) = 0.0_JPRB
!   PV10M(JL) = 0.0_JPRB
!   PT2M (JL) = 0.0_JPRB
!   PD2M (JL) = 0.0_JPRB
!   PQ2M (JL) = 0.0_JPRB
! ENDDO
!xxx

!amk dummy (unused!!)
! DO JL=KIDIA,KFDIA
!   PAHFLEV(JL) = 0.0_JPRB
!   PAHFLSB(JL) = 0.0_JPRB
!   PFWSB(JL)   = 0.0_JPRB
! ENDDO
!xxx

PDIFTL  (KIDIA:KFDIA,KLEV) = 0.0_JPRB
PDIFTI  (KIDIA:KFDIA,KLEV) = 0.0_JPRB
ZDIFTQT (KIDIA:KFDIA,KLEV) = 0.0_JPRB
ZDIFTSLG(KIDIA:KFDIA,KLEV) = 0.0_JPRB
ZDIFTQT (KIDIA:KFDIA,KLEV) = PDIFTQ(KIDIA:KFDIA,KLEV)
ZDIFTSLG(KIDIA:KFDIA,KLEV) = PDIFTS(KIDIA:KFDIA,KLEV)



!     ------------------------------------------------------------------

!*         8.     SLG, QT, U, V FLUX COMPUTATIONS AND T,SKIN TENDENCY
!                 ---------------------------------------------------

DO JL=KIDIA,KFDIA
  ZDIFTQT (JL,0) = 0.0_JPRB
  PDIFTQ  (JL,0) = 0.0_JPRB
  PDIFTL  (JL,0) = 0.0_JPRB
  PDIFTI  (JL,0) = 0.0_JPRB
  PDIFTS  (JL,0) = 0.0_JPRB
  ZDIFTSLG(JL,0) = 0.0_JPRB
  PSTRTU  (JL,0) = 0.0_JPRB
  PSTRTV  (JL,0) = 0.0_JPRB
ENDDO

DO JK=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
    ZGDPH = - (PAPHM1(JL,JK)-PAPHM1(JL,JK+1)) * ZRG
!...change in slg,qt,u,v tendencies are converted to fluxes
    ZDIFTSLG(JL,JK) = ( ZSLGEWODIS(JL,JK+1) - ZSLGEA(JL,JK+1) ) * ZGDPH &
                    & + ZDIFTSLG(JL,JK+1)  
    ZDIFTQT(JL,JK)  = (ZQTE (JL,JK+1)-ZQTEA(JL,JK+1))*ZGDPH + ZDIFTQT(JL,JK+1)
  ENDDO
ENDDO



!     ------------------------------------------------------------------
! 
!*         9.     OLD SATURATION SPECIFIC HUMIDITY
!                 --------------------------------
!
!          This part will be redundant soon.
!          Only ZALFAW is still used to separate ice and liquid water.
!


  !-- state after dynamics and radiation --
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
      
      ZAUPD(JL,JK)  = PAM1(JL,JK) + PAE(JL,JK) * ZTMST
      ZLUPD(JL,JK)  = PLM1(JL,JK) + PLE(JL,JK) * ZTMST
      ZIUPD(JL,JK)  = PIM1(JL,JK) + PIE(JL,JK) * ZTMST
      ZQUPD(JL,JK)  = PQM1(JL,JK) + PQE(JL,JK) * ZTMST
      ZTUPD(JL,JK)  = PTM1(JL,JK) + PTE(JL,JK) * ZTMST
!          total condensate (water + ice)
      ZLIUPD(JL,JK) = ZLUPD(JL,JK) + ZIUPD(JL,JK)
      
  ENDDO
ENDDO

  !-- calculate alfa --
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
      ZALFAW(JL,JK)=FOEALFA(ZTUPD(JL,JK))
  ENDDO
ENDDO


  
!     ------------------------------------------------------------------

!*         10.    QT VARIANCE BUDGET
!                 ------------------
!
!          Flux-gradient production, transport and dissipation
!

!-- reset variance ----
!xmk: prognostic total water variance - turn off reset
!DO JK=1,KLEV
! DO JL=KIDIA,KFDIA
!    PVAR(JL,JK) = 0._JPRB
!  ENDDO
!ENDDO
!xxx

DO JL=KIDIA,KFDIA
  PVAR(JL,1)    = 0._JPRB
  PVAR(JL,KLEV) = 0._JPRB
  ZTAU(JL,1)    = 0._JPRB
  ZTAU(JL,KLEV) = 0._JPRB
ENDDO


ZDQTDZ(:,:)  = 0._JPRB

DO JK=KLEV-1,2,-1
  DO JL=KIDIA,KFDIA
    
      
      !-- variance dissipation timescale ZTAU --
!xmk  ZTAU(JL,JK) =  PZINV(JL) / MAX(PWUAVG(JL),0.01_JPRB)

!amk ... optional 2x time-scale of sigma(qt) decay: more cloud
!     ZTAU(JL,JK) = ZTAU(JL,JK) * 2.0_JPRB
!xxx
!amk ... take local cloud height instead of zi:
!     ZTAU(JL,JK) = PGEOM1(JL,JK)*ZRG / MAX(PWUAVG(JL),0.01_JPRB)
!xxx
!amk ... take cloud top height instead of zi (if no cloud use zi):
      IF (KPBLTYPE(JL)==2 .OR. KPBLTYPE(JL)==3) THEN
        ZTAU(JL,JK) = ZCLDTOP(JL) / MAX(PWUAVG(JL),0.01_JPRB)
      ELSE
        ZTAU(JL,JK) = PZINV(JL)   / MAX(PWUAVG(JL),0.01_JPRB)
      ENDIF
!xxx
      
      !-- do the individual variance budget terms --
      ZVARDISS   = 0._JPRB
      ZVARGEN    = 0._JPRB
      ZVARTRANS  = 0._JPRB
      ZWQTF   = 0._JPRB
      
      IF (LLPBL(JL,JK)) THEN
      
!dmk    !--- I  flux-gradient variance production at full level ---
!       IF (KPBLTYPE(JL)==2 .AND. JK<=KVARTOP(JL)+1 ) THEN
!         !for stratocumulus, protect variance production in top PBL layer against strong capping gradient
!         ZDQTDZ = (ZQTUPD(JL,JK)-ZQTUPD(JL,JK+1)) * RG / (PGEOM1(JL,JK)-PGEOM1(JL,JK+1))
!       ELSE
!         !otherwise, do it truely centered
!         ZDQTDZ = (ZQTUPD(JL,JK-1)-ZQTUPD(JL,JK+1)) * RG / (PGEOM1(JL,JK-1)-PGEOM1(JL,JK+1))
!    ENDIF
!
!       ZWQTF = -(RD * PTM1(JL,JK) / PAPM1(JL,JK)) * &       ! flux/rho = w'qt'  
!              & (ZDIFTQT(JL,JK-1) + ZDIFTQT(JL,JK))/2._JPRB
!
!xxx    ZVARGEN = -2._JPRB * ZWQTF * ZDQTDZ

        ZDQTDZ(JL,JK-1) = (ZQTM1(JL,JK-1)-ZQTM1(JL,JK  )) / (PGEOM1(JL,JK-1)-PGEOM1(JL,JK  )) * RG 
        ZDQTDZ(JL,JK)   = (ZQTM1(JL,JK  )-ZQTM1(JL,JK+1)) / (PGEOM1(JL,JK  )-PGEOM1(JL,JK+1)) * RG 

!       - 2 * ( w'qt'|k-1 * dqt/dz|k-1 + w'qt'|k * dqt/dz|k ) / 2
!       note that: w'qt' = flux/rho and flux is downward
        ZVARGEN = (RD * PTM1(JL,JK) / PAPM1(JL,JK))    &
              & * ( ZDIFTQT(JL,JK-1) * ZDQTDZ(JL,JK-1) &
              &   + ZDIFTQT(JL,JK  ) * ZDQTDZ(JL,JK  ) ) 
         
        ZVARGEN = MAX(ZVARGEN,0._JPRB)   ! exclude countergradient flow
          
          
        !--- II   variance transport at full level ---
        ZVARTRANS = - ( ZWQT2(JL,JK-1) - ZWQT2(JL,JK) ) / (PGEOH(JL,JK-1)-PGEOH(JL,JK)) * RG
             
             
        !--- III  variance dissipation at full level (for output only) ---
        ZVARDISS = MIN(0._JPRB, -ZVARGEN -ZVARTRANS)  


!amk: alternative qt variance equation (as suggested by Bechtold and Chaboureaux):
!!      d(var)/dt = 2 K (dqt/dz)^2 + (qt,up-qt,mean)^2/tau - var^2/tau
!       ZVARGEN = PKH(JL,JK-1) * ZDQTDZ(JL,JK-1)**2 &
!             & + PKH(JL,JK  ) * ZDQTDZ(JL,JK  )**2
!
!       IF ( KPBLTYPE(JL)==3 .AND. JK>=IPTOP(JL,3) ) THEN
!!        ZQTENH(JL,JK) = ( ZQTM1(JL,JK+1) *(PGEOH(JL,JK-1)-PGEOH(JL,JK  )) &
!!                    & +   ZQTM1(JL,JK)   *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
!!                    &   )                /(PGEOH(JL,JK-1)-PGEOH(JL,JK+1))
!!        ZVARGEN = ZVARGEN + (ZQTUH(JL,JK,3) - ZQTENH(JL,JK))**2 / ZTAU(JL,JK)
!         ! upwind as in solver
!         ZVARGEN = ZVARGEN + (ZQTUH(JL,JK,3) - ZQTM1(JL,JK))**2 / ZTAU(JL,JK)
!       ENDIF
!
!       ZVARTRANS = 0.0_JPRB
!xxx    
        
      ENDIF

  
      !--- update the variance ---

!xmk: prognostic total water variance 
!     PVAR(JL,JK) = (ZVARGEN+ZVARTRANS) * ZTAU(JL,JK) 

!     ...only prognostic when within PBL now and last time-step
!     IF ( LLPBL(JL,JK) .AND. PVAR(JL,JK) .GT. 0.0_JPRB ) THEN
!       PVAR(JL,JK) = ( PVAR(JL,JK) + (ZVARGEN + ZVARTRANS) * ZTMST ) &
!                 & * ZTAU(JL,JK) / (ZTMST + ZTAU(JL,JK))
!     ELSE
!     ...diagnostic if variance was 0 last step (outside PBL last step)
!       PVAR(JL,JK) = (ZVARGEN+ZVARTRANS) * ZTAU(JL,JK)
!     ENDIF

!     ...prognostic within PBL
      IF ( LLPBL(JL,JK) ) THEN
       !exact solution:
        PVAR(JL,JK) = ( PVAR(JL,JK) - ZTAU(JL,JK) * (ZVARGEN + ZVARTRANS) ) *  &
          &               EXP(-ZTMST/ZTAU(JL,JK))                              &
          &             + ZTAU(JL,JK) * (ZVARGEN + ZVARTRANS)
       !implicit solution:
       !PVAR(JL,JK) = ( PVAR(JL,JK) + (ZVARGEN + ZVARTRANS) * ZTMST ) &
       !          & * ZTAU(JL,JK) / (ZTMST + ZTAU(JL,JK))
      ELSE
!     ...prognostic outside PBL (generation and transport are zero)
        ZTAU(JL,JK) = 3*3600.0_JPRB  ! decay time-scale = 3 hours
        PVAR(JL,JK) = PVAR(JL,JK) * exp(-ZTMST/ZTAU(JL,JK))              !exact
!       PVAR(JL,JK) = PVAR(JL,JK) * ZTAU(JL,JK) / (ZTMST + ZTAU(JL,JK))  !implicit
      ENDIF
!xxx

      PVAR(JL,JK) = MAX(PVAR(JL,JK),0._JPRB)

  ENDDO
ENDDO


      
!     ------------------------------------------------------------------
!
!*         11.    VECTORIZED BIMODAL CLOUD SCHEME
!                 -------------------------------
!          
!          The EDMF decomposition is extended into the cloud scheme, by doing a bimodal
!          PDF: one diffusive, one updraft. Each PDF is Gaussian, their 1st and 2nd moments are
!          parameterized. The moments of all PDFs are related, see Lewellen and Yoh (JAS, 1993).
!          The scheme is formulated in {thl,qt} space, using vector calculus.
! 
!          Variance closures:
!            * The overall variance is done through the full budget (see section 10).
!            * The updraft PDF variance is done using properties of multiple updrafts.
!
!          See Neggers et al. (In preparation for JAS, 2006)
!

  !  use new {QT,SLG} state
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      
      ZQTTEST(JL,JK) = 0._JPRB

      SELECT CASE (KPBLTYPE(JL))
      
        CASE (0,1,4)
          ZSLGUP(JL,JK)  = ZSLGUPD(JL,JK)
          ZQTUP(JL,JK)   = ZQTUPD(JL,JK)
          ZQTTEST(JL,JK) = ZQTUPD(JL,JK)
          
        CASE (2,3)
          IF (LLPBL(JL,JK) .AND. ZSLGUH(JL,JK,3)>0._JPRB .AND. ZSLGUH(JL,JK-1,3)>0._JPRB ) THEN
            !-- interpolate the updraft fields to full levels --
!xmk: fixed below: *ZTMST!!!
!           ZSLGUP(JL,JK)  = ( ZSLGUH(JL,JK,3)+ ZSLGUH(JL,JK-1,3) ) / 2._JPRB + ZSLGE(JL,JK)
!           ZQTUP(JL,JK)   = ( ZQTUH(JL,JK,3) + ZQTUH(JL,JK-1,3)  ) / 2._JPRB + ZQTE(JL,JK)
!           ZQTTEST(JL,JK) = ( ZQTUH(JL,JK,1) + ZQTUH(JL,JK-1,1)  ) / 2._JPRB + ZQTE(JL,JK)
            ZSLGUP(JL,JK)  = ( ZSLGUH(JL,JK,3)+ ZSLGUH(JL,JK-1,3) ) / 2._JPRB + ZSLGE(JL,JK) *ZTMST
            ZQTUP(JL,JK)   = ( ZQTUH(JL,JK,3) + ZQTUH(JL,JK-1,3)  ) / 2._JPRB + ZQTE(JL,JK)  *ZTMST
            ZQTTEST(JL,JK) = ( ZQTUH(JL,JK,1) + ZQTUH(JL,JK-1,1)  ) / 2._JPRB + ZQTE(JL,JK)  *ZTMST
!xxx
          ELSE
            ZSLGUP(JL,JK)  = ZSLGUPD(JL,JK)
            ZQTUP(JL,JK)   = ZQTUPD(JL,JK)
            ZQTTEST(JL,JK) = ZQTUPD(JL,JK)
          ENDIF
          
      END SELECT
      
      ZRHO = PAPHM1(JL,JK)/( RD*PTM1(JL,JK)*(1.0_JPRB+RETV*PQM1(JL,JK)) )
      ZAWAKE(JL,JK) = MAX(0._JPRB, ZDETR(JL,JK) * ZTAU(JL,JK) / ZRHO )
      !ZAWAKE(JL,JK) = MAX(0._JPRB, ZDETR(JL,JK) * ZTAU(JL,JK)*2._JPRB / ZRHO )
        
    ENDDO
  ENDDO


  CALL VDFCLOUD ( KIDIA     , KFDIA   , KLON    , KLEV   , IDRAFT , &
                & PAPM1     , PGEOM1  , PGEOH   , &
                & ZQTUPD    , ZSLGUPD , &
                & ZFRACB    , KVARTOP , &
                & ZQTUP     , ZSLGUP  , &
                & ZQTTEST   , &
                & PVAR      , ZAWAKE  , &
! DIAGNOSTIC OUTPUT
                & PEXTR2    , KFLDX2  , PEXTRA  , KLEVX  , KFLDX  , &
!              
                & ZCLDFRAC  , ZQLAV   , PLDIFF)



!     ------------------------------------------------------------------
!
!*         12.    NET TENDENCIES
!                 --------------
!
!          Calculate net tendencies of the prognostic model variables QL, QI, A, QV and T


  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA


      !RN --- temporarily switch off cloudiness, for testing ---
      !ZCLDFRAC(JL,JK) = 0._JPRB
      !ZQLAV(JL,JK)    = 0._JPRB
      

      !--- Convert back from conserved variables QT and SLG to non-conserved QV and T ---
      !
      !  New cloud variables (liquid, ice, fraction) for tendency calculation:
      !
      !    Within PBL: use new values from VDFCLOUD
      !    Above PBL:  use t-1 input values + tendencies
      !
      !    Result: * above PBL, VDFMAIN does not change cloud variables:
      !                tendencies only contain contributions by rad+dyn.
      !
      IF (LLPBL(JL,JK)) THEN    !reset cloudiness within PBL
      !IF (.TRUE.) THEN          !reset cloudiness everywhere! (also above PBL)

        !-- safety: total condensate can not be larger than total specific humidity ---
        ZQLAV(JL,JK) = MIN( ZQTUPD(JL,JK) ,ZQLAV(JL,JK) )
        ZAUPD(JL,JK) = ZCLDFRAC(JL,JK)

        !-- decomposition of total condensate into ice and liquid ---
        !ZALFAW(JL,JK) = FOEALFA(PTM1(JL,JK))   !new alpha?
        ZLUPD(JL,JK) = ZQLAV(JL,JK) * ZALFAW(JL,JK) 
        ZIUPD(JL,JK) = ZQLAV(JL,JK) * ( 1.0_JPRB - ZALFAW(JL,JK))

      ELSE
      
        !-- outside PBL, maintain tendencies from rad+dyn ---
        ZAUPD(JL,JK) = PAM1(JL,JK) + PAE(JL,JK)*ZTMST
        ZLUPD(JL,JK) = PLM1(JL,JK) + PLE(JL,JK)*ZTMST
        ZIUPD(JL,JK) = PIM1(JL,JK) + PIE(JL,JK)*ZTMST
        
      ENDIF  
      
      
      !--- Derive non-conserved properties QV and T ---
      ZQUPD(JL,JK)  = ZQTUPD(JL,JK) - ZLUPD(JL,JK) - ZIUPD(JL,JK)
      ZTUPD(JL,JK)  = ( ZSLGUPD(JL,JK) - PGEOM1(JL,JK) &
        &     + RLVTT * ZLUPD(JL,JK) + RLSTT * ZIUPD(JL,JK) &
        &   ) / ( RCPD * ( 1.0_JPRB + RVTMP2 * ZQUPD(JL,JK) ) )   !compare to T->SLG conversion in section 2.2


      !--- Calculate the final tendencies between state at t-1 and state after rad + dyn + pbl ---
      PQE(JL,JK) = ( ZQUPD(JL,JK) - PQM1(JL,JK) ) * ZRTMST
      PTE(JL,JK) = ( ZTUPD(JL,JK) - PTM1(JL,JK) ) * ZRTMST
      PLE(JL,JK) = ( ZLUPD(JL,JK) - PLM1(JL,JK) ) * ZRTMST
      PIE(JL,JK) = ( ZIUPD(JL,JK) - PIM1(JL,JK) ) * ZRTMST
      PAE(JL,JK) = ( ZAUPD(JL,JK) - PAM1(JL,JK) ) * ZRTMST
      
!amk: debug
!IF ( PTE(JL,JK)  > 50.0/3600 ) THEN
!  WRITE(*,*) 'PTE>50K/h PTE(JL,JK), JK',     PTE(JL,JK),  JK
!ENDIF
!IF ( PVOM(JL,JK) > 50.0/3600 ) THEN
!  WRITE(*,*) 'PVOM>50m/s/h PVOM(JL,JK), JK', PVOM(JL,JK), JK
!ENDIF
!IF ( PVOL(JL,JK) > 50.0/3600 ) THEN
!  WRITE(*,*) 'PVOL>50m/s/h PVOL(JL,JK), JK', PVOL(JL,JK), JK
!ENDIF
!xxx

    ENDDO
  ENDDO

!amk: debug
! DO JK=1,KLEV
!   DO JL=KIDIA,KFDIA
!      IF ( ZTUPD(JL,JK) < 100.0_JPRB .OR. ZTUPD(JL,JK) > 400.0_JPRB ) THEN
!        WRITE(*,*) 'vdfmain T<100 or T>400, kstep, JL,JK,T:', KSTEP, JL, JK, ZTUPD(JL,JK)
!      ENDIF
!      IF ( ZQUPD(JL,JK) < -0.01_JPRB .OR. ZQUPD(JL,JK) > 0.1_JPRB ) THEN
!        WRITE(*,*) 'vdfmain q<-10g/kg or q>100g/kg, kstep, JL,JK,Q:', KSTEP, JL, JK, ZQUPD(JL,JK)
!      ENDIF
!   ENDDO
! ENDDO
!xxx

! DO JL=KIDIA,KFDIA
!   IF ( (SUM(PAHFSTI(JL,:)) == 0.0) .or. (SUM(PEVAPTI(JL,:)) == 0.0) .or. &
!        (PDIFTS(JL,KLEV)    == 0.0) .or. (PDIFTQ(JL,KLEV)    == 0.0) ) THEN
!     write(*,*) 'vdfmain5: ', PDIFTS(JL,KLEV), PDIFTQ(JL,KLEV), PAHFSTI(JL,:), PEVAPTI(JL,:)*RLVTT, &
!       & PFRTI(JL,:)
!   ENDIF
! ENDDO

!     ------------------------------------------------------------------

!*         13.    FLUX COMPUTATIONS
!                 -----------------
!


!*         13.1    Q, QL, QI AND S FLUX

DO JK=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
    ZGDPH = - (PAPHM1(JL,JK)-PAPHM1(JL,JK+1)) * ZRG
!...changes in q,l,i tendencies are converted to fluxes
    PDIFTQ(JL,JK) = (PQE (JL,JK+1) - ZQEA(JL,JK+1)) * ZGDPH + PDIFTQ(JL,JK+1)
    PDIFTL(JL,JK) = (PLE (JL,JK+1) - ZLEA(JL,JK+1)) * ZGDPH + PDIFTL(JL,JK+1)
    PDIFTI(JL,JK) = (PIE (JL,JK+1) - ZIEA(JL,JK+1)) * ZGDPH + PDIFTI(JL,JK+1)
!...slg=s-Lc*ql-Ld*qi (same for fluxes)
    PDIFTS(JL,JK) = ZDIFTSLG(JL,JK) &
                & + RLVTT * PDIFTL(JL,JK) + RLSTT * PDIFTI(JL,JK)  
  ENDDO
ENDDO

!*         13.2  PBL PRECIPITATION ENTHALPY FLUXES

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PFHPVL(JL,JK) = -RLVTT*PFPLVL(JL,JK)
    PFHPVN(JL,JK) = -RLSTT*PFPLVN(JL,JK)
  ENDDO    
ENDDO


IF ( LLDIAG ) THEN
  DO JL=KIDIA,KFDIA
    PEXTRA(JL,23,41) = PDIFTQ(JL,IPLCL(JL,3)) * RLVTT
    PEXTRA(JL,24,41) = PDIFTS(JL,IPLCL(JL,3))
    PEXTRA(JL,25,41) = PSTRTU(JL,IPLCL(JL,3))
    PEXTRA(JL,26,41) = PSTRTV(JL,IPLCL(JL,3))
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('VDFMAIN',1,ZHOOK_HANDLE)
END SUBROUTINE VDFMAIN


END MODULE mo_vdfmain
