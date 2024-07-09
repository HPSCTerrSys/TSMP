SUBROUTINE receive_fld_2cos(kt, lcoupled)

!---------------------------------------------------------------------
! Description:
!  This routine receives the atmospheric fields from COSMO and updates
!  the CLM3.53.53.5 variable x(:,:,:)
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
! 2.1        2012/09/26 Markus Uebel, Prabhakar Shrestha 
!   CO2 coupling included
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

USE oas_clm_vardef
USE shr_kind_mod ,            ONLY : r8 => shr_kind_r8
USE clm_time_manager ,        ONLY : dtime,                            &! timestep in second
                                     nelapse,                          &!
                                     get_nstep,                        &! return timestep number
                                     get_curr_date                      ! get current date

USE spmdMod ,                 ONLY : masterproc, dummy => MPI_LOGICAL   ! CPS
USE netcdf
use clm_atmlnd  ,             only :  atm_a2l
USE decompMod   ,          ONLY :  get_proc_bounds, adecomp
use clm_varcon, only : sb,rair, cpair, co2_ppmv_const, o2_molar_const, tcrit, c13ratio
use shr_const_mod, only : SHR_CONST_TKFRZ, SHR_CONST_PSTD
USE domainMod   , ONLY :  amask
!==============================================================================

IMPLICIT NONE

!==============================================================================

! Parameter List
INTEGER, INTENT(IN)                    ::   kt                           ! time step

LOGICAL, INTENT(OUT)                   ::  lcoupled                      ! CPS flag for coupling 
!
INTEGER, PARAMETER ::   jps_t   =  1            ! temperature
INTEGER, PARAMETER ::   jps_u   =  2            ! u wind
INTEGER, PARAMETER ::   jps_v   =  3            ! v wind
INTEGER, PARAMETER ::   jps_q   =  4            ! specific water vapor content
INTEGER, PARAMETER ::   jps_th  =  5            ! thickness of lowest level (m)
INTEGER, PARAMETER ::   jps_pr  =  6            ! surface pressure (Pa)
INTEGER, PARAMETER ::   jps_rs  =  7            ! direct shortwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_fs  =  8            ! diffuse shortwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_lw  =  9            ! longwave downward radiation (W/m2) 
INTEGER, PARAMETER ::   jps_cr  = 10            ! convective rain precipitation      (kg/m2*s)
INTEGER, PARAMETER ::   jps_cs  = 11            ! convective snow precipitation      (kg/m2*s)
INTEGER, PARAMETER ::   jps_gr  = 12            ! gridscale rain precipitation
INTEGER, PARAMETER ::   jps_gs  = 13            ! gridscale snow precipitation
INTEGER, PARAMETER ::   jps_gg  = 14            ! gridscale graupel precipitation
INTEGER, PARAMETER ::   jps_cp  = 15            ! total convective precipitation
INTEGER, PARAMETER ::   jps_gp  = 16            ! total gridscale precipitation
!MU (13.09.2012)
INTEGER, PARAMETER ::   jps_co2 = 17            ! CO2 partial pressure (Pa)
!MU (13.09.2012)

! Local Variables

INTEGER                                :: jn, isec, ier , begg,endg,begl,endl,begc,endc,begp,endp
INTEGER , DIMENSION(3)              ::   nrcvinfo           ! OASIS info argument
INTEGER                                :: ji, jj, g


real(r8) ea                    !atmospheric emissivity
real(r8) esatw                 !saturation vapor pressure over water (Pa)
real(r8) esati                 !saturation vapor pressure over ice (Pa)
real(r8) e                     !vapor pressure (Pa)
real(r8) qsat                  !saturation specific humidity (kg/kg)
real(r8) a0,a1,a2,a3,a4,a5,a6  !coefficients for esat over water
real(r8) b0,b1,b2,b3,b4,b5,b6  !coefficients for esat over ice
real(r8) tdc, t                !Kelvins to Celcius function and its input

parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
                a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
                a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
                a6=6.136820929e-11_r8)
 
parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
                b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
                b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
                b6=1.838826904e-10_r8)


!
! function declarations
!
tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))



REAL(KIND=r8)                          :: tmp1,tmp2

REAL(KIND=r8) ,allocatable ,  DIMENSION(:,:) ::   ztmp1   

CHARACTER(8)                           :: dateoas
CHARACTER(10)                          :: timeoas
CHARACTER(5)                           :: zoneoas
INTEGER ,DIMENSION(8)                  :: valuesoas
REAL(KIND=r8)                          :: millisec
REAL(KIND=r8)                          :: sec_date
INTEGER ,DIMENSION(3)                  :: dimids                ! CPS 2 to 3
INTEGER ,DIMENSION(krcv)               :: il_var_id
INTEGER                                :: il_file_id, status
REAL(KIND=r8), SAVE                    :: tic, toc
INTEGER                                :: year, mon, day, sec   ! curent date
REAL(KIND=r8)                          :: co2_ppmv_ssp(151), co2_ppmv_ssp126(151), co2_ppmv_ssp370(151),co2_ppmv_ssp585(151)

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine receive_fld2cos 
!------------------------------------------------------------------------------

 isec = dtime * kt
 call get_curr_date(year, mon, day, sec) ! for ssp scenarios

 ! Default value (changed when coupling)
 lcoupled = .FALSE.
 nrcvinfo = OASIS_idle

 CALL get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

 allocate(ztmp1(begg:endg,3))


! Update CLM variable with the received input from COSMO
! FORC_TXY
ztmp1(:,:) = -1._r8
nrcvinfo(:) = OASIS_idle
if(srcv(jps_t)%laction) CALL oas_clm_rcv( jps_t, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) )
IF( nrcvinfo(1)==OASIS_Rcv) then     ! air temperature (K)
        do g=begg,endg
  !      if( amask(g).ne.1) ztmp1(g,1) = 60.
        if (nint(ztmp1(g,1)) == -1) then
              if(masterproc) write(6,*)'ATM error: TBOT has not been read by atmrd'
        else if (ztmp1(g,1) < 50._r8) then
              if(masterproc) write(6,*)'ATM error: TBOT appears to be in deg C'
              if(masterproc) write(6,*)'Converting to Kelvins now'
              atm_a2l%forc_t(g) = ztmp1(g,1) + SHR_CONST_TKFRZ
        else
              atm_a2l%forc_t(g) = ztmp1(g,1)
        end if
        enddo
        lcoupled = .TRUE.
endif




!########################################

 ! FORC_UXY, FORC_VXY
ztmp1(:,:) = -999._r8
nrcvinfo(:) = OASIS_idle
if(srcv(jps_u)%laction) CALL oas_clm_rcv( jps_u, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) )
IF( nrcvinfo(1)==OASIS_Rcv )  then 
        do g=begg,endg
           atm_a2l%forc_u(g) = abs(ztmp1(g,1))
        enddo
        lcoupled = .TRUE.
endif    

if(srcv(jps_v)%laction) CALL oas_clm_rcv( jps_v, isec, ztmp1(:,2),begg,endg, nrcvinfo(2) )
IF( nrcvinfo(2)==OASIS_Rcv )  then                   
        do g=begg,endg
           atm_a2l%forc_v(g) = abs(ztmp1(g,2))
        enddo
        lcoupled = .TRUE.
endif


IF( srcv(jps_u)%laction .AND. srcv(jps_v)%laction .and. nrcvinfo(1)==OASIS_Rcv .and. nrcvinfo(2)==OASIS_Rcv) then         
        do g=begg,endg                                          ! wind (m/s)
          if (nint(ztmp1(g,1)).eq.-999 .or.  nint( ztmp1(g,2)).eq. -999) then
              if(masterproc) write(6,*)'ATM error: WIND has not been read by atmrd'
          endif
          atm_a2l%forc_wind(g) = sqrt(ztmp1(g,1)**2 + ztmp1(g,2)**2)
        enddo
endif


!########################################

! FORC_SOLSXY, FORC_SOLLXY, FORC_SOLSDXY, FORC_SOLLDXY
ztmp1(:,:) = -1._r8
nrcvinfo(:) = OASIS_idle
if (srcv(jps_rs)%laction) CALL oas_clm_rcv( jps_rs, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) )
if (srcv(jps_fs)%laction) CALL oas_clm_rcv( jps_fs, isec, ztmp1(:,2),begg,endg, nrcvinfo(2) )
!ztmp(:,3) (former x(:,:,8) is never coupled)
IF( (nrcvinfo(1)==OASIS_Rcv .or. nrcvinfo(2)==OASIS_Rcv )   ) then
lcoupled = .TRUE.
do g=begg,endg
if (nint(ztmp1(g,1))==-1.or.nint(ztmp1(g,2))==-1) then
              if (nint(ztmp1(g,3)) /= -1) then
                 atm_a2l%forc_solad(g,1)  = 0.7_r8 * (0.5_r8 * ztmp1(g,3))
                 atm_a2l%forc_solad(g,2)  = atm_a2l%forc_solad(g,1)
                 atm_a2l%forc_solai(g,1) = 0.3_r8 * (0.5_r8 * ztmp1(g,3))
                 atm_a2l%forc_solai(g,2) = atm_a2l%forc_solai(g,1)
              else
                 if(masterproc) write(6,*)'ATM error: neither FSDSdir/dif nor'
                 if(masterproc) write(6,*)'       FSDS have been read in by atmrd'
              end if
else
              atm_a2l%forc_solad(g,1)  = 0.5_r8 * ztmp1(g,1)
              atm_a2l%forc_solad(g,2)  = atm_a2l%forc_solad(g,1)
              atm_a2l%forc_solai(g,1) = 0.5_r8 * ztmp1(g,2)
              atm_a2l%forc_solai(g,2) = atm_a2l%forc_solai(g,1)
end if
atm_a2l%forc_solar(g) = atm_a2l%forc_solad(g,1) + atm_a2l%forc_solad(g,2)  + atm_a2l%forc_solai(g,1) + atm_a2l%forc_solai(g,2)

enddo 
endif

!########################################

! FORC_PSRFXY, FORC_PBOTXY
ztmp1(:,:) = -1._r8
nrcvinfo(:) = OASIS_idle
if(srcv(jps_pr)%laction) CALL oas_clm_rcv( jps_pr, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) )
IF( nrcvinfo(1)==OASIS_Rcv) then
lcoupled = .TRUE.
do g=begg,endg
      if (nint(ztmp1(g,1)) == -1) then
              atm_a2l%forc_psrf(g) = SHR_CONST_PSTD
      else
              atm_a2l%forc_psrf(g) = ztmp1(g,1)
      end if
      atm_a2l%forc_pbot(g)  = atm_a2l%forc_psrf(g)
enddo
endif

!########################################

!FORC_QXY
ztmp1(:,:) = -1._r8
nrcvinfo(:) = OASIS_idle
if(srcv(jps_q)%laction) CALL oas_clm_rcv( jps_q, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) ) ! former x(:,:,3)
!former x(:,:,4) and x(:,:,5) are never coupled

if( nrcvinfo(1)==OASIS_Rcv) then
lcoupled = .TRUE.
do g=begg,endg
if (nint(ztmp1(g,1)) == -1) then
    if (nint(ztmp1(g,2)) == -1) then
        if (nint(ztmp1(g,3)) == -1) then
            if(masterproc) write(6,*)'ATM error: Humidity has not been'
            if(masterproc) write(6,*)'read by atmrd'
        else          !using RH as %
            if (atm_a2l%forc_t(g) > SHR_CONST_TKFRZ) then
               e = ztmp1(g,3)/100._r8 * esatw(tdc(atm_a2l%forc_t(g)))
            else
               e = ztmp1(g,3)/100._r8 *esati(tdc(atm_a2l%forc_t(g)))
            end if
        end if
        atm_a2l%forc_q(g) = 0.622_r8*e /(atm_a2l%forc_pbot(g) - 0.378_r8*e)
    else             !using Tdew
        if (ztmp1(g,2) < 50._r8) then
             if(masterproc) write(6,*)'ATM warning: Tdew appears to be in'
             if(masterproc) write(6,*)'deg C, so converting to Kelvin'
             ztmp1(g,2) = ztmp1(g,2) + SHR_CONST_TKFRZ
        end if
        if (ztmp1(g,2) > atm_a2l%forc_t(g)) then
             if(masterproc) write(6,*)'ATM warning: Dewpt temp > temp!'
        end if
        if (ztmp1(g,2) > SHR_CONST_TKFRZ) then
             e = esatw(tdc(ztmp1(g,2)))
        else
             e = esati(tdc(ztmp1(g,2)))
        end if
        atm_a2l%forc_q(g) = 0.622_r8*e /(atm_a2l%forc_pbot(g) - 0.378_r8*e)
    end if
else                !using QBOT in kg/kg
    if (atm_a2l%forc_t(g) > SHR_CONST_TKFRZ) then
       e = esatw(tdc(atm_a2l%forc_t(g)))
    else
       e = esati(tdc(atm_a2l%forc_t(g)))
    end if
    qsat = 0.622_r8*e / (atm_a2l%forc_pbot(g) -0.378_r8*e)
    if (qsat < ztmp1(g,1)) then
      atm_a2l%forc_q(g) = qsat
!                 write(6,*)'ATM warning: qsat < q!'
    else
      atm_a2l%forc_q(g) = ztmp1(g,1)
    end if
end if

        
  atm_a2l%forc_th(g)  = atm_a2l%forc_t(g) * (atm_a2l%forc_psrf(g) / atm_a2l%forc_pbot(g))**(rair/cpair)
  atm_a2l%forc_vp(g)  = atm_a2l%forc_q(g) * atm_a2l%forc_pbot(g) / (0.622_r8 + 0.378_r8 * atm_a2l%forc_q(g))
  atm_a2l%forc_rho(g) = (atm_a2l%forc_pbot(g) - 0.378_r8 * atm_a2l%forc_vp(g)) / (rair * atm_a2l%forc_t(g))

enddo
endif

!########################################

! ZGCMXY

ztmp1(:,:) = -1._r8
nrcvinfo(:) = OASIS_idle
IF( srcv(jps_th)%laction) CALL oas_clm_rcv( jps_th, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) )
IF( nrcvinfo(1)==OASIS_Rcv) then
lcoupled = .TRUE.
do g=begg,endg
          if (nint(ztmp1(g,1)) == -1) then
             atm_a2l%forc_hgt(g) = 30._r8
             atm_a2l%forc_hgt_u(g) = 30._r8
             atm_a2l%forc_hgt_t(g) = 30._r8
             atm_a2l%forc_hgt_q(g) = 30._r8   
          else
              atm_a2l%forc_hgt(g) = ztmp1(g,1)
              atm_a2l%forc_hgt_u(g) = ztmp1(g,1)
              atm_a2l%forc_hgt_t(g) = ztmp1(g,1)
              atm_a2l%forc_hgt_q(g) = ztmp1(g,1)
          end if
enddo
endif

!########################################

! PRCXY, PRLXY
ztmp1(:,:) = -1._r8
nrcvinfo(:) = OASIS_idle
IF( srcv(jps_gg)%laction )CALL oas_clm_rcv( jps_gg, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) ) 
IF( srcv(jps_gs)%laction) CALL oas_clm_rcv( jps_gs, isec, ztmp1(:,2),begg,endg, nrcvinfo(2) )
IF( srcv(jps_gr)%laction) CALL oas_clm_rcv( jps_gr, isec, ztmp1(:,3),begg,endg, nrcvinfo(3) ) !former x(:,:,14) 

IF( srcv(jps_gr)%laction .AND. srcv(jps_gs)%laction ) THEN
          ztmp1(:,3) = ztmp1(:,3) + ztmp1(:,2)    ! gridscale precip (mm    /s) 
  IF( srcv(jps_gg)%laction )                                 &
          ztmp1(:,3) = ztmp1(:,3) +  ztmp1(:,1)
 
ENDIF     ! if coupled graupel


if(nrcvinfo(1)==OASIS_Rcv .or. nrcvinfo(2)==OASIS_Rcv .or. nrcvinfo(3)==OASIS_Rcv) lcoupled = .TRUE.

ztmp1(:,1) = -1._r8
nrcvinfo(1) = OASIS_idle
IF( srcv(jps_gp)%laction) CALL oas_clm_rcv( jps_gp, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) )
IF( srcv(jps_gp)%laction )  ztmp1(:,3) = ztmp1(:,1)

if(nrcvinfo(1)==OASIS_Rcv) lcoupled = .TRUE.


ztmp1(:,1:2) = -1._r8
nrcvinfo(1:2) = OASIS_idle
IF( srcv(jps_cr)%laction) CALL oas_clm_rcv( jps_cr, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) )
IF( srcv(jps_cs)%laction) CALL oas_clm_rcv( jps_cs, isec, ztmp1(:,2),begg,endg, nrcvinfo(2) ) !former x(:,:,13    )

IF( srcv(jps_cr)%laction .AND. srcv(jps_cs)%laction ) ztmp1(:,2) = ztmp1(:,1) + ztmp1(:,2)

if(nrcvinfo(1)==OASIS_Rcv .or. nrcvinfo(2)==OASIS_Rcv) lcoupled = .TRUE.


ztmp1(:,1) = -1._r8
nrcvinfo(1) = OASIS_idle
IF( srcv(jps_cp)%laction) CALL oas_clm_rcv( jps_cp, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) )
IF( srcv(jps_cp)%laction )  ztmp1(:,2) = ztmp1(:,1)

if(nrcvinfo(1)==OASIS_Rcv) lcoupled = .TRUE.

!ztmp1(:,1) (former x(:,:,12)) is never coupled
ztmp1(:,1) = -1._r8
nrcvinfo(1) = OASIS_idle


do g=begg,endg
 
           if (nint(ztmp1(g,2))==-1.or.nint(ztmp1(g,3))==-1) then
              if (nint(ztmp1(g,1)).ne.-1) then
                 tmp1 = 0.1_r8 * ztmp1(g,1)
                 tmp2 = 0.9_r8 * ztmp1(g,1)
              else
                 if(masterproc) write(6,*)'ATM error: neither PRECC/L nor PRECT'
                 if(masterproc) write(6,*)'           have been read in by atmrd'
              end if
           else
              tmp1 = ztmp1(g,2)
              tmp2 = ztmp1(g,3)
           end if

! Snow and Rain
! Set upper limit of air temperature for snowfall at 275.65K.
! This cut-off was selected based on Fig. 1, Plate 3-1, of Snow
! Hydrology (1956).

        if (tmp1 + tmp2 > 0._r8) then
              if (atm_a2l%forc_t(g) > (SHR_CONST_TKFRZ + tcrit)) then
                 atm_a2l%forc_rain(g) = tmp1 + tmp2
                 atm_a2l%forc_snow(g) = 0._r8
                 atm_a2l%flfall(g) = 1._r8
              else
                 atm_a2l%forc_rain(g) = 0._r8
                 atm_a2l%forc_snow(g) = tmp1 + tmp2
                 if (atm_a2l%forc_t(g) <= SHR_CONST_TKFRZ) then
                    atm_a2l%flfall(g) = 0._r8
                 else if (atm_a2l%forc_t(g) <= SHR_CONST_TKFRZ+2._r8) then
                    atm_a2l%flfall(g) = -0.2_r8*SHR_CONST_TKFRZ + 0.2_r8*atm_a2l%forc_t(g)
                 else
                    atm_a2l%flfall(g) = 0.4_r8
                 endif
              endif
           else
              atm_a2l%forc_rain(g) = 0._r8
              atm_a2l%forc_snow(g) = 0._r8
              atm_a2l%flfall(g) = 1._r8
           endif
enddo


!########################################

! FLWDSXY
ztmp1(:,:) = -1._r8
nrcvinfo(:) = OASIS_idle
if( srcv(jps_lw)%laction) CALL oas_clm_rcv( jps_lw, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) )

if(nrcvinfo(1)==OASIS_Rcv) then
lcoupled = .TRUE.
do g=begg,endg
          if (nint(ztmp1(g,1)) == -1) then
             e = atm_a2l%forc_psrf(g) * atm_a2l%forc_q(g) / (0.622_r8 + 0.378_r8 * atm_a2l%forc_q(g))
             ea = 0.70_r8 + 5.95e-05_r8 * 0.01_r8*e *exp(1500.0_r8/atm_a2l%forc_t(g))
             atm_a2l%forc_lwrad(g) = ea * sb *atm_a2l%forc_t(g)**4
          else
             atm_a2l%forc_lwrad(g) = ztmp1(g,1)
          end if
enddo
endif


!########################################


! FORC_PCO2
ztmp1(:,:) = -1._r8
nrcvinfo(:) = OASIS_idle
IF( srcv(jps_co2)%laction) CALL oas_clm_rcv( jps_co2, isec, ztmp1(:,1),begg,endg, nrcvinfo(1) )
IF( nrcvinfo(1)==OASIS_Rcv) then   ! CO2 partial pressu    re (Pa)
lcoupled = .TRUE.
! CO2 values from 1950 to 2100, origin https://github.com/wrf-model/WRF/blob/master/run/CAMtr_volume_mixing_ratio.SSP???
co2_ppmv_ssp126= (/312.821_r8, 313.014_r8, 313.342_r8, 313.730_r8, 314.095_r8, 314.415_r8, 314.698_r8, 314.992_r8, 315.345_r8,&
315.807_r8, 316.625_r8, 317.299_r8, 318.044_r8, 318.650_r8, 319.333_r8, 319.816_r8, 320.880_r8, 321.480_r8, 322.389_r8, 323.251_r8,&
324.783_r8, 325.400_r8, 327.349_r8, 329.909_r8, 330.756_r8, 330.827_r8, 331.545_r8, 333.353_r8, 335.010_r8, 336.605_r8, 338.705_r8,&
340.059_r8, 340.644_r8, 342.266_r8, 344.008_r8, 345.459_r8, 346.903_r8, 348.775_r8, 351.276_r8, 352.894_r8, 354.073_r8, 355.353_r8,&
356.229_r8, 356.925_r8, 358.254_r8, 360.239_r8, 362.005_r8, 363.252_r8, 365.933_r8, 367.845_r8, 369.125_r8, 370.673_r8, 372.835_r8,&
375.411_r8, 376.987_r8, 378.907_r8, 381.010_r8, 382.603_r8, 384.739_r8, 386.280_r8, 388.717_r8, 390.944_r8, 393.016_r8, 395.725_r8,&
397.547_r8, 399.949_r8, 403.117_r8, 405.752_r8, 408.592_r8, 411.417_r8, 414.232_r8, 417.042_r8, 419.809_r8, 422.499_r8, 425.120_r8,&
427.674_r8, 430.166_r8, 432.598_r8, 434.974_r8, 437.295_r8, 439.564_r8, 441.781_r8, 443.929_r8, 445.990_r8, 447.969_r8, 449.866_r8,&
451.685_r8, 453.426_r8, 455.092_r8, 456.684_r8, 458.201_r8, 459.647_r8, 461.019_r8, 462.315_r8, 463.536_r8, 464.683_r8, 465.755_r8,&
466.753_r8, 467.678_r8, 468.530_r8, 469.310_r8, 470.019_r8, 470.663_r8, 471.250_r8, 471.779_r8, 472.251_r8, 472.665_r8, 473.022_r8,&
473.322_r8, 473.564_r8, 473.750_r8, 473.878_r8, 473.959_r8, 473.999_r8, 473.999_r8, 473.957_r8, 473.873_r8, 473.747_r8, 473.578_r8,&
473.365_r8, 473.110_r8, 472.811_r8, 472.458_r8, 472.041_r8, 471.561_r8, 471.018_r8, 470.414_r8, 469.749_r8, 469.022_r8, 468.236_r8,&
467.390_r8, 466.484_r8, 465.536_r8, 464.561_r8, 463.558_r8, 462.527_r8, 461.466_r8, 460.375_r8, 459.255_r8, 458.104_r8, 456.924_r8,&
455.713_r8, 454.502_r8, 453.319_r8, 452.160_r8, 451.024_r8, 449.909_r8, 448.812_r8, 447.734_r8, 446.672_r8, 445.625_r8 /)
co2_ppmv_ssp370= (/312.821_r8, 313.014_r8, 313.342_r8, 313.730_r8, 314.095_r8, 314.415_r8, 314.698_r8, 314.992_r8, 315.345_r8,&
315.807_r8, 316.625_r8, 317.299_r8, 318.044_r8, 318.650_r8, 319.333_r8, 319.816_r8, 320.880_r8, 321.480_r8, 322.389_r8, 323.251_r8,&
324.783_r8, 325.400_r8, 327.349_r8, 329.909_r8, 330.756_r8, 330.827_r8, 331.545_r8, 333.353_r8, 335.010_r8, 336.605_r8, 338.705_r8,&
340.059_r8, 340.644_r8, 342.266_r8, 344.008_r8, 345.459_r8, 346.903_r8, 348.775_r8, 351.276_r8, 352.894_r8, 354.073_r8, 355.353_r8,&
356.229_r8, 356.925_r8, 358.254_r8, 360.239_r8, 362.005_r8, 363.252_r8, 365.933_r8, 367.845_r8, 369.125_r8, 370.673_r8, 372.835_r8,&
375.411_r8, 376.987_r8, 378.907_r8, 381.010_r8, 382.603_r8, 384.739_r8, 386.280_r8, 388.717_r8, 390.944_r8, 393.016_r8, 395.725_r8,&
397.547_r8, 399.948_r8, 403.117_r8, 405.813_r8, 408.832_r8, 411.946_r8, 415.155_r8, 418.461_r8, 421.843_r8, 425.283_r8, 428.784_r8,&
432.347_r8, 435.975_r8, 439.670_r8, 443.432_r8, 447.263_r8, 451.164_r8, 455.136_r8, 459.166_r8, 463.242_r8, 467.365_r8, 471.539_r8,&
475.765_r8, 480.045_r8, 484.380_r8, 488.772_r8, 493.221_r8, 497.731_r8, 502.292_r8, 506.901_r8, 511.557_r8, 516.264_r8, 521.021_r8,&
525.830_r8, 530.692_r8, 535.609_r8, 540.581_r8, 545.610_r8, 550.691_r8, 555.822_r8, 561.004_r8, 566.239_r8, 571.526_r8, 576.868_r8,&
582.266_r8, 587.719_r8, 593.230_r8, 598.798_r8, 604.424_r8, 610.107_r8, 615.847_r8, 621.647_r8, 627.505_r8, 633.423_r8, 639.402_r8,&
645.442_r8, 651.544_r8, 657.708_r8, 663.934_r8, 670.222_r8, 676.572_r8, 682.985_r8, 689.461_r8, 696.002_r8, 702.606_r8, 709.275_r8,&
716.008_r8, 722.808_r8, 729.679_r8, 736.629_r8, 743.657_r8, 750.764_r8, 757.951_r8, 765.216_r8, 772.561_r8, 779.984_r8, 787.488_r8,&
795.070_r8, 802.736_r8, 810.489_r8, 818.328_r8, 826.254_r8, 834.267_r8, 842.367_r8, 850.554_r8, 858.829_r8, 867.191_r8 /)
co2_ppmv_ssp585= (/312.821_r8, 313.014_r8, 313.342_r8, 313.730_r8, 314.095_r8, 314.415_r8, 314.698_r8, 314.992_r8, 315.345_r8,&
315.807_r8, 316.625_r8, 317.299_r8, 318.044_r8, 318.650_r8, 319.333_r8, 319.816_r8, 320.880_r8, 321.480_r8, 322.389_r8, 323.251_r8,&
324.783_r8, 325.400_r8, 327.349_r8, 329.909_r8, 330.756_r8, 330.827_r8, 331.545_r8, 333.353_r8, 335.010_r8, 336.605_r8, 338.705_r8,&
340.059_r8, 340.644_r8, 342.266_r8, 344.008_r8, 345.459_r8, 346.903_r8, 348.775_r8, 351.276_r8, 352.894_r8, 354.073_r8, 355.353_r8,&
356.229_r8, 356.925_r8, 358.254_r8, 360.239_r8, 362.005_r8, 363.252_r8, 365.933_r8, 367.845_r8, 369.125_r8, 370.673_r8, 372.835_r8,&
375.411_r8, 376.987_r8, 378.907_r8, 381.010_r8, 382.603_r8, 384.739_r8, 386.280_r8, 388.717_r8, 390.944_r8, 393.016_r8, 395.725_r8,&
397.547_r8, 399.949_r8, 403.117_r8, 405.794_r8, 408.760_r8, 411.790_r8, 414.891_r8, 418.065_r8, 421.335_r8, 424.721_r8, 428.220_r8,&
431.830_r8, 435.550_r8, 439.378_r8, 443.315_r8, 447.361_r8, 451.515_r8, 455.779_r8, 460.164_r8, 464.682_r8, 469.331_r8, 474.112_r8,&
479.025_r8, 484.072_r8, 489.252_r8, 494.567_r8, 500.019_r8, 505.608_r8, 511.342_r8, 517.228_r8, 523.267_r8, 529.460_r8, 535.809_r8,&
542.314_r8, 548.977_r8, 555.800_r8, 562.784_r8, 569.930_r8, 577.258_r8, 584.784_r8, 592.507_r8, 600.429_r8, 608.550_r8, 616.870_r8,&
625.391_r8, 634.114_r8, 643.041_r8, 652.172_r8, 661.507_r8, 671.045_r8, 680.789_r8, 690.739_r8, 700.899_r8, 711.269_r8, 721.851_r8,&
732.648_r8, 743.662_r8, 754.893_r8, 766.323_r8, 777.930_r8, 789.719_r8, 801.691_r8, 813.848_r8, 826.193_r8, 838.727_r8, 851.452_r8,&
864.369_r8, 877.480_r8, 890.709_r8, 903.981_r8, 917.301_r8, 930.671_r8, 944.092_r8, 957.567_r8, 971.096_r8, 984.681_r8, 998.321_r8,&
1012.020_r8, 1025.740_r8, 1039.450_r8, 1053.150_r8, 1066.850_r8, 1080.530_r8, 1094.210_r8, 1107.890_r8, 1121.550_r8, 1135.210_r8 /)
co2_ppmv_ssp=co2_ppmv_ssp370
!    if(masterproc) write(6,*)'ATM DEBUG CO2=',co2_ppmv_ssp370(year-1950+1)
do g=begg,endg   
! Ignore the values send by COSMO and use pre-defined CO2 values
! if (nint(ztmp1(g,1)) == -1) then
!             if(masterproc) write(6,*)'ATM error: PCO2 has not been read by atmrd'
!    else
!        IF ( nint(ztmp1(g,1)) == -999 ) THEN
!             atm_a2l%forc_pco2(g) = co2_ppmv_const * 1.e-6_r8 * atm_a2l%forc_pbot(g)
!           ELSE
!             atm_a2l%forc_pco2(g) = ztmp1(g,1)
!           ENDIF
!    end if
!    atm_a2l%forc_po2(g)  = o2_molar_const * atm_a2l%forc_pbot(g)
!        ! 4/14/05: PET
!        ! Adding isotope code
!    atm_a2l%forc_pc13o2(g) = co2_ppmv_const * c13ratio * 1.e-6_r8 *atm_a2l%forc_pbot(g)
    atm_a2l%forc_pco2(g) = co2_ppmv_ssp(year-1950+1) * 1.e-6_r8 * atm_a2l%forc_pbot(g)
    atm_a2l%forc_po2(g)  = o2_molar_const * atm_a2l%forc_pbot(g)
    atm_a2l%forc_pc13o2(g) = co2_ppmv_ssp(year-1950+1) * c13ratio * 1.e-6_r8 *atm_a2l%forc_pbot(g)
enddo
endif

!########################################

deallocate(ztmp1)

 !  Wait coupling result (involved PE and non involved PE)
 CALL MPI_Bcast( lcoupled, 1, dummy, 0, kl_comm, ier )

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE receive_fld_2cos
