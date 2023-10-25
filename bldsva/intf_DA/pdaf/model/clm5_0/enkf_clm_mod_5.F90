!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
!
!This file is part of TerrSysMP-PDAF
!
!TerrSysMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TerrSysMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!enkf_clm_mod.F90: Module for CLM
!-------------------------------------------------------------------------------------------

module enkf_clm_mod

  use iso_c_binding

  use shr_kind_mod    , only : r8 => shr_kind_r8, SHR_KIND_CL

#if (defined CLMFIVE)
  integer :: da_comm_clm
  integer :: clm_statevecsize
  integer :: clm_varsize
  integer :: clm_begg,clm_endg,clm_begc,clm_endc
  real(r8),allocatable :: clm_statevec(:)
  real(r8),allocatable :: clm_paramarr(:)  !hcp LAI
  integer(c_int),bind(C,name="clmupdate_swc")     :: clmupdate_swc
  integer(c_int),bind(C,name="clmupdate_T")     :: clmupdate_T  ! by hcp
  integer(c_int),bind(C,name="clmupdate_texture") :: clmupdate_texture
  integer(c_int),bind(C,name="clmupdate_snow")    :: clmupdate_snow
  integer(c_int),bind(C,name="clmupdate_snow_repartitioning") :: clmupdate_snow_repartitioning
  integer(c_int),bind(C,name="clmprint_swc")      :: clmprint_swc
#endif
  integer(c_int),bind(C,name="clmprint_et")       :: clmprint_et

  integer  :: nstep     ! time step index
  real(r8) :: dtime     ! time step increment (sec)
  integer  :: ier       ! error code

  character(c_char),bind(C,name="outdir"),target :: outdir

  logical  :: log_print    ! true=> print diagnostics
  real(r8) :: eccf         ! earth orbit eccentricity factor
  logical  :: mpi_running  ! true => MPI is initialized 
  integer  :: mpicom_glob  ! MPI communicator

  character(len=SHR_KIND_CL) :: nlfilename = " "
  integer :: ierror, lengths_of_types, i
  logical :: flag
  integer(c_int),bind(C,name="clmprefixlen") :: clmprefixlen
  integer :: statcomm

  logical :: newgridcell
  contains

#if defined CLMFIVE 
  subroutine define_clm_statevec()
    use shr_kind_mod, only: r8 => shr_kind_r8
    use decompMod , only : get_proc_bounds
    use clm_varpar   , only : nlevsoi

    implicit none

    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices


    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    !write(*,*) "----",begg,",",endg,",",begl,",",endl,",",begc,",",endc,",",begp,",",endp," -------"
    clm_begg     = begg
    clm_endg     = endg
    clm_begc     = begc
    clm_endc     = endc

    if(clmupdate_swc.eq.1) then
      clm_varsize      =  (clm_endg-clm_begg+1) * nlevsoi
      clm_statevecsize =  (clm_endg-clm_begg+1) * nlevsoi
    endif

    if(clmupdate_swc.eq.2) then
      error stop "Not implemented swc update 2"
    endif

    if(clmupdate_texture.eq.1) then
        clm_statevecsize = clm_statevecsize + 2*((endg-begg+1)*nlevsoi)
    endif

    if(clmupdate_texture.eq.2) then
        clm_statevecsize = clm_statevecsize + 3*((endg-begg+1)*nlevsoi)
    endif

    ! Snow assimilation
    ! Case 1: Assimilation of snow depth : allocated 1 per column in CLM5
    ! But observations and history file 1 per grid cell and therefore statevecsize 1 per grid cell
    if(clmupdate_snow.eq.1) then 
        clm_varsize      =  (clm_endg-clm_begg+1) ! Currently no combination of SWC and snow DA
        clm_statevecsize =  (clm_endg-clm_begg+1) ! So like this if snow is set it takes priority
    endif

    IF (allocated(clm_statevec)) deallocate(clm_statevec)
    allocate(clm_statevec(clm_statevecsize))

    !write(*,*) 'clm_paramsize is ',clm_paramsize
    IF (allocated(clm_paramarr)) deallocate(clm_paramarr)         !hcp
    ! if ((clmupdate_T.NE.0)) allocate(clm_paramarr(clm_paramsize))  !hcp
    if ((clmupdate_T.NE.0)) error stop "Not implemented clmupdate_T.NE.0"

  end subroutine

  subroutine set_clm_statevec()
    use clm_instMod, only : soilstate_inst, waterstate_inst
    use clm_varpar   , only : nlevsoi
    use ColumnType , only : col
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    real(r8), pointer :: snow_depth(:)

    integer :: i,j,jj,g,cc=1,offset=0

    swc   => waterstate_inst%h2osoi_vol_col
    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm => soilstate_inst%cellorg_col

    snow_depth => waterstate_inst%snow_depth_col ! snow height of snow covered area (m) 

    if(clmupdate_swc.ne.0) then
        ! write swc values to state vector
        cc = 1
        do i=1,nlevsoi
          do j=clm_begg,clm_endg
            ! Only get the SWC from the first column of each gridcell
            ! and add it to the clm_statevec at the position of the gridcell (cc)
            newgridcell = .true.
            do jj=clm_begc,clm_endc
              g = col%gridcell(jj)
              if (g .eq. j) then
                if (newgridcell) then
                  newgridcell = .false.
                  clm_statevec(cc+offset) = swc(jj,i)
                endif
              endif 
            end do
            cc = cc + 1
          end do
        end do
    endif

    ! write texture values to state vector (if desired)
    if(clmupdate_texture.eq.1) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
            ! Only get the parameters from the first column of each gridcell
            ! and add it to the clm_statevec at the position of the gridcell (cc)
            newgridcell = .true.
            do jj=clm_begc,clm_endc
              g = col%gridcell(jj)
              if (g .eq. j) then
                if (newgridcell) then
                  newgridcell = .false.
                  clm_statevec(cc+1*clm_varsize+offset) = psand(jj,i)
                  clm_statevec(cc+2*clm_varsize+offset) = pclay(jj,i)
                endif
              endif 
            end do
            cc = cc + 1
        end do
      end do
    endif

    ! write texture incl. organic matter values to state vector (if desired)
    if(clmupdate_texture.eq.2) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
            ! Only get the parameters from the first column of each gridcell
            ! and add it to the clm_statevec at the position of the gridcell (cc)
            newgridcell = .true.
            do jj=clm_begc,clm_endc
              g = col%gridcell(jj)
              if (g .eq. j) then
                if (newgridcell) then
                  newgridcell = .false.
                  clm_statevec(cc+1*clm_varsize+offset) = psand(jj,i)
                  clm_statevec(cc+2*clm_varsize+offset) = pclay(jj,i)
                  clm_statevec(cc+3*clm_varsize+offset) = porgm(jj,i)
                endif
              endif 
            end do
            cc = cc + 1
        end do
      end do
    endif

    ! Snow assimilation
    ! Case 1: Snow depth
    if(clmupdate_snow.eq.1) then
        cc = 1
        do j=clm_begg,clm_endg
        ! Only get the snow_depth from the first column of each gridcell
        ! and add it to the clm_statevec at the position of the gridcell (cc)
        newgridcell = .true.
        do jj=clm_begc,clm_endc
          g = col%gridcell(jj)
          if (g .eq. j) then
            if (newgridcell) then
              newgridcell = .false.
              clm_statevec(cc+offset) = snow_depth(jj)
            endif
          endif 
        end do
        cc = cc + 1
        end do
    endif

  end subroutine 

  subroutine update_clm() bind(C,name="update_clm")
    use clm_varpar   , only : nlevsoi
    use clm_time_manager  , only : update_DA_nstep
    use shr_kind_mod , only : r8 => shr_kind_r8
    use ColumnType , only : col
    use clm_instMod, only : soilstate_inst, waterstate_inst
    use clm_varcon      , only : denh2o, denice, watmin

    implicit none
    
    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: watsat(:,:)
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    real(r8), pointer :: snow_depth(:)

    real(r8), pointer :: dz(:,:)          ! layer thickness depth (m)
    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)
    real(r8)  :: rliq,rice
    real(r8) :: rsnow(clm_statevecsize)
    integer :: i,j,jj,g,cc=1,offset=0

    swc   => waterstate_inst%h2osoi_vol_col
    watsat => soilstate_inst%watsat_col
    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm => soilstate_inst%cellorg_col

    snow_depth => waterstate_inst%snow_depth_col ! snow height of snow covered area (m) 

    dz            => col%dz
    h2osoi_liq    => waterstate_inst%h2osoi_liq_col
    h2osoi_ice    => waterstate_inst%h2osoi_ice_col

    call update_DA_nstep()

    ! write updated swc back to CLM
    if(clmupdate_swc.ne.0) then
        cc = 1
        do i=1,nlevsoi
          do j=clm_begg,clm_endg
            ! iterate through the columns and copy from the same gridcell
            ! i.e. statevec position (cc) for each column
            do jj=clm_begc,clm_endc
              rliq = h2osoi_liq(jj,i)/(dz(jj,i)*denh2o*swc(jj,i))
              rice = h2osoi_ice(jj,i)/(dz(jj,i)*denice*swc(jj,i))
     
              if(clm_statevec(cc+offset).le.watmin) then
                swc(jj,i)   = watmin
              else if(clm_statevec(cc+offset).ge.watsat(jj,i)) then
                swc(jj,i) = watsat(jj,i)
              else
                swc(jj,i)   = clm_statevec(cc+offset)
              endif

              if (isnan(swc(jj,i))) then
                      swc(jj,i) = watmin
                      print *, "WARNING: swc at j,i is nan: ", jj, i
              endif

              ! update liquid water content
              h2osoi_liq(jj,i) = swc(jj,i) * dz(jj,i)*denh2o*rliq
              ! update ice content
              h2osoi_ice(j,i) = swc(j,i) * dz(j,i)*denice*rice
            end do
            cc = cc + 1
          end do
        end do
    endif 

    ! write updated texture back to CLM
    if(clmupdate_texture.eq.1) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
            ! iterate through the columns and copy from the same gridcell
            ! i.e. statevec position (cc) for each column
            do jj=clm_begc,clm_endc
                psand(jj,i) = clm_statevec(cc+1*clm_varsize+offset)
                pclay(jj,i) = clm_statevec(cc+2*clm_varsize+offset)
            end do
          cc = cc + 1
        end do
      end do
      call clm_correct_texture
      call clm_texture_to_parameters
    endif

    ! write updated texture incl. organic matter back to CLM
    if(clmupdate_texture.eq.2) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
            ! iterate through the columns and copy from the same gridcell
            ! i.e. statevec position (cc) for each column
            do jj=clm_begc,clm_endc
                psand(jj,i) = clm_statevec(cc+1*clm_varsize+offset)
                pclay(jj,i) = clm_statevec(cc+2*clm_varsize+offset)
                porgm(jj,i) = clm_statevec(cc+3*clm_varsize+offset)
            end do
          cc = cc + 1
        end do
      end do
      call clm_correct_texture
      call clm_texture_to_parameters
    endif

    ! Snow assimilation:
    ! Case 1: Snow depth
    ! Write updated snow depth back to CLM and then repartition snow and adjust related variables
    if(clmupdate_snow.eq.1) then
        cc = 1
        do j=clm_begg,clm_endg
        ! iterate through the columns and copy from the same gridcell
        ! i.e. statevec position (cc) for each column
          do jj=clm_begc,clm_endc
              ! Catch negative or 0 values from DA
              if (clm_statevec(cc+offset).lt.0.0) then
                print *, "WARNING: snow depth at g,c is negative: ", j, jj, clm_statevec(cc+offset)
              else
                rsnow(jj) = snow_depth(jj) - clm_statevec(cc+offset)
                snow_depth(jj)   = clm_statevec(cc+offset)
              endif
          end do
        cc = cc + 1
        end do
        if ( clmupdate_snow_repartitioning.ne.0) then
          if ( ABS(SUM(rsnow(:))).gt.0.000001) then
            call clm_repartition_snow()
          end if
        end if
    endif

  end subroutine

  subroutine clm_repartition_snow()
    use ColumnType,    only : col
    use clm_instMod,   only : waterstate_inst
    use clm_varpar,    only : nlevsno, nlevsoi
    use clm_varcon,    only : bdsno, denice
    use shr_kind_mod,  only : r8 => shr_kind_r8

    implicit none

    real(r8), pointer :: snow_depth(:)
    real(r8), pointer :: h2osno(:)
    real(r8), pointer :: h2oliq(:,:)
    real(r8), pointer :: h2oice(:,:)
    real(r8), pointer :: frac_sno(:)
    real(r8), pointer :: snowdp(:)
    integer, pointer :: snlsno(:)

    real(r8) :: dzsno(clm_begc:clm_endc,nlevsno)
    real(r8) :: h2osno_po(clm_begc:clm_endc)
    real(r8) :: snowden, frac_swe, frac_liq, frac_ice
    real(r8) :: gain_h2osno, gain_h2oliq, gain_h2oice, gain_dzsno
    real(r8) :: rho_avg, z_avg
    integer :: i,ii,j,jj,g,cc=1,offset=0

    snow_depth => waterstate_inst%snow_depth_col ! snow height of snow covered area (m)
    snowdp     => waterstate_inst%snowdp_col     ! area-averaged snow height (m) 
    h2osno     => waterstate_inst%h2osno_col     ! col snow water (mm H2O)
    h2oliq     => waterstate_inst%h2osoi_liq_col ! col liquid water (kg/m2) (-nlevsno+1:nlevgrnd) 
    h2oice     => waterstate_inst%h2osoi_ice_col ! col ice lens (kg/m2) (-nlevsno+1:nlevgrnd)

    snlsno     => col%snl                        ! number of snow layers (negative)

    frac_sno   => waterstate_inst%frac_sno_eff_col ! fraction of ground covered by snow 
    ! dz for snow layers is defined like in the history output as col%dz for the snow layers
    dzsno(clm_begc:clm_endc, -nlevsno+1:0) = col%dz(clm_begc:clm_endc,-nlevsno+1:0)

    ! Iterate through all columns
    do jj=clm_begc,clm_endc
      if (h2osno(jj).lt.0.0) then ! No snow in column
        print *, "WARNING: negative snow in col: ", jj, h2osno
!        ! Set existing layers to near zero and let CLM do the layer aggregation
        do i=0,snlsno(jj),-1
            h2oliq(jj,i) = 0.0_r8
            h2oice(jj,i) = 0.00000001_r8
            dzsno(jj,i)  = 0.00000001_r8
        end do
        snow_depth(jj) = sum(dzsno(jj,:))
        h2osno(jj) = sum(h2oice(jj,:))
      else ! snow (h2osno) present
        if (snlsno(jj).lt.0) then ! snow layers in the column
          ! Formulas below from DART use h2osno_po / h2osno_pr for after / before DA SWE
          ! Here we have snow_depth after and h2osno before DA snow depth
          ! Therefore need to have a transform to get h2osno_po
          ! v1 init
          ! h2osno_po(jj) = (snow_depth(jj) * bdsno) ! calculations from Init using constant SBD
          ! v2 SoilTemperatureMod
          if (snowdp(jj).gt.0.0_r8) then
            rho_avg = min(800.0_r8, h2osno(jj)/snowdp(jj))
          else
            rho_avg = 200.0_r8
          end if
          if (frac_sno(jj).gt.0.0_r8 .and. snlsno(jj).lt.0.0_r8) then
            h2osno_po(jj) = snow_depth(jj) * (rho_avg*frac_sno(jj))
          else
            h2osno_po(jj) = snow_depth(jj) * denice
          end if

          do ii=0,-snlsno(jj),-1 ! iterate through the snow layers
            ! ii = nlevsoi - i + 1            ! DART VERSION: ii = nlevsoi - i + 1
            ! snow density prior for each layer
            if (dzsno(jj,ii).gt.0.0_r8) then
              snowden = (h2oliq(jj,ii) + h2oice(jj,ii) / dzsno(jj,ii))
            else
              snowden = 0.0_r8
            endif
            ! fraction of SWE in each active layers
            if(h2osno(jj).gt.0.0_r8) then
              ! repartition Option 1: Adjust bottom layer only (set weight to 1 for bottom 0 else)
              if(clmupdate_snow_repartitioning.eq.1) then
                if (ii .eq. 0) then ! DART version indexing check against nlevsno but for us 0
                  frac_swe = 1.0_r8
                else
                  frac_swe = 0.0_r8
                end if
              ! repartition Option 2: Adjust all active layers
              else if (clmupdate_snow_repartitioning.eq.2) then
                frac_swe = (h2oliq(jj,ii) + h2oice(jj,ii)) / h2osno(jj)
              end if
            else
              frac_swe = 0.0_r8 ! no fraction SWE if no snow is present in column
            end if ! end SWE fraction if

            ! fraction of liquid and ice
            if ((h2oliq(jj,ii) + h2oice(jj,ii)).gt.0.0_r8) then
                frac_liq = h2oliq(jj,ii) / (h2oliq(jj,ii) + h2oice(jj,ii))
                frac_ice = 1.0_r8 - frac_liq
            else
                frac_liq = 0.0_r8
                frac_ice = 0.0_r8
            end if

            ! SWE adjustment per layer 
            ! assumes identical layer distribution of liq and ice than before DA (frac_*)
            if (abs(h2osno_po(jj) - h2osno(jj)).gt.0.0_r8) then
              gain_h2osno = (h2osno_po(jj) - h2osno(jj)) * frac_swe
              gain_h2oliq = gain_h2osno * frac_liq
              gain_h2oice = gain_h2osno * frac_ice
            else
              gain_h2osno = 0.0_r8
              gain_h2oliq = 0.0_r8
              gain_h2oice = 0.0_r8
            end if
            ! layer level adjustments
            if (snowden.gt.0.0_r8) then
              gain_dzsno = gain_h2osno / snowden
            else
              gain_dzsno = 0.0_r8
            end if
            h2oliq(jj,ii) = h2oliq(jj,ii) + gain_h2oliq
            h2oice(jj,ii) = h2oice(jj,ii) + gain_h2oice

            ! Adjust snow layer dimensions so that CLM5 can calculate compaction / aggregation
            ! in the DART code dzsno is adjusted directly but in CLM5 dzsno is local and diagnostic
            ! i.e. calculated / assigned from frac_sno and dz(:, snow_layer) in SnowHydrologyMod
            ! therefore we adjust dz(:, snow_layer) here
            if (abs(h2osno_po(jj) - h2osno(jj)).gt.0.0_r8) then
              col%dz(jj, ii) = col%dz(jj, ii) + gain_dzsno
              ! mid point and interface adjustments
              ! i.e. zsno (col%z(:, snow_layers)) and zisno (col%zi(:, snow_layers))
              ! DART version the sum goes from ilevel:nlevsno to fit with our indexing:
              col%zi(jj, ii) = sum(col%dz(jj,ii:0))*-1.0_r8
              ! In DART the check is ilevel == nlevsno but here
              if (ii.eq.0) then
                col%z(jj,ii) = col%zi(jj,ii) / 2.0_r8
              else
                col%z(jj,ii) = sum(col%zi(jj, ii:ii+1)) / 2.0_r8
              end if
            end if

          end do ! End iteration of snow layers
        endif ! End of snow layers present check
        ! Finally adjust SWE (h2osno) since the prior value is no longer needed
        ! column level variables so we can adjust it outside the layer loop
        h2osno(jj) = h2osno_po(jj) 
      end if ! End of snow present check
    end do ! End of column iteration

  end subroutine

  subroutine clm_correct_texture()

    use clm_varpar   , only : nlevsoi
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_instMod, only : soilstate_inst
    use CNSharedParamsMod, only : CNParamsShareInst

    implicit none

    integer :: c,lev
    real(r8) :: clay,sand,orgm,ttot
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm => soilstate_inst%cellorg_col

    do c = clm_begg,clm_endg
      do lev = 1,nlevsoi
         clay = pclay(c,lev)
         sand = psand(c,lev)

         if(sand.le.0.0) sand = 1.0
         if(clay.le.0.0) clay = 1.0

         ttot = sand + clay
         if(ttot.gt.100) then
             sand = sand/ttot * 100.0
             clay = clay/ttot * 100.0
         end if

         pclay(c,lev) = clay
         psand(c,lev) = sand

         if(clmupdate_texture.eq.2) then
           orgm = (porgm(c,lev) / CNParamsShareInst%organic_max) * 100.0
           if(orgm.le.0.0) orgm = 0.0
           if(orgm.ge.100.0) orgm = 100.0
           porgm(c,lev) = (orgm / 100.0)* CNParamsShareInst%organic_max
         end if

      end do
    end do
  end subroutine

  subroutine clm_texture_to_parameters()
    use clm_varpar   , only : nlevsoi
    use clm_varcon   , only : zsoi, secspday
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_instMod, only : soilstate_inst
    use FuncPedotransferMod , only : pedotransf, get_ipedof
    use CNSharedParamsMod, only : CNParamsShareInst
    implicit none

    real(r8)           :: om_tkm         = 0.25_r8      
    ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(r8)           :: om_watsat_lake = 0.9_r8       
    ! porosity of organic soil
    real(r8)           :: om_hksat_lake  = 0.1_r8       
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat_lake = 10.3_r8      
    ! saturated suction for organic matter (Letts, 2000)
    real(r8)           :: om_b_lake      = 2.7_r8       
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000) (lake)
    real(r8)           :: om_watsat                     
    ! porosity of organic soil
    real(r8)           :: om_hksat                      
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat                     
    ! saturated suction for organic matter (mm)(Letts, 2000)
    real(r8)           :: om_csol        = 2.5_r8       
    ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(r8)           :: om_tkd         = 0.05_r8      
    ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(r8)           :: om_b                          
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8)           :: zsapric        = 0.5_r8
    ! depth (m) that organic matter takes on characteristics of sapric peat
    real(r8)           :: pcalpha        = 0.5_r8       ! percolation threshold
    real(r8)           :: pcbeta         = 0.139_r8     ! percolation exponent
    real(r8)           :: perc_frac    ! "percolating" fraction of organic soil
    real(r8)           :: perc_norm    ! normalize to 1 when 100% organic soil
    real(r8)           :: uncon_hksat  ! series conductivity of mineral/organic soil
    real(r8)           :: uncon_frac   ! fraction of "unconnected" soil
    real(r8)           :: bd           ! bulk density of dry soil material [kg/m^3]
    real(r8)           :: tkm          ! mineral conductivity 
    real(r8)           :: xksat        ! maximum hydraulic conductivity of soil [mm/s]
    
    integer  :: ipedof,c,lev
    real(r8) :: clay,sand,om_frac
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm   => soilstate_inst%cellorg_col

    do c = clm_begg,clm_endg
      do lev = 1,nlevsoi
         clay = pclay(c,lev)
         sand = psand(c,lev)
         om_frac = porgm(c,lev) / CNParamsShareInst%organic_max

         ipedof=get_ipedof(0)
         call pedotransf(ipedof, sand, clay, &
                         soilstate_inst%watsat_col(c,lev), &
                         soilstate_inst%bsw_col(c,lev), &
                         soilstate_inst%sucsat_col(c,lev), &
                         xksat)

         om_watsat         = max(0.93_r8 - 0.1_r8   *(zsoi(lev)/zsapric), 0.83_r8)
         om_b              = min(2.7_r8  + 9.3_r8   *(zsoi(lev)/zsapric), 12.0_r8)
         om_sucsat         = min(10.3_r8 - 0.2_r8   *(zsoi(lev)/zsapric), 10.1_r8)
         om_hksat          = max(0.28_r8 - 0.2799_r8*(zsoi(lev)/zsapric), xksat)

         soilstate_inst%bd_col(c,lev) = &
              (1._r8 - soilstate_inst%watsat_col(c,lev))*2.7e3_r8
         
         soilstate_inst%watsat_col(c,lev) = &
              (1._r8 - om_frac) * soilstate_inst%watsat_col(c,lev) + om_watsat*om_frac
         
         tkm = (1._r8 - om_frac) * (8.80_r8*sand+2.92_r8*clay)/(sand+clay)+om_tkm*om_frac 
         ! W/(m K)
         soilstate_inst%bsw_col(c,lev) = &
              (1._r8-om_frac) * (2.91_r8 + 0.159_r8*clay) + om_frac*om_b   
         
         soilstate_inst%sucsat_col(c,lev) = &
              (1._r8-om_frac) * soilstate_inst%sucsat_col(c,lev) + om_sucsat*om_frac  
         
         soilstate_inst%hksat_min_col(c,lev) = xksat

         ! perc_frac is zero unless perf_frac greater than percolation threshold
         if (om_frac > pcalpha) then
            perc_norm = (1._r8 - pcalpha)**(-pcbeta)
            perc_frac = perc_norm*(om_frac - pcalpha)**pcbeta
         else
            perc_frac = 0._r8
         endif

         ! uncon_frac is fraction of mineral soil plus fraction of
         ! "nonpercolating" organic soil
         uncon_frac = (1._r8-om_frac)+(1._r8-perc_frac)*om_frac

         ! uncon_hksat is series addition of mineral/organic
         ! conductivites
         if (om_frac < 1._r8) then
            uncon_hksat = uncon_frac/((1._r8-om_frac)/xksat &
                 +((1._r8-perc_frac)*om_frac)/om_hksat)
         else
            uncon_hksat = 0._r8
         end if

         soilstate_inst%hksat_col(c,lev)  = uncon_frac*uncon_hksat + &
             (perc_frac*om_frac)*om_hksat

         soilstate_inst%tkmg_col(c,lev)   = tkm ** (1._r8 - &
             soilstate_inst%watsat_col(c,lev))           

         soilstate_inst%tksatu_col(c,lev) = &
             soilstate_inst%tkmg_col(c,lev)*0.57_r8**soilstate_inst%watsat_col(c,lev)

         soilstate_inst%tkdry_col(c,lev)  = &
             ((0.135_r8*soilstate_inst%bd_col(c,lev) + 64.7_r8) / &
             (2.7e3_r8 - 0.947_r8*soilstate_inst%bd_col(c,lev)) &
             )*(1._r8-om_frac) + om_tkd*om_frac  

         soilstate_inst%csol_col(c,lev)   = &
             ((1._r8-om_frac)*(2.128_r8*sand+2.385_r8*clay) / (sand+clay) + &
             om_csol*om_frac)*1.e6_r8

         soilstate_inst%watdry_col(c,lev) = &
             soilstate_inst%watsat_col(c,lev) * &
             (316230._r8/soilstate_inst%sucsat_col(c,lev)  &
             ) ** (-1._r8/soilstate_inst%bsw_col(c,lev)) 

         soilstate_inst%watopt_col(c,lev) = & 
             soilstate_inst%watsat_col(c,lev) * &
             (158490._r8/soilstate_inst%sucsat_col(c,lev)  &
             ) ** (-1._r8/soilstate_inst%bsw_col(c,lev)) 

         ! secspday (day/sec)
         soilstate_inst%watfc_col(c,lev) = &
             soilstate_inst%watsat_col(c,lev) * &
             (0.1_r8 / (soilstate_inst%hksat_col(c,lev)*secspday) &
             )**(1._r8/(2._r8*soilstate_inst%bsw_col(c,lev)+3._r8))

      end do
    end do
  end subroutine
 
  subroutine  average_swc_crp(profdat,profave)
    error stop "Not implemented average_swc_crp"
  end subroutine average_swc_crp
#endif

  subroutine domain_def_clm(lon_clmobs, lat_clmobs, dim_obs, &
                            longxy, latixy, longxy_obs, latixy_obs)
    use spmdMod,   only : npes, iam
    use domainMod, only : ldomain
    use decompMod, only : get_proc_total, get_proc_bounds, ldecomp
    real, intent(in) :: lon_clmobs(:)
    real, intent(in) :: lat_clmobs(:)
    integer, intent(in) :: dim_obs
    integer, allocatable, intent(inout) :: longxy(:)
    integer, allocatable, intent(inout) :: latixy(:)
    integer, allocatable, intent(inout) :: longxy_obs(:)
    integer, allocatable, intent(inout) :: latixy_obs(:)
    integer :: ni, nj, ii, jj, kk, cid, ier, ncells, nlunits, &
               ncols, npatches, ncohorts, counter
    real :: minlon, minlat, maxlon, maxlat
    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices


    lon => ldomain%lonc
    lat => ldomain%latc
    ni = ldomain%ni
    nj = ldomain%nj

    ! get total number of gridcells, landunits,
    ! columns, patches and cohorts on processor

    call get_proc_total(iam, ncells, nlunits, ncols, npatches, ncohorts)

    ! beg and end gridcell
    call get_proc_bounds(begg=begg, endg=endg)
    
    allocate(longxy(ncells), stat=ier)
    allocate(latixy(ncells), stat=ier)

    longxy(:) = 0
    latixy(:) = 0
  
    counter = 1
    do ii = 1, nj
      do jj = 1, ni
        cid = (ii-1)*ni + jj
        do kk = begg, endg
          if(cid == ldecomp%gdc2glo(kk)) then
            latixy(counter) = ii
            longxy(counter) = jj
            counter = counter + 1
          end if
        end do
      end do
    end do

    ! cell centers to min/max longitude and latitude
    minlon = MINVAL(lon(:) + 180)
    maxlon = MAXVAL(lon(:) + 180)
    minlat = MINVAL(lat(:) + 90)
    maxlat = MAXVAL(lat(:) + 90)

    allocate(longxy_obs(dim_obs), stat=ier)
    allocate(latixy_obs(dim_obs), stat=ier)

    do i = 1, dim_obs
       if(((lon_clmobs(i) + 180) - minlon) /= 0 .and. &
          ((lat_clmobs(i) + 90) - minlat) /= 0) then
          longxy_obs(i) = ceiling(((lon_clmobs(i)+180) - minlon) * ni / (maxlon-minlon))
          latixy_obs(i) = ceiling(((lat_clmobs(i)+90) - minlat) * nj / (maxlat-minlat))
       else if(((lon_clmobs(i) + 180) - minlon) == 0 .and. &
               ((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs(i) = 1
          latixy_obs(i) = 1
       else if(((lon_clmobs(i) + 180) - minlon) == 0) then
          longxy_obs(i) = 1
          latixy_obs(i) = ceiling(((lat_clmobs(i)+90) - minlat) * nj / (maxlat-minlat))
       else if(((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs(i) = ceiling(((lon_clmobs(i)+180) - minlon) * ni / (maxlon-minlon))
          latixy_obs(i) = 1
       endif
    end do
  end subroutine  

  !> @author  Mukund Pondkule, Johannes Keller
  !> @date    27.03.2023
  !> @brief   Set indices of grid cells with lon/lat smaller than observation locations
  !> @details
  !>    This routine sets the indices of grid cells with lon/lat
  !>    smaller than observation locations.
  subroutine get_interp_idx(lon_clmobs, lat_clmobs, dim_obs, longxy_obs_floor, latixy_obs_floor)

    USE domainMod, ONLY: ldomain
    ! USE decompMod, ONLY: get_proc_total, get_proc_bounds_atm, adecomp
    ! USE spmdMod,   ONLY: npes, iam ! number of processors for clm and processor number
   ! USE mod_read_obs, ONLY: lon_clmobs, lat_clmobs
    ! USE clmtype,    ONLY : clm3
   ! USE mod_assimilation, ONLY: dim_obs
!#endif

    implicit none
    real, intent(in) :: lon_clmobs(:)
    real, intent(in) :: lat_clmobs(:)
    integer, intent(in) :: dim_obs
    integer, allocatable, intent(inout) :: longxy_obs_floor(:)
    integer, allocatable, intent(inout) :: latixy_obs_floor(:)
    integer :: i

    integer :: ni, nj
    ! integer :: ii, jj, kk, cid
    integer :: ier
    ! integer :: ncells, nlunits,ncols, npfts
    ! integer :: counter

    real :: minlon, minlat, maxlon, maxlat
    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)
    ! integer :: begg, endg   ! per-proc gridcell ending gridcell indices

    lon => ldomain%lonc
    lat => ldomain%latc
    ni = ldomain%ni
    nj = ldomain%nj

    ! ! get total number of gridcells, landunits,
    ! ! columns and pfts for any processor
    ! call get_proc_total(iam, ncells, nlunits, ncols, npfts)

    ! ! beg and end gridcell for atm
    ! call get_proc_bounds_atm(begg, endg)

    ! set intial values for max/min of lon/lat
    minlon = 999
    minlat = 999
    maxlon = -999
    maxlat = -999

    ! looping over all cell centers to get min/max longitude and latitude
    minlon = MINVAL(lon(:) + 180)
    maxlon = MAXVAL(lon(:) + 180)
    minlat = MINVAL(lat(:) + 90)
    maxlat = MAXVAL(lat(:) + 90)

    if(allocated(longxy_obs_floor)) deallocate(longxy_obs_floor)
    allocate(longxy_obs_floor(dim_obs), stat=ier)
    if(allocated(latixy_obs_floor)) deallocate(latixy_obs_floor)
    allocate(latixy_obs_floor(dim_obs), stat=ier)
    do i = 1, dim_obs
       if(((lon_clmobs(i) + 180) - minlon) /= 0 .and. ((lat_clmobs(i) + 90) - minlat) /= 0) then
          longxy_obs_floor(i) = floor(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon)) !+ 1
          latixy_obs_floor(i) = floor(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat)) !+ 1
          !print *,'longxy_obs(i) , latixy_obs(i) ', longxy_obs(i) , latixy_obs(i)
       else if(((lon_clmobs(i) + 180) - minlon) == 0 .and. ((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs_floor(i) = 1
          latixy_obs_floor(i) = 1
       else if(((lon_clmobs(i) + 180) - minlon) == 0) then
          longxy_obs_floor(i) = 1
          latixy_obs_floor(i) = floor(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat))
       else if(((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs_floor(i) = floor(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon))
          latixy_obs_floor(i) = 1
       endif
    end do

  end subroutine get_interp_idx

#if defined CLMFIVE
  subroutine init_clm_l_size(dim_l)
    use clm_varpar   , only : nlevsoi

    implicit none

    integer, intent(out) :: dim_l
    integer              :: nshift

    if(clmupdate_swc.eq.1) then
      dim_l = nlevsoi
      nshift = nlevsoi
    endif

    if(clmupdate_swc.eq.2) then
      dim_l = nlevsoi + 1
      nshift = nlevsoi + 1
      error stop "Not implemented swc update 2"
    endif

    if(clmupdate_texture.eq.1) then
      dim_l = 2*nlevsoi + nshift
    endif

    if(clmupdate_texture.eq.2) then
      dim_l = 3*nlevsoi + nshift
    endif

  end subroutine    
#endif

end module enkf_clm_mod

