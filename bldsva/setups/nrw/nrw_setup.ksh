#! /bin/ksh

initSetup(){
  forcingdir_clm="/homea/slts/slts06/forcings/testdata_NRW_std"
#  forcingdir_cos="/work/slts/slts06/tsmp/TSMPForecastNRW2015-11-09-00/forcing/cosmoinput"
  forcingdir_cos="/work/slts/slts00/tsmp/TerrSysMPdb/testdata_NRW_std/cosmo/int2lm_output"
  forcingdir_oas="/homea/slts/slts06/forcings/testdata_NRW_std/oasis3/"
  forcingdir_pfl="/homea/slts/slts06/forcings/testdata_NRW_std/ParFlow/Rur_NRW/"

  namelist_clm=$rootdir/bldsva/setups/nrw/lnd.stdin
  namelist_cos=$rootdir/bldsva/setups/nrw/lmrun_uc
#  namelist_cos=~/debug/tsmp_full_old_jq/cosmo/lmrun_uc
  namelist_pfl=$rootdir/bldsva/setups/nrw/coup_oas.tcl

  gx_clm=300
  gy_clm=300
  dt_clm=900

  gx_cos=150
  gy_cos=150
  dt_cos=10

  gx_pfl=300
  gy_pfl=300
  res="0300x0300"
  dt_pfl=0.25

  cplfreq1=900
  cplfreq2=900

}

finalizeSetup(){

  cp $rootdir/bldsva/machines/$platform/loadenvs $rundir
  cp $forcingdir_clm/clm3.5/Rur_NRW/grid* $rundir/clmgrid.nc
  
  if [[ $withOAS == "true" ]] then
    cp $forcingdir_oas/* $rundir
    if [[ $withOASMCT == "true" ]] then
      for x in $rundir/*BILINEA* ;do 
        mv $x $(echo $x | sed "s/BILINEA/BILINEAR/") 
      done
    fi  
  fi  

  if [[ $withPFL == "true" ]] then
        cd $rundir
        echo "preprocessing of ParFlow static fields (distribution)"

        cp $forcingdir_pfl/slopes/ascii2pfb.tcl.template_new $rundir/ascii2pfb.tcl
        cp $forcingdir_pfl/slopes/*slope.pfb* $rundir
        chmod u+w $rundir/*slope*  $rundir/ascii2pfb.tcl
        sed "s,__svaroot__.*,$pfldir/bin," -i $rundir/ascii2pfb.tcl
        sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb.tcl
        sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb.tcl
        tclsh ./ascii2pfb.tcl

        cp $forcingdir_pfl/soilInd/ascii2pfb.tcl.template_new $rundir/ascii2pfb.tcl
        cp $forcingdir_pfl/soilInd/*Soil* $rundir
        chmod u+w $rundir/*Soil* $rundir/ascii2pfb.tcl
        sed "s,__svaroot__.*,$pfldir/bin," -i $rundir/ascii2pfb.tcl
        sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb.tcl
        sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb.tcl
        tclsh ./ascii2pfb.tcl
  fi 
}
