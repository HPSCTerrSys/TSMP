#! /bin/ksh

initSetup(){
  if [[ $forcingdir_clm == "" ]] ; then ;forcingdir_clm="/work/slts/slts00/tsmp/TerrSysMPdb/testdata_NRW_std";fi
  if [[ $forcingdir_cos == "" ]] ; then ;forcingdir_cos="/work/slts/slts00/tsmp/TerrSysMPdb/testdata_NRW_std/cosmo/int2lm_output.20080508";fi
  if [[ $forcingdir_oas == "" ]] ; then ;forcingdir_oas="/work/slts/slts00/tsmp/TerrSysMPdb/testdata_NRW_std/oasis3";fi
  if [[ $forcingdir_pfl == "" ]] ; then ;forcingdir_pfl="/work/slts/slts00/tsmp/TerrSysMPdb/testdata_NRW_std/ParFlow/Rur_NRW";fi


  if [[ $namelist_clm == "" ]] ; then ; namelist_clm=$rootdir/bldsva/setups/nrw/lnd.stdin ; fi
  if [[ $namelist_cos == "" ]] ; then ; namelist_cos=$rootdir/bldsva/setups/nrw/lmrun_uc ; fi
  if [[ $namelist_pfl == "" ]] ; then ; namelist_pfl=$rootdir/bldsva/setups/nrw/coup_oas.tcl ; fi


  defaultNppn=48
  defaultCLMProcX=2
  defaultCLMProcY=2
  defaultCOSProcX=4
  defaultCOSProcY=8
  defaultPFLProcX=3
  defaultPFLProcY=4


  gx_clm=300
  gy_clm=300
  dt_clm=900
  res="0300x0300"

  gx_cos=150
  gy_cos=150
  dt_cos=10
  nbndlines=4

  gx_pfl=300
  gy_pfl=300
  dt_pfl=0.25
  pflrunname="rurlaf"

  cplfreq1=900
  cplfreq2=900


  if [[ $namelist_oas == "" ]] ; then  
    if [[ $withPFL == "false" && $withCOS == "true" ]]; then
      namelist_oas=$rootdir/bldsva/data_oas3/namcouple_cos_clm
    fi	
    if [[ $withPFL == "true" && $withCOS == "false" ]]; then
      namelist_oas=$rootdir/bldsva/data_oas3/namcouple_pfl_clm
    fi
    if [[ $withPFL == "true" && $withCOS == "true" ]]; then
      namelist_oas=$rootdir/bldsva/data_oas3/namcouple_cos_clm_pfl
    fi
  fi

  restDir="/work/slts/slts15/tsmp/TSMPForecastNRW$(date '+%Y-%m-%d-%H' -d "$restDate")/run"
  fn_finidat="$restDir/clmoas.clm2.r.$(date '+%Y-%m-%d' -d "$startDate")-00000.nc"
  pfbfilename="$restDir/rurlaf.out.press.00024.pfb"

}

finalizeSetup(){
route "${cblue}>> finalizeSetup${cnormal}"
  comment "   copy clmgrid into rundir"
    cp $forcingdir_clm/clm3.5/Rur_NRW/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
  check  
  if [[ $withOAS == "true" ]] then
    comment "   copy oasis remappingfiles into rundir"
      cp $forcingdir_oas/* $rundir >> $log_file 2>> $err_file
    check
    if [[ $withOASMCT == "true" ]] then
      for x in $rundir/*BILINEA* ;do 
        comment "   rename oasis3 remapping files" 
          mv $x $(echo $x | sed "s/BILINEA/BILINEAR/") >> $log_file 2>> $err_file
        check 
      done
    fi  
  fi  

  if [[ $withPFL == "true" ]] then
        comment "   cd to rundir"
          cd $rundir >> $log_file 2>> $err_file
        check

        comment "   copy slopes and slope script into rundir"
          cp $forcingdir_pfl/slopes/ascii2pfb.tcl.template_new $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
          cp $forcingdir_pfl/slopes/*slope.pfb* $rundir >> $log_file 2>> $err_file
	check
          chmod u+w $rundir/*slope*  $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
	comment "   sed procs into slopescript"
          sed "s,__svaroot__.*,$pfldir/bin," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
	comment "   create sloap pfb with tclsh"
          tclsh ./ascii2pfb.tcl >> $log_file 2>> $err_file
	check		
	
	comment "   copy soilind and soilind script into rundir"
          cp $forcingdir_pfl/soilInd/ascii2pfb.tcl.template_new $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
          cp $forcingdir_pfl/soilInd/*Soil* $rundir >> $log_file 2>> $err_file
	check
          chmod u+w $rundir/*Soil* $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
	comment "   sed procs into soilindscript"
          sed "s,__svaroot__.*,$pfldir/bin," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
	comment "   create soilInd pfb with tclsh"
        tclsh ./ascii2pfb.tcl >> $log_file 2>> $err_file
	check
  fi 
route "${cblue}<< finalizeSetup${cnormal}"
}
