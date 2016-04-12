#! /bin/ksh

initSetup(){
  if [[ $forcingdir_clm == "" ]] ; then ;forcingdir_clm="/homea/slts/slts06/forcings/testdata_EU_std";fi
  if [[ $forcingdir_cos == "" ]] ; then ;forcingdir_cos="";fi
  if [[ $forcingdir_oas == "" ]] ; then ;forcingdir_oas="/homea/slts/slts06/forcings/testdata_EU_std/oasis3";fi
  if [[ $forcingdir_pfl == "" ]] ; then ;forcingdir_pfl="/homea/slts/slts06/forcings/testdata_EU_std/ParFlow";fi


  if [[ $namelist_clm == "" ]] ; then ; namelist_clm=$rootdir/bldsva/setups/cordex/lnd.stdin ; fi
  if [[ $namelist_cos == "" ]] ; then ; namelist_cos=$rootdir/bldsva/setups/cordex/lmrun_uc ; fi
  if [[ $namelist_pfl == "" ]] ; then ; namelist_pfl=$rootdir/bldsva/setups/cordex/coup_oas.tcl ; fi


  defaultNppn=16
  defaultCLMProcX=8
  defaultCLMProcY=8
  defaultCOSProcX=24
  defaultCOSProcY=16
  defaultPFLProcX=8
  defaultPFLProcY=8


  gx_clm=436
  gy_clm=424
  dt_clm=3600
  res="436x424"

  gx_cos=444
  gy_cos=432
  dt_cos=60
  nbndlines=4

  gx_pfl=436
  gy_pfl=424
  dt_pfl=1.0
  pflrunname="cordex0.11"

  cplfreq1=3600
  cplfreq2=3600


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

  restDir="/work/slts/slts15/tsmp/TSMPForecastEU$(date '+%Y-%m-%d-%H' -d "$restDat")/run"
  fn_finidat="$restDir/clmoas.clm2.r.$(date '+%Y-%m-%d' -d "$startDate")-43200.nc"
  pfbfilename="$restDir/${pflrunname}.out.press.00024.pfb"

}

finalizeSetup(){
route "${cblue}>> finalizeSetup${cnormal}"
  comment "   copy clmgrid into rundir"
    cp $forcingdir_clm/clm3.5/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
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
          cp $forcingdir_pfl/slopes/ascii2pfb_slopes.tcl.template $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
          cp $forcingdir_pfl/slopes/slope*.sa $rundir >> $log_file 2>> $err_file
	check
          chmod u+w $rundir/slope*.sa  $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
	comment "   sed procs into slopescript"
          sed "s,lappend auto_path.*,lappend auto_path $pfldir/bin," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
	comment "   create sloap pfb with tclsh"
          tclsh ./ascii2pfb.tcl >> $log_file 2>> $err_file
	check		
	
	comment "   copy soilind and soilind script into rundir"
          cp $forcingdir_pfl/soilInd/ascii2pfb_ind.tcl.template $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	check
          cp $forcingdir_pfl/soilInd/parflow_436x424x15_cosmomask_indicator_FAOonly.sa $rundir/pfl_ind.sa >> $log_file 2>> $err_file
	check
          chmod u+w $rundir/pfl_ind.sa $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
	comment "   sed procs into soilindscript"
          sed "s,lappend auto_path.*,lappend auto_path $pfldir/bin," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
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
