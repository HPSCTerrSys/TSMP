#! /bin/ksh

initSetup(){
  defaultFDCLM="/homea/slts/slts06/forcings/testdata_EU_std"
  defaultFDCOS="/work/slts/slts15/tsmp/TSMPForecastEU2016-04-10-12/forcing/cosmoinput"
  defaultFDOAS="/homea/slts/slts06/forcings/testdata_EU_std/oasis3"
  defaultFDPFL="/homea/slts/slts06/forcings/testdata_EU_std/ParFlow"


  defaultNLCLM=$rootdir/bldsva/setups/cordex/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/cordex/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/cordex/coup_oas.tcl 


  defaultNppn=48
  defaultCLMProcX=3
  defaultCLMProcY=8
  defaultCOSProcX=12
  defaultCOSProcY=16
  defaultPFLProcX=9
  defaultPFLProcY=8

  defaultStartDate="2016-04-10 12"
  defaultRestDate=""
  defaultRunhours=3

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


  if [[ $withPFL == "false" && $withCOS == "true" ]]; then
    if [[ $cplscheme == "false" ]]; then
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm_a1
    else
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm
    fi
  fi
  if [[ $withPFL == "true" && $withCOS == "false" ]]; then
    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_pfl_clm
  fi
  if [[ $withPFL == "true" && $withCOS == "true" ]]; then
    if [[ $cplscheme == "false" ]]; then
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm_pfl_a1
    else
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm_pfl
    fi
  fi

 

  restDir="/work/slts/slts15/tsmp/TSMPForecastEU$(date '+%Y-%m-%d-%H' -d "$restDate")/run"
  fn_finidat="$restDir/clmoas.clm2.r.$(date '+%Y-%m-%d' -d "$startDate")-43200.nc"
  pfbfilename="$restDir/${pflrunname}.out.press.00024.pfb"


}

finalizeSetup(){
route "${cblue}>> finalizeSetup${cnormal}"
  if [[ $withOAS == "true" ]] then
    comment "   copy clmgrid into rundir"
      cp $forcingdir_clm/clm3.5/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
    check

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
