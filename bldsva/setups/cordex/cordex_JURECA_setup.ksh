#! /bin/ksh

initSetup(){
  defaultFDCLM="$rootdir/tsmp_eur11_eraint_eval/input/clm"
  defaultFDCOS="$rootdir/tsmp_eur11_eraint_eval/input/cosmo"
  defaultFDOAS="$rootdir/tsmp_eur11_eraint_eval/input/oasis3"
  defaultFDPFL="$rootdir/tsmp_eur11_eraint_eval/input/parflow"


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

  defaultStartDate="2016-05-01 12"
  defaultInitDate="2016-05-01 12"
  
  defaultDumpCLM=1
  defaultDumpCOS=1
  defaultDumpPFL=1
  
  defaultRunhours=3

  defaultDumpCLM=1
  defaultDumpCOS=1
  defaultDumpPFL=1 

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
  base_pfl=0.0025

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

}

finalizeSetup(){
route "${cblue}>> finalizeSetup${cnormal}"
  if [[ $withOAS == "true" ]] then
    comment "   copy clmgrid into rundir"
      cp $forcingdir_clm/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
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
          cp $forcingdir_pfl/ascii2pfb_slopes.tcl $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
          cp $forcingdir_pfl/slope*.sa $rundir >> $log_file 2>> $err_file
	check
          chmod u+w $rundir/slope*.sa  $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
        check
	comment "   sed procs into slopescript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
        if [ ! -f "$rundir/coup_oas.tcl" ]; then
          comment "  parflow namelist coup_oas.tcl is copied to rundir. See cordex_JURECA_setup.ksh"
          cp $namelist_pfl $rundir/coup_oas.tcl >> $log_file 2>> $err_file
          check
        fi
          sed "s,__pfl_solidinput_filename__,$defaultFDPFL/geom_cordex0.11_436x424.pfsol," -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
        comment "   create parflow db with tclsh (coup_oas.tcl)"
          tclsh $rundir/coup_oas.tcl >> $log_file 2>> $err_file
        check
	comment "   create sloap pfb with tclsh"
          tclsh ./ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check		
	
	comment "   copy soilind and soilind script into rundir"
          cp $forcingdir_pfl/ascii2pfb_SoilInd.tcl $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	check
          cp $forcingdir_pfl/parflow_436x424x15_cosmomask_indicator_FAOonly.sa $rundir/pfl_ind.sa >> $log_file 2>> $err_file
	check
          chmod u+w $rundir/pfl_ind.sa $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
        check
	comment "   sed procs into soilindscript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	check
	comment "   create soilInd pfb with tclsh"
        tclsh ./ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	check
  fi 
route "${cblue}<< finalizeSetup${cnormal}"
}
