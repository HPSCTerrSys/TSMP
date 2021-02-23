#! /bin/ksh

initSetup(){
  defaultFDCLM="$rootdir/tsmp_nrw/input/clm"
  defaultFDCOS="$rootdir/tsmp_nrw/input/cosmo"
  defaultFDOAS="$rootdir/tsmp_nrw/input/oasis3"
  defaultFDPFL="$rootdir/tsmp_nrw/input/parflow"

  defaultNLCLM=$rootdir/bldsva/setups/nrw/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/nrw/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/nrw/coup_oas.tcl
  defaultNLDA=$rootdir/bldsva/setups/nrw/DA-nl

  defaultNppn=48
  defaultCLMProcX=2
  defaultCLMProcY=2
  defaultCOSProcX=4
  defaultCOSProcY=8
  defaultPFLProcX=3
  defaultPFLProcY=4
 
  defaultStartDate="2008-05-08 00"
  defaultInitDate="2008-05-08 00"
  defaultRunhours=24

  defaultDumpCLM=1
  defaultDumpCOS=1
  defaultDumpPFL=1	

  gx_clm=300
  gy_clm=300
  dt_clm=90
  res="0300x0300"

  gx_cos=150
  gy_cos=150
  dt_cos=10
  nbndlines=3

  gx_pfl=300
  gy_pfl=300
  dt_pfl=0.025
  pflrunname="rurlaf"
  base_pfl=0.0025

  cplfreq1=90
  cplfreq2=90

  delta_obs=1
 
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

        comment "   copy slopes and ascii2pfb into rundir"
          cp $forcingdir_pfl/ascii2pfb_slopes.tcl $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
          cp $forcingdir_pfl/*slope.pfb* $rundir >> $log_file 2>> $err_file
	check
          chmod u+w $rundir/*slope*  $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
        check
	comment "   sed procs into slopescript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
	comment "   create sloap pfb with tclsh"
          tclsh ./ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check		
	
	comment "   copy soilind, ics and  asc2pfb into rundir"
          cp $forcingdir_pfl/ascii2pfb_SoilInd.tcl $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	check
          cp $forcingdir_pfl/*Soil* $rundir >> $log_file 2>> $err_file
        check
          cp $forcingdir_pfl/*.ics.pfb* $rundir >> $log_file 2>> $err_file
	check
          chmod u+w $rundir/*Soil* $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
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
