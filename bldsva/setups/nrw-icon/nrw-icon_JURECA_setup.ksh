#! /bin/ksh

initSetup(){
  defaultFDCLM="/work/slts/slts00/tsmp/TestCases/nrw/clm"
  defaultFDOAS="/work/slts/slts00/tsmp/TestCases/nrw/oasis3"
  #defaultFDPFL="/work/slts/slts00/tsmp/TestCases/nrw/parflow"

  defaultNLICON="$rootdir/bldsva/setups/nrw-icon/icon_master.namelist $rootdir/bldsva/setups/nrw-icon/NAMELIST_icon"
  defaultNLCLM=$rootdir/bldsva/setups/nrw/lnd.stdin 
  #defaultNLPFL=$rootdir/bldsva/setups/nrw/coup_oas.tcl

  defaultNppn=24
  defaultICONProc=20
  defaultCLMProcX=2
  defaultCLMProcY=2
  #defaultPFLProcX=3
  #defaultPFLProcY=4
  
  defaultStartDate="2013-04-24 00"
  defaultInitDate="2013-04-24 00"
  defaultRunhours=3

  defaultDumpCLM=1
  #defaultDumpPFL=1	

  gx_icon=15472

  gx_clm=300
  gy_clm=300
  dt_clm=900
  res="0300x0300"

  #gx_pfl=300
  #gy_pfl=300
  #dt_pfl=0.25
  #pflrunname="rurlaf"
  #base_pfl=0.0025

  cplfreq1=9000
  #cplfreq2=900

  delta_obs=1
 
  if [[ $withPFL == "false" && $withICON == "true" ]]; then
    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_icon_clm
  fi
  if [[ $withPFL == "true" && $withICON == "false" ]]; then
    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_pfl_clm
  fi
  if [[ $withPFL == "true" && $withICON == "true" ]]; then
    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_icon_clm_pfl
  fi


}

finalizeSetup(){
route "${cyellow}>> finalizeSetup${cnormal}"
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

  comment "  copy grid file for icon"
    ln -s /work/jibg34/jibg3401/input_nrw/input_nrw_150km_R1247/icon_nrw_150km_R1247m.nc $rundir/ >> $log_file 2>> $err_file
  check

  comment "  copy initial data for icon"
    ln -s /work/jibg34/jibg3401/input_nrw/input_nrw_150km_R1247/init_icon_nrw_150km_R1247m_2013042400.nc $rundir/ >> $log_file 2>> $err_file
    ln -s /work/jibg34/jibg3401/input_nrw/input_nrw_150km_R1247/init_icon_nrw_150km_R1247m_2013042400.nc $rundir/ifs2icon_R2B11_DOM01.nc >> $log_file 2>> $err_file
    ln -s /work/jibg34/jibg3401/input_nrw/input_nrw_150km_R1247/init_icon_nrw_150km_R1247m_2013042400.nc $rundir/dwdFG_R2B11_DOM01.nc >> $log_file 2>> $err_file
  check

  comment "  copy extpart for icon"
    ln -s /work/jibg34/jibg3401/input_nrw/input_nrw_150km_R1247/extpar_icon_nrw_150km_R1247m.nc $rundir/ >> $log_file 2>> $err_file
  check

  comment "  copy remap data for icon"
    ln -s /work/jibg34/jibg3401/input_nrw/input_nrw_150km_R1247/rmp_* $rundir/ >> $log_file 2>> $err_file
  check

  comment "  copy 2mom lookup table"
    ln -s $icondir/data/dmin_wetgrowth_lookup.dat $rundir >> $log_file 2>> $err_file
  check

  if [[ $withPFL == "true" ]] then
        comment "   cd to rundir"
          cd $rundir >> $log_file 2>> $err_file
        check

        comment "   copy slopes and slope script into rundir"
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
	
	comment "   copy soilind and soilind script into rundir"
          cp $forcingdir_pfl/ascii2pfb_SoilInd.tcl $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	check
          cp $forcingdir_pfl/*Soil* $rundir >> $log_file 2>> $err_file
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
route "${cyellow}<< finalizeSetup${cnormal}"
}
