#! /bin/ksh

initSetup(){
  defaultFDCLM="$static_files/clm"
  defaultFDCOS="$static_files/cosmo"
  defaultFDOAS="$static_files/oasis3"
  defaultFDPFL="$static_files/parflow"
  defaultFDICON="$static_files/icon"

  defaultNLCLM=$rootdir/bldsva/setups/$refSetup/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/$refSetup/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/$refSetup/coup_oas.tcl
  defaultNLDA=$rootdir/bldsva/setups/$refSetup/DA-nl




  defaultNppn=$Npp
  defaultNgpn=$Ngp

  if [[ $processor == "GPU" || $processor == "MSA" ]]; then
    defaultPFLProcX=$PFLProcXg
    defaultPFLProcY=$PFLProcYg
    defaultCLMProcX=$CLMProcXg
	defaultCLMProcY=$CLMProcYg
	defaultCOSProcX=$COSProcXg
	defaultCOSProcY=$COSProcYg
  else
    defaultPFLProcX=$PFLProcX
    defaultPFLProcY=$PFLProcY
    defaultCLMProcX=$CLMProcX
	defaultCLMProcY=$CLMProcY
	defaultCOSProcX=$COSProcX
	defaultCOSProcY=$COSProcY
  fi
  defaultStartDate=$StartDate
  defaultInitDate=$InitDate
  
  defaultDumpCLM=$DumpCLM
  defaultDumpCOS=$DumpCOS
  defaultDumpPFL=$DumpPFL
  
  defaultRunhours=$Runhours

  gx_clm=$gxCLM
  gy_clm=$gyCLM
  dt_clm=$dtCLM
  res=$resCLM

  gx_cos=$gxCOS
  gy_cos=$gyCOS
  dt_cos=$dtCOS
  nbndlines=$nboundlinesCOS

  gx_pfl=$gxPFL
  gy_pfl=$gyPFL
  dt_pfl=$dtPFL
  pflrunname=$runnamePFL
  base_pfl=$basePFL

  cplfreq1=$freq1OAS
  cplfreq2=$freq2OAS
  
  delta_obs=$deltaobs


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
  if [[ $withPFL == "false" && $withICON == "true" ]]; then
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_icon_clm
  fi
  if [[ $withPFL == "true" && $withICON == "true" ]]; then
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_icon_clm_pfl
  fi
}


finalizeSetup(){
route "${cyellow}>> finalizeSetup${cnormal}"
if [[ $withICON == "true" ]]; then
  comment "  copy icon grid"
    ln -s $gicon $rundir/ >> $log_file 2>> $err_file
  check
  comment "   copy clmgrid into rundir"
      cp $forcingdir_clm/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
    check_pfl

  comment "  copy remap-files"
    cp $rmp_folder/rmp* $rundir/ >> $log_file 2>> $err_file
    rm -r $rmp_folder
  check
else
  if [[ $withOAS == "true" ]] then
    comment "   copy clmgrid into rundir"
      cp $forcingdir_clm/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
    check_pfl

    comment "   copy oasis remappingfiles into rundir"
      cp $forcingdir_oas/* $rundir >> $log_file 2>> $err_file
    check
    if [[ $withOASMCT == "true" ]] then
      if [ -f "$rundir"/*BILINEAR* ]; then
        comment "  no renaming of oasis3 remapping files needed" 
      else
        for x in $rundir/*BILINEA* ;do 
          comment "   rename oasis3 remapping files" 
            mv $x $(echo $x | sed "s/BILINEA/BILINEAR/") >> $log_file 2>> $err_file
          check 
        done
      fi
    fi  
  fi  

  if [[ $withCOS == "true" ]] then
    comment "  sed gribapi definitions and samples to namelist"
     p_samp=$gribPath/share/eccodes/samples/
     p_def=$gribPath/share/eccodes/definitions/
     sed "s,__definitions__,$p_def," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
     sed "s,__samples__,$p_samp," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
    check
  fi

  if [[ $withPFL == "true" ]] then
        comment "   cd to rundir"
          cd $rundir >> $log_file 2>> $err_file
        check
		if [ -f $forcingdir_pfl/ascii2pfb_slopes.tcl ]; then
        comment "   copy slopes and slope script into rundir"
          cp $forcingdir_pfl/ascii2*.tcl $rundir/ >> $log_file 2>> $err_file
		check_pfl
          cp $forcingdir_pfl/$slope $rundir >> $log_file 2>> $err_file
		check_pfl
          chmod u+w $rundir/$slope  $rundir/ascii2*.tcl >> $log_file 2>> $err_file
		comment "   sed procs into slopescript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	
          sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	
          sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
          tclsh ./ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
		check_pfl
		fi
	if [ -f $forcingdir_pfl/$inipress ]; then
    comment "   copy initial pressure and script into rundir"
		  cp $forcingdir_pfl/$inipress $rundir/ >> $log_file 2>> $err_file
		  if [ -f $forcingdir_pfl/ascii2pfb.tcl ]; then
		  comment "   distribute initial pressure"
		  cp $forcingdir_pfl/ascii2pfb.tcl $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
		  check_pfl
		  chmod u+w $rundir/rur_ic_press.pfb  $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
		  check_pfl
		  sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
		  check_pfl
	      sed "s,pfset Process\.Topology\.P.*,pfset Process\.Topology\.P $px_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
	      check_pfl
          sed "s,pfset Process\.Topology\.Q.*,pfset Process\.Topology\.Q $py_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
          check_pfl
          tclsh ./ascii2pfb.tcl >> $log_file 2>> $err_file
          fi
	check_pfl
	fi
	if [ -f $forcingdir_pfl/ascii2pfb_SoilInd.tcl ]; then
	comment "   copy soilind and soilind script into rundir"
          cp $forcingdir_pfl/ascii2pfb_SoilInd.tcl $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	
          cp $forcingdir_pfl/$indPFL $rundir/$indPFL2 >> $log_file 2>> $err_file
    check_pfl	
          chmod u+w $rundir/$indPFL2  $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
        
	comment "   sed procs into soilindscript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	
          sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	
          sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file
	check_pfl
	comment "   create soilInd pfb with tclsh"
        tclsh ./ascii2pfb_SoilInd.tcl >> $log_file 2>> $err_file		
          sed "s,__pfl_solidinput_filename__,$defaultFDPFL/$pfsolPFL," -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
	fi
	
        tclsh ./coup_oas.tcl >> $log_file 2>> $err_file
    check_pfl
  fi 
fi


route "${cyellow}<< finalizeSetup${cnormal}"
}
