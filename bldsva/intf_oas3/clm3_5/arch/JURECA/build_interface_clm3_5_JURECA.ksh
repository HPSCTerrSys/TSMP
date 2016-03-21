#! /bin/ksh

always_clm(){
rout "${cblue}>> always_clm${cnormal}"
rout "${cblue}<< always_clm${cnormal}"
}

configure_clm(){
rout "${cblue}>> configure_clm${cnormal}"
  cplLib="-lnetcdff "
  flags=""
  flags+="-cc $mpiPath/bin/mpicc "
  flags+="-fc $mpiPath/bin/mpif90 "
  c_configure_clm
rout "${cblue}<< configure_clm${cnormal}"
}

make_clm(){
rout "${cblue}>> make_clm${cnormal}"
  c_make_clm
rout "${cblue}<< make_clm${cnormal}"
}


substitutions_clm(){
rout "${cblue}>> substitutions_clm${cnormal}"
   c_substitutions_clm
  comment "   cp m_FileResolve.F90 and shr_sys_mod.F90 to usr.src folder"
    cp $rootdir/bldsva/intf_oas3/clm3_5/arch/JURECA/src/m_FileResolv.F90 $clmdir/bld/usr.src	
  check
    cp $rootdir/bldsva/intf_oas3/clm3_5/arch/JURECA/src/shr_sys_mod.F90 $clmdir/bld/usr.src
  check
  if [[ $withOASMCT == "true" ]] ; then
    comment "   replace files for oasis3-mct and parallel clm coupling"
        cp $rootdir/bldsva/intf_oas3/clm3_5/arch/JURECA/src/mct/atmdrvMod.F90 $clmdir/bld/usr.src/ >> $log_file 2>> $err_file
    check
        cp $rootdir/bldsva/intf_oas3/clm3_5/arch/JURECA/src/mct/decompMod.F90 $clmdir/bld/usr.src/ >> $log_file 2>> $err_file
    check
        cp $rootdir/bldsva/intf_oas3/clm3_5/arch/JURECA/src/mct/oas* $clmdir/src/oas3/ >> $log_file 2>> $err_file
    check
        cp $rootdir/bldsva/intf_oas3/clm3_5/arch/JURECA/src/mct/receive* $clmdir/src/oas3/ >> $log_file 2>> $err_file
    check
        cp $rootdir/bldsva/intf_oas3/clm3_5/arch/JURECA/src/mct/send* $clmdir/src/oas3/ >> $log_file 2>> $err_file
    check
  fi

  comment "   cp new clm configure & Makefile.in to clm/bld/"
    cp $rootdir/bldsva/intf_oas3/clm3_5/arch/JURECA/config/configure $clmdir/bld >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/clm3_5/arch/JURECA/config/Makefile.in $clmdir/bld >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/clm3_5/arch/JURECA/config/config_clm_defaults.xml $clmdir/bld >> $log_file 2>> $err_file
  check
rout "${cblue}<< substitutions_clm${cnormal}"
}


setup_clm(){
rout "${cblue}>> setupClm${cnormal}"
  seconds_clm=$(($hh*3600))
  runstep_clm=$(($runhours*3600/$dt_clm))
  rpointer=$rundir/lnd.clmoas.rpointer

comment "  cp namelist to rundir"
  cp $namelist_clm $rundir/lnd.stdin >> $log_file 2>> $err_file
check



comment "  sed starttime to namelist"
  sed "s,__seconds_clm_bldsva__,$seconds_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed dt to namelist"
  sed "s,__dt_clm_bldsva__,$dt_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed forcingdir to namelist"
  sed "s,__forcingdir__,$forcingdir_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed gridsize to namelist"
  sed "s,__gridsize__,$res," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  create rpointer dummy file"
  touch $rpointer >> $log_file 2>> $err_file
check
comment "  sed rpointer path to namelist"
  sed "s,__rundir_rpointerdir__,$rpointer," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed date to namelist"
  sed "s,__yyyymmdd_bldsva__,${yyyy}${mm}${dd}," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed runtime to namelist"
  sed "s,__runstep_clm_bldsva__,$runstep_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
  if [[ $restart == 0 ]] then
comment "  sed no restart file path to namelist"
    sed "s,__finidat__,," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
  else
comment "  sed restart file path to namelist"
    sed "s,__finidat__,${fn_finidat}," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check

fi

    echo "${fn_finidat}" > $rpointer

rout "${cblue}<< setupClm${cnormal}"
}
