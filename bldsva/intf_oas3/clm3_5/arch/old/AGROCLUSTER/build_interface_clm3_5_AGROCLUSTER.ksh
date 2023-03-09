#! /bin/ksh

always_clm(){
route "${cyellow}>> always_clm${cnormal}"
route "${cyellow}<< always_clm${cnormal}"
}

configure_clm(){
route "${cyellow}>> configure_clm${cnormal}"
  cplLib=" -lnetcdf -lnetcdff "
  flags=""
  ccc="$mpiPath/bin/mpicc "
  cfc="mpif90 "
  flags+="-mpi_lib $mpiPath/lib "
  c_configure_clm
route "${cyellow}<< configure_clm${cnormal}"
}

make_clm(){
route "${cyellow}>> make_clm${cnormal}"
  c_make_clm
route "${cyellow}<< make_clm${cnormal}"
}


substitutions_clm(){
route "${cyellow}>> substitutions_clm${cnormal}"
   c_substitutions_clm
  comment "   cp m_FileResolve.F90 and shr_sys_mod.F90 to usr.src folder"
    patch $rootdir/bldsva/intf_oas3/clm3_5/arch/$platform/src/m_FileResolv.F90 $clmdir/bld/usr.src	
  check
    patch $rootdir/bldsva/intf_oas3/clm3_5/arch/$platform/src/shr_sys_mod.F90 $clmdir/bld/usr.src
  check
  if [[ $withOASMCT == "true" ]] ; then
    comment "   replace files for oasis3-mct and parallel clm coupling"
        patch $rootdir/bldsva/intf_oas3/clm3_5/mct/atmdrvMod.F90 $clmdir/bld/usr.src/
    check
        patch $rootdir/bldsva/intf_oas3/clm3_5/mct/decompMod.F90 $clmdir/bld/usr.src/ 
    check
        patch "$rootdir/bldsva/intf_oas3/clm3_5/mct/oas*" $clmdir/src/oas3/ 
    check
        patch "$rootdir/bldsva/intf_oas3/clm3_5/mct/receive*" $clmdir/src/oas3/ 
    check
        patch "$rootdir/bldsva/intf_oas3/clm3_5/mct/send*" $clmdir/src/oas3/ 
    check
  fi

  comment "   cp new clm configure & Makefile.in to clm/bld/"
    patch $rootdir/bldsva/intf_oas3/clm3_5/arch/$platform/config/configure $clmdir/bld 
  check
    patch $rootdir/bldsva/intf_oas3/clm3_5/arch/$platform/config/Makefile.in $clmdir/bld 
  check
    patch $rootdir/bldsva/intf_oas3/clm3_5/arch/$platform/config/config_clm_defaults.xml $clmdir/bld 
  check
route "${cyellow}<< substitutions_clm${cnormal}"
}


setup_clm(){
route "${cyellow}>> setupClm${cnormal}"
  seconds_clm=$(($hh*3600))
  runstep_clm=$(($runhours*3600/$dt_clm))
  rpointer=$rundir/lnd.clmoas.rpointer

comment "  cp namelist to rundir"
  cp $namelist_clm $rundir/lnd.stdin >> $log_file 2>> $err_file
check

  c_setup_clm

route "${cyellow}<< setupClm${cnormal}"
}
