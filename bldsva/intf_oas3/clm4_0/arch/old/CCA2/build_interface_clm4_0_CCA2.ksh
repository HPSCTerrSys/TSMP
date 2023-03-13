#! /bin/ksh

always_clm(){
route "${cyellow}>> always_clm${cnormal}"
route "${cyellow}<< always_clm${cnormal}"
}

configure_clm(){
route "${cyellow}>> configure_clm${cnormal}"
  cplFlag=""
  cplLib=""
  cplInc=""
  USEOASLIB="FALSE"
  if [[ $withOAS == "true" ]]; then
       cplLib+="$liboas $libpsmile"
       cplInc=$incpsmile
       cplFlag="-DCOUP_OAS_COS"
       USEOASLIB="TRUE"
  fi
  LMASK="navy"
  if [[ $cplscheme == "true" ]] ; then ; cplFlag+=" -DCPL_SCHEME_F " ; fi

  comment "   cp fresh makefile to $clmdir/scripts"
    cp -R $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/scripts/ccsm_utils/Machines/Makefile $clmdir/scripts/ccsm_utils/Machines/ >> $log_file 2>> $err_file
  check
    cp -R $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/Makefile $clmdir/models/utils/oas3/ >> $log_file 2>> $err_file
  check
    cp -R $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/scripts/ccsm_utils/Machines/config_compilers.xml $clmdir/scripts/ccsm_utils/Machines/ >> $log_file 2>> $err_file
  check
  file=$clmdir/scripts/ccsm_utils/Machines/Makefile
  file2=$clmdir/models/utils/oas3/Makefile
  file3=$clmdir/scripts/ccsm_utils/Machines/config_compilers.xml
  comment "   sed comflg to clm Makefiles"
    sed -i "s@__comflg__@$optComp $cplInc -Duse_libMPI -Duse_netCDF -Duse_comm_MPI1 -DVERBOSE -DTREAT_OVERLAY $cplFlag@" $file $file2 >> $log_file 2>> $err_file
  check
  comment "   sed ldflg to clm Makefiles"
    sed -i "s@__ldflg__@$cplLib@" $file $file2 >> $log_file 2>> $err_file
  check

    comment "   sed -loas3 to clm Makefile if withOAS"
      if [[ $withOAS == "true" ]]; then
        sed -i 's@__oaslib__@-L$(OAS_LIBDIR) -loas3@' $file >> $log_file 2>> $err_file
      else		
        sed -i 's@__oaslib__@@' $file >> $log_file 2>> $err_file
      fi		
    check


  export NETCDF_PATH=$ncdfPath
  comment "   sed ldflg to machine definition"
    sed -i "s@__slibs__@ -L$lapackPath -llapack -lblas -L$liboas -L$libpsmile -L$ncdfPath/lib/ -lnetcdf -lnetcdff @" $file3 >> $log_file 2>> $err_file
  check
  comment "   sed pnetcdf to machine definition"
    sed -i "s@__pnetcdf__@ $pncdfPath @" $file3 >> $log_file 2>> $err_file
  check  
  comment "   sed netcdf to machine definition"
    sed -i "s@__netcdf__@ $ncdfPath @" $file3 >> $log_file 2>> $err_file
  check
  comment "   sed mpi to machine definition"
    sed -i "s@__mpi__@ $mpiPath @" $file3 >> $log_file 2>> $err_file
  check


  comment "   rm $clmdir/clmoas dir to make clean"
    rm -rf $clmdir/clmoas >> $log_file 2>> $err_file
  check
  comment "   cd to $clmdir/scripts/"
    cd $clmdir/scripts/ >> $log_file 2>> $err_file
  check
  comment "   create new cesm case" 
    ./create_newcase -v -case $clmdir/clmoas -res CLM_USRDAT -compset I -mach cluma2  >> $log_file 2>> $err_file
  check
  comment "   cd to $clmdir/clmoas"
    cd $clmdir/clmoas >> $log_file 2>> $err_file
  check
  comment "   xmlchange USE_OAS_LIB"
    ./xmlchange USE_OAS_LIB=$USEOASLIB >> $log_file 2>> $err_file
  check
  comment "   xmlchange MASK_GRID"
    ./xmlchange MASK_GRID=$LMASK >> $log_file 2>> $err_file
  check
  comment "   cesm_setup to create macros"
    ./cesm_setup >> $log_file 2>> $err_file
  check

  comment "   mkdir $clmdir/clmoas/bld tree"
    mkdir -p $clmdir/clmoas/bld/lib/include/ >> $log_file 2>> $err_file
  check  
    for model in cpl atm lnd ice ocn glc wav rof
    do
      mkdir -p "$clmdir/clmoas/bld/$model/obj" >> $log_file 2>> $err_file
    check
    done
  comment "   configure clm"
    ./Buildconf/clm.buildnml.csh  >> $log_file 2>> $err_file
  check

route "${cyellow}<< configure_clm${cnormal}"
}

make_clm(){
route "${cyellow}>> make_clm${cnormal}"
  comment "   cd to $clmdir/clmoas"
    cd $clmdir/clmoas >> $log_file 2>> $err_file
  check
  comment "   make clm"
    ./clmoas.build >> $log_file 2>> $err_file
  check
  comment "   copy clm exe to $bindir"
    cp $clmdir/clmoas/bld/cesm.exe $bindir/clm >> $log_file 2>> $err_file
  check
route "${cyellow}<< make_clm${cnormal}"
}


substitutions_clm(){
route "${cyellow}>> substitutions_clm${cnormal}"

  comment "   cp new clm configure, Makefile XMLs and scripts to $clmdir"
    patch "$rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/scripts/ccsm_utils/Machines/*" $clmdir/scripts/ccsm_utils/Machines/
  check
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/scripts/ccsm_utils/Case.template/config_definition.xml $clmdir/scripts/ccsm_utils/Case.template
  check
    patch "$rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/scripts/ccsm_utils/Tools/*" $clmdir/scripts/ccsm_utils/Tools/ 
  check
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/scripts/create_newcase $clmdir/scripts/
  check
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/models/lnd/clm/bld/clm.buildnml.csh $clmdir/models/lnd/clm/bld/ 
  check
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/models/lnd/clm/bld/configure $clmdir/models/lnd/clm/bld/ 
  check
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/models/lnd/clm/bld/config_files/config_definition.xml $clmdir/models/lnd/clm/bld/config_files/ 
  check
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/models/atm/datm/bld/datm.buildexe.csh $clmdir/models/atm/datm/bld/ 
  check
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/models/atm/datm/bld/namelist_defaults_datm.xml $clmdir/models/atm/datm/bld/namelist_files/ 
  check
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/models/drv/bld/cesm.buildexe.csh $clmdir/models/drv/bld/ 
  check


  comment "   cp oasis interface to $clmdir"
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/oas3/utils/oas3 $clmdir/models/utils/
  check
    patch  $rootdir/bldsva/intf_oas3/${mList[1]}/oas3/clm/cpl_oas3 $clmdir/models/lnd/clm/src/ 
  check
    patch  $rootdir/bldsva/intf_oas3/${mList[1]}/oas3/atm/cpl_oas3 $clmdir/models/atm/datm/ 
  check



  if [[ $withOASMCT == "true" ]] ; then
     comment "   sed replace old mod_prism includes from clm oas files"
       sed -i "s/mod_prism_proto/mod_prism/" $clmdir/models/utils/oas3/oas_clm_vardef.F90 >> $log_file 2>> $err_file
     check
       sed -i "s/USE mod_prism.*//" $clmdir/models/utils/oas3/oas_clm_rcv.F90 >> $log_file 2>> $err_file
     check
       sed -i "s/USE mod_prism.*//" $clmdir/models/utils/oas3/oas_clm_snd.F90 >> $log_file 2>> $err_file
     check
       sed -i "s/USE mod_prism.*//" $clmdir/models/utils/oas3/oas_clm_define.F90 >> $log_file 2>> $err_file
     check
  fi




  comment "   mkdir folder with source modifications:  $clmdir/SourceMods/"
    mkdir -p $clmdir/SourceMods/ >> $log_file 2>> $err_file
  check
  comment "   cp source modifications to $clmdir"
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/src.drv $clmdir/SourceMods/ 
  check
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/src.datm $clmdir/SourceMods/ 
  check
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/src.clm $clmdir/SourceMods/ 
  check
    #FG: I'm not sure which consequences this can have. I replaced ISO_C_BINDING:C_SIZE_OF by something intrinsic.
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/src/iompi_mod.F90 $clmdir/models/utils/pio/ 
  check

route "${cyellow}<< substitutions_clm${cnormal}"
}


setup_clm(){
route "${cyellow}>> setupClm${cnormal}"

  withCESM="true"
  seconds_clm=$(($hh*3600))
  runstep_clm=$runhours
  rpointer=$rundir/lnd.clmoas.rpointer

comment "  cp namelist to rundir"
  cp ${namelist_clm} $rundir/lnd.stdin >> $log_file 2>> $err_file
check

comment "  sed num procs to namelist"
  sed "s,__nprocs__,$(($px_clm * $py_clm))," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check

  c_setup_clm  

route "${cyellow}<< setupClm${cnormal}"
}
