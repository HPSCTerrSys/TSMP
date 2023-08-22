#! /bin/ksh

always_clm(){
route "${cyellow}>> always_clm${cnormal}"
route "${cyellow}<< always_clm${cnormal}"
}

configure_clm(){
route "${cyellow}>> configure_clm${cnormal}"

  comment " source software and cesm input paths"
    source $rootdir/bldsva/intf_oas3/clm5_0/arch/config/softwarepaths.ksh
  check

#  cplFlag=""
#  cplLib=""
#  cplInc=""
#  USEOASLIB="FALSE"
  if [[ $withOAS == "true" ]]; then
       cplLib+="$libpsmile"
       cplInc=$incpsmile
       cplFlag="-DCOUP_OAS_COS"
#       USEOASLIB="TRUE"
  fi
#  LMASK="navy"
#  if [[ $cplscheme == "true" ]] ; then ; cplFlag+=" -DCPL_SCHEME_F " ; fi

#  comment "   cp fresh makefile to $clmdir/scripts"
#    cp -R $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/scripts/cime/scripts/Tools/Makefile $clmdir/cime/scripts/Tools/ >> $log_file 2>> $err_file
#  check
#    mkdir -p $clmdir/components/utils/oas3/
#    cp -R $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/Makefile $clmdir/components/utils/oas3/ >> $log_file 2>> $err_file
#  check
#    cp -R $rootdir/bldsva/intf_oas3/${mList[1]}/arch/$platform/config/scripts/ccsm_utils/Machines/config_compilers.xml $clmdir/scripts/ccsm_utils/Machines/ >> $log_file 2>> $err_file
#  check
#  file=$clmdir/cime/scripts/Tools/Makefile
#  file2=$clmdir/components/utils/oas3/Makefile
#  file3=$clmdir/cime/config/cesm/machines/config_compilers.xml
#  comment "   sed comflg to clm Makefiles"
#    sed -i "s@__comflg__@$optComp $cplInc -Duse_libMPI -Duse_netCDF -Duse_comm_MPI1 -DVERBOSE -DTREAT_OVERLAY $cplFlag@" $file $file2 >> $log_file 2>> $err_file
#  check
#  comment "   sed ldflg to clm Makefiles"
#    sed -i "s@__ldflg__@$cplLib@" $file $file2 >> $log_file 2>> $err_file
#  check

#    comment "   sed -loas3 to clm Makefile if withOAS"
#      if [[ $withOAS == "true" ]]; then
#        sed -i 's@__oaslib__@-L$(OAS_LIBDIR) -loas3@' $file >> $log_file 2>> $err_file
#      else		
#        sed -i 's@__oaslib__@@' $file >> $log_file 2>> $err_file
#      fi		
#    check


#  export NETCDF_PATH=$ncdfPath

#  comment "   sed ldflg to machine definition"
#    sed -i "s@__slibs__@ -L$lapackPath -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lblacs_intelmpi_lp64 -L$liboas -L$libpsmile -L$ncdfPath/lib/ -lnetcdff -lnetcdf @" $file3 >> $log_file 2>> $err_file
#  check
#  comment "   sed pnetcdf to machine definition"
#    sed -i "s@__pnetcdf__@ $pncdfPath @" $file3 >> $log_file 2>> $err_file
#  check  
#  comment "   sed netcdf to machine definition"
#    sed -i "s@__netcdf__@ $ncdfPath @" $file3 >> $log_file 2>> $err_file
#  check
#  comment "   sed mpi to machine definition"
#    sed -i "s@__mpi__@ $mpiPath @" $file3 >> $log_file 2>> $err_file
#  check

  comment "   rm $clmdir/clmoas dir to make clean"
    rm -rf $clmdir/clmoas >> $log_file 2>> $err_file
  check
  comment "   cd to $clmdir/cime/scripts/"
    cd $clmdir/cime/scripts/ >> $log_file 2>> $err_file
  check
  comment "   create new cesm case" 
    ./create_newcase --case $clmdir/clmoas --res CLM_USRDAT --compset I2000Clm50BgcCropGs --run-unsupported >> $log_file 2>> $err_file
  check

  comment "   cd to $clmdir/clmoas"
    cd $clmdir/clmoas >> $log_file 2>> $err_file
  check
#  comment "   xmlchange USE_OAS_LIB"
#    ./xmlchange USE_OAS_LIB=$USEOASLIB >> $log_file 2>> $err_file
#  check
#  comment "   xmlchange MASK_GRID"
#    ./xmlchange MASK_GRID=$LMASK >> $log_file 2>> $err_file
#  check


# PATH TO BLD
  comment "   xmlchange BLD Paths"
    ./xmlchange CIME_OUTPUT_ROOT=$clmdir >> $log_file 2>> $err_file
  check
    ./xmlchange EXEROOT=$clmdir/build >> $log_file 2>> $err_file
  check
    ./xmlchange INCROOT=$clmdir/build/lib/include >> $log_file 2>> $err_file
  check
   ./xmlchange LIBROOT=$clmdir/build/lib >> $log_file 2>> $err_file
  check
    ./xmlchange OBJROOT=$clmdir/build >> $log_file 2>> $err_file
  check
    ./xmlchange SHAREDLIBROOT=$clmdir/build >> $log_file 2>> $err_file
  check
    ./xmlchange RUNDIR=$clmdir/clmoas/run >> $log_file 2>> $err_file
  check
    ./xmlchange DOUT_S_ROOT=$clmdir >> $log_file 2>> $err_file
  check



  comment "   case.setup to create macros"
    ./case.setup >> $log_file 2>> $err_file
  check

#  comment "   mkdir $clmdir/clmoas/bld tree"
#    mkdir -p $clmdir/clmoas/bld/lib/include/ >> $log_file 2>> $err_file
#  check  
#    for model in cpl atm lnd ice ocn glc wav rof
#    do
#      mkdir -p "$clmdir/clmoas/bld/$model/obj" >> $log_file 2>> $err_file
#    check
#    done


# Need to add empty fsurdat otherwise 
# preview_namelist / case.build scripts will look for default values and not find em
  comment " add fsurdat to user_nl_clm"
    sed -i '$ a fsurdat=""' user_nl_clm >> $log_file 2>> $err_file
  check
# preview_namelists generates some default namelists.
  comment "   configure clm"
    ./preview_namelists >> $log_file 2>> $err_file
  check

route "${cyellow}<< configure_clm${cnormal}"
}

make_clm(){
route "${cyellow}>> make_clm${cnormal}"
  comment "   cd to $clmdir/clmoas"
    cd $clmdir/clmoas >> $log_file 2>> $err_file
  check
  comment "   make clm"
    ./case.build >> $log_file 2>> $err_file
  check

  comment "  create libcpl.a and move to lib folder"
    ar rcs libcpl.a $clmdir/build/cpl/obj/*.o 
    ranlib libcpl.a
    mv libcpl.a $clmdir/build/lib/
    echo "cpl library compilation suceeded"
  check
  comment "   copy clm exe to $bindir"
    cp $clmdir/build/cesm.exe $bindir/clm >> $log_file 2>> $err_file
  check
route "${cyellow}<< make_clm${cnormal}"
}


substitutions_clm(){
route "${cyellow}>> substitutions_clm${cnormal}"

  comment "   cp clm config_ XMLs to $clmdir"
    cp $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/cime/config/cesm/machines/*.xml $clmdir/cime/config/cesm/machines/
  check

  if [[ $withPDAF == "true" ]] ; then
    comment " substitutions with pdaf"
    patch $rootdir/bldsva/intf_DA/${mList[4]}/tsmp/${mList[1]}/cime_comp_mod.F90 $clmdir/cime/src/drivers/mct/main
    check
    patch $rootdir/bldsva/intf_DA/${mList[4]}/tsmp/${mList[1]}/seq_comm_mct.F90 $clmdir/cime/src/drivers/mct/shr
    check
  fi
  
  if [[ ${mList[2]} == "icon-lem" ]] ; then
    # start of ICON coupling by Slavko

    # Following lines commented, conflict with necessary cime_comp_mod.F90
    # replacements from PDAF (see directly above). 
    # Need to resolve when ICON coupling is implemented.
    # comment "   cp cime_comp_mod.F90 to $clmdir/drivers/mct/main"
    #  cp $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/cime/cime_comp_mod.F90 $clmdir/cime/src/drivers/mct/main/cime_comp_mod.F90

    comment "   cp ESMF_Stubs.F90 to $clmdir/cime/src/share/esmf_wrf_timemgr"
      cp $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/cime/ESMF_Stubs.F90 $clmdir/cime/src/share/esmf_wrf_timemgr

    comment "   cp CMakeLists.txt to $clmdir/cime/src/share/util"
      cp $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/cime/CMakeLists.txt $clmdir/cime/src/share/util

    comment "   cp oas_clm_define.F90 to $clmdir/cime/src/share/util"
      cp $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/cime/oas_clm_define.F90 $clmdir/cime/src/share/util
  fi

#  comment "   cp oasis interface to $clmdir"
#    patch  $rootdir/bldsva/intf_oas3/${mList[1]}/oas3/utils/oas3 $clmdir/models/utils/ 
#  check
#    patch  $rootdir/bldsva/intf_oas3/${mList[1]}/oas3/clm/cpl_oas3 $clmdir/models/lnd/clm/src/ 
#  check
#    patch  $rootdir/bldsva/intf_oas3/${mList[1]}/oas3/atm/cpl_oas3 $clmdir/models/atm/datm/ 
#  check

#  if [[ $withOASMCT == "true" ]] ; then
#     comment "   sed replace old mod_prism includes from clm oas files"
#       sed -i "s/mod_prism_proto/mod_prism/" $clmdir/models/utils/oas3/oas_clm_vardef.F90 >> $log_file 2>> $err_file
#     check
#       sed -i "s/USE mod_prism.*//" $clmdir/models/utils/oas3/oas_clm_rcv.F90 >> $log_file 2>> $err_file
#     check
#       sed -i "s/USE mod_prism.*//" $clmdir/models/utils/oas3/oas_clm_snd.F90 >> $log_file 2>> $err_file
#     check
#       sed -i "s/USE mod_prism.*//" $clmdir/models/utils/oas3/oas_clm_define.F90 >> $log_file 2>> $err_file
#     check
#  fi




#  comment "   mkdir folder with source modifications:  $clmdir/SourceMods/"
#    mkdir -p $clmdir/SourceMods/ >> $log_file 2>> $err_file
#  check
#  comment "   cp source modifications to $clmdir"
#    patch $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/src.drv $clmdir/SourceMods/
#  check
#    patch $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/src.datm $clmdir/SourceMods/ 
#  check
#    patch $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/src.clm $clmdir/SourceMods/ 
#  check


route "${cyellow}<< substitutions_clm${cnormal}"
}


setup_clm(){
route "${cyellow}>> setupClm${cnormal}"

#  withCESM="true"
  seconds_clm=$(($hh*3600))
#  rpointer=$rundir/lnd.clmoas.rpointer
  runstep_clm=$runhours


comment "  cp namelist to rundir"
  cp ${namelist_clm} $rundir/lnd.stdin >> $log_file 2>> $err_file
check


comment "  sed num procs to namelist"
  sed "s,__nprocs__,$(($px_clm * $py_clm))," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check

  c_setup_clm

route "${cyellow}<< setupClm${cnormal}"
}
