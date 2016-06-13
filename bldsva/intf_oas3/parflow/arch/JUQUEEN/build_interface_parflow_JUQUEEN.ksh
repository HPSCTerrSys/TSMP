#! /bin/ksh

always_pfl(){
route "${cblue}>> always_pfl${cnormal}"
route "${cblue}<< always_pfl${cnormal}"
}

configure_pfl(){
route "${cblue}>> configure_pfl${cnormal}"
  comment "   cp new Makefile.in to /pfsimulator/parflow_exe/"
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/config/Makefile.in $pfldir/pfsimulator/parflow_exe/ >> $log_file 2>> $err_file
  check

    if [[ $withOAS == "true" ]]; then
      cplLib="$liboas $libpsmile"
      cplInc="$incpsmile"
    fi  

    flagsSim+="CC=$mpiPath/bin/mpixlc_r  CXX=$mpiPath/bin/mpixlcxx_r FC=$mpiPath/bin/mpixlf90_r F77=$mpiPath/bin/mpixlf77_r "
    flagsTools+="CC=gcc FC=gfortran F77=gfortran "
    libsSim="$cplLib -L$ncdfPath/lib -lnetcdf"
    fcflagsSim="-qfree=f90 -qsuffix=cpp=F90 -qnoextname $cplInc -WF,-Duse_libMPI -WF,-Duse_netCDF -WF,-Duse_comm_MPI1 -WF,-DVERBOSE -WF,-DDEBUG -WF,-DTREAT_OVERLAY -I$ncdfPath/include "
    c_configure_pfl

  comment "   sed correct linker command in pfsimulator"
    sed "s/ gfortran /  -lgfortran  /g" -i $pfldir/pfsimulator/config/Makefile.config >> $log_file 2>> $err_file
  check
    sed "s/ m /  -lm  /g" -i $pfldir/pfsimulator/config/Makefile.config >> $log_file 2>> $err_file
  check
    sed "s/-l /  /g" -i $pfldir/pfsimulator/config/Makefile.config >> $log_file 2>> $err_file
  check
  print -n "   sed correct linker command in pftools"
    sed "s/ gfortran /  -lgfortran  /g" -i $pfldir/pftools/config/Makefile.config >> $log_file 2>> $err_file
  check
    sed "s/ m /  -lm  /g" -i $pfldir/pftools/config/Makefile.config >> $log_file 2>> $err_file
  check
    sed "s/-l /  /g" -i $pfldir/pftools/config/Makefile.config >> $log_file 2>> $err_file
  check

route "${cblue}<< configure_pfl${cnormal}"
}

make_pfl(){
route "${cblue}>> make_pfl${cnormal}"
  c_make_pfl
route "${cblue}<< make_pfl${cnormal}"
}


substitutions_pfl(){
route "${cblue}>> substitutions_pfl${cnormal}"
  c_substitutions_pfl
 

  print -n "   cp new files with Fortran underscore fix"
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/amps_init.c         $pfldir/pfsimulator/amps/oas3/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/oas3_coupler.h      $pfldir/pfsimulator/amps/oas3/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/oas3_external.h     $pfldir/pfsimulator/amps/oas3/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/oas_pfl_vardef.F90  $pfldir/pfsimulator/amps/oas3/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/parflow_proto_f90.h $pfldir/pfsimulator/parflow_lib/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/parflow_proto_f.h   $pfldir/pfsimulator/parflow_lib/ >> $log_file 2>> $err_file
  check
  print -n "   cp new files with little endian fix"
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/amps_proto.h        $pfldir/pfsimulator/amps/mpi1/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/amps_proto.h        $pfldir/pfsimulator/amps/oas3/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/amps_io.c           $pfldir/pfsimulator/amps/common/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/parflow_config.h.in $pfldir/pfsimulator/config/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/tools_io.h          $pfldir/pftools/ >> $log_file 2>> $err_file
  check

    if [[ $withOASMCT == "true" ]] ; then 
      comment "   sed replace old mod_prism includes from pfl oas files"
        sed -i "s/mod_prism_proto/mod_prism/" $pfldir/pfsimulator/amps/oas3/oas_pfl_vardef.F90 >> $log_file 2>> $err_file
      check
        sed -i "s/USE mod_prism.*//" $pfldir/pfsimulator/amps/oas3/oas_pfl_define.F90 >> $log_file 2>> $err_file
      check
        sed -i "s/USE mod_prism.*//" $pfldir/pfsimulator/amps/oas3/oas_pfl_snd.F90 >> $log_file 2>> $err_file
      check
        sed -i "s/USE mod_prism.*//" $pfldir/pfsimulator/amps/oas3/oas_pfl_rcv.F90 >> $log_file 2>> $err_file
      check
    fi

route "${cblue}<< substitutions_pfl${cnormal}"
}


setup_pfl(){
route "${cblue}>> setup_pfl${cnormal}"
 
  c_setup_pfl 

route "${cblue}<< setup_pfl${cnormal}"
}


