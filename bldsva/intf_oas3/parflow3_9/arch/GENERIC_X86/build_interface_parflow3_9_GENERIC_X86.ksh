#! /bin/ksh

always_pfl(){
route "${cblue}>> always_pfl${cnormal}"
route "${cblue}<< always_pfl${cnormal}"
}

configure_pfl(){
route "${cblue}>> configure_pfl${cnormal}"
  comment "   cp new Makefile.in to /pfsimulator/parflow_exe/"
    cp $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/config/Makefile.in $pfldir/pfsimulator/parflow_exe/ >> $log_file 2>> $err_file
  check

    if [[ $withOAS == "true" ]]; then
      cplLib="$liboas $libpsmile"
      cplInc="$incpsmile"
    fi  

    if [[ $readCLM == "true" ]] ; then ; cplInc+=" -DREADCLM " ; fi
    flagsSim=" "
    pcc="${profComp} $mpiPath/bin/mpicc"
    pfc="${profComp} $mpiPath/bin/mpif90"
    pf77="${profComp} $mpiPath/bin/mpif77"
    pcxx="${profComp} $mpiPath/bin/mpic++"
#    flagsTools+="CC=$mpiPath/bin/mpicc FC=$mpiPath/bin/mpif90 F77=$mpiPath/bin/mpif77i FCFLAGS=\"-I/home/ag/lib/libhdf5/include -L/home/ag/lib/libhdf5/lib -lhdf5 -lhdf5_hl -lhdf5_fortran\" CFLAGS=\"-I/home/ag/lib/libhdf5/include -L/home/ag/lib/libhdf5/lib -lhdf5 -lhdf5_hl -lhdf5_fortran\" "
    flagsTools+="CC=$mpiPath/bin/mpicc FC=$mpiPath/bin/mpif90 F77=$mpiPath/bin/mpif77 "
   # LDFLAGS=\"-L/home/ag/lib/libhdf5/lib -lhdf5_hl\" "

    comment "pcc : $pcc , pfc : $pfc , pf77 : $pf77   "
    check
    if [[ $profiling == "scalasca" ]]; then
      pcc="scorep-mpicc"
      pfc="scorep-mpif90"
      pf77="scorep-mpif77"
      pcxx="scorep-mpic++"
      flagsTools+="CC=scorep-mpicc FC=scorep-mpif90 F77=scorep-mpif77 "
    fi
    libsSim="$cplLib -L$ncdfPath/lib -lnetcdff -lnetcdf"
    comment "libsSim: $libsSim "
    check
    #fcflagsSim="$cplInc -Duse_libMPI -Duse_netCDF -Duse_comm_MPI1 -DVERBOSE -DDEBUG -DTREAT_OVERLAY -I$ncdfPath/include "
    fcflagsSim="$cplInc -Duse_libMPI -Duse_netCDF -Duse_comm_MPI1 -DVERBOSE -DDEBUG -DTREAT_OVERLAY -I$ncdfPath/include -L$ncdfPath/lib -lnetcdff -lnetcdf -I$hyprePath/include -L$hyprePath/lib -lHYPRE -I$siloPath/include -L$siloPath/lib -lsiloh5"
    cflagsSim=" -fopenmp "
    comment "fcflagsSim : $fcflagsSim "
    check
    comment "flagsTools flagsTools flagsTools flagsTools : $flagsTools "
    check
    if [[ $freeDrain == "true" ]] ; then ; cflagsSim+="-DFREEDRAINAGE" ; fi

    c_configure_pfl


  comment "   sed correct linker command in pfsimulator"
    sed -i 's@-lmpi"@@' $pfldir/pfsimulator/config/Makefile.config >> $log_file 2>> $err_file
  check
  comment "   sed CXX to CC in Makefile.config in  pfsimulator"
    sed -i 's@CXX@CC@' $pfldir/pfsimulator/config/Makefile.config >> $log_file 2>> $err_file
  check
  comment "   sed correct linker command in pftools"
    sed -i 's@-lmpi"@@' $pfldir/pftools/config/Makefile.config >> $log_file 2>> $err_file
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
    comment "   cp amps_init.c and oas3_external.h to amps/oas3 folder"
    patch $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/src.$compiler/amps_init.c $pfldir/pfsimulator/amps/oas3
  check
    patch $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/src.$compiler/oas3_external.h $pfldir/pfsimulator/amps/oas3
  check
 
  comment "   cp new pf_pfmg_octree.c to /parflow_lib/"
    patch $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/src.$compiler/pf_pfmg_octree.c  $pfldir/pfsimulator/parflow_lib/ 
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
  echo "bindir setup_pfl setup_pfl setup_pfl:   $bindir"
  #if [[ $bindir == "" ]] then
  #   bindir="$rootdir/bin/${platform}_${version}_${combination}"
  #fi
  #check
  echo "************* c_setup_pfl ************** : $c_setup_pfl"
  c_setup_pfl

route "${cblue}<< setup_pfl${cnormal}"
}


