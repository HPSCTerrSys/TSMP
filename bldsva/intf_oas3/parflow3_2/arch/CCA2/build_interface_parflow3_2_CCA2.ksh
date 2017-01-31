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

    flagsSim+="CC=cc  CXX=cc FC=ftn F77=ftn "
    flagsTools+="CC=cc FC=ftn F77=ftn "
    libsSim="$cplLib -L$ncdfPath/lib -lnetcdff"
    fcflagsSim="$cplInc -Duse_libMPI -Duse_netCDF -Duse_comm_MPI1 -DVERBOSE -DDEBUG -DTREAT_OVERLAY -I$ncdfPath/include "
    if [[ $freeDrain == "true" ]] ; then ; cflagsSim+="-DFREEDRAINAGE" ; fi
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
    comment "   cp amps_init.c and oas3_external.h to amps/oas3 folder"
    cp $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/src/amps_init.c $pfldir/pfsimulator/amps/oas3
  check
    cp $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/src/oas3_external.h $pfldir/pfsimulator/amps/oas3
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


