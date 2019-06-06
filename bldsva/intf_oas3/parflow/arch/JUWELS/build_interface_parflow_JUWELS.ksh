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

    if [[ $readCLM == "true" ]] ; then ; cplInc+=" -DREADCLM " ; fi

    flagsSim=" "
    pcc="${profComp} $mpiPath/bin/mpicc"
    pfc="${profComp} $mpiPath/bin/mpif90"
    pf77="${profComp} $mpiPath/bin/mpif77"
    pcxx="${profComp} $mpiPath/bin/mpic++"
    if [[ $profiling == "scalasca" ]]; then
      pcc="scorep-mpicc"
      pfc="scorep-mpif90"
      pf77="scorep-mpif77"
      pcxx="scorep-mpic++"
      flagsTools+="CC=scorep-mpicc FC=scorep-mpif90 F77=scorep-mpif77 "
    fi
    flagsTools+="CC=$mpiPath/bin/mpicc FC=$mpiPath/bin/mpif90 F77=$mpiPath/bin/mpif77 "
    libsSim="$cplLib -L$ncdfPath/lib -lnetcdff"
    fcflagsSim="$cplInc -Duse_libMPI -Duse_netCDF -Duse_comm_MPI1 -DVERBOSE -DDEBUG -DTREAT_OVERLAY -I$ncdfPath/include "
    cflagsSim=" -fopenmp "
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
  comment "   cp amps_init.c and oas3_external.h to amps/oas3 folder"
    patch $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/amps_init.c $pfldir/pfsimulator/amps/oas3
  check
    patch $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/oas3_external.h $pfldir/pfsimulator/amps/oas3
  check
 
  c_substitutions_pfl

#  comment "   sed hypre Boxcreate fix into /parflow_lib/pf_pfmg_octree.c"
#    sed -i "s/hypre_BoxCreate()/hypre_BoxCreate(ndim)/" $pfldir/pfsimulator/parflow_lib/pf_pfmg_octree.c >> $log_file 2>> $err_file
#  check

  comment "    copy nl_function_eval.c with free drainage feature to parflow/pfsimulator/parflow_lib "
    patch $rootdir/bldsva/intf_oas3/${mList[3]}/tsmp/nl_function_eval.c $pfldir/pfsimulator/parflow_lib/nl_function_eval.c 
  check 
  comment "   cp new pf_pfmg_octree.c to /parflow_lib/"
    patch $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/pf_pfmg_octree.c  $pfldir/pfsimulator/parflow_lib/ 
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


