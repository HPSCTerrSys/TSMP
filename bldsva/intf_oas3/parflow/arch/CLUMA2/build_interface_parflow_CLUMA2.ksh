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

    flagsSim+="CC=$mpiPath/bin/mpicc  CXX=$mpiPath/bin/mpic++ FC=$mpiPath/bin/mpif90 F77=$mpiPath/bin/mpif77 "
    flagsTools+="CC=$mpiPath/bin/mpicc FC=$mpiPath/bin/mpif90 F77=$mpiPath/bin/mpif77 "
    libsSim="$cplLib -L$ncdfPath/lib -lnetcdff"
    fcflagsSim="$cplInc -Duse_libMPI -Duse_netCDF -Duse_comm_MPI1 -DVERBOSE -DDEBUG -DTREAT_OVERLAY -I$ncdfPath/include "
    c_configure_pfl


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
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/amps_init.c $pfldir/pfsimulator/amps/oas3
  check
    cp $rootdir/bldsva/intf_oas3/parflow/arch/$platform/src/oas3_external.h $pfldir/pfsimulator/amps/oas3
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
  comment "   copy parflow namelist to rundir."
    cp $namelist_pfl $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed nproc x to pfl namelist."
    sed "s/__nprocx_pfl_bldsva__/$px_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed nproc y to pfl namelist."
    sed "s/__nprocy_pfl_bldsva__/$py_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed gridpoints x to pfl namelist."
    sed "s/__ngpflx_bldsva__/$gx_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed gridpoints y to pfl namelist."
    sed "s/__ngpfly_bldsva__/$gy_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed forcingdir to pfl namelist."
    sed "s,__forcingdir__,$rundir," -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed dt to pfl namelist."
    sed "s/__dt_pfl_bldsva__/$dt_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed end time to pfl namelist."
    sed "s/__stop_pfl_bldsva__/$runhours/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed start counter to pfl namelist."
    sed "s/__start_cnt_pfl__/0/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
    if [[ $restDate == "" ]] then
  comment "   sed initial condition to pfl namelist."
      sed "s/__pfl_ICPpressureType__/HydroStaticPatch/" -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file      # HydrostaticPatch > PFBFile
  check
      sed "s/__pfl_ICPpressureValue__/-5.0/" -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file      # comment this during restart run
  check
    else
  comment "   sed initial condition to pfl namelist."
      sed "s/__pfl_ICPpressureType__/PFBFile/" -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file      # HydrostaticPatch > PFBFile
  check
      sed "s,__pfl_ICPpressureFileName__,$pfbfilename," -i $rundir/coup_oas.tcl  >> $log_file 2>> $err_file       # comment this during restart run
  check
    fi


  comment "   cd to rundir."
    cd $rundir >> $log_file 2>> $err_file
  check
  comment "   sed pfl_dir into coup_oas.tcl"
    sed "s,lappend auto_path.*,lappend auto_path $pfldir/bin," -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   create parflow db with tclsh from namelist."
    tclsh $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
route "${cblue}<< setup_pfl${cnormal}"
}


