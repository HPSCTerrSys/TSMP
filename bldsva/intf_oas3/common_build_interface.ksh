#! /bin/ksh

###
# This are the default implementations for the interfaces
# Interface methods are:
#
# - always
# - configure
# - make
# - substitutions



# only platform and version independent stuff should be done here



############################ 
#Cosmo interface methods
############################


c_configure_cos(){
route "${cblue}>>> c_configure_cos${cnormal}"
  comment "    make clean cosmo"
    make -f $cosdir/Makefile clean >> $log_file 2>> $err_file
  check

    file=$cosdir/Makefile
    cplFlag=""
    cplLib=""
    cplInc=""
    if [[ $withOAS == "true" ]] ; then
      cplLib="$liboas $libpsmile"
      cplInc="$incpsmile"
    fi
  comment "    sed cosmo rootdir to Makefile"
    sed -i "s@__cosmoroot__@$cosdir@" $file >> $log_file 2>> $err_file
  check
  comment "    sed OAS flag to Makefile"
    sed -i "s@__withoas__@$withOAS@" $file >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_configure_cos${cnormal}"
}

c_make_cos(){
route "${cblue}>>> c_make_cos${cnormal}"
  comment "    make cosmo"
    make -f $cosdir/Makefile >> $log_file 2>> $err_file
  check
  comment "    cp cosmo binary to $bindir"
    cp $cosdir/lmparbin_pur $bindir >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_make_cos${cnormal}"
}


c_substitutions_cos(){
route "${cblue}>>> c_substitutions_cos${cnormal}"
  comment "    copy oas3 interface to cosmo/src "
    cp -R $rootdir/bldsva/intf_oas3/${mList[2]}/oas3 $cosdir/src >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_substitutions_cos${cnormal}"
}



############################ 
#OASIS interface methods
############################

c_configure_oas(){
route "${cblue}>>> c_configure_oas${cnormal}"
  comment "    sed oasis rootdir to Makefile"
    sed -i "s@__oasisroot__@$oasdir@" $file >> $log_file 2>> $err_file
  check
  comment "    sed platform to Makefile"
    sed -i "s@__platform__@$platform@" $file >> $log_file 2>> $err_file
  check
  comment "    make clean oasis"
    make -f $oasdir/util/make_dir/TopMakefileOasis3 realclean >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_configure_oas${cnormal}"
}

c_make_oas(){
route "${cblue}>>> c_make_oas${cnormal}"
    export SKIN_MODE=none
  comment "    make oasis" 
    make -f $oasdir/util/make_dir/TopMakefileOasis3 oasis3_psmile >> $log_file 2>> $err_file
  check
    export SKIN_MODE=mpi
route "${cblue}<<< c_make_oas${cnormal}"
}


c_substitutions_oas(){
route "${cblue}>>> c_substitutions_oas${cnormal}"
  comment "    sed absolut include paths to Makefile"
    sed -i "s@include make.inc@include $oasdir/util/make_dir/make.inc@" ${oasdir}/util/make_dir/TopMakefileOasis3 >> $log_file 2>> $err_file
  check
  comment "    sed absolut makefile path to Makefile"
    sed -i "s@ TopMakefileOasis3@ $oasdir/util/make_dir/TopMakefileOasis3@" ${oasdir}/util/make_dir/TopMakefileOasis3 >> $log_file 2>> $err_file
  check
  comment "    sed usermakefile to make.inc"
    sed -i "s@include.*@include $oasdir/util/make_dir/make.oas3@" ${oasdir}/util/make_dir/make.inc >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_substitutions_oas${cnormal}"
}




############################ 
#CLM interface methods
############################


c_configure_clm(){
route "${cblue}>>> c_configure_clm${cnormal}"
  comment "    clean clm by removing build dir"
    rm -rf $clmdir/build >> $log_file 2>> $err_file
  check
  comment "    create new build dir"
    mkdir -p $clmdir/build >> $log_file 2>> $err_file
  check

    spmd="on"       # settings are [on   | off       ] (default is off)
    maxpft="1"        # settings are 4->17               (default is 4)
    rtm="off"      # settings are [on   | off       ] (default is off) 
    cps_catch="off"       # settings are [on   | off       ] (default is off)
    usr_src="$clmdir/bld/usr.src "
    if [[ $spmd == "on" ]] ; then ; flags+="-spmd " ; fi
    if [[ $spmd == "off" ]] ; then ; flags+="-nospmd " ; fi
    flags+="-maxpft $maxpft -rtm $rtm -usr_src $usr_src "
    if [[ $withCOS == "true" ]] ; then ; flags+="-oas3_cos " ; fi
    if [[ $withPFL == "true" ]] ; then ; flags+="-oas3_pfl " ; fi
    if [[ $withPFL == "false" && $withCOS == "false" ]] ; then ; flags+="-cps_catch $cps_catch " ; fi
    flags+="-nc_inc $ncdfPath/include "
    flags+="-nc_lib $ncdfPath/lib "
    flags+="-nc_mod $ncdfPath/include "
    flags+="-mpi_inc $mpiPath/include "
    flags+="-mpi_lib $mpiPath/lib "
    flags+="-clm_bld $clmdir/build "
    flags+="-clm_exedir $clmdir/build "
    cplInc=""
    if [[ $withOAS == "true" ]]; then
      cplLib+="$liboas $libpsmile"
      cplInc=$incpsmile
    fi
  comment "    cd to clm build"
    cd $clmdir/build >> $log_file 2>> $err_file
  check
  comment "    configure clm"
    $clmdir/bld/configure $flags -fflags "$cplInc" -ldflags "$cplLib" -fopt "$optComp" >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_configure_clm${cnormal}"
}

c_make_clm(){
route "${cblue}>>> c_make_clm${cnormal}"
  comment "    cd to clm build"
    cd $clmdir/build >> $log_file 2>> $err_file
  check
  comment "    make clm"
    gmake -f $clmdir/build/Makefile >> $log_file 2>> $err_file
  check
  comment "    cp clm binary to $bindir"
    cp $clmdir/build/clm $bindir >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_make_clm${cnormal}"
}


c_substitutions_clm(){
route "${cblue}>>> c_substitutions_clm${cnormal}"
  comment "    copy oas3 interface to clm/src "
    cp -R $rootdir/bldsva/intf_oas3/${mList[1]}/oas3 $clmdir/src >> $log_file 2>> $err_file
  check
  comment "    replace hydrology. Add files to clm/bld/usr.src "
    cp $rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/* $clmdir/bld/usr.src >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_substitutions_clm${cnormal}"
}


############################ 
#Parflow interface methods
############################


c_configure_pfl(){


route "${cblue}>>> c_configure_pfl${cnormal}"


    if [[ $withOAS == "true" ]] ; then 
      flagsSim+="--with-amps=oas3 --with-oas3 "  
      flagsTools+="--with-amps=oas3 --with-oas3 "
    else 
      flagsSim+="--with-amps=mpi1 " 
      flagsTools+="--with-amps=mpi1 "
    fi

    flagsSim+="--prefix=$pfldir --with-hypre=$hyprePath --with-silo=$siloPath --with-amps-sequential-io --enable-timing"
    flagsTools+="--prefix=$pfldir --with-hypre=$hyprePath --with-silo=$siloPath --with-tcl=$tclPath --with-amps-sequential-io"

    export SKIN_MODE=none
  comment "    cd to pfsimulator"
    cd $pfldir/pfsimulator >> $log_file 2>> $err_file
  check

    if [[ -e "$pfldir/pfsimulator/Makefile" ]] ; then
      comment "    make pfsimulator very clean"
        make -f $pfldir/pfsimulator/Makefile veryclean >> $log_file 2>> $err_file
      check
    fi 

  comment "    configure pfsimulator"
    $pfldir/pfsimulator/configure $flagsSim --enable-opt="$optComp" FCFLAGS="$fcflagsSim" >> $log_file 2>> $err_file
  check
  comment "    cd to pftools"
    cd $pfldir/pftools >> $log_file 2>> $err_file
  check

    if [[ -e "$pfldir/pftools/Makefile" ]] ; then
      comment "    make pftools very clean"
        make -f $pfldir/pftools/Makefile veryclean >> $log_file 2>> $err_file
      check
    fi

  comment "    configure pftools"
    $pfldir/pftools/configure $flagsTools >> $log_file 2>> $err_file
  check
  export SKIN_MODE=mpi

  comment "    sed libs to /parflow_exe/Makefile"
    sed -i "s@__libs__@$libsSim@" $pfldir/pfsimulator/parflow_exe/Makefile >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_configure_pfl${cnormal}"
}

c_make_pfl(){
route "${cblue}>>> c_make_pfl${cnormal}"
comment "    cd to pfsimulator" 
  cd $pfldir/pfsimulator >> $log_file 2>> $err_file
check
comment "    make pfsimulator"
  make -f $pfldir/pfsimulator/Makefile >> $log_file 2>> $err_file
check
comment "    make install pfsimulator"
  make -f $pfldir/pfsimulator/Makefile install >> $log_file 2>> $err_file
check


comment "    cd to pftools"
  cd $pfldir/pftools >> $log_file 2>> $err_file
check
comment "    make pftools"
  make -f $pfldir/pftools/Makefile >> $log_file 2>> $err_file
check
comment "    make install pftools"
  make -f $pfldir/pftools/Makefile install >> $log_file 2>> $err_file
check

comment "    cp binary to $bindir"
  cp $pfldir/bin/parflow $bindir >> $log_file 2>> $err_file
check

route "${cblue}<<< c_make_pfl${cnormal}"
}

c_substitutions_pfl(){
route "${cblue}>>> c_substitutions_pfl${cnormal}"
  comment "    copy oas3 interface to parflow/pfsimulator/amps "
    cp -R $rootdir/bldsva/intf_oas3/${mList[3]}/oas3 $pfldir/pfsimulator/amps >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_substitutions_pfl${cnormal}"
}


############################ 
#setup interface methods
############################

