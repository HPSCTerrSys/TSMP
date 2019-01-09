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
# ICON interface methods
############################


c_configure_icon(){
route "${cblue}>>> c_configure_icon${cnormal}"
  file=$icondir/Makefile
  cplFlag=""
  cplLib=""
  cplInc=""
  if [[ $withOAS == "true" ]] ; then
    cplLib="$liboas $libpsmile"
    cplInc="$incpsmile"
    comment "    sed OAS flag to Makefile"
    sed -i "s@__withoas__@COUP_OAS_ICON@" $file >> $log_file 2>> $err_file
    check
    comment "    sed make.oas3 to Makefile"
    sed -i "s@__oasismakefile__@\$(oasisdir)/util/make_dir/make.oas3@" $file >> $log_file 2>> $err_file
    check
  else
    comment "    remove make.oas3 from Makefile"
    sed -i "/__oasismakefile__/d" $file >> $log_file 2>> $err_file
    check
  fi
  comment "    make clean icon"
    make clean >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_configure_icon${cnormal}"
}

c_make_icon(){
route "${cblue}>>> c_make_icon${cnormal}"
  comment "    cd to icon dir"
    cd $icondir >> $log_file 2>> $err_file
  check
  comment "    make icon"
    make -j8 -f $icondir/Makefile >> $log_file 2>> $err_file
  check

  comment "    cp icon binary to $bindir"
  cp $icondir/build/x86_64-unknown-linux-gnu/bin/icon $bindir >> $log_file 2>> $err_file
  check

route "${cblue}<<< c_make_icon${cnormal}"
}


c_substitutions_icon(){
route "${cblue}>>> c_substitutions_icon${cnormal}"
if [[ $withOAS == "true" ]]; then
  comment "    copy oas3 interface to icon/src "
    cp -R $rootdir/bldsva/intf_oas3/${mList[2]}/oas3 $icondir/src >> $log_file 2>> $err_file
  check
  comment "    replace files with coupling. Add files to icon/src "
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/mo_mpi.f90 $icondir/src/parallel_infrastructure/ >> $log_file 2>> $err_file
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon.f90 $icondir/src/drivers/ >> $log_file 2>> $err_file
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/mo_atmo_model.f90 $icondir/src/drivers/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/mo_nh_stepping.f90 $icondir/src/atm_dyn_iconam/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/mo_nwp_sfc_interface.f90 $icondir/src/lnd_phy_nwp/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/mo_nwp_rad_interface.f90 $icondir/src/atm_phy_nwp/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/mo_nwp_rrtm_interface.f90 $icondir/src/atm_phy_nwp/ >> $log_file 2>> $err_file
  check
fi
route "${cblue}<<< c_substitutions_icon${cnormal}"
}

c_setup_icon(){
route "${cblue}>>> c_setup_icon${cnormal}"

comment "  cp add_run_routines to rundir"
  cp $rootdir/bldsva/setups/common/add_run_routines $rundir >> $log_file 2>> $err_file
check

comment "  cp namelist to rundir"
  cp ${namelist_icon} $rundir >> $log_file 2>> $err_file
check

route "${cblue}<<< c_setup_icon${cnormal}"
}



############################ 
#Cosmo interface methods
############################


c_configure_cos(){
route "${cblue}>>> c_configure_cos${cnormal}"
  comment "    cd to cosmo dir"
    cd $cosdir >> $log_file 2>> $err_file
  check
  file=$cosdir/Makefile
  comment "    sed cosmo rootdir to Makefile"
    sed -i "s@__cosmoroot__@$cosdir@" $file >> $log_file 2>> $err_file
  check
  comment "    sed OAS flag to Makefile"
    sed -i "s@__withoas__@$withOAS@" $file >> $log_file 2>> $err_file
  check
  comment "    make clean cosmo"
    make clean >> $log_file 2>> $err_file
  check

    cplFlag=""
    cplLib=""
    cplInc=""
    if [[ $withOAS == "true" ]] ; then
      cplLib="$liboas $libpsmile"
      cplInc="$incpsmile"
    fi
route "${cblue}<<< c_configure_cos${cnormal}"
}

c_make_cos(){
route "${cblue}>>> c_make_cos${cnormal}"
  comment "    cd to cosmo dir"
    cd $cosdir >> $log_file 2>> $err_file
  check
  comment "    make cosmo"
    make -f $cosdir/Makefile >> $log_file 2>> $err_file
  check

#DA
  if [[ $withPDAF == "true" ]]; then
    comment "    cd to cos build dir"
      cd $cosdir/obj >> $log_file 2>> $err_file
    check
    comment "    ar cos libs"
      ar rc libcosmo.a *.o >> $log_file 2>> $err_file
    check
    comment "    cp libs to $bindir/libs"
      cp $cosdir/obj/libcosmo.a $bindir/libs >> $log_file 2>> $err_file
    check
  else
    comment "    cp cosmo binary to $bindir"
      cp $cosdir/lmparbin_pur $bindir >> $log_file 2>> $err_file
    check
  fi

route "${cblue}<<< c_make_cos${cnormal}"
}


c_substitutions_cos(){
route "${cblue}>>> c_substitutions_cos${cnormal}"
  comment "    copy oas3 interface to cosmo/src "
    patch $rootdir/bldsva/intf_oas3/${mList[2]}/oas3 $cosdir/src 
  check
  comment "    replace files with coupling. Add files to cosmo/src "
    patch "$rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/*" $cosdir/src 
  check

#DA
  if [[ $withPDAF == "true" ]]  then
    comment "    sed PDAF fix into cosmo files "  
	patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[2]}/data_parallel.f90 $cosdir/src/ 
    check
        patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[2]}/organize_data.f90 $cosdir/src/ 
    check
        patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[2]}/src_meanvalues.f90 $cosdir/src/ 
    check
        patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[2]}/src_setup.f90 $cosdir/src/ 
    check
  fi
route "${cblue}<<< c_substitutions_cos${cnormal}"
}

c_setup_cos(){
route "${cblue}>>> c_setup_cos${cnormal}"

comment "  cp namelist to rundir"
  cp ${namelist_cos} $rundir/lmrun_uc >> $log_file 2>> $err_file
check

nstop_cos=$((  ($runhours*3600 + ($(date -u '+%s' -d "${startDate}") - $(date -u '+%s' -d "${initDate}")) )  /$dt_cos  ))
#if [[ $withCESM == "false" ]] ; then ; nstop_cos=$(($nstop_cos-($cplfreq1/$dt_cos))) ; fi

comment "  sed dt to namelist"
  sed "s,dt_cos_bldsva,$dt_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check

comment "  sed number of procs to namelist"
  sed "s,nprocx_cos_bldsva,$px_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s,nprocy_cos_bldsva,$py_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check

comment "  sed gridpoints to namelist"
  sed "s,ie_tot_bldsva,$gx_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s,je_tot_bldsva,$gy_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check

comment "  sed gridpoints to namelist"
  sed "s,nbdl_cos_bldsva,$nbndlines," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check


comment "  create input dir for cosmo"
  mkdir -p $rundir/cosmo_in >> $log_file 2>> $err_file
check
comment "  fill cosmo input dir with softlinks from cosmo forcing dir"
  ln -sf $forcingdir_cos/* $rundir/cosmo_in >> $log_file 2>> $err_file
check

comment "  sed forcingdir to namelist"
  sed "s,__forcingdir__,$rundir/cosmo_in," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed rundir to namelist"
  sed "s,__rundir__,$rundir," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed stop time to namelist"
  sed "s/nstop_cos_bldsva/$nstop_cos/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed date to namelist"
  sed "s/init_y_bldsva/$(date '+%Y' -d "$initDate")/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_m_bldsva/$(date '+%m' -d "$initDate")/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_d_bldsva/$(date '+%d' -d "$initDate")/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_h_bldsva/$(date '+%H' -d "$initDate")/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check

cnt=$(( ($(date -u '+%s' -d "${startDate}") - $(date -u '+%s' -d "${initDate}"))/3600))
comment "  sed start hour to namelist"
sed "s/__hstart__/$cnt/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed restart interval to namelist"
sed "s/__nhour_restart_start__/$(($cnt+$runhours))/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check
sed "s/__nhour_restart_stop__/$(($cnt+$runhours))/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check
sed "s/__nhour_restart_incr__/1/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check

cnts=$(( ( $(date -u '+%s' -d "${startDate}") - $(date -u '+%s' -d "${initDate}")) / ${dt_cos} ))
comment "  sed output interval to namelist"
sed "s/__ncomb_start__/$cnts/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check
sed "s/__dump_cos_interval__/ $(python -c "print $dump_cos*(3600/$dt_cos)")/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check

if [[ $restfile_cos != "" ]] then
comment "  softlink restart file to input dir"
ln -s $restfile_cos $rundir/cosmo_in  >> $log_file 2>> $err_file
check
fi



comment "  cd to rundir"
  cd $rundir >> $log_file 2>> $err_file
check
comment "  run lmrun_uc clean"
  $rundir/lmrun_uc cleancluma >> $log_file 2>> $err_file
check
comment "  run lmrun_uc exe"
  $rundir/lmrun_uc execluma >> $log_file 2>> $err_file
check

if [[ $withPDAF == "true" ]] ; then
  cp $rundir/INPUT_IO $rundir/INPUT_IO_$(printf "%05d" $(($instance-$startInst)))     
fi

route "${cblue}<<< c_setup_cos${cnormal}"
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
  comment "    make oasis j16" 
    make -j16 -f $oasdir/util/make_dir/TopMakefileOasis3 oasis3_psmile >> $log_file 2>> $err_file
  check
  export SKIN_MODE=mpi
#DA
  if [[ $withPDAF == "true" ]]; then
  comment "    cp oas libs to $bindir/libs" 
    cp $libpsmile $bindir/libs >> $log_file 2>> $err_file
  check
  fi
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

#DA
  if [[ $withPDAF == "true" ]] ; then
     comment "    cp PDAF fix to ${oasdir}/lib/psmile/src"
       patch "$rootdir/bldsva/intf_DA/pdaf1_1/tsmp/mod_oasis*"  ${oasdir}/lib/psmile/src
     check
  fi
route "${cblue}<<< c_substitutions_oas${cnormal}"
}

c_setup_oas(){
route "${cblue}>>> c_setup_oas${cnormal}"

  comment "   copy cf_name_table to rundir"
    cp $rootdir/bldsva/data_oas3/cf_name_table.txt $rundir >> $log_file 2>> $err_file
  check
  comment "   copy oas namelist to rundir"
    cp $namelist_oas $rundir/namcouple >> $log_file 2>> $err_file
  check
  comment "   sed procs, gridsize & coupling freq into namcouple"

  ncpl_exe1=$nproc_cos
  ncpl_exe2=$nproc_pfl
  ncpl_exe3=1
  if [[ $withCESM == "true" || $withOASMCT == "true" ]] ; then ; ncpl_exe3=$nproc_clm ; fi


  if [[ $withICON == "true" ]]; then
    sed "s/ngiconx/$gx_icon/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq1/$cplfreq1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

    sed "s/ngclmx/$(($gx_clm*$gy_clm))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  fi


  if [[ $withPFL == "true" && $withCOS == "true" ]] then

    sed "s/nproc_exe1/$nproc_cos/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe3/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe3/$ncpl_exe3/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq1/$cplfreq1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq2/$cplfreq2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    
    sed "s/ngpflx/$gx_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngpfly/$gy_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmx/$(($gx_clm*$gy_clm))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngcosx/$(($gx_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngcosy/$(($gy_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

  fi  
if [[ $withPFL == "true" && $withCOS == "false" ]] then

    sed "s/nproc_exe1/$nproc_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe3/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq2/$cplfreq2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

    sed "s/ngpflx/$gx_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngpfly/$gy_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmx/$gx_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/$gy_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  fi
  if [[ $withPFL == "false" && $withCOS == "true" ]] then
    sed "s/nproc_exe1/$nproc_cos/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe3/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq1/$cplfreq1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

    sed "s/ngcosx/$(($gx_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngcosy/$(($gy_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmx/$gx_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/$gy_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

  fi
  rtime=$(($runhours*3600 + $cplfreq1))
  if [[ $withCESM == "true" ]] ; then ; rtime=$(($rtime+$cplfreq1)) ; fi
  comment "   sed sim time into namcouple"
    sed "s/totalruntime/$rtime/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  comment "   sed startdate into namcouple"
    sed "s/yyyymmdd/${yyyy}${mm}${dd}/" -i $rundir/namcouple  >> $log_file 2>> $err_file
  check


route "${cblue}<<< c_setup_oas${cnormal}"
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
    if [[ $withICON == "true" ]] ; then ; flags+="-oas3_icon " ; fi
    if [[ $withPFL == "true" ]] ; then ; flags+="-oas3_pfl " ; fi
    if [[ $withPFL == "false" && $withCOS == "false" && $withICON == "false" ]] ; then ; flags+="-cps_catch $cps_catch " ; fi
    flags+="-nc_inc $ncdfPath/include "
    flags+="-nc_lib $ncdfPath/lib "
    flags+="-nc_mod $ncdfPath/include "
    flags+="-mpi_inc $mpiPath/include "
    flags+="-clm_bld $clmdir/build "
    flags+="-clm_exedir $clmdir/build "
    cplInc=""

    cplInc+="-g -traceback -heap-arrays " # Mukund

      comment "adding OAS libs"
    if [[ $withOAS == "true" ]]; then
      comment "adding OAS libs"
      cplLib+="$liboas $libpsmile"
      cplInc=$incpsmile
    fi
      comment " OAS libs cplLib: $cplLib"
      comment " OAS libs cplInc: $cplInc"
      comment " OAS libs cplInc: $incpsmile"
  comment "    cd to clm build"
    cd $clmdir/build >> $log_file 2>> $err_file
  check
  cppdef=""
  if [ $cplscheme == "true" ] && [ $withICON == "false" ] ; then ; cppdef+=" -DCPL_SCHEME_F " ; fi
  comment "    configure clm"
  comment "    $clmdir/bld/configure -fc $cfc -cc $ccc $flags -fflags $cplInc -ldflags $cplLib -fopt $optComp -cppdefs $cppdef"
    $clmdir/bld/configure -fc "$cfc" -cc "$ccc" $flags -fflags "$cplInc" -ldflags "$cplLib" -fopt "$optComp" -cppdefs "$cppdef"  >> $log_file 2>> $err_file
  check
route "${cblue}<<< c_configure_clm${cnormal}"
}

c_make_clm(){
route "${cblue}>>> c_make_clm${cnormal}"
  comment "    cd to clm build"
    cd $clmdir/build >> $log_file 2>> $err_file
  check
  comment "    make clm"
    gmake -j16 -f $clmdir/build/Makefile >> $log_file 2>> $err_file
  check

#DA
  if [[ $withPDAF == "true" ]]; then
    comment "    cd to clm build dir"
      cd $clmdir/build >> $log_file 2>> $err_file
    check
    comment "    ar clm libs"
      ar rc libclm.a *.o >> $log_file 2>> $err_file
    check
    comment "    cp libs to $bindir/libs"
      cp $clmdir/build/libclm.a $bindir/libs >> $log_file 2>> $err_file
    check
  else
    comment "    cp clm binary to $bindir"
      cp $clmdir/build/clm $bindir >> $log_file 2>> $err_file
    check
  fi 
route "${cblue}<<< c_make_clm${cnormal}"
}


c_substitutions_clm(){
route "${cblue}>>> c_substitutions_clm${cnormal}"
  comment "    copy oas3 interface to clm/src "
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/mct $clmdir/src
  check
  comment "    replace hydrology. Add files to clm/bld/usr.src "
    patch "$rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/*" $clmdir/bld/usr.src 
  check
#DA
  if [[ $withPDAF == "true" ]] ; then
  comment "    copy PDAF fix to ${mList[1]}/bld/usr.src "
    patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/clmtype.F90 $clmdir/bld/usr.src
  check
    patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[1]}/clmtypeInitMod.F90 $clmdir/bld/usr.src
  check
    patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[1]}/iniTimeConst.F90 $clmdir/bld/usr.src	
  check
  fi	
route "${cblue}<<< c_substitutions_clm${cnormal}"
}


c_setup_clm(){
route "${cblue}>>> c_setup_clm${cnormal}"

comment "  sed rundir to namelist"
  sed "s,__rundir__,$rundir," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
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
comment "  sed dump interval namelist"
  sed "s,__dump_clm_interval__,$dump_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed runtime to namelist"
  runstep_clm=$((($runhours*3600 + $cplfreq1)/$dt_clm))
  sed "s,__runstep_clm_bldsva__,$runstep_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed restart file path to namelist"
    sed "s,__finidat__,${restfile_clm}," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check


comment "  add axe rights to clm namelist"
    chmod 755 $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  run clm namelist"
    $rundir/lnd.stdin >> $log_file 2>> $err_file
check

if [[ $withPDAF == "true" ]] ; then
  cp $rundir/lnd.stdin $rundir/lnd.stdin_$(printf "%05d" $(($instance-$startInst)))     
fi

route "${cblue}<<< c_setup_clm${cnormal}"
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
#DA
      if [[ $withPDAF == "true" ]] ; then
        flagsSim+="--with-amps=da " 
        flagsTools+="--with-amps=da "
      else
        flagsSim+="--with-amps=mpi1 " 
        flagsTools+="--with-amps=mpi1 "
      fi
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
    echo "$pfldir/pfsimulator/configure CC="$pcc" FC="$pfc" F77="$pf77" CXX="$pcxx" $flagsSim --enable-opt="$optComp" FCFLAGS="$fcflagsSim" CFLAGS="$cflagsSim" >> $log_file 2>> $err_file"
    $pfldir/pfsimulator/configure CC="$pcc" FC="$pfc" F77="$pf77" CXX="$pcxx" $flagsSim --enable-opt="$optComp" FCFLAGS="$fcflagsSim" CFLAGS="$cflagsSim" >> $log_file 2>> $err_file
  check
  comment "    patch pfsimulator/parflow_lib/problem_phase_rel_perm.c "
    sed -i "s@inline double VanGLookupSpline@double VanGLookupSpline@" $pfldir/pfsimulator/parflow_lib/problem_phase_rel_perm.c >> $log_file 2>> $err_file
    sed -i "s@inline double VanGLookupLinear@double VanGLookupLinear@" $pfldir/pfsimulator/parflow_lib/problem_phase_rel_perm.c >> $log_file 2>> $err_file
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
comment "    cp pfl bin to $bindir"
  cp -R $pfldir/bin $bindir >> $log_file 2>> $err_file
check
#DA
  if [[ $withPDAF == "true" ]]; then
    comment "    cp libs to $bindir/libs"
      cp $pfldir/pfsimulator/lib/* $bindir/libs >> $log_file 2>> $err_file
    check
  else
    comment "    cp binary to $bindir"
      cp $pfldir/bin/parflow $bindir >> $log_file 2>> $err_file
    check
  fi
route "${cblue}<<< c_make_pfl${cnormal}"
}

c_substitutions_pfl(){
route "${cblue}>>> c_substitutions_pfl${cnormal}"
  comment "    copy oas3 interface to parflow/pfsimulator/amps "
    patch $rootdir/bldsva/intf_oas3/${mList[3]}/oas3 $pfldir/pfsimulator/amps 
  check

  comment "    copy fix for hardwired MPI_COMM_WORLD in amps "
    patch "$rootdir/bldsva/intf_oas3/${mList[3]}/tsmp/amps*" $pfldir/pfsimulator/amps/mpi1
  check
    patch "$rootdir/bldsva/intf_oas3/${mList[3]}/tsmp/pf_pfmg*" $pfldir/pfsimulator/parflow_lib 
  check
#DA
  if [[ $withPDAF == "true" ]]; then
    comment "    sed DA amps into configure"
      sed "/\"\$with_amps\" in\s*/ a\
  da\)\\
    AMPS=da\\
  ;;
      " -i $pfldir/pfsimulator/configure $pfldir/pftools/configure >> $log_file 2>> $err_file
    check
    comment "    copy fix for PDAF into $pfldir"
      patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[3]}/parflow_proto.h $pfldir/pfsimulator/parflow_lib 
    check
      patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[3]}/solver_richards.c $pfldir/pfsimulator/parflow_lib 
    check
      patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[3]}/da $pfldir/pfsimulator/amps
    check
      sed "s/MPI_COMM_WORLD/amps_CommWorld/g" -i $pfldir/pfsimulator/parflow_lib/pf_pfmg.c
    check
      sed "s/MPI_COMM_WORLD/amps_CommWorld/g" -i $pfldir/pfsimulator/parflow_lib/pf_pfmg_octree.c
    check
      sed "s/MPI_COMM_WORLD/amps_CommWorld/g" -i $pfldir/pfsimulator/parflow_lib/pf_smg.c
    check
  fi

route "${cblue}<<< c_substitutions_pfl${cnormal}"
}

c_setup_pfl(){
route "${cblue}>>> c_setup_pfl${cnormal}"

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
    sed "s/__stop_pfl_bldsva__/$(python -c "print ${runhours} + ${base_pfl}")/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed dump interval to pfl namelist."
    sed "s/__dump_pfl_interval__/$dump_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check

  comment "   sed timing base to pfl namelist."
    sed "s/__base_pfl__/$base_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check

  comment "   sed start counter to pfl namelist."
      cnt=$(( ($(date -u '+%s' -d "${startDate}") - $(date -u '+%s' -d "${initDate}"))))
      cnt=$(python -c "print $cnt/($dump_pfl*3600.)")
      sed "s/__start_cnt_pfl__/$cnt/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check


    if [[ $restfile_pfl == "" ]] then
  comment "   sed initial condition to pfl namelist."
      sed "s/__pfl_ICPpressureType__/HydroStaticPatch/" -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file      # HydrostaticPatch > PFBFile
  check
      sed "s/__pfl_ICPpressureValue__/-5.0/" -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file      # delete this during restart run
  check
  comment "   sed delete restart file name from pfl namelist."
      sed '/__pfl_ICPpressureFileName__/d' -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file
  check
    else
  comment "   sed delete IC_Pressure value  from pfl namelist." 
      sed '/__pfl_ICPpressureValue__/d' -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed initial condition to pfl namelist."
      sed "s/__pfl_ICPpressureType__/PFBFile/" -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file      # HydrostaticPatch > PFBFile
  check
      sed "s,__pfl_ICPpressureFileName__,$restfile_pfl," -i $rundir/coup_oas.tcl  >> $log_file 2>> $err_file
  check
    fi

  export PARFLOW_DIR=$bindir
  comment "   cd to rundir."
    cd $rundir >> $log_file 2>> $err_file
  check

  comment "   create parflow db with tclsh from namelist."
    tclsh $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check


route "${cblue}<<< c_setup_pfl${cnormal}"
}


############################ 
#PDAF interface methods
############################


c_setup_pdaf(){
route "${cblue}>>> c_setup_da${cnormal}"
  comment "   copy pdaf namelist to rundir."
    cp $namelist_da $rundir/enkfpf.par >> $log_file 2>> $err_file
  check 
  comment "   sed num instances into pdaf namelist."
    sed "s/__ninst__/$(($numInst-$startInst))/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed pflname into pdaf namelist."
    sed "s/__pflname__/$pflrunname/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed pflproc into pdaf namelist."
    sed "s/__pflproc__/$nproc_pfl/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed dt into pdaf namelist."
    sed "s/__dt__/$dt_pfl/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed endtime into pdaf namelist."
    sed "s/__endtime__/$(python -c "print ${runhours} + ${base_pfl}")/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed clmproc into pdaf namelist."
    sed "s/__clmproc__/$nproc_clm/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed cosproc into pdaf namelist."
    sed "s/__cosproc__/$nproc_cos/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check 
  comment "   sed dtmult into pdaf namelist."
    sed "s/__dtmult__/$(python -c "print ${dt_pfl} * 3600 / ${dt_cos}")/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check 

route "${cblue}<<< c_setup_da${cnormal}"
}
