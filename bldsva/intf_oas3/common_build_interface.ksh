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
route "${cyellow}>>> c_configure_icon${cnormal}"
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
route "${cyellow}<<< c_configure_icon${cnormal}"
}

c_make_icon(){
route "${cyellow}>>> c_make_icon${cnormal}"
  comment "    cd to icon dir"
    cd $icondir >> $log_file 2>> $err_file
  check
  comment "    make icon"
    export SCOREP_WRAPPER=on
    make -j16 -f $icondir/Makefile >> $log_file 2>> $err_file
  check

  comment "    cp icon binary to $bindir"
#  comment " SPo binary refSetup $refSetup mList ${mList[2]} "
  if [[ ${mList[2]} == "icon2-622" ]] ; then
    cp $icondir/bin/icon $bindir >> $log_file 2>> $err_file
  else
    cp $icondir/build/x86_64-unknown-linux-gnu/bin/icon $bindir >> $log_file 2>> $err_file  
  fi
  check

route "${cyellow}<<< c_make_icon${cnormal}"
}


c_substitutions_icon(){
route "${cyellow}>>> c_substitutions_icon${cnormal}"
if [[ $withOAS == "true" ]]; then
  comment "    copy oas3 interface to icon/src "
    cp -R $rootdir/bldsva/intf_oas3/${mList[2]}/oas3 $icondir/src >> $log_file 2>> $err_file
  check
  comment "    replace coupling files in ICON code. Add files to icon/src "
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
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/mo_nwp_rg_interface.f90 $icondir/src/atm_phy_nwp/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/mo_nwp_turbtrans_interface.f90 $icondir/src/atm_phy_nwp/ >> $log_file 2>> $err_file
  check
fi # SPo
if [[ ${mList[2]} == "icon2-1" ]] ; then
  comment "    replace icon-ccs files in ICON code. Add files to icon/src "
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_nonhydro_types.f90 $icondir/src/atm_dyn_iconam/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_nonhydro_state.f90 $icondir/src/atm_dyn_iconam/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_vdfmain.f90 $icondir/src/atm_phy_edmf/ >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_les_turb_interface.f90 $icondir/src/atm_phy_les >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_sgs_turbulence.f90 $icondir/src/atm_phy_les >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_turbulent_diagnostic.f90 $icondir/src/atm_phy_les >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_interface_les.f90 $icondir/src/atm_phy_les >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_surface_les.f90 $icondir/src/atm_phy_les >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_sgs_turbmetric.f90 $icondir/src/atm_phy_les >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_nwp_turbtrans_interface.f90 $icondir/src/atm_phy_nwp >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_radiation.f90 $icondir/src/atm_phy_schemes >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_radiation_config.f90 $icondir/src/configure_model >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_les_config.f90 $icondir/src/configure_model >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_nwp_lnd_state.f90 $icondir/src/lnd_phy_nwp >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_soil_ml.f90 $icondir/src/lnd_phy_schemes >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_nwp_soil_init.f90 $icondir/src/lnd_phy_schemes >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_nh_testcases_nml.f90 $icondir/src/namelists >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_les_nml.f90 $icondir/src/namelists >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_radiation_nml.f90 $icondir/src/namelists >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_ext_data_init.f90 $icondir/src/shr_horizontal >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_nh_torus_exp.f90 $icondir/src/testcases >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/icon-ccs/mo_nh_testcases.f90 $icondir/src/testcases >> $log_file 2>> $err_file
  check
fi
route "${cyellow}<<< c_substitutions_icon${cnormal}"
}

c_setup_icon(){
route "${cyellow}>>> c_setup_icon${cnormal}"

comment "  cp add_run_routines to rundir"
  cp $rootdir/bldsva/setups/common/add_run_routines $rundir >> $log_file 2>> $err_file
check

comment "  cp namelist to rundir"
  cp ${namelist_icon} $rundir >> $log_file 2>> $err_file
check

comment "  sed dt to namelist"
  sed "s,__dt_icon_bldsva__,$dt_icon," -i $rundir/NAMELIST_icon >> $log_file 2>> $err_file
check

comment "  sed start time to namelist"
  dSD=($defaultStartDate)
  sed "s,__starttime_icon_bldsva__,${dSD[0]}T${dSD[1]}:00:00Z," -i $rundir/icon_master.namelist >> $log_file 2>> $err_file
  sed "s,__starttime_icon_bldsva__,${dSD[0]}T${dSD[1]}:00:00Z," -i $rundir/NAMELIST_icon >> $log_file 2>> $err_file
check

comment "  sed end time to namelist"
  dED=($defaultEndDate)
  sed "s,__endtime_icon_bldsva__,${dED[0]}T${dED[1]}:00:00Z," -i $rundir/icon_master.namelist >> $log_file 2>> $err_file
  sed "s,__endtime_icon_bldsva__,${dED[0]}T${dED[1]}:00:00Z," -i $rundir/NAMELIST_icon >> $log_file 2>> $err_file
check

route "${cyellow}<<< c_setup_icon${cnormal}"
}



############################ 
#Cosmo interface methods
############################


c_configure_cos(){
route "${cyellow}>>> c_configure_cos${cnormal}"
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
route "${cyellow}<<< c_configure_cos${cnormal}"
}

c_make_cos(){
route "${cyellow}>>> c_make_cos${cnormal}"
  comment "    cd to cosmo dir"
    cd $cosdir >> $log_file 2>> $err_file
  check
  comment "    make cosmo"
    export SCOREP_WRAPPER=on
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

route "${cyellow}<<< c_make_cos${cnormal}"
}


c_substitutions_cos(){

route "${cyellow}>>> c_substitutions_cos${cnormal}"

  comment "    copy oas3 interface to cosmo/src "
    patch $rootdir/bldsva/intf_oas3/${mList[2]}/oas3 $cosdir/src 
  check

  if [[ ${mList[2]} == cosmo5_1 ]] ; then
    cp  $rootdir/cosmo5_1/LOCAL/TWOMOM/src_twomom_sb* $cosdir/src
    check
  fi

  if [[ ${mList[2]} == cosmo4_21 ]] ; then
    comment "    replace files with coupling. Add files to cosmo/src "
      patch "$rootdir/bldsva/intf_oas3/${mList[2]}/tsmp/*" $cosdir/src 
    check
  fi

  if [[ ${mList[2]} == cosmo5_1 ]] ; then
    comment "	copy the diff files to cosmo src : from $rootdir/bldsva/intf_oas3/${mList[2]}/pfile"
      cp $rootdir/bldsva/intf_oas3/${mList[2]}/pfile/patch* $cosdir/src 
    check
    comment "     apply diff files on the original files using patch command in $cosdir/src "
     /usr/bin/patch  -d $cosdir/src -i patch_src_radiation.f90.diff -o src_radiation1.f90
      cp $cosdir/src/src_radiation1.f90 $cosdir/src/src_radiation.f90
    check
      /usr/bin/patch  -d $cosdir/src -i patch_data_fields.f90.diff -o data_fields1.f90  
      cp $cosdir/src/data_fields1.f90  $cosdir/src/data_fields.f90
    check
    #    patch  -d $cosdir/src -i patch_phillips_nucleation.incf.diff  -o phillips_nucleation1.incf >> $log_pfile 2>> $err_pfile
    #    cp $cosdir/src/phillips_nucleation1.incf $cosdir/src/phillips_nucleation.incf
    check
      /usr/bin/patch  -d $cosdir/src -i patch_lmorg.f90.diff -o lmorg1.f90 
      cp $cosdir/src/lmorg1.f90 $cosdir/src/lmorg.f90
    check
      /usr/bin/patch  -d $cosdir/src -i patch_src_artifdata.f90.diff -o src_artifdata1.f90
      cp $cosdir/src/src_artifdata1.f90 $cosdir/src/src_artifdata.f90
    check
      /usr/bin/patch  -d $cosdir/src -i patch_src_setup_vartab.f90.diff -o src_setup_vartab1.f90
      cp $cosdir/src/src_setup_vartab1.f90 $cosdir/src/src_setup_vartab.f90
    check
      /usr/bin/patch -d $cosdir/src -i patch_src_twomom_sb.f90.diff -o src_twomom_sb1.f90
      cp $cosdir/src/src_twomom_sb1.f90  $cosdir/src/src_twomom_sb.f90
    check
      /usr/bin/patch  -d $cosdir/src -i patch_environment.f90.diff -o environment1.f90
      cp $cosdir/src/environment1.f90 $cosdir/src/environment.f90
    check
      /usr/bin/patch  -d $cosdir/src -i patch_organize_physics.f90.diff -o organize_physics1.f90
      cp $cosdir/src/organize_physics1.f90 $cosdir/src/organize_physics.f90
    check
      /usr/bin/patch  -d $cosdir/src -i patch_src_allocation.f90.diff -o src_allocation1.f90
      cp $cosdir/src/src_allocation1.f90 $cosdir/src/src_allocation.f90
    check
     /usr/bin/patch  -d $cosdir/src -i patch_src_gridpoints.f90.diff -o src_gridpoints1.f90
      cp $cosdir/src/src_gridpoints1.f90 $cosdir/src/src_gridpoints.f90
    check
      /usr/bin/patch -d $cosdir/src -i patch_src_runge_kutta.f90.diff -o src_runge_kutta1.f90
      cp $cosdir/src/src_runge_kutta1.f90 $cosdir/src/src_runge_kutta.f90
    check
     /usr/bin/patch -d $cosdir/src -i patch_src_slow_tendencies_rk.f90.diff -o src_slow_tendencies_rk1.f90
      cp $cosdir/src/src_slow_tendencies_rk1.f90 $cosdir/src/src_slow_tendencies_rk.f90
    check
    #    patch -d $cosdir/src -i patch_src_twomom_sb_interface.f90.diff  -o src_twomom_sb_interface1.f90
    #    cp $cosdir/src/src_twomom_sb_interface1.f90 $cosdir/src/src_twomom_sb_interface.f90
    #  check

      rm -rf $cosdir/src/*1.f90
      rm -rf $cosdir/src/*1.incf

    # copy the changed file to $rootdir/bldsva/cosmo5_1/tsmp
    comment "    copy the cosmo changed file for coupling to $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp "
      cp $cosdir/src/src_slow_tendencies_rk.f90 $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/src_runge_kutta.f90 $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/src_gridpoints.f90  $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/src_allocation.f90  $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/organize_physics.f90 $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/environment.f90 $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/src_twomom_sb.f90 $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/src_setup_vartab.f90 $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/src_artifdata.f90 $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/lmorg.f90 $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/data_fields.f90 $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
      cp $cosdir/src/src_radiation.f90 $rootdir/bldsva/intf_oas3/cosmo5_1/tsmp
    check
  fi

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
route "${cyellow}<<< c_substitutions_cos${cnormal}"
}

c_setup_cos(){
route "${cyellow}>>> c_setup_cos${cnormal}"

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
comment "  sed init_m_bldsva to namelist"  
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
#sed "s/__nhour_restart_start__/$(($cnt+$runhours))/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
sed "s/__nhour_restart_start__/$cnt/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check
sed "s/__nhour_restart_stop__/$(($cnt+$runhours))/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check
sed "s/__nhour_restart_incr__/1/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check

cnts=$(( ( $(date -u '+%s' -d "${startDate}") - $(date -u '+%s' -d "${initDate}")) / ${dt_cos} ))
comment "  sed output interval to namelist"
sed "s/__ncomb_start__/$cnts/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check
sed "s/__dump_cos_interval__/ $(python -c "print ($dump_cos*(3600/$dt_cos))")/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
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
  cp $rundir/INPUT_ORG $rundir/INPUT_ORG_$(printf "%05d" $(($instance-$startInst)))
fi

route "${cyellow}<<< c_setup_cos${cnormal}"
}

############################ 
#OASIS interface methods
############################

c_configure_oas(){
route "${cyellow}>>> c_configure_oas${cnormal}"
  comment "    sed oasis rootdir to Makefile"
    sed -i "s@__oasisroot__@$oasdir@" $file >> $log_file 2>> $err_file
  check
  comment "    sed platform to Makefile"
    sed -i "s@__platform__@$platform@" $file >> $log_file 2>> $err_file
  check
  comment "    make clean oasis"
    make -f $oasdir/util/make_dir/TopMakefileOasis3 realclean >> $log_file 2>> $err_file
  check
route "${cyellow}<<< c_configure_oas${cnormal}"
}

c_make_oas(){
route "${cyellow}>>> c_make_oas${cnormal}"
  comment "    make oasis"
    export SCOREP_WRAPPER=on
    make -j16 -f $oasdir/util/make_dir/TopMakefileOasis3 oasis3_psmile >> $log_file 2>> $err_file
  check
#DA
  if [[ $withPDAF == "true" ]]; then
  comment "    cp oas libs to $bindir/libs" 
    cp $libpsmile $bindir/libs >> $log_file 2>> $err_file
  check
  fi
route "${cyellow}<<< c_make_oas${cnormal}"
}


c_substitutions_oas(){
route "${cyellow}>>> c_substitutions_oas${cnormal}"
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
       patch "$rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[0]}/mod_oasis*"  ${oasdir}/lib/psmile/src
     check
  fi
route "${cyellow}<<< c_substitutions_oas${cnormal}"
}

c_setup_oas(){
route "${cyellow}>>> c_setup_oas${cnormal}"

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
  if [[ $withICON == "true" ]]; then
  rtime=$(( ($runhours*3600 + $cplfreq1/10)*10 ))  # with icon in tenths of second
  else
  rtime=$(($runhours*3600 + $cplfreq1))
  fi
  if [[ $withCESM == "true" ]] ; then ; rtime=$(($rtime+$cplfreq1)) ; fi
  comment "   sed sim time into namcouple"
    sed "s/totalruntime/$rtime/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  comment "   sed startdate into namcouple"
    sed "s/yyyymmdd/${yyyy}${mm}${dd}/" -i $rundir/namcouple  >> $log_file 2>> $err_file
  check


route "${cyellow}<<< c_setup_oas${cnormal}"
}


############################ 
#CLM interface methods
############################


c_configure_clm(){
route "${cyellow}>>> c_configure_clm${cnormal}"
  comment "    clean clm by removing build dir"
    rm -rf $clmdir/build >> $log_file 2>> $err_file
  check
  comment "    create new build dir"
    mkdir -p $clmdir/build >> $log_file 2>> $err_file
  check
  if [[ ${mList[1]} == clm3_5 ]] ; then
  comment "    copy oas_clm_init.F90 to  $clmdir/src/oas3"
    cp $rootdir/bldsva/intf_oas3/${mList[1]}/oas3/oas_clm_init.F90 $clmdir/src/oas3
  fi

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

    # cplInc+="-g -traceback -heap-arrays " # Mukund

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
  cppdef=""			# add "-DWATSAT3D" for input of "watsat3d" in "iniTimeConst.F90"
  if [ $cplscheme == "true" ] && [ $withICON == "false" ] ; then ; cppdef+=" -DCPL_SCHEME_F " ; fi
  comment "    configure clm"
  comment "    $clmdir/bld/configure -fc $cfc -cc $ccc $flags -fflags $cplInc -ldflags $cplLib -fopt $optComp -cppdefs $cppdef"
    export SCOREP_WRAPPER=off
    $clmdir/bld/configure -fc "$cfc" -cc "$ccc" $flags -fflags "$cplInc" -ldflags "$cplLib" -fopt "$optComp" -cppdefs "$cppdef"  >> $log_file 2>> $err_file
  check
route "${cyellow}<<< c_configure_clm${cnormal}"
}

c_make_clm(){
route "${cyellow}>>> c_make_clm${cnormal}"
  comment "    cd to clm build"
    cd $clmdir/build >> $log_file 2>> $err_file
  check
  comment "    make clm"
    export SCOREP_WRAPPER=on
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
route "${cyellow}<<< c_make_clm${cnormal}"
}


c_substitutions_clm(){
route "${cyellow}>>> c_substitutions_clm${cnormal}"
# CLM 3.5 substitutions
  comment "    create oas3 dir in $clmdir/src"
    mkdir -p $clmdir/src/oas3 >> $log_file 2>> $err_file
  check
  comment "    copy oas3 interface to clm/src "
    patch $rootdir/bldsva/intf_oas3/${mList[1]}/mct $clmdir/src
  check
  comment "    replace hydrology. Add files to clm/bld/usr.src "
    patch "$rootdir/bldsva/intf_oas3/${mList[1]}/tsmp/*" $clmdir/bld/usr.src 
  check
#DA
  if [[ $withPDAF == "true" ]] ; then
  comment "    copy PDAF fix to ${mList[1]}/bld/usr.src "
    patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[1]}/clmtype.F90 $clmdir/bld/usr.src
  check
    patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[1]}/clmtypeInitMod.F90 $clmdir/bld/usr.src
  check
    patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[1]}/iniTimeConst.F90 $clmdir/bld/usr.src	
  check
  fi	
route "${cyellow}<<< c_substitutions_clm${cnormal}"
}


c_setup_clm(){
route "${cyellow}>>> c_setup_clm${cnormal}"
if [[ ${mList[1]} == "clm5_0" ]]; then
 comment "  CLM5.0 setup"

 comment "  sed rundir to namelist"
   sed "s,__rundir__,$rundir," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
 check

 comment "  sed dt to namelist"
   sed "s,__dt_clm_bldsva__,$dt_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
 check

 comment "  sed forcingdir to namelist"
   sed "s,__forcingdir__,$forcingdir_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
 check
 comment "  sed dump interval namelist"
   sed "s,__dump_clm_interval__,$dump_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
 check
 comment "  sed runtime to namelist"
   if [[ $withICON == "true" ]]; then
     runstep_clm=$((($runhours*3600 + $cplfreq1/10)/$dt_clm))
   else
     runstep_clm=$((($runhours*3600 + $cplfreq1)/$dt_clm))
   fi
   sed "s,__runstep_clm_bldsva__,$runstep_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
 check
 comment "  add axe rights to clm namelist"
     chmod 755 $rundir/lnd.stdin >> $log_file 2>> $err_file
 check
 comment "  run clm namelist"
     $rundir/lnd.stdin >> $log_file 2>> $err_file
 check
# Below default clm 3.5 setup
else

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
  if [[ $withICON == "true" ]]; then
    runstep_clm=$((($runhours*3600 + $cplfreq1/10)/$dt_clm))
  else
    runstep_clm=$((($runhours*3600 + $cplfreq1)/$dt_clm))
  fi
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

fi 
route "${cyellow}<<< c_setup_clm${cnormal}"
}

############################ 
# eCLM interface methods
############################


c_configure_eclm(){
route "${cyellow}>>> c_configure_eclm${cnormal}"
  comment "    Using land component model eCLM (experimental) \n"
  comment "    Checking if eCLM repo is valid"
  cd ${clmdir}
  git status >> $log_file 2>> $err_file
  check

  if [[ -z $ECLM_CC || "$ECLM_CC" == " " ]]; then
    ECLM_CC=mpicc
  fi
  if [[ -z $ECLM_FC || "$ECLM_FC" == " " ]]; then
    ECLM_FC=mpifort
  fi

  if [[ $withOASMCT == "true" ]]; then
    ECLM_CMAKE_VARS+=" -DCMAKE_PREFIX_PATH="$oasdir/$platform""
  fi

  comment "    Running CMake configure step..."
  ECLM_BUILD_DIR="$clmdir/build"
  cmake -S src -B "$ECLM_BUILD_DIR" \
    -DCMAKE_INSTALL_PREFIX="$bindir" \
    -DCMAKE_C_COMPILER=$ECLM_CC \
    -DCMAKE_Fortran_COMPILER=$ECLM_FC \
     $ECLM_CMAKE_VARS >> $log_file 2>> $err_file
  check

route "${cyellow}<<< c_configure_eclm${cnormal}"
}

c_make_eclm(){
route "${cyellow}>>> c_make_eclm${cnormal}"
  comment "    Building eCLM (this will take approximately 30 mins)..."
  timer_start=$(date +%s)
  cmake --build "$ECLM_BUILD_DIR" >> $log_file 2>> $err_file
  check
  timer_end=$(date +%s)
  comment "    Build duration: $(date -u -d "0 $timer_end sec - $timer_start sec" +"%H:%M:%S")\n"
  comment "    Installing eCLM"
  cmake --install "$ECLM_BUILD_DIR" >> $log_file 2>> $err_file
  check
  comment "    Installing clm5nl-gen"
  pip3 install --user $clmdir/namelist_generator >> $log_file 2>> $err_file
  check
route "${cyellow}<<< c_make_eclm${cnormal}"
}

c_substitutions_eclm(){
route "${cyellow}>>> c_substitutions_eclm${cnormal}"
route "${cyellow}<<< c_substitutions_eclm${cnormal}"
}

c_setup_eclm(){
  route "${cyellow}>>> c_setup_eclm${cnormal}"
  route "${cyellow}<<< c_setup_eclm${cnormal}"
}

############################ 
#Parflow interface methods
############################


c_configure_pfl(){


route "${cyellow}>>> c_configure_pfl${cnormal}"

  if [[ ${mList[3]} == parflow ]] ; then

  comment "    cd to pfl build directory "
  cd $PARFLOW_BLD >> $log_file 2>> $err_file
  check
  export CC=$pcc
  export FC=$pfc
  export F77=$pf77
  export CXX=$pcxx

  comment "    configure pfsimulator and pftools"
  export SCOREP_WRAPPER=off
  cmake ../ $flagsSim >> $log_file 2>> $err_file
  check

  fi

  if [[ ${mList[3]} == parflow3_2 || ${mList[3]} == parflow3_0 ]] ; then
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

  comment "    cd to pfsimulator"
    cd $pfldir/pfsimulator >> $log_file 2>> $err_file
  check

    if [[ -e "$pfldir/pfsimulator/Makefile" ]] ; then
      comment "    make pfsimulator very clean"
        make -f $pfldir/pfsimulator/Makefile veryclean >> $log_file 2>> $err_file
      check
    fi 

  comment "    configure pfsimulator"
    export SCOREP_WRAPPER=off
    if [[ $withICON == "true" ]]; then
    $pfldir/pfsimulator/configure CC="$pcc" FC="$pfc" F77="$pf77" CXX="$pcxx" $flagsSim --enable-opt="$optComp" FCFLAGS="$fcflagsSim -DCOUP_OAS_ICON" CFLAGS="$cflagsSim -DCOUP_OAS_ICON" >> $log_file 2>> $err_file
    else
    $pfldir/pfsimulator/configure CC="$pcc" FC="$pfc" F77="$pf77" CXX="$pcxx" $flagsSim --enable-opt="$optComp" FCFLAGS="$fcflagsSim" CFLAGS="$cflagsSim" >> $log_file 2>> $err_file
    fi
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

  fi

route "${cyellow}<<< c_configure_pfl${cnormal}"
}

c_make_pfl(){
route "${cyellow}>>> c_make_pfl${cnormal}"

  if [[ ${mList[3]} == parflow ]] ; then
comment "    cd to pfl build directory "
  cd $PARFLOW_BLD >> $log_file 2>> $err_file
check
if [[ $profiling == "scalasca" ]]; then
  comment "    fix link.txt files for scalasca"
    export cpp_compiler=$(echo `which mpicc` | sed 's_/_\\/_g')
    find ${PARFLOW_BLD} -name 'link.txt' -exec sed -i "s/${cpp_compiler}/scorep-mpicc/g" {} \;
    export cpp_compiler=$(echo `which mpic++` | sed 's_/_\\/_g')
    find ${PARFLOW_BLD} -name 'link.txt' -exec sed -i "s/${cpp_compiler}/scorep-mpicxx/g" {} \;
    export cpp_compiler=$(echo `which g++` | sed 's_/_\\/_g')
    find ${PARFLOW_BLD} -name 'link.txt' -exec sed -i "s/${cpp_compiler}/scorep-mpicxx/g" {} \;
  check
  comment "    make pfsimulator and pftools"
    SCOREP_WRAPPER=off make pftools >> $log_file 2>> $err_file
    SCOREP_WRAPPER=on make -j8 >> $log_file 2>> $err_file
  check
else
  comment "    make pfsimulator and pftools"
    make -j8 >> $log_file 2>> $err_file
  check
fi
comment "    make install pfsimulator and pftools"
  make install >> $log_file 2>> $err_file
check
comment "    cp pfl bin to $bindir"
  cp -R $pfldir/bin/bin $bindir >> $log_file 2>> $err_file
check

comment "    cp binary to $bindir"
 cp $pfldir/bin/bin/parflow $bindir >> $log_file 2>> $err_file
 check

 if [[ $withPDAF == "true" ]]; then
    comment "    cp libs to $bindir/libs"
      cp $pfldir/bin/lib/* $bindir/libs >> $log_file 2>> $err_file
    check
    if [[ $processor == "GPU" ]]; then
      comment "    GPU: cp rmm libs to $bindir/libs"
        cp $pfldir/rmm/lib/* $bindir/libs >> $log_file 2>> $err_file
      check
    fi
 fi
  fi

  if [[ ${mList[3]} == parflow3_2 || ${mList[3]} == parflow3_0 ]] ; then
comment "    cd to pfsimulator" 
  cd $pfldir/pfsimulator >> $log_file 2>> $err_file
check
comment "    make pfsimulator"
  export SCOREP_WRAPPER=on
  make -f $pfldir/pfsimulator/Makefile >> $log_file 2>> $err_file
check
comment "    make install pfsimulator"
  make -f $pfldir/pfsimulator/Makefile install >> $log_file 2>> $err_file
check

  export SCOREP_WRAPPER=off
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
  export SCOREP_WRAPPER=on
  fi

route "${cyellow}<<< c_make_pfl${cnormal}"
}

c_substitutions_pfl(){
route "${cyellow}>>> c_substitutions_pfl${cnormal}"

  if [[ ${mList[3]} == parflow ]] ; then
    if [[ $withPDAF == "true" ]]; then

      comment "    sed DA amps into CMakeLists.txt"
        sed "s/PARFLOW_AMPS_LAYER PROPERTY STRINGS seq/PARFLOW_AMPS_LAYER PROPERTY STRINGS da seq/g" -i $pfldir/CMakeLists.txt >> $log_file 2>> $err_file
      check
      comment "    copy fix for PDAF into $pfldir"
        patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[3]}/parflow_proto.h $pfldir/pfsimulator/parflow_lib 
      check
        patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[3]}/solver_richards.c $pfldir/pfsimulator/parflow_lib 
      check
        patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[3]}/problem_saturation.c $pfldir/pfsimulator/parflow_lib
      check
        patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[3]}/problem_phase_rel_perm.c $pfldir/pfsimulator/parflow_lib
      check
        patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[3]}/da $pfldir/pfsimulator/amps
      check
    fi
  fi

  if [[ ${mList[3]} == parflow3_2 || ${mList[3]} == parflow3_0 ]] ; then
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
      patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[3]}/problem_saturation.c $pfldir/pfsimulator/parflow_lib
    check
      patch $rootdir/bldsva/intf_DA/pdaf1_1/tsmp/${mList[3]}/problem_phase_rel_perm.c $pfldir/pfsimulator/parflow_lib
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

  fi

route "${cyellow}<<< c_substitutions_pfl${cnormal}"
}

c_setup_pfl(){
route "${cyellow}>>> c_setup_pfl${cnormal}"

  if [ $numInst > 1 ] || [! -f "$rundir/coup_oas.tcl" ]; then
    comment "  $rundir/coup_oas.tcl does not exist (or for ensemble runs: numInst > 1), is copied, see c_setup_pfl()"
    cp $namelist_pfl $rundir/coup_oas.tcl >> $log_file 2>> $err_file
    check
  fi
 
  if [[ $withICON == "true" ]]; then
   comment "  sed ccs_ic_press file to pfl namelist."
      sed "s,__ccs_ic_press__,$defaultFDPFL/ccs_ic_press.pfb," -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
   check
  fi
 
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
#    sed "s/__stop_pfl_bldsva__/$runstep_clm/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
    sed "s/__stop_pfl_bldsva__/$(python -c "print (${runhours} + ${base_pfl})")/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed dump interval to pfl namelist."
    sed "s/__dump_pfl_interval__/$dump_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check

  comment "   sed timing base to pfl namelist."
    sed "s/__base_pfl__/$base_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check

  comment "   sed start counter to pfl namelist."
      cnt=$(( ($(date -u '+%s' -d "${startDate}") - $(date -u '+%s' -d "${initDate}"))))
      cnt=$(python -c "print ($cnt/($dump_pfl*3600.))")
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
  check
  tclsh $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check


route "${cyellow}<<< c_setup_pfl${cnormal}"
}


############################ 
#PDAF interface methods
############################


c_substitutions_pdaf(){
route "${cyellow}>>> c_substitutions_pdaf${cnormal}"

  comment "   mkdir  $dadir/interface"
    mkdir -p $dadir/interface  >> $log_file 2>> $err_file
  check

  comment "   cp pdaf interface model to $dadir/interface"
    patch $rootdir/bldsva/intf_DA/pdaf1_1/model $dadir/interface
  check

  comment "   cp pdaf interface framework to $dadir/interface"
    patch $rootdir/bldsva/intf_DA/pdaf1_1/framework $dadir/interface
  check

  comment "   mkdir $dadir/lib"
    mkdir -p $dadir/lib >> $log_file 2>> $err_file
  check

route "${cyellow}<<< c_substitutions_pdaf${cnormal}"
}

c_configure_pdaf_arch(){
route "${cyellow}>>> c_configure_pdaf_arch${cnormal}"

#PDAF arch part
  file=$dadir/make.arch/${PDAF_ARCH}.h

  comment "   cp pdaf config to $dadir"
    cp $rootdir/bldsva/intf_DA/pdaf1_1/arch/$platform/config/${PDAF_ARCH}.h $file >> $log_file 2>> $err_file
  check

  comment "   sed comFC dir to $file"
  sed -i "s@__comFC__@${comFC}@" $file >> $log_file 2>> $err_file
  check

  comment "   sed comCC dir to $file"
  sed -i "s@__comCC__@${comCC}@" $file >> $log_file 2>> $err_file
  check

  comment "   sed MPI dir to $file"
    sed -i "s@__MPI_INC__@-I${mpiPath}/include@" $file >> $log_file 2>> $err_file
  check

  comment "   sed LIBS to $file"
    sed -i "s@__LIBS__@${libs_src}@" $file >> $log_file 2>> $err_file
  check

  comment "   sed optimizations to $file"
    sed -i "s@__OPT__@${optComp}@" $file >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/src"
    cd $dadir/src >> $log_file 2>> $err_file
  check
  comment "   make clean pdaf"
    make clean >> $log_file 2>> $err_file
  check

route "${cyellow}<<< c_configure_pdaf_arch${cnormal}"
}

c_configure_pdaf(){
route "${cyellow}>>> c_configure_pdaf${cnormal}"

#PDAF interface part
  file1=$dadir/interface/model/Makefile
  file2=$dadir/interface/framework/Makefile
  comment "   cp pdaf interface Makefiles to $dadir"
    cp $rootdir/bldsva/intf_DA/pdaf1_1/model/Makefile  $file1 >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_DA/pdaf1_1/framework/Makefile  $file2 >> $log_file 2>> $err_file
  check

  comment "   sed bindir to Makefiles"
    sed -i "s,__bindir__,$bindir," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed comp flags to Makefiles"
    sed -i "s,__fflags__,-cpp -I$dadir/interface/model -I$ncdfPath/include $importFlags," $file1 $file2 >> $log_file 2>> $err_file
  check
    sed -i "s,__ccflags__,-I$dadir/interface/model -I$ncdfPath/include $importFlags," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed preproc flags to Makefiles"
    sed -i "s,__cpp_defs__,$cppdefs," $file1 $file2 >> $log_file 2>> $err_file
  check
    sed -i "s,__fcpp_defs__,$cppdefs," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed libs to Makefiles"
    sed -i "s,__libs__,$libs," $file2 >> $log_file 2>> $err_file
  check
  comment "   sed obj to Makefiles"
    sed -i "s,__obj__,$obj," $file1 >> $log_file 2>> $err_file
  check
  comment "   sed -D prefix to Makefiles"
    sed -i "s,__pf__,$pf," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed clm directory to Makefiles"
    sed -i "s,__clmdir__,${mList[1]}," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed cosmo directory to Makefiles"
    sed -i "s,__cosdir__,${mList[2]}," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed parflow directory to Makefiles"
    sed -i "s,__pfldir__,${mList[3]}," $file1 $file2 >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/interface/model"
    cd $dadir/interface/model >> $log_file 2>> $err_file
  check
  comment "   make clean model"
    make clean >> $log_file 2>> $err_file
  check
  comment "   cd to $dadir/src/interface/framework"
    cd $dadir/interface/framework >> $log_file 2>> $err_file
  check
  comment "   make clean framework"
    make clean >> $log_file 2>> $err_file
  check


route "${cyellow}<<< c_configure_pdaf${cnormal}"
}

c_make_pdaf(){
route "${cyellow}>>> c_make_pdaf${cnormal}"

  comment "   cd to $dadir/src"
    cd $dadir/src >> $log_file 2>> $err_file
  check
  comment "   make pdaf"
    make >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/interface/model"
    cd $dadir/interface/model >> $log_file 2>> $err_file
  check
  comment "   make pdaf model"
    make >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/interface/framework"
    cd $dadir/interface/framework >> $log_file 2>> $err_file
  check
  comment "   make pdaf framework"
    make >> $log_file 2>> $err_file
  check

route "${cyellow}<<< c_make_pdaf${cnormal}"
}

c_setup_pdaf(){
route "${cyellow}>>> c_setup_pdaf${cnormal}"
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
    sed "s/__endtime__/$(python -c "print (${runhours} + ${base_pfl})")/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed clmproc into pdaf namelist."
    sed "s/__clmproc__/$nproc_clm/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed cosproc into pdaf namelist."
    sed "s/__cosproc__/$nproc_cos/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check 
  comment "   sed dtmult into pdaf namelist."
    sed "s/__dtmult__/$(python -c "print (${dt_pfl} * 3600 / ${dt_cos})")/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check 

route "${cyellow}<<< c_setup_pdaf${cnormal}"
}

c_setup_rst(){

 comment " copy $restart_script to $rundir"
   cp $restart_script $rundir >> $log_file 2>> $err_file
 check

 comment "   sed startDate into restart template."
    sed 's/__startDate_bldsva__/"'"$startDate"'"/' -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
  check

 comment "   sed initDate into restart template."
    sed 's/__initDate_bldsva__/"'"$initDate"'"/' -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
  check

 comment "   sed dt_clm into restart template."
    sed "s/__dt_clm_bldsva__/$dt_clm/" -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
  check

 comment "   sed dt_cosmo into restart template."
    sed "s/__dt_cos_bldsva__/$dt_cos/" -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
  check

 comment "   sed PARFLOW_DIR into restart template $bindir."
#    sed "/__PARFLOW_DIR__/ \$bindir" -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
    sed -i "s|__PARFLOW_DIR__|$bindir|" $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
#    sed "s/__PARFLOW_DIR__/$bindir/" -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
  check

}

