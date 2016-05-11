#! /bin/ksh

always_cos(){
route "${cblue}>> always_cos${cnormal}"
route "${cblue}<< always_cos${cnormal}"
}

configure_cos(){
route "${cblue}>> configure_cos${cnormal}"
comment "   cp Makefile to cosmo dir"
   cp $rootdir/bldsva/intf_oas3/${mList[2]}/arch/$platform/config/Makefile $cosdir >> $log_file 2>> $err_file
check
   cp $rootdir/bldsva/intf_oas3/${mList[2]}/arch/$platform/config/Fopts $cosdir >> $log_file 2>> $err_file
check
  c_configure_cos
  if [[ $withOAS == "true" ]]; then
    cplFlag="-DCOUP_OAS_COS " 
  fi
  if [[ $cplscheme == "true" ]] ; then ; cplFlag+=" -DCPL_SCHEME_F " ; fi 
  file=$cosdir/Fopts 
comment "   sed comflg to cos Makefile"
  sed -i "s@__comflg__@$optComp -ffree-line-length-0 -I$ncdfPath/include $cplInc -cpp -DGRIBDWD -DNETCDF -D__COSMO__ $cplFlag -DHYMACS@" $file >> $log_file 2>> $err_file
check
comment "   sed ldflg to cos Makefile"
  sed -i "s@__ldflg__@@" $file >> $log_file 2>> $err_file
check
comment "   sed comF90 to cos Makefile"
  sed -i "s@__comF90__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
check
comment "   sed ld to cos Makefile"
  sed -i "s@__ld__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
check
comment "   sed libs to cos Makefile"
  sed -i "s@__lib__@$grib1Path/libgrib1.a $cplLib -L$ncdfPath/lib/ -lnetcdff -lnetcdf@" $file >> $log_file 2>> $err_file
check
route "${cblue}<< configure_cos${cnormal}"
}

make_cos(){
route "${cblue}>> make_cos${cnormal}"
  c_make_cos
route "${cblue}<< make_cos${cnormal}"
}


substitutions_cos(){
route "${cblue}>> substitutions_cos${cnormal}"
 c_substitutions_cos
   comment "   cp ObjFiles & ObjDependencies in $cosdir"
     cp $rootdir/bldsva/intf_oas3/${mList[2]}/arch/$platform/config/Obj* $cosdir >> $log_file 2>> $err_file
   check
   if [[ $withOASMCT == "true" ]] ; then
     comment "   sed replace old mod_prism includes from cos oas files"
       sed -i "s/mod_prism_proto/mod_prism/" $cosdir/src/oas3/oas_cos_vardef.F90 >> $log_file 2>> $err_file
     check
       sed -i "s/USE mod_prism.*//" $cosdir/src/oas3/oas_cos_rcv.F90 >> $log_file 2>> $err_file
     check
       sed -i "s/USE mod_prism.*//" $cosdir/src/oas3/oas_cos_snd.F90 >> $log_file 2>> $err_file
     check
       sed -i "s/USE mod_prism.*//" $cosdir/src/oas3/oas_cos_define.F90 >> $log_file 2>> $err_file
     check
   fi
route "${cblue}<< substitutions_cos${cnormal}"
}

setup_cos(){
route "${cblue}>> setupCos${cnormal}"
  nstop_cos=$((($runhours*3600-$cplfreq1)/$dt_cos))
comment "  cp namelist to rundir"
  cp ${namelist_cos}5_1 $rundir/lmrun_uc >> $log_file 2>> $err_file
check


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

comment "  sed gridpoints to namelist"
  sed "s,nbdl_cos_bldsva,$nbndlines," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check

check
comment "  sed forcingdir to namelist"
  sed "s,__forcingdir__,$forcingdir_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed rundir to namelist"
  sed "s,__rundir__,$rundir," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed stop time to namelist"
  sed "s/nstop_cos_bldsva/$nstop_cos/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed date to namelist"
  sed "s/init_y_bldsva/$yyyy/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_m_bldsva/$mm/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_d_bldsva/$dd/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_h_bldsva/$hh/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check


comment "  cd to rundir"
  cd $rundir >> $log_file 2>> $err_file
check
comment "  run lmrun_uc clean"
  $rundir/lmrun_uc cleancluma >> $log_file 2>> $err_file
check
comment "  run lmrun_uc exe"
  $rundir/lmrun_uc execluma >> $log_file 2>> $err_file
check
route "${cblue}<< setupCos${cnormal}" 
}
