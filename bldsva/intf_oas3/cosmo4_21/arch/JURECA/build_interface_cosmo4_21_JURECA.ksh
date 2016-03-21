#! /bin/ksh

always_cos(){
print "${cblue}>> always_cos${cnormal}"
print "${cblue}<< always_cos${cnormal}"
}

configure_cos(){
print "${cblue}>> configure_cos${cnormal}"
  c_configure_cos
  if [[ $withOAS == "true" ]]; then
    cplFlag="-DCOUP_OAS_COS" 
  fi  
print -n "   sed comflg to cos Makefile"
  sed -i "s@__comflg__@$optComp -I$ncdfPath/include $cplInc -cpp -DGRIBDWD -DNETCDF $cplFlag -DHYMACS@" $file >> $log_file 2>> $err_file
check
print -n "   sed ldflg to cos Makefile"
  sed -i "s@__ldflg__@@" $file >> $log_file 2>> $err_file
check
print -n "   sed comF90 to cos Makefile"
  sed -i "s@__comF90__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
check
print -n "   sed ld to cos Makefile"
  sed -i "s@__ld__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
check
print -n "   sed libs to cos Makefile"
  sed -i "s@__lib__@$grib1Path/libgrib1.a $cplLib -L$ncdfPath/lib/ -lnetcdff@" $file >> $log_file 2>> $err_file
check
print "${cblue}<< configure_cos${cnormal}"
}

make_cos(){
print "${cblue}>> make_cos${cnormal}"
  c_make_cos
print "${cblue}<< make_cos${cnormal}"
}


substitutions_cos(){
print "${cblue}>> substitutions_cos${cnormal}"
 c_substitutions_cos

 print -n "   currently a fixed receive_fld2clm.F90 is necessary"
   cp $rootdir/bldsva/intf_oas3/cosmo4_21/arch/JURECA/src/receive_fld_2clm.F90 $cosdir/src/oas3 >> $log_file 2>> $err_file
 check
 print -n "   cp Makefile to cosmo dir"
   cp $rootdir/bldsva/intf_oas3/cosmo4_21/arch/JURECA/config/Makefile $cosdir >> $log_file 2>> $err_file
 check
   if [[ $withOASMCT == "true" ]] ; then
     print -n "   sed replace old mod_prism includes from cos oas files"
       sed -i "s/mod_prism_proto/mod_prism/" $cosdir/src/oas3/oas_cos_vardef.F90 >> $log_file 2>> $err_file
     check
       sed -i "s/USE mod_prism.*//" $cosdir/src/oas3/oas_cos_rcv.F90 >> $log_file 2>> $err_file
     check
       sed -i "s/USE mod_prism.*//" $cosdir/src/oas3/oas_cos_snd.F90 >> $log_file 2>> $err_file
     check
       sed -i "s/USE mod_prism.*//" $cosdir/src/oas3/oas_cos_define.F90 >> $log_file 2>> $err_file
     check
   fi
print "${cblue}<< substitutions_cos${cnormal}"
}

setup_cos(){
print "${cblue}>> setupCos${cnormal}"
  nstop_cos=$((($runhours*3600-$cplfreq1)/$dt_cos))
print -n "  cp namelist to rundir"
  cp $namelist_cos $rundir/lmrun_uc >> $log_file 2>> $err_file
check


print -n "  sed dt to namelist"
  sed "s,dt_cos_bldsva,$dt_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
print -n "  sed number of procs to namelist"
  sed "s,nprocx_cos_bldsva,$px_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s,nprocy_cos_bldsva,$py_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
print -n "  sed gridpoints to namelist"
  sed "s,ie_tot_bldsva,$gx_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s,je_tot_bldsva,$gy_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file

print -n "  sed gridpoints to namelist"
  sed "s,nbdl_cos_bldsva,$nbndlines," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check

check
print -n "  sed forcingdir to namelist"
  sed "s,__forcingdir__,$forcingdir_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
print -n "  sed rundir to namelist"
  sed "s,__rundir__,$rundir," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
print -n "  sed stop time to namelist"
  sed "s/nstop_cos_bldsva/$nstop_cos/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
print -n "  sed date to namelist"
  sed "s/init_y_bldsva/$yyyy/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_m_bldsva/$mm/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_d_bldsva/$dd/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_h_bldsva/$hh/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check


print -n "  cd to rundir"
  cd $rundir >> $log_file 2>> $err_file
check
print -n "  run lmrun_uc clean"
  $rundir/lmrun_uc cleancluma >> $log_file 2>> $err_file
check
print -n "  run lmrun_uc exe"
  $rundir/lmrun_uc execluma >> $log_file 2>> $err_file
check
print "${cblue}<< setupCos${cnormal}" 
}
