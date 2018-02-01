#! /bin/ksh

always_icon(){
route "${cblue}>> always_icon${cnormal}"
route "${cblue}<< always_icon${cnormal}"
}

configure_icon(){
route "${cblue}>> configure_icon${cnormal}"
comment "    cd to icon dir"
  cd $icondir >> $log_file 2>> $err_file
check
comment "    fix JURECA specific library path"
git checkout -- configure
sed -i.bak s@NETCDFROOT/lib@NETCDFROOT/lib64@g configure
check
NETCDFROOT=$EBROOTNETCDF NETCDFFROOT=$EBROOTNETCDFMINFORTRAN ./configure --with-fortran=intel --with-mpi --disable-ocean --disable-jsbach --without-yac --with-grib-api=/homea/slts/slts23/sw/grib_api-1.25.0 >> $log_file 2>> $err_file
comment "   cp Makefile to icon dir"
   cp $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/config/Makefile $icondir >> $log_file 2>> $err_file
check
route "cp $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/config/Makefile $icondir"
check
comment "   sed oasisdir to icon Makefile"
  sed -i "s@__oasisdir__@$oasdir@" $icondir/Makefile >> $log_file 2>> $err_file
check
  c_configure_icon
  if [[ $withOAS == "true" ]]; then
    cplFlag="-DCOUP_OAS_ICON " 
  fi
  file=$icondir/Makefile
echo "SLAVKO: ---- $file"
comment "   sed cplFlag to icon Makefile"
  sed -i "s@__cplflg__@$cplFlag@" $file >> $log_file 2>> $err_file
check
comment "   sed cplLib to icon Makefile"
  sed -i "s@__cpllib__@$cplLib@" $file >> $log_file 2>> $err_file
check
comment "   sed cplInc to icon Makefile"
  sed -i "s@__cplinc__@$cplInc@" $file >> $log_file 2>> $err_file
check
route "${cblue}<< configure_icon${cnormal}"
}

make_icon(){
route "${cblue}>> make_icon${cnormal}"
  c_make_icon
route "${cblue}<< make_icon${cnormal}"
}


substitutions_icon(){
route "${cblue}>> substitutions_icon${cnormal}"
 c_substitutions_icon
route "${cblue}<< substitutions_icon${cnormal}"
}

setup_icon(){
route "${cblue}>> setupIcon${cnormal}"

  c_setup_icon

route "${cblue}<< setupIcon${cnormal}" 
}
