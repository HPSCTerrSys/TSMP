#! /bin/ksh

always_icon(){
route "${cyellow}>> always_icon${cnormal}"
route "${cyellow}<< always_icon${cnormal}"
}

configure_icon(){
route "${cyellow}>> configure_icon${cnormal}"
comment "    cd to icon dir"
  cd $icondir >> $log_file 2>> $err_file
check
comment "    fix JURECA specific library path"
git checkout -- configure
sed -i.bak s@NETCDFROOT/lib@NETCDFROOT/lib64@g configure
check
comment "   sed libeccode to icon config file  "
  sed -i "s@libgrib_api.a@libeccodes.so@g" configure >> $log_file 2>> $err_file
check

if [[ $profiling == "scalasca" ]]; then
  export SCOREP_WRAPPER=off
  CC=scorep-mpicc F90=scorep-mpif90 F77=scorep-mpif77 NETCDFROOT=$EBROOTNETCDF NETCDFFROOT=$EBROOTNETCDFMINFORTRAN ./configure --with-fortran=intel --with-mpi --disable-ocean --disable-jsbach --without-yac --with-grib-api=/p/project/cslts/brdar1/sw/grib_api-1.25.0 >> $log_file 2>> $err_file
else
  NETCDFROOT=$EBROOTNETCDF NETCDFFROOT=$EBROOTNETCDFMINFORTRAN ./configure --with-fortran=intel --with-mpi --disable-ocean --disable-jsbach --without-yac --with-grib-api=$gribPath >> $log_file 2>> $err_file
fi
  export SCOREP_WRAPPER=on
comment "   cp Makefile to icon dir"
   cp $rootdir/bldsva/intf_oas3/${mList[2]}/arch/$platform/config/Makefile $icondir >> $log_file 2>> $err_file
check
route "cp $rootdir/bldsva/intf_oas3/${mList[2]}/arch/$platform/config/Makefile $icondir"
check
comment "   sed oasisdir to icon Makefile"
  sed -i "s@__oasisdir__@$oasdir@" $icondir/Makefile >> $log_file 2>> $err_file
check
  c_configure_icon
  if [[ $withOAS == "true" ]]; then
    cplFlag="-DCOUP_OAS_ICON " 
  fi
  file=$icondir/Makefile
comment "   sed netcdffroot to icon Makefile"
  sed -i "s@__netcdffroot__@$EBROOTNETCDFMINFORTRAN@" $file >> $log_file 2>> $err_file
check
comment "   sed netcdfroot to icon Makefile"
  sed -i "s@__netcdfroot__@$EBROOTNETCDF@" $file >> $log_file 2>> $err_file
check
comment "   sed sziproot to icon Makefile"
  sed -i "s@__sziproot__@$EBROOTSZIP@" $file >> $log_file 2>> $err_file
check
#
comment "   sed eccode root to icon Makefile"
  sed -i "s@__eccoderoot__@$EBROOTECCODES@" $file >> $log_file 2>> $err_file
check
comment "   sed zlibroot to icon Makefile"
  sed -i "s@__zlibroot__@$EBROOTZLIB@" $file >> $log_file 2>> $err_file
check
comment "   sed hdf5root to icon Makefile"
  sed -i "s@__hdf5root__@$EBROOTHDF5@" $file >> $log_file 2>> $err_file
check
comment "   sed cplFlag to icon Makefile"
  sed -i "s@__cplflg__@$cplFlag@" $file >> $log_file 2>> $err_file
check
comment "   sed cplLib to icon Makefile"
  sed -i "s@__cpllib__@$cplLib@" $file >> $log_file 2>> $err_file
check
comment "   sed cplInc to icon Makefile"
  sed -i "s@__cplinc__@$cplInc@" $file >> $log_file 2>> $err_file
check
comment "   sed comCC to icon Makefile"
  if [[ $profiling == "scalasca" ]]; then
    sed -i "s@__comCC__@scorep-mpicc@" $file >> $log_file 2>> $err_file
  else
    sed -i "s@__comCC__@$mpiPath/bin/mpicc@" $file >> $log_file 2>> $err_file
  fi
check
comment "   sed comF90 to icon Makefile"
  if [[ $profiling == "scalasca" ]]; then
    sed -i "s@__comF90__@scorep-mpif90@" $file >> $log_file 2>> $err_file
  else
    sed -i "s@__comF90__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
  fi
check
comment "   sed comF77 to icon Makefile"
  if [[ $profiling == "scalasca" ]]; then
    sed -i "s@__comF77__@scorep-mpif77@" $file >> $log_file 2>> $err_file
  else
    sed -i "s@__comF77__@$mpiPath/bin/mpif77@" $file >> $log_file 2>> $err_file
  fi
check
route "${cyellow}<< configure_icon${cnormal}"
}

make_icon(){
route "${cyellow}>> make_icon${cnormal}"
  c_make_icon
route "${cyellow}<< make_icon${cnormal}"
}


substitutions_icon(){
route "${cyellow}>> substitutions_icon${cnormal}"
 c_substitutions_icon
route "${cyellow}<< substitutions_icon${cnormal}"
}

setup_icon(){
route "${cyellow}>> setupIcon${cnormal}"

  c_setup_icon

route "${cyellow}<< setupIcon${cnormal}" 
}
