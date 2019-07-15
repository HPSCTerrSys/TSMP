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
comment "    patch icon configure"
  sed -i "s@CDFROOT/lib@CDFROOT/lib64@" configure >> $log_file 2>> $err_file
  sed -i "s@CDFROOT)/lib@CDFROOT)/lib64@" configure >> $log_file 2>> $err_file
check
if [[ $profiling == "scalasca" ]]; then
  export SCOREP_WRAPPER=off
  CC=scorep-mpicc F90=scorep-mpif90 F77=scorep-mpif77 NETCDFROOT=$EBROOTNETCDF NETCDFFROOT=$EBROOTNETCDFMINFORTRAN ./configure --with-mpi --disable-ocean --disable-jsbach --without-yac --with-grib-api=/p/project/cslts/brdar1/sw/grib_api-1.25.0 >> $log_file 2>> $err_file
else
  NETCDFROOT=$EBROOTNETCDF NETCDFFROOT=$EBROOTNETCDFMINFORTRAN ./configure --with-mpi --disable-ocean --disable-jsbach --without-yac --with-grib-api=/p/project/cslts/brdar1/sw/grib_api-1.25.0 >> $log_file 2>> $err_file
fi
comment "   cp Makefile to icon dir"
   cp $rootdir/bldsva/intf_oas3/${mList[2]}/arch/$platform/config/Makefile $icondir >> $log_file 2>> $err_file
check
comment "   sed oasisdir to icon Makefile"
  sed -i "s@__oasisdir__@$oasdir@" $icondir/Makefile >> $log_file 2>> $err_file
check
  c_configure_icon
  if [[ $withOAS == "true" ]]; then
    cplFlag="-DCOUP_OAS_ICON " 
  fi
  file=$icondir/Makefile
comment "   sed cplFlag to icon Makefile"
  sed -i "s@__cplflg__@$cplFlag@" $file >> $log_file 2>> $err_file
check
if [[ $profiling == "scalasca" ]]; then
  comment "   sed compilers + profiling to icon Makefile"
    sed -i "s@__comCC__@scorep-mpicc@" $file >> $log_file 2>> $err_file
  check
    sed -i "s@__comF77__@scorep-mpif77@" $file >> $log_file 2>> $err_file
  check
    sed -i "s@__comF90__@scorep-mpif90@" $file >> $log_file 2>> $err_file
  check
else
  comment "   sed compilers to icon Makefile"
    sed -i "s@__comCC__@$mpiPath/bin/mpicc@" $file >> $log_file 2>> $err_file
  check
    sed -i "s@__comF77__@$mpiPath/bin/mpif77@" $file >> $log_file 2>> $err_file
  check
    sed -i "s@__comF90__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
  check
fi
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
