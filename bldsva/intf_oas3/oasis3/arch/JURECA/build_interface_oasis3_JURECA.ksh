#! /bin/ksh
#

always_oas(){
print "${cblue}>> always_oas${cnormal}"
  liboas="$oasdir/$platform/lib/oasis3/liboasis3.MPI1.a"
  libpsmile="$oasdir/$platform/lib/libanaisg.a $oasdir/$platform/lib/libanaism.a $oasdir/$platform/lib/libclim.MPI1.a $oasdir/$platform/lib/libpsmile.MPI1.a $oasdir/$platform/lib/libfscint.a $oasdir/$platform/lib/libmpp_io.a $oasdir/$platform/lib/libscrip.a $oasdir/$platform/lib/libdownscal.a"
  incpsmile="-I$oasdir/$platform/build/lib/psmile.MPI1 -I$oasdir/$platform/build/lib/clim.MPI1 -I$oasdir/$platform/build/lib/mpp_io"
print "${cblue}<< always_oas${cnormal}"
}

substitutions_oas(){
print "${cblue}>> substitutions_oas${cnormal}"
    c_substitutions_oas    
print "${cblue}<< substitutions_oas${cnormal}"
}

configure_oas(){
print "${cblue}>> configure_oas${cnormal}"
  file=${oasdir}/util/make_dir/make.oas3
  print -n "   cp jureca oasis3 makefile to /util/make_dir/"
    cp $rootdir/bldsva/intf_oas3/oasis3/arch/$platform/config/make.intel_jureca_oa3 $file >> $log_file 2>> $err_file
  check
  c_configure_oas
  print -n "   sed new psmile includes to Makefile"
    sed -i 's@__inc__@-I$(LIBBUILD)/psmile.$(CHAN) -I$(LIBBUILD)/clim.$(CHAN) -I$(LIBBUILD)/mpp_io'" -I$ncdfPath/include@" $file >> $log_file 2>> $err_file
  check
  print -n "   sed ldflg to oas Makefile"
    sed -i "s@__ldflg__@@" $file >> $log_file 2>> $err_file
  check
  print -n "   sed comF90 to oas Makefile"
    sed -i "s@__comF90__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
  check
  print -n "   sed comCC to oas Makefile"
    sed -i "s@__comCC__@$mpiPath/bin/mpicc@" $file >> $log_file 2>> $err_file
  check
  print -n "   sed ld to oas Makefile"
    sed -i "s@__ld__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
  check
  print -n "   sed libs to oas Makefile"
    sed -i "s@__lib__@-L$ncdfPath/lib/ -lnetcdff@" $file >> $log_file 2>> $err_file
  check
  print -n "   sed precision to oas Makefile"
    sed -i "s@__precision__@-i4 -r8@" $file >> $log_file 2>> $err_file
  check

print "${cblue}<< configure_oas${cnormal}"
}

make_oas(){
print "${cblue}>> make_oas${cnormal}"
    c_make_oas
  print -n "    cp oasis binary to $bindir"
    cp $oasdir/$platform/bin/oasis3.MPI1.x $bindir >> $log_file 2>> $err_file
  check
print "${cblue}<< make_oas${cnormal}"
}


setup_oas(){
print "${cblue}>> setupOas${cnormal}"
  print -n "   copy cf_name_table to rundir"
    cp $rootdir/bldsva/data_oas3/cf_name_table.txt $rundir
  check
  print -n "   copy oas namelist to rundir"
    cp $namelist_oas $rundir/namcouple
  check
  print -n "   sed procs, gridsize & coupling freq into namcouple"
  if [[ $withPFL == "true" && $withCOS == "true" ]] then
    ncpl_exe1=$nproc_cos
    ncpl_exe2=$nproc_pfl
    ncpl_exe3=1

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
    ncpl_exe1=$nproc_pfl
    ncpl_exe2=$nproc_clm

    sed "s/nproc_exe1/$nproc_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe2/" -i $rundir/namcouple >> $log_file 2>> $err_file
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
    ncpl_exe1=$nproc_cos
    ncpl_exe2=$nproc_clm
    sed "s/nproc_exe1/$nproc_cos/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe2/" -i $rundir/namcouple >> $log_file 2>> $err_file
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
  print -n "   sed sim time into namcouple"
    sed "s/totalruntime/$(($runhours*3600))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  print -n "   sed startdate into namcouple"
    sed "s/yyyymmdd/${yyyy}${mm}${dd}/" -i $rundir/namcouple  >> $log_file 2>> $err_file
  check


print "${cblue}<< setupOas${cnormal}"
}

