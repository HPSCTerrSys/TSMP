#! /bin/ksh
#

always_oas(){
route "${cblue}>> always_oas${cnormal}"
  liboas="$oasdir/$platform/lib/oasis3/liboasis3.MPI1.a"
  libpsmile="$oasdir/$platform/lib/libanaisg.a $oasdir/$platform/lib/libanaism.a $oasdir/$platform/lib/libclim.MPI1.a $oasdir/$platform/lib/libpsmile.MPI1.a $oasdir/$platform/lib/libfscint.a $oasdir/$platform/lib/libmpp_io.a $oasdir/$platform/lib/libscrip.a $oasdir/$platform/lib/libdownscal.a"
  incpsmile="-I$oasdir/$platform/build/lib/psmile.MPI1 -I$oasdir/$platform/build/lib/clim.MPI1 -I$oasdir/$platform/build/lib/mpp_io"
route "${cblue}<< always_oas${cnormal}"
}

substitutions_oas(){
route "${cblue}>> substitutions_oas${cnormal}"
    c_substitutions_oas    
route "${cblue}<< substitutions_oas${cnormal}"
}

configure_oas(){
route "${cblue}>> configure_oas${cnormal}"
  file=${oasdir}/util/make_dir/make.oas3
  comment "   cp jureca oasis3 makefile to /util/make_dir/"
    cp $rootdir/bldsva/intf_oas3/oasis3/arch/$platform/config/make.gnu_ecmwf_oa3 $file >> $log_file 2>> $err_file
  check
  c_configure_oas
  comment "   sed new psmile includes to Makefile"
    sed -i 's@__inc__@-I$(LIBBUILD)/psmile.$(CHAN) -I$(LIBBUILD)/clim.$(CHAN) -I$(LIBBUILD)/mpp_io'" -I$ncdfPath/include@" $file >> $log_file 2>> $err_file
  check
  comment "   sed ldflg to oas Makefile"
    sed -i "s@__ldflg__@@" $file >> $log_file 2>> $err_file
  check
  comment "   sed comF90 to oas Makefile"
    sed -i "s@__comF90__@ftn $optComp@" $file >> $log_file 2>> $err_file
  check
  comment "   sed comCC to oas Makefile"
    sed -i "s@__comCC__@cc $optComp@" $file >> $log_file 2>> $err_file
  check
  comment "   sed ld to oas Makefile"
    sed -i "s@__ld__@ftn@" $file >> $log_file 2>> $err_file
  check
  comment "   sed libs to oas Makefile"
    sed -i "s@__lib__@-L$ncdfPath/lib/ -lnetcdff@" $file >> $log_file 2>> $err_file
  check
  comment "   sed precision to oas Makefile"
    sed -i "s@__precision__@-fdefault-real-8@" $file >> $log_file 2>> $err_file
  check

route "${cblue}<< configure_oas${cnormal}"
}

make_oas(){
route "${cblue}>> make_oas${cnormal}"
    c_make_oas
  comment "    cp oasis binary to $bindir"
    cp $oasdir/$platform/bin/oasis3.MPI1.x $bindir >> $log_file 2>> $err_file
  check
route "${cblue}<< make_oas${cnormal}"
}


setup_oas(){
route "${cblue}>> setupOas${cnormal}"
  comment "   copy cf_name_table to rundir"
    cp $rootdir/bldsva/data_oas3/cf_name_table.txt $rundir
  check
  comment "   copy oas namelist to rundir"
    cp $namelist_oas $rundir/namcouple
  check
  comment "   sed procs, gridsize & coupling freq into namcouple"
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
    ncpl_exe2=1

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
    ncpl_exe2=1
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
  comment "   sed sim time into namcouple"
    sed "s/totalruntime/$(($runhours*3600))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  comment "   sed startdate into namcouple"
    sed "s/yyyymmdd/${yyyy}${mm}${dd}/" -i $rundir/namcouple  >> $log_file 2>> $err_file
  check


route "${cblue}<< setupOas${cnormal}"
}

