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
    cp $rootdir/bldsva/intf_oas3/oasis3/arch/$platform/config/make.intel_jureca_oa3 $file >> $log_file 2>> $err_file
  check
  c_configure_oas
  comment "   sed new psmile includes to Makefile"
    sed -i 's@__inc__@-I$(LIBBUILD)/psmile.$(CHAN) -I$(LIBBUILD)/clim.$(CHAN) -I$(LIBBUILD)/mpp_io'" -I$ncdfPath/include@" $file >> $log_file 2>> $err_file
  check
  comment "   sed ldflg to oas Makefile"
    sed -i "s@__ldflg__@@" $file >> $log_file 2>> $err_file
  check
  comment "   sed comF90 to oas Makefile"
    sed -i "s@__comF90__@$mpiPath/bin/mpif90 $optComp@" $file >> $log_file 2>> $err_file
  check
  comment "   sed comCC to oas Makefile"
    sed -i "s@__comCC__@$mpiPath/bin/mpicc $optComp@" $file >> $log_file 2>> $err_file
  check
  comment "   sed ld to oas Makefile"
    sed -i "s@__ld__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
  check
  comment "   sed libs to oas Makefile"
    sed -i "s@__lib__@-L$ncdfPath/lib/ -lnetcdff@" $file >> $log_file 2>> $err_file
  check
  comment "   sed precision to oas Makefile"
    sed -i "s@__precision__@-i4 -r8@" $file >> $log_file 2>> $err_file
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

  c_setup_oas

route "${cblue}<< setupOas${cnormal}"
}

