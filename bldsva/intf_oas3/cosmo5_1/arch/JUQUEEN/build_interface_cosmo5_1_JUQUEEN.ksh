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
    cplFlag="-WF,-DCOUP_OAS_COS " 
  fi
  if [[ $cplscheme == "true" ]] ; then ; cplFlag+=" -WF,-DCPL_SCHEME_F " ; fi
  if [[ $readCLM == "true" ]] ; then ; cplFlag+=" -DREADCLM " ; fi 
  file=$cosdir/Fopts  
comment "   sed comflg to cos Makefile"
  sed -i "s@__comflg__@$optComp -I$ncdfPath/include $cplInc -qextname=flush -qsuffix=cpp=f90 -WF,-DGRIBDWD -WF,-DNETCDF -WF,-D__COSMO__ $cplFlag -WF,-DHYMACS@" $file >> $log_file 2>> $err_file
check
comment "   sed ldflg to cos Makefile"
  sed -i "s@__ldflg__@@" $file >> $log_file 2>> $err_file
check
comment "   sed comF90 to cos Makefile"
  sed -i "s@__comF90__@$mpiPath/bin/mpixlf90_r@" $file >> $log_file 2>> $err_file
check
comment "   sed ld to cos Makefile"
  sed -i "s@__ld__@$mpiPath/bin/mpixlf90_r@" $file >> $log_file 2>> $err_file
check
comment "   sed libs to cos Makefile"
  sed -i "s@__lib__@$grib1Path/libgrib1.a $cplLib $ncdfPath/lib/libnetcdf.a@" $file >> $log_file 2>> $err_file
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
 comment "   currently a fixed receive_fld2clm.F90 is necessary"
   cp $rootdir/bldsva/intf_oas3/${mList[2]}/arch/$platform/src/receive_fld_2clm.F90 $cosdir/src/oas3 >> $log_file 2>> $err_file
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

  c_setup_cos

route "${cblue}<< setupCos${cnormal}" 
}
