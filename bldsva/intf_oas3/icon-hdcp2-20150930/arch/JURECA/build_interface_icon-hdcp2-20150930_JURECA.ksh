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
./configure --with-fortran=intel >> $log_file 2>> $err_file
comment "   cp Makefile to icon dir"
   cp $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/config/Makefile $icondir >> $log_file 2>> $err_file
route "cp $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/config/Makefile $icondir"
check
  c_configure_icon
  if [[ $withOAS == "true" ]]; then
    cplFlag="-DCOUP_OAS_ICON " 
  fi
  if [[ $cplscheme == "true" ]] ; then ; cplFlag+=" -DCPL_SCHEME_F " ; fi
  if [[ $readCLM == "true" ]] ; then ; cplFlag+=" -DREADCLM " ; fi 
  file=$icondir/Makefile
comment "   sed cplFlag to icon Makefile"
  sed -i "s@__cplflg__@$cplFlag@" $file >> $log_file 2>> $err_file
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
 #comment "   cp ObjFiles & ObjDependencies in $icondir"
 #  cp $rootdir/bldsva/intf_oas3/${mList[3]}/arch/$platform/config/Obj* $icondir >> $log_file 2>> $err_file
 #check
   #if [[ $withOASMCT == "true" ]] ; then
     #comment "   sed replace old mod_prism includes from icon oas files"
     #comment " == not implemented == "
       #sed -i "s/mod_prism_proto/mod_prism/" $icondir/src/oas3/oas_icon_vardef.F90 >> $log_file 2>> $err_file
     #check
       #sed -i "s/USE mod_prism.*//" $icondir/src/oas3/oas_icon_rcv.F90 >> $log_file 2>> $err_file
     #check
       #sed -i "s/USE mod_prism.*//" $icondir/src/oas3/oas_icon_snd.F90 >> $log_file 2>> $err_file
     #check
       #sed -i "s/USE mod_prism.*//" $icondir/src/oas3/oas_icon_define.F90 >> $log_file 2>> $err_file
     #check
   #fi
route "${cblue}<< substitutions_icon${cnormal}"
}

setup_icon(){
route "${cblue}>> setupIcon${cnormal}"

  c_setup_icon

route "${cblue}<< setupIcon${cnormal}" 
}
