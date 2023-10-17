#! /bin/ksh

always_cos(){
route "${cyellow}>> always_cos${cnormal}"
route "${cyellow}<< always_cos${cnormal}"
}

configure_cos(){
route "${cyellow}>> configure_cos${cnormal}"
comment "   cp Makefile to cosmo dir"
   cp $rootdir/bldsva/intf_oas3/${mList[2]}/arch/config/Makefile $cosdir >> $log_file 2>> $err_file
check
   cp $rootdir/bldsva/intf_oas3/${mList[2]}/arch/config/Fopts $cosdir >> $log_file 2>> $err_file
check
  c_configure_cos
  if [[ $withOAS == "true" ]]; then
    cplFlag="-DCOUP_OAS_COS " 
  fi
  if [[ $cplscheme == "true" ]] ; then ; cplFlag+=" -DCPL_SCHEME_F " ; fi
  if [[ $readCLM == "true" ]] ; then ; cplFlag+=" -DREADCLM " ; fi 
  file=$cosdir/Fopts 
comment "   sed ldflg to cos Makefile"
  sed -i "s@__ldflg__@@" $file >> $log_file 2>> $err_file
check
  if  echo "$compiler" | grep -qE 'Gnu' ; then
     comment "   sed comF90 based on Gnu to cos Makefile"
          if [[ $profiling == "scalasca" ]]; then
        sed -i "s@__comF90__@scorep-mpif90 -cpp -c -fallow-invalid-boz -fallow-argument-mismatch -ffree-line-length-0 -fstack-protector-all -finit-real=nan -finit-integer=-2147483648 -finit-character=127 -ffpe-trap=invalid,zero,overflow@" $file >> $log_file 2>> $err_file
     else
        sed -i "s@__comF90__@$profComp $mpiPath/bin/mpif90 -cpp -c -fallow-invalid-boz -fallow-argument-mismatch -ffree-line-length-0 -fstack-protector-all -finit-real=nan -finit-integer=-2147483648 -finit-character=127 -ffpe-trap=invalid,zero,overflow@" $file >> $log_file 2>> $err_file
     fi
     check
     comment "   sed comflg to cos Makefile"
     sed -i "s@__comflg__@$optComp -I$ncdfPath/include -I$gribPath/include $cplInc -cpp -DGRIBAPI -DNETCDF -D__COSMO__ $cplFlag -DHYMACS@" $file >> $log_file 2>> $err_file
     check
  else
     comment "   sed comF90 based on Intel compiler to cos Makefile"
    if [[ $profiling == "scalasca" ]]; then
        sed -i "s@__comF90__@scorep-mpif90  -fpp -O2 -fp-model source@" $file >> $log_file 2>> $err_file
     else
        sed -i "s@__comF90__@$profComp $mpiPath/bin/mpif90  -fpp -O2 -fp-model source@" $file >> $log_file 2>> $err_file
     fi
     check
     comment "   sed comflg to cos Makefile"
     sed -i "s@__comflg__@$optComp -I$ncdfPath/include -I$gribPath/include $cplInc -fpp -DGRIBAPI -DNETCDF -D__COSMO__ $cplFlag -DHYMACS@" $file >> $log_file 2>> $err_file
     check
  fi
comment "   sed ld to cos Makefile"
  if [[ $profiling == "scalasca" ]]; then
     sed -i "s@__ld__@scorep-mpif90@" $file >> $log_file 2>> $err_file
  else
     sed -i "s@__ld__@$profComp $mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
  fi
check
comment "   sed libs to cos Makefile"
  sed -i "s@__lib__@-L$gribPath/lib/ $cplLib -L$ncdfPath/lib/ -leccodes_f90 -leccodes -lnetcdff@" $file >> $log_file 2>> $err_file
check
route "${cyellow}<< configure_cos${cnormal}"
}

make_cos(){
route "${cyellow}>> make_cos${cnormal}"
  c_make_cos
route "${cyellow}<< make_cos${cnormal}"
}


substitutions_cos(){
route "${cyellow}>> substitutions_cos${cnormal}"
 c_substitutions_cos
 comment "   cp ObjFiles & ObjDependencies in $cosdir"
   patch "$rootdir/bldsva/intf_oas3/${mList[2]}/arch/config/Obj*" $cosdir 
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
route "${cyellow}<< substitutions_cos${cnormal}"
}


