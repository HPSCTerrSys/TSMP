#! /bin/ksh
#

always_oas(){
route "${cblue}>> always_oas${cnormal}"
  liboas=""
  libpsmile="$oasdir/$platform/lib/libpsmile.MPI1.a $oasdir/$platform/lib/libmct.a $oasdir/$platform/lib/libmpeu.a $oasdir/$platform/lib/libscrip.a"
  incpsmile="-I$oasdir/$platform/build/lib/psmile.MPI1"
route "${cblue}<< always_oas${cnormal}"
}

substitutions_oas(){
route "${cblue}>> substitutions_oas${cnormal}"
  comment "    cp TopMakefileOasis3 to make_dir"
    patch $rootdir/bldsva/intf_oas3/oasis3-mct/tsmp/TopMakefileOasis3 $oasdir/util/make_dir/
  check
  comment "   cp new  mod_oasis_method.F90 mod_oasis_grid.F90 to psmile/src"
    patch "$rootdir/bldsva/intf_oas3/oasis3-mct/tsmp/mod_oasis_*" ${oasdir}/lib/psmile/src 
  check
    c_substitutions_oas
  comment "   sed prism_get_freq functionality to mod_prism.F90"
    sed -i "/oasis_get_debug/a   use mod_oasis_method ,only: prism_get_freq            => oasis_get_freq" ${oasdir}/lib/psmile/src/mod_prism.F90  >> $log_file 2>> $err_file   # critical anchor
  check
    # set prefix for all mct files to run with other mct
    prefix="oas_"
  comment "   cd to ${oasdir}/lib/mct/mct"
    cd ${oasdir}/lib/mct/mct >> $log_file 2>> $err_file
  check
  comment "   mv all m_* to ${prefix}m_*"
    for i in m_* ; do mv $i ${prefix}${i}  >> $log_file 2>> $err_file ; check ; done
  comment "   mv mct_mod.F90 to ${prefix}mct_mod.F90"
    mv mct_mod.F90 ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
  comment "   rename all interfaces in mct_mod.F90 and Makefile"
    sed -i "s/\(${prefix}\)*mct_/${prefix}mct_/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/\(${prefix}\)*mct_/${prefix}mct_/g" Makefile >> $log_file 2>> $err_file
  check
  comment "   rename all sources in mct_mod.F90 and Makefile"
    sed -i "s/use m_/use ${prefix}m_/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/\(${prefix}\)*m_/${prefix}m_/g" Makefile >> $log_file 2>> $err_file
  check
  comment "   except those which are redirecting mpeu sources"
    sed -i "s/${prefix}m_List/m_List/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/${prefix}m_string/m_string/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/${prefix}m_die/m_die/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/${prefix}m_MergeSort/m_MergeSort/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/${prefix}m_inpak90/m_inpak90/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/${prefix}m_Permuter/m_Permuter/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
  comment "   rename all module names"
    sed -i "s/module m_/module ${prefix}m_/g" * >> $log_file 2>> $err_file
  check
  comment "   rename all dependencies"
    sed -i "s/use m_MCTW/use ${prefix}m_MCTW/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Glob/use ${prefix}m_Glob/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_AttrV/use ${prefix}m_AttrV/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Accum/use ${prefix}m_Accum/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Gener/use ${prefix}m_Gener/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Spa/use ${prefix}m_Spa/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Nav/use ${prefix}m_Nav/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Con/use ${prefix}m_Con/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Ex/use ${prefix}m_Ex/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Rou/use ${prefix}m_Rou/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Rear/use ${prefix}m_Rear/g" * >> $log_file 2>> $err_file
  check
  comment "   cd in psmile-source dir"
    cd ${oasdir}/lib/psmile/src >> $log_file 2>> $err_file
  check
  comment "   rename all mct references"
    sed -i "s/\(${prefix}\)*mct_/${prefix}mct_/g" * >> $log_file 2>> $err_file
  check
    
route "${cblue}<< substitutions_oas${cnormal}"
}

configure_oas(){
route "${cblue}>> configure_oas${cnormal}"
  file=${oasdir}/util/make_dir/make.oas3
  comment "   cp jureca oasis3-mct makefile to /util/make_dir/"
    cp $rootdir/bldsva/intf_oas3/oasis3-mct/arch/$platform/config/make.intel_jureca_oa3 $file >> $log_file 2>> $err_file
  check
  c_configure_oas
  comment "   sed new psmile includes to Makefile"
    sed -i 's@__inc__@-I$(LIBBUILD)/psmile.$(CHAN) -I$(LIBBUILD)/scrip  -I$(LIBBUILD)/mct'" -I$ncdfPath/include@" $file >> $log_file 2>> $err_file
  check
  comment "   sed ldflg to oas Makefile"
    sed -i "s@__ldflg__@@" $file >> $log_file 2>> $err_file
  check
  comment "   sed comF90 to oas Makefile"
    if [[ $profiling == "scalasca" ]]; then
      sed -i "s@__comF90__@scorep-mpif90 $optComp@" $file >> $log_file 2>> $err_file
    else
      sed -i "s@__comF90__@${profComp} $mpiPath/bin/mpif90 $optComp@" $file >> $log_file 2>> $err_file
    fi
  check
  comment "   sed comCC to oas Makefile"
    if [[ $profiling == "scalasca" ]]; then
      sed -i "s@__comCC__@scorep-mpicc $optComp@" $file >> $log_file 2>> $err_file
    else
      sed -i "s@__comCC__@${profComp} $mpiPath/bin/mpicc $optComp@" $file >> $log_file 2>> $err_file
    fi
  check
  comment "   sed ld to oas Makefile"
    if [[ $profiling == "scalasca" ]]; then
      sed -i "s@__ld__@scorep-mpif90@" $file >> $log_file 2>> $err_file
    else
      sed -i "s@__ld__@${profComp} $mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
    fi
  check
  comment "   sed libs to oas Makefile"
    sed -i "s@__lib__@-L$ncdfPath/lib/ -lnetcdff@" $file >> $log_file 2>> $err_file
  check
  comment "   sed precision to oas Makefile"
    if [[ $compiler == "Intel" ]]; then
      sed -i "s@__precision__@-i4 -r8@" $file >> $log_file 2>> $err_file
    else
      sed -i "s@__precision__@@" $file >> $log_file 2>> $err_file
    fi
  check

route "${cblue}<< configure_oas${cnormal}"
}

make_oas(){
route "${cblue}>> make_oas${cnormal}"
  c_make_oas
route "${cblue}<< make_oas${cnormal}"
}


setup_oas(){
route "${cblue}>> setupOas${cnormal}"

  c_setup_oas

route "${cblue}<< setupOas${cnormal}"
}

