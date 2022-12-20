#! /bin/ksh
#

always_da(){
route "${cyellow}>> always_da${cnormal}"
route "${cyellow}<< always_da${cnormal}"
}

substitutions_da(){
route "${cyellow}>> substitutions_da${cnormal}"
  c_substitutions_pdaf
route "${cyellow}<< substitutions_da${cnormal}"
}

configure_da(){
route "${cyellow}>> configure_da${cnormal}"

#PDAF part configuration variables
  export PDAF_DIR=$dadir
  export PDAF_ARCH=linux_gnu_generic_x86

  comFC="/opt/mpich-3.2.1/bin/mpif90"
  comCC="/opt/mpich-3.2.1/bin/mpicc"

  libs_src=" -L${lapackPath}/lib64 -llapack  -lblas -L${mpiPath}/lib64 "

  c_configure_pdaf_arch

#PDAF interface part
  file1=$dadir/interface/model/Makefile
  file2=$dadir/interface/framework/Makefile
  comment "   cp pdaf interface Makefiles to $dadir"
    cp $rootdir/bldsva/intf_DA/pdaf1_1/model/Makefile  $file1 >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_DA/pdaf1_1/framework/Makefile  $file2 >> $log_file 2>> $err_file
  check

  importFlags=" "
  importFlagsOAS=" "
  importFlagsPFL=" "
  importFlagsCLM=" "
  importFlagsCOS=" "
  importFlagsDA=" "
  cppdefs=" "
  obj=' '
  libs=" -L$mpiPath -lmpich -L$netcdfPath/lib/ -lnetcdff -lnetcdf "
  libsOAS=" "
  libsPFL=" "
  libsCLM=" "
  libsCOS=" "
  pf=""

  # Oasis include dirs
  importFlagsOAS+="-I$oasdir/$platform/build/lib/psmile.MPI1 "
  importFlagsOAS+="-I$oasdir/$platform/build/lib/scrip "

  # CLM include dirs
  importFlagsCLM+="-I$clmdir/build/ "

  # COSMO include dirs
  importFlagsCOS+="-I$cosdir/obj "

  # ParFlow include dirs
  importFlagsPFL+="-I$pfldir/pfsimulator/parflow_lib "
  importFlagsPFL+="-I$pfldir/pfsimulator/amps/oas3 "
  # importFlagsPFL+="-I$pfldir/pfsimulator/amps/common "
  if [[ ${mList[3]} == parflow ]] ; then
    importFlagsPFL+="-I$pfldir/build/include "
    if [[ $processor == "GPU" ]]; then
      importFlagsPFL+="-I$pfldir/rmm/include/rmm "
    fi
  else
    importFlagsPFL+="-I$pfldir/pfsimulator/include "
  fi

  # DA include dirs
  importFlagsDA+="-I$dadir/interface/model/common "
  if [[ $withPFL == "true" ]] ; then
    importFlagsDA+="-I$dadir/interface/model/${mList[3]} "
  fi

  # Oasis libraries
  libsOAS+="-lpsmile.MPI1 "
  libsOAS+="-lmct "
  libsOAS+="-lmpeu "
  libsOAS+="-lscrip "

  # CLM libraries
  libsCLM+="-lclm "

  # COSMO libraries
  libsCOS+="-lcosmo "
  if [[ ${mList[2]} == "cosmo5_1" ]] ; then
    libsCOS+="-L$gribPath/lib/ "
    libsCOS+="-leccodes_f90 "
    libsCOS+="-leccodes "
  else
    libsCOS+="$grib1Path/libgrib1.a "
  fi

  # ParFlow library paths and libraries
  if [[ ${mList[3]} == parflow ]] ; then
    libsPFL+="-lpfsimulator "
  else
    libsPFL+="-lparflow "
  fi
  libsPFL+="-lamps "
  if [[ ${mList[3]} == parflow3_2 || ${mList[3]} == parflow3_0 ]] ; then
    libsPFL+="-lamps_common "
    libsPFL+="-lamps "
    libsPFL+="-lamps_common "
  fi
  if [[ ${mList[3]} == parflow ]] ; then
    libsPFL+="-lpfkinsol "
  else
    libsPFL+="-lkinsol "
  fi
  libsPFL+="-lgfortran "
  if [[ ${mList[3]} == parflow ]] ; then
    libsPFL+="-lcjson "
    if [[ $processor == "GPU" ]]; then
      libsPFL+="-lstdc++ "
      libsPFL+="-lcudart "
      libsPFL+="-lrmm "
      libsPFL+="-lnvToolsExt "
    fi
  fi
  libsPFL+="-L$hyprePath/lib -lHYPRE "
  libsPFL+="-L$siloPath/lib -lsilo "

  if [[ $withOAS == "false" && $withPFL == "true" ]] ; then
     importFlags+=$importFlagsPFL
     importFlags+=$importFlagsDA
     cppdefs+=" ${pf}-DPARFLOW_STAND_ALONE "
     libs+=$libsPFL
     obj+=' $(OBJPF) '
  fi

  if [[ $withOAS == "false" && $withCLM == "true" ]] ; then
     importFlags+=$importFlagsCLM
     importFlags+=$importFlagsDA
     cppdefs+=" ${pf}-DCLMSA "
     libs+=$libsCLM
     obj+=' $(OBJCLM) print_update_clm.o'
  fi

  if [[ $withCLM == "true" && $withCOS == "true" && $withPFL == "false" ]] ; then
     importFlags+=$importFlagsCLM
     importFlags+=$importFlagsOAS
     importFlags+=$importFlagsCOS
     importFlags+=$importFlagsDA
     cppdefs+=" ${pf}-Duse_comm_da ${pf}-DCOUP_OAS_COS ${pf}-DGRIBDWD ${pf}-DNETCDF ${pf}-DHYMACS ${pf}-DMAXPATCH_PFT=1 "
     if [[ $cplscheme == "true" ]] ; then ; cppdefs+=" ${pf}-DCPL_SCHEME_F " ; fi
     if [[ $readCLM == "true" ]] ; then ; cppdefs+=" ${pf}-DREADCLM " ; fi
     libs+=$libsCLM
     libs+=$libsCOS
     libs+=$libsOAS
     obj+=' $(OBJCLM) $(OBJCOSMO) '
  fi

  if [[ $withCLM == "true" && $withCOS == "false" && $withPFL == "true" ]] ; then
     importFlags+=$importFlagsCLM
     importFlags+=$importFlagsOAS
     importFlags+=$importFlagsPFL
     importFlags+=$importFlagsDA
     cppdefs+=" ${pf}-Duse_comm_da ${pf}-DCOUP_OAS_PFL ${pf}-DMAXPATCH_PFT=1 "
     cppdefs+=" ${pf}-DOBS_ONLY_PARFLOW " # Remove for observations from both ParFlow + CLM
     if [[ $readCLM == "true" ]] ; then ; cppdefs+=" ${pf}-DREADCLM " ; fi
     if [[ $freeDrain == "true" ]] ; then ; cppdefs+=" ${pf}-DFREEDRAINAGE " ; fi
     libs+=$libsCLM
     libs+=$libsOAS
     libs+=$libsPFL
     obj+=' $(OBJCLM) $(OBJPF) '
  fi
  if [[ $withCLM == "true" && $withCOS == "true" && $withPFL == "true" ]] ; then
     importFlags+=$importFlagsCLM
     importFlags+=$importFlagsOAS
     importFlags+=$importFlagsPFL
     importFlags+=$importFlagsCOS
     importFlags+=$importFlagsDA
     cppdefs+=" ${pf}-Duse_comm_da ${pf}-DCOUP_OAS_COS ${pf}-DGRIBDWD ${pf}-DNETCDF ${pf}-DHYMACS ${pf}-DMAXPATCH_PFT=1 ${pf}-DCOUP_OAS_PFL "
     if [[ $cplscheme == "true" ]] ; then ; cppdefs+=" ${pf}-DCPL_SCHEME_F " ; fi
     if [[ $readCLM == "true" ]] ; then ; cppdefs+=" ${pf}-DREADCLM " ; fi
     if [[ $freeDrain == "true" ]] ; then ; cppdefs+=" ${pf}-DFREEDRAINAGE " ; fi
     libs+=$libsCLM
     libs+=$libsOAS
     libs+=$libsCOS
     libs+=$libsPFL
     obj+=' $(OBJCLM) $(OBJCOSMO) $(OBJPF) '
  fi

  comment "   sed bindir to Makefiles"
    sed -i "s,__bindir__,$bindir," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed comp flags to Makefiles"
    sed -i "s,__fflags__,-cpp -I$dadir/interface/model -I$ncdfPath/include $importFlags," $file1 $file2 >> $log_file 2>> $err_file
  check
    sed -i "s,__ccflags__,-I$dadir/interface/model -I$ncdfPath/include $importFlags," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed preproc flags to Makefiles"
    sed -i "s,__cpp_defs__,$cppdefs," $file1 $file2 >> $log_file 2>> $err_file
  check
    sed -i "s,__fcpp_defs__,$cppdefs," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed libs to Makefiles"
    sed -i "s,__libs__,$libs," $file2 >> $log_file 2>> $err_file
  check
  comment "   sed obj to Makefiles"
    sed -i "s,__obj__,$obj," $file1 >> $log_file 2>> $err_file
  check
  comment "   sed -D prefix to Makefiles"
    sed -i "s,__pf__,$pf," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed clm directory to Makefiles"
    sed -i "s,__clmdir__,${mList[1]}," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed cosmo directory to Makefiles"
    sed -i "s,__cosdir__,${mList[2]}," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed parflow directory to Makefiles"
    sed -i "s,__pfldir__,${mList[3]}," $file1 $file2 >> $log_file 2>> $err_file
  check


  comment "   cd to $dadir/interface/model"
    cd $dadir/interface/model >> $log_file 2>> $err_file
  check
  comment "   make clean model"
    make clean >> $log_file 2>> $err_file
  check
  comment "   cd to $dadir/src/interface/framework"
    cd $dadir/interface/framework >> $log_file 2>> $err_file
  check
  comment "   make clean framework"
    make clean >> $log_file 2>> $err_file
  check


route "${cyellow}<< configure_da${cnormal}"
}

make_da(){
route "${cyellow}>> make_da${cnormal}"
  c_make_pdaf
route "${cyellow}<< make_da${cnormal}"
}

setup_da(){
route "${cyellow}>> setup_da${cnormal}"
  c_setup_pdaf
route "${cyellow}<< setup_da${cnormal}"
}
