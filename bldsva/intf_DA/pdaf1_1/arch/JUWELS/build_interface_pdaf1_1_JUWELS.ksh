#! /bin/ksh
#

always_da(){
route "${cblue}>> always_da${cnormal}"
route "${cblue}<< always_da${cnormal}"
}

substitutions_da(){
route "${cblue}>> substitutions_da${cnormal}"
if [[ ${mList[1]} == clm5_0  ]] ; then
  comment "   mkdir  $dadir/interface"
    mkdir -p $dadir/interface  >> $log_file 2>> $err_file
  check

  comment "   cp pdaf interface model to $dadir/interface"
    patch $rootdir/bldsva/intf_DA/pdaf1_1/model $dadir/interface
  check

  comment "   cp pdaf interface framework to $dadir/interface"
    patch $rootdir/bldsva/intf_DA/pdaf1_1/framework $dadir/interface
  check

# CLM not 5.0
else
  comment "   mkdir  $dadir/interface"
    mkdir -p $dadir/interface  >> $log_file 2>> $err_file
  check

  comment "   cp pdaf interface model to $dadir/interface"
    patch $rootdir/bldsva/intf_DA/pdaf1_1/model $dadir/interface 
  check

  comment "   cp pdaf interface framework to $dadir/interface"
    patch $rootdir/bldsva/intf_DA/pdaf1_1/framework $dadir/interface 
  check
fi
route "${cblue}<< substitutions_da${cnormal}"
}

configure_da(){
route "${cblue}>> configure_da${cnormal}"
if [[ ${mList[1]} == clm5_0  ]] ; then
  export PDAF_DIR=$dadir
  export PDAF_ARCH=linux_ifort_juwels

#PDAF part
  file=$dadir/make.arch/linux_ifort_juwels.h
  
  comment "   cp pdaf config to $dadir"
    cp $rootdir/bldsva/intf_DA/pdaf1_1/arch/$platform/config/linux_ifort_juwels.h $file >> $log_file 2>> $err_file
  check

  comment "   sed comFC dir to $file"
    if [[ $profiling == "scalasca" ]]; then
      sed -i "s@__comFC__@scorep-mpif90@" $file >> $log_file 2>> $err_file
    else
      sed -i "s@__comFC__@${mpiPath}/bin/mpif90@" $file >> $log_file 2>> $err_file
    fi
  check

  comment "   sed comCC dir to $file"
    if [[ $profiling == "scalasca" ]]; then
      sed -i "s@__comCC__@scorep-mpicc@" $file >> $log_file 2>> $err_file
    else
      sed -i "s@__comCC__@${mpiPath}/bin/mpicc@" $file >> $log_file 2>> $err_file
    fi
  check

  comment "   sed MPI dir to $file"
    sed -i "s@__MPI_INC__@-I${mpiPath}/include@" $file >> $log_file 2>> $err_file
  check

  comment "   sed LIBS to $file"
    sed -i "s@__LIBS__@ $lapackPath/mkl/lib/intel64/libmkl_intel_lp64.a $lapackPath/mkl/lib/intel64/libmkl_intel_thread.a $lapackPath/mkl/lib/intel64/libmkl_core.a -L${mpiPath}/lib64@" $file >> $log_file 2>> $err_file
  check

  comment "   sed optimizations to $file"
    sed -i "s@__OPT__@${optComp}@" $file >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/src"
    cd $dadir/src >> $log_file 2>> $err_file
  check
  comment "   make clean pdaf"
    make clean >> $log_file 2>> $err_file
  check

#PDAF interface part
  file1=$dadir/interface/model/Makefile
  file2=$dadir/interface/framework/Makefile
  comment "   cp pdaf interface Makefiles to $dadir"
    cp $rootdir/bldsva/intf_DA/pdaf1_1/model/Makefile  $file1 >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_DA/pdaf1_1/framework/Makefile  $file2 >> $log_file 2>> $err_file
  check

  importFlags=" "
  cppdefs=" "
  obj=' '
  libs=" -L$mpiPath -lmpich -L$netcdfPath/lib/ -lnetcdff -lnetcdf "
  pf=""

  if [[ $withOAS == "false" && $withPFL == "true" ]] ; then
     comment "Not yet implemented combination of pdaf and models"
  fi

  if [[ $withOAS == "false" && $withCLM == "true" ]] ; then
     importFlags+=" -I$clmdir/build/intel/mpi/nodebug/nothreads/include "
     importFlags+=" -I$clmdir/build/intel/mpi/nodebug/nothreads/mct/noesmf/c1a1l1i1o1r1g1w1e1/include "
     importFlags+=" -I$clmdir/build/intel/mpi/nodebug/nothreads/mct/noesmf/include "
     importFlags+=" -I$clmdir/build/intel/mpi/nodebug/nothreads/mct/noesmf/clm/obj "
     importFlags+=" -I$clmdir/build/atm/obj "
     importFlags+=" -I$clmdir/build/ice/obj "
     importFlags+=" -I$clmdir/build/ocn/obj "
     importFlags+=" -I$clmdir/build/glc/obj "
     importFlags+=" -I$clmdir/build/rof/obj "
     importFlags+=" -I$clmdir/build/wav/obj "
     importFlags+=" -I$clmdir/build/esp/obj "
     importFlags+=" -I$clmdir/build/cpl/obj "
     importFlags+=" -I$clmdir/build/lib/include "
     importFlags+=" -I$clmdir/build/ "
     cppdefs+=" ${pf} -DCLMFIVE"
     libs+=" -L$clmdir/build/lib/ -lcpl "
     libs+=" -L$clmdir/build/lib/ -latm -lice "
     libs+=" -L$clmdir/build/intel/mpi/nodebug/nothreads/mct/noesmf/lib/ -lclm "
     libs+=" -L$clmdir/build/lib/ -locn -lrof -lglc -lwav -lesp "
     libs+=" -L$clmdir/build/intel/mpi/nodebug/nothreads/mct/noesmf/c1a1l1i1o1r1g1w1e1/lib -lcsm_share "
     libs+=" -L$clmdir/build/intel/mpi/nodebug/nothreads/lib -lpio -lgptl -lmct -lmpeu  "
     libs+=" -lpnetcdf  -mkl -lnetcdff -lnetcdf "
     obj+=' $(OBJCLM5)'
  fi

  if [[ $withCLM == "true" && $withCOS == "true" && $withPFL == "false" ]] ; then
     comment "Not yet implemented combination of pdaf and models"
  fi

  if [[ $withCLM == "true" && $withCOS == "false" && $withPFL == "true" ]] ; then
     comment "Not yet implemented combination of pdaf and models"
  fi
  
  if [[ $withCLM == "true" && $withCOS == "true" && $withPFL == "true" ]] ; then
     comment "Not yet implemented combination of pdaf and models"
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

# CLM not 5.0
else
  export PDAF_DIR=$dadir
  export PDAF_ARCH=linux_ifort_juwels

#PDAF part
  file=$dadir/make.arch/linux_ifort_juwels.h
  
  comment "   cp pdaf config to $dadir"
    cp $rootdir/bldsva/intf_DA/pdaf1_1/arch/$platform/config/linux_ifort_juwels.h $file >> $log_file 2>> $err_file
  check

  comment "   sed comFC dir to $file" 
    if [[ $profiling == "scalasca" ]]; then
      sed -i "s@__comFC__@scorep-mpif90@" $file >> $log_file 2>> $err_file  
    else
      sed -i "s@__comFC__@${mpiPath}/bin/mpif90@" $file >> $log_file 2>> $err_file  
    fi
  check

  comment "   sed comCC dir to $file" 
    if [[ $profiling == "scalasca" ]]; then
      sed -i "s@__comCC__@scorep-mpicc@" $file >> $log_file 2>> $err_file  
    else
      sed -i "s@__comCC__@${mpiPath}/bin/mpicc@" $file >> $log_file 2>> $err_file  
    fi
  check

  comment "   sed MPI dir to $file" 
    sed -i "s@__MPI_INC__@-I${mpiPath}/include@" $file >> $log_file 2>> $err_file  
  check

  comment "   sed LIBS to $file"
    sed -i "s@__LIBS__@ -L$lapackPath -lopenblas -L${mpiPath}/lib64@" $file >> $log_file 2>> $err_file
    #sed -i "s@__LIBS__@ -L$lapackPath -llapack -lblas -L${mpiPath}/lib64@" $file >> $log_file 2>> $err_file
    #sed -i "s@__LIBS__@ -L$lapackPath/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_scalapack_ilp64 -lmkl_cdft_core -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_ilp64 -lm -ldl -L${mpiPath}/lib64 -lirc -lintlc@" $file >> $log_file 2>> $err_file
    #sed -i "s@__LIBS__@ $lapackPath/mkl/lib/intel64/libmkl_intel_lp64.a $lapackPath/mkl/lib/intel64/libmkl_intel_thread.a $lapackPath/mkl/lib/intel64/libmkl_core.a -L${mpiPath}/lib64@" $file >> $log_file 2>> $err_file
  check

  comment "   sed optimizations to $file"
    sed -i "s@__OPT__@${optComp}@" $file >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/src"
    cd $dadir/src >> $log_file 2>> $err_file
  check
  comment "   make clean pdaf"
    make clean >> $log_file 2>> $err_file
  check

#PDAF interface part
  file1=$dadir/interface/model/Makefile
  file2=$dadir/interface/framework/Makefile
  comment "   cp pdaf interface Makefiles to $dadir"
    cp $rootdir/bldsva/intf_DA/pdaf1_1/model/Makefile  $file1 >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_DA/pdaf1_1/framework/Makefile  $file2 >> $log_file 2>> $err_file
  check

  importFlags="-g "
  cppdefs=" "
  obj=' ' 
  libs=" -L$mpiPath -lmpich -L$netcdfPath/lib/ -lnetcdff -lnetcdf " 
  pf=""
 
  if [[ $withOAS == "false" && $withPFL == "true" ]] ; then
     importFlags+="-I$pfldir/pfsimulator/parflow_lib -I$pfldir/pfsimulator/amps/common -I$pfldir/pfsimulator/amps/oas3 -I$pfldir/pfsimulator/include "
     cppdefs+=" ${pf}-DPARFLOW_STAND_ALONE "
     libs+=" -L$hyprePath/lib -L$siloPath/lib -lparflow -lamps -lamps_common -lamps -lamps_common -lkinsol -lgfortran -lHYPRE -lsilo "
     obj+=' $(OBJPF) '
  fi

  if [[ $withOAS == "false" && $withCLM == "true" ]] ; then
     importFlags+=" -I$clmdir/build/ "
     cppdefs+=" ${pf}-DCLMSA "
     libs+=" -lclm "
     obj+=' $(OBJCLM) print_update_clm.o'
  fi

  if [[ $withCLM == "true" && $withCOS == "true" && $withPFL == "false" ]] ; then
     importFlags+=" -I$clmdir/build/ -I$oasdir/$platform/build/lib/psmile.MPI1 -I$oasdir/$platform/build/lib/scrip -I$cosdir/obj "
     cppdefs+=" ${pf}-Duse_comm_da ${pf}-DCOUP_OAS_COS ${pf}-DGRIBDWD ${pf}-DNETCDF ${pf}-DHYMACS ${pf}-DMAXPATCH_PFT=1 "
     if [[ $cplscheme == "true" ]] ; then ; cppdefs+=" ${pf}-DCPL_SCHEME_F " ; fi
     if [[ $readCLM == "true" ]] ; then ; cppdefs+=" ${pf}-DREADCLM " ; fi
     libs+=" -lclm -lcosmo -lpsmile.MPI1 -lmct -lmpeu -lscrip $grib1Path/libgrib1.a "
     obj+=' $(OBJCLM) $(OBJCOSMO) '
  fi 

  if [[ $withCLM == "true" && $withCOS == "false" && $withPFL == "true" ]] ; then
     importFlags+=" -I$clmdir/build/ -I$oasdir/$platform/build/lib/psmile.MPI1 -I$oasdir/$platform/build/lib/scrip -I$pfldir/pfsimulator/parflow_lib -I$pfldir/pfsimulator/amps/oas3 -I$pfldir/pfsimulator/include -I$pfldir/pfsimulator/amps/common "
     cppdefs+=" ${pf}-Duse_comm_da ${pf}-DCOUP_OAS_PFL ${pf}-DMAXPATCH_PFT=1 "
     if [[ $readCLM == "true" ]] ; then ; cppdefs+=" ${pf}-DREADCLM " ; fi
     if [[ $freeDrain == "true" ]] ; then ; cppdefs+=" ${pf}-DFREEDRAINAGE " ; fi
     libs+=" -lclm -lpsmile.MPI1 -lmct -lmpeu -lscrip -L$hyprePath/lib -L$siloPath/lib -lparflow -lamps -lamps_common -lamps -lamps_common -lkinsol -lgfortran -lHYPRE -lsilo "
     obj+=' $(OBJCLM) $(OBJPF) '
  fi
  if [[ $withCLM == "true" && $withCOS == "true" && $withPFL == "true" ]] ; then
     importFlags+=" -I$clmdir/build/ -I$oasdir/$platform/build/lib/psmile.MPI1 -I$oasdir/$platform/build/lib/scrip -I$pfldir/pfsimulator/parflow_lib -I$pfldir/pfsimulator/amps/oas3 -I$pfldir/pfsimulator/include -I$pfldir/pfsimulator/amps/common -I$cosdir/obj "
     cppdefs+=" ${pf}-Duse_comm_da ${pf}-DCOUP_OAS_COS ${pf}-DGRIBDWD ${pf}-DNETCDF ${pf}-DHYMACS ${pf}-DMAXPATCH_PFT=1 ${pf}-DCOUP_OAS_PFL "
     if [[ $cplscheme == "true" ]] ; then ; cppdefs+=" ${pf}-DCPL_SCHEME_F " ; fi
     if [[ $readCLM == "true" ]] ; then ; cppdefs+=" ${pf}-DREADCLM " ; fi
     if [[ $freeDrain == "true" ]] ; then ; cppdefs+=" ${pf}-DFREEDRAINAGE " ; fi
     libs+=" -lclm -lpsmile.MPI1 -lmct -lmpeu -lscrip $grib1Path/libgrib1.a -L$hyprePath/lib -L$siloPath/lib -lparflow -lamps -lamps_common -lamps -lamps_common -lkinsol -lgfortran -lHYPRE -lsilo -lcosmo $grib1Path/libgrib1.a"
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

fi
route "${cblue}<< configure_da${cnormal}"
}

make_da(){
route "${cblue}>> make_da${cnormal}"
if [[ ${mList[1]} == clm5_0  ]] ; then
  export PDAF_DIR=$dadir
  export PDAF_ARCH=linux_ifort_juwels

  comment "   cd to $dadir/src"
    cd $dadir/src >> $log_file 2>> $err_file
  check
  comment "   make pdaf"
    make >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/interface/model"
    cd $dadir/interface/model >> $log_file 2>> $err_file
  check
  comment "   make pdaf model"
    make >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/interface/framework"
    cd $dadir/interface/framework >> $log_file 2>> $err_file
  check
  comment "   make pdaf framework"
    make >> $log_file 2>> $err_file
  check

# CLM not 5.0
else
  export PDAF_DIR=$dadir
  export PDAF_ARCH=linux_ifort_juwels

  comment "   cd to $dadir/src"
    cd $dadir/src >> $log_file 2>> $err_file
  check
  comment "   make pdaf"
    make >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/interface/model"
    cd $dadir/interface/model >> $log_file 2>> $err_file
  check
  comment "   make pdaf model"
    make >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/interface/framework"
    cd $dadir/interface/framework >> $log_file 2>> $err_file
  check
  comment "   make pdaf framework"
    make >> $log_file 2>> $err_file
  check

fi
route "${cblue}<< make_da${cnormal}"
}


setup_da(){
route "${cblue}>> setup_da${cnormal}"
  c_setup_pdaf
route "${cblue}<< setup_da${cnormal}"
}

