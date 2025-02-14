#! /bin/ksh

always_pfl(){
route "${cyellow}>> always_pfl${cnormal}"
route "${cyellow}<< always_pfl${cnormal}"
}

configure_pfl(){
route "${cyellow}>> configure_pfl${cnormal}"
    export PARFLOW_INS="$pfldir/bin"
    export PARFLOW_BLD="$pfldir/build"
#    export PFV="oas-gpu"
    export RMM_ROOT=$pfldir/rmm
#
    C_FLAGS="-fopenmp -Wall -Werror"
    flagsSim="  -DMPIEXEC_EXECUTABLE=$(which srun)"
    if [[ $withOAS == "true" ]]; then
        flagsSim+=" -DPARFLOW_AMPS_LAYER=oas3"
    else
      if [[ $withPDAF == "true" ]] ; then
        flagsSim+=" -DPARFLOW_AMPS_LAYER=da"
      else
        flagsSim+=" -DPARFLOW_AMPS_LAYER=mpi1"
      fi
    fi
    flagsSim+=" -DOAS3_ROOT=$oasdir/$platform"
    flagsSim+=" -DSILO_ROOT=$siloPath"
    flagsSim+=" -DHYPRE_ROOT=$hyprePath"
    flagsSim+=" -DCMAKE_C_FLAGS=$C_FLAGS"
    flagsSim+=" -DCMAKE_BUILD_TYPE=Release"
    flagsSim+=" -DPARFLOW_ENABLE_TIMING=TRUE"
    flagsSim+=" -DCMAKE_INSTALL_PREFIX=$PARFLOW_INS"
    flagsSim+=" -DNETCDF_DIR=$ncdfPath"
    flagsSim+=" -DNETCDF_Fortran_ROOT=$ncdfPath"
    flagsSim+=" -DTCL_TCLSH=$tclPath/bin/tclsh8.6"
    flagsSim+=" -DPARFLOW_AMPS_SEQUENTIAL_IO=on"
    if [[ $withPDAF == "true" ]] ; then
      # PDAF:
      # Turn off SLURM for PDAF due to error
      # Only used for finishing job close to SLURM time limit
      # Could be added in future effort
      flagsSim+=" -DPARFLOW_ENABLE_SLURM=FALSE"
    else
      flagsSim+=" -DPARFLOW_ENABLE_SLURM=TRUE"
    fi
	# Define compilers
    if [[ $profiling == "scalasca" ]]; then
      pcc="scorep-mpicc"
      pfc="scorep-mpif90"
      pf77="scorep-mpif77"
      pcxx="scorep-mpicxx"
      flagsTools+="CC=scorep-mpicc FC=scorep-mpif90 F77=scorep-mpif77 "
    else
	  pcc="$mpiPath/bin/mpicc"
      pfc="$mpiPath/bin/mpif90"
      pf77="$mpiPath/bin/mpif77"
      pcxx="$mpiPath/bin/mpic++"
    fi
#
    comment "    add parflow paths $PARFLOW_INS, $PARFLOW_BLD "
     mkdir -p $PARFLOW_INS
     mkdir -p $PARFLOW_BLD
    check

    comment " parflow is configured for $processor "
    check
    if [[ $processor == "GPU"|| $processor == "MSA" ]]; then
       cd $pfldir
       comment "module load CUDA  mpi-settings/CUDA "
        module load CUDA  $gpuMpiSettings >> $log_file 2>> $err_file
       check
       comment "    additional configuration options for GPU are set "
        flagsSim+=" -DPARFLOW_ACCELERATOR_BACKEND=cuda"
        flagsSim+=" -DRMM_ROOT=$RMM_ROOT/install"
        flagsSim+=" -DCMAKE_CUDA_RUNTIME_LIBRARY=Shared"
        flagsSim+=" -DCMAKE_CUDA_ARCHITECTURES=90"
        #flagsSim+=" -DCMAKE_EXE_LINKER_FLAGS=\"-lcurand -lcusparse -lcublas\"" 
        flagsSim+=" -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON"   
       check
       comment "    git clone  RAPIDS Memory Manager "
       if [ -d $RMM_ROOT ] ; then 
        comment "  remove $RMM_ROOT "
        rm -rf $RMM_ROOT >> $log_file 2>> $err_file
        check
       fi
       #git clone -b branch-0.10 --single-branch --recurse-submodules https://github.com/hokkanen/rmm.git >> $log_file 2>> $err_file
       git clone -b branch-0.10 --single-branch --recurse-submodules https://github.com/mfahdaz/rmm.git >> $log_file 2>> $err_file
       check
        mkdir -p $RMM_ROOT/build
        cd $RMM_ROOT/build
       comment "    configure RMM: RAPIDS Memory Manager "
        cmake ../ -DCMAKE_INSTALL_PREFIX=$RMM_ROOT/install ${cuda_architectures} >> $log_file 2>> $err_file
       check
       comment "    make RMM "
        make -j  >> $log_file 2>> $err_file
       check
       comment "    make install RMM "
        make install >> $log_file 2>> $err_file
       check
    fi

    c_configure_pfl

route "${cyellow}<< configure_pfl${cnormal}"
}

make_pfl(){
route "${cyellow}>> make_pfl${cnormal}"
  c_make_pfl
route "${cyellow}<< make_pfl${cnormal}"
}


substitutions_pfl(){
route "${cyellow}>> substitutions_pfl${cnormal}"

  c_substitutions_pfl

route "${cyellow}<< substitutions_pfl${cnormal}"
}




