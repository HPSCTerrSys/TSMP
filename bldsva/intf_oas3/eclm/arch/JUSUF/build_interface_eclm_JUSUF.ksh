#!/bin/ksh

always_clm(){
route "${cyellow}>> always_clm${cnormal}"
route "${cyellow}<< always_clm${cnormal}"
}

configure_clm(){
route "${cyellow}>> configure_clm${cnormal}"
  #
  # Uncomment these lines to customize build options.
  #
  # ECLM_CC=/path/to/mpicc
  # ECLM_FC=/path/to/mpifort
  # ECLM_CMAKE_VARS+=" -DBUILD_VAR1=VALUE1"
  # ECLM_CMAKE_VARS+=" -DBUILD_VAR2=VALUE2"
  #
  c_configure_eclm
route "${cyellow}<< configure_clm${cnormal}"
}

make_clm(){
route "${cyellow}>> make_clm${cnormal}"
    c_make_eclm
route "${cyellow}<< make_clm${cnormal}"
}

substitutions_clm(){
route "${cyellow}>> substitutions_clm${cnormal}"
route "${cyellow}<< substitutions_clm${cnormal}"
}

setup_clm(){
route "${cyellow}>> setupClm${cnormal}"
    c_setup_eclm
route "${cyellow}<< setupClm${cnormal}"
}
