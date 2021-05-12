#!/bin/ksh

always_clm(){
route "${cyellow}>> always_clm${cnormal}"
route "${cyellow}<< always_clm${cnormal}"
}

configure_clm(){
route "${cyellow}>> configure_clm${cnormal}"
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
