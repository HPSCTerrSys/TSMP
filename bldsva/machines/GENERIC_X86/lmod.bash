#!/bin/bash

LMOD_PKG="/usr/share/lmod/lmod"
LMOD_DIR="/usr/local/software/lmod/lmod/libexec/"
LMOD_CMD="/usr/local/software/lmod/lmod/libexec/lmod"
MODULESHOME="/usr/local/software/lmod/lmod"
export LMOD_CMD
export LMOD_PKG
export LMOD_DIR
export MODULESHOME

module()
{
  eval `$LMOD_CMD sh "$@"`
}
export LMOD_DIR

clearMT()
{
  eval $($LMOD_DIR/clearMT_cmd bash)
}

########################################################################
#  ml is a shorthand tool for people who can't type moduel, err, module
#  It is also a combination command:
#     ml            -> module list
#     ml gcc        -> module load gcc
#     ml -gcc intel -> module unload gcc; module load intel
#  It does much more do: "ml --help" for more information.
ml()
{
  eval $($LMOD_DIR/ml_cmd "$@")
}

# Local Variables:
# mode: shell-script
# indent-tabs-mode: nil
# End:
