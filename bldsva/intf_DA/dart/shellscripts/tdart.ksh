#!/bin/ksh
#-----------------------------------------------------------------------
SPATH="$HOME/terrsysmp/bldsva/intf_DA/dart/shellscripts/"
DPATH="$HOME/DART/lanai/models/terrsysmp/"
NENS=49
ICYCLE=2
NRST=0
MACHINE="JUWELS"
#
MAP_FN="$SPATH/map_fn.txt"
#-----------------------------------------------------------------------
# Manually setup TerrSysMP-DART run
#-----------------------------------------------------------------------
# Update DART source codes including namelists and  clean compile 
if [[ $ICYCLE == 1 ]] then
  ./dart_setup.csh $NENS $DPATH
fi

# Setup TerrSysMP testcase and perturb model states and parameters
./tsmp_setup.csh $ICYCLE $NRST $NENS $SPATH/map_fn.txt $MACHINE

exit 0
