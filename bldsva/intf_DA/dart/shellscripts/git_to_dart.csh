#!/bin/csh -f
#
# This script copies all of the TerrSysMP git-controlled source files and scripts to the
# rest of the DART framework (which is all under svn).
# If an environment variable TERRSYSMP points to the TerrSysMP installation, it will be used.
# If an environment variable DART points to the DART installation, it will be used.
# If not, the default locations 
#
# The clm, parflow, and cosmo directories are all installed into a DART directory called
# 'models/terrsysmp/'.
#
# After the files are copied, it is required to build all the executables by kjj

cat << EndOfText >! exclude-list.txt
*.mod
*.swp
Makefile
advance_time
clm_to_dart
cosmo_to_dart
parflow_to_dart
create_obs_sequence
dart_log.*
dart_to_clm
dart_to_cosmo
dart_to_parflow
fill_inflation_restart
filter
input*default
model_mod_check
obs_sequence_tool
perfect_model_obs
preprocess
EndOfText

if ( ! $?TERRSYSMP ) then
   set   TERRSYSMP = $HOME/terrsysmp
endif

if ( ! $?DART ) then
   set  DART = $HOME/DART/lanai
endif

echo "Installing $TERRSYSMP into $DART"

foreach FILE ( clm cosmo parflow )

   echo ""
   echo "========================================"
   echo "Synching $FILE directory"

   set SOURCE = ${TERRSYSMP}/bldsva/intf_DA/dart/$FILE/
   set   DEST = ${DART}/models/terrsysmp/$FILE/

   if (! -d $DEST) then
      mkdir -p $DEST
   endif

   rsync -Cavzr --exclude-from=exclude-list.txt ${SOURCE} ${DEST}
   echo ""

end

echo "to finish up ... "
echo "cd ${DART}/models/terrsysmp  and check/build each of the models."

exit 0

