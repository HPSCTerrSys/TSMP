#!/bin/csh
echo "Usage: check_spinup.csh  rundir N"

set rundir  = $1
set numInst = $2

cd ${rundir}
foreach instance (`seq 0 $numInst`)
  cd tsmp_instance_$instance
  if ( -f timing_all ) then
    echo "Spinup Complete : " $instance
  else
    echo "XXXXX"
  endif
  cd ..
end
exit 0

