#!/bin/csh
# Copy it to the archive spinup directory
# Execute it from within the directory
Usage: ./archive_restart_hydroSpinup.csh N
set Nens = $2
#
mkdir restart
cd restart
foreach inst (`seq 0 ${Nens}`)
mkdir tsmp_instance_${inst}
cp ../tsmp_instance_${inst}/pflout/restart/* tsmp_instance_${inst}/ 
cp ../tsmp_instance_${inst}/clmout/restart/* tsmp_instance_${inst}/
end
exit 0
