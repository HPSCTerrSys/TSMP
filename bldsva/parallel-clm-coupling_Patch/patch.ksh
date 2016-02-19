#! /usr/bin/ksh

################  User Setup section  ###################

export SVAROOT=/home/pshrestha/coup_oas
patchroot=${SVAROOT}/bldsva/parallel-clm-coupling_Patch

#########################################################
main(){
echo "########################################"
echo "start patching parallel clm coupling" 
echo "########################################"
echo " "

echo "cd in patchroot"
cd ${patchroot}
echo "copy new atmdrvMod in clm/bld/usr.src/"
cp atmdrvMod.F90 ${SVAROOT}/clm/bld/usr.src/
echo "copy new decompMod in clm/bld/usr.src/"
cp decompMod.F90 ${SVAROOT}/clm/bld/usr.src/
echo "copy all new clm coupling files to clm/src/oas3"
cp oas* ${SVAROOT}/clm/src/oas3/
cp receive* ${SVAROOT}/clm/src/oas3/
cp send* ${SVAROOT}/clm/src/oas3/
echo "set coupling processes accordingly in namcouple"
sed 's@set ncpl_exe3 = 1@set ncpl_exe3 = $nproc_clm@' -i $SVAROOT/bldsva/build_oas3
sed 's@set ncpl_exe2 = 1@set ncpl_exe2 = $nproc_clm@' -i $SVAROOT/bldsva/build_oas3
echo " "
echo "########################################"
echo "finished patching parallel clm coupling" 
echo "########################################"

}


logDate=$(date +"%Y%m%d%H%M")
main 2>&1 | tee -a patch${logDate}.log

