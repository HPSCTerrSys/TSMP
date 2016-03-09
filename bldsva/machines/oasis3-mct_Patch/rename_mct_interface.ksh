#! /usr/bin/ksh

################  User Setup section  ###################

#SVAROOT=$HOME/coup_oas  #activate for running standalone
mctSourceDir=${SVAROOT}/oasis3/lib/mct/mct
psmileSourceDir=${SVAROOT}/oasis3/lib/psmile/src
prefix="oas_"

#########################################################

main(){

echo "########################################"
echo "start patching mct and oasis3-mct" 
echo "for not interacting with other mct"
echo "   -> all files and interfaces" 
echo "      are prefixed with: ${prefix}"
echo "########################################"
echo " "

echo "cd in mct-source Dir: ${mctSourceDir}"
cd ${mctSourceDir}

echo "rename all mct files"	
for i in m_* ; do mv $i ${prefix}${i} ; done 

echo "also rename mct_mod.F90"
mv mct_mod.F90 ${prefix}mct_mod.F90

echo "rename all interfaces in mct_mod.F90 and Makefile"
sed -i "s/mct_/${prefix}mct_/g" ${prefix}mct_mod.F90
sed -i "s/mct_/${prefix}mct_/g" Makefile

echo "rename all sources in mct_mod.F90 and Makefile"
sed -i "s/m_/${prefix}m_/g" ${prefix}mct_mod.F90 
sed -i "s/m_/${prefix}m_/g" Makefile

echo "except those which are redirecting mpeu sources"
sed -i "s/${prefix}m_List/m_List/g" ${prefix}mct_mod.F90
sed -i "s/${prefix}m_string/m_string/g" ${prefix}mct_mod.F90
sed -i "s/${prefix}m_die/m_die/g" ${prefix}mct_mod.F90
sed -i "s/${prefix}m_MergeSort/m_MergeSort/g" ${prefix}mct_mod.F90
sed -i "s/${prefix}m_inpak90/m_inpak90/g" ${prefix}mct_mod.F90
sed -i "s/${prefix}m_Permuter/m_Permuter/g" ${prefix}mct_mod.F90

echo "rename all module names"
sed -i "s/module m_/module ${prefix}m_/g" *

echo "rename all dependencies"
sed -i "s/use m_MCTW/use ${prefix}m_MCTW/g" *
sed -i "s/use m_Glob/use ${prefix}m_Glob/g" *
sed -i "s/use m_AttrV/use ${prefix}m_AttrV/g" *
sed -i "s/use m_Accum/use ${prefix}m_Accum/g" *
sed -i "s/use m_Gener/use ${prefix}m_Gener/g" *
sed -i "s/use m_Spa/use ${prefix}m_Spa/g" *
sed -i "s/use m_Nav/use ${prefix}m_Nav/g" *
sed -i "s/use m_Con/use ${prefix}m_Con/g" *
sed -i "s/use m_Ex/use ${prefix}m_Ex/g" *
sed -i "s/use m_Rou/use ${prefix}m_Rou/g" *
sed -i "s/use m_Rear/use ${prefix}m_Rear/g" *

echo "cd in psmile-source dir: ${psmileSourceDir}" 
cd ${psmileSourceDir}

echo "rename all mct references"
sed -i "s/mct_/${prefix}mct_/g" *

echo " "
echo "########################################"
echo "Finished"
echo "########################################"
}

logDate=$(date +"%Y%m%d%H%M")
main 2>&1 | tee -a patch_mct${logDate}.log

