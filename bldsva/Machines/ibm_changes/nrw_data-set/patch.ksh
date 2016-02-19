#! /usr/bin/ksh

################  User Setup section  ###################

export SVAROOT=$HOME/coup_oas
patchroot=${SVAROOT}/bldsva/ibm_changes/nrw_data-set

#########################################################
main(){
echo "########################################"
echo "start patching NRW-domain" 
echo "########################################"
echo " "

echo "cd in patchroot"
cd ${patchroot}
echo "copy namelists from clm, cosmo and parflow"
cp ${patchroot}/lmrun_uc $SVAROOT/cosmo/
cp ${patchroot}/lnd.stdin $SVAROOT/clm/bld/
cp ${patchroot}/coup_oas.tcl $SVAROOT/parflow/coup_oasrun/coup_oas.tcl
echo "change path of forcing-files in build_oas3"
sed "s@set forcingdir = .*@set forcingdir = /work/slts/slts00/tsmp/TerrSysMPdb/testdata_NRW_std@" -i $SVAROOT/bldsva/build_oas3
sed 's@set clm_data = .*@set clm_data = $forcingdir/clm3.5/Rur_NRW@' -i $SVAROOT/bldsva/build_oas3

echo "uncomment ascii2pfb.tcl-stuff in build_oas3"
sed "s/set nrw = 0/set nrw = 1/" -i $SVAROOT/bldsva/build_oas3

echo "set readclm = 0 in cosmo/src/oas3/oas_cos_define.F90"
sed "s/readclm = 1/readclm = 0/" -i $SVAROOT/cosmo/src/oas3/oas_cos_define.F90

echo " "
echo "########################################"
echo "finished patching NRW-domain" 
echo "########################################"

}


logDate=$(date +"%Y%m%d%H%M")
main 2>&1 | tee -a patch${logDate}.log

