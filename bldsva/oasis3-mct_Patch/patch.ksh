#! /usr/bin/ksh

################  User Setup section  ###################

host=gnu
export SVAROOT=/home/pshrestha/coup_oas
patchroot=${SVAROOT}/bldsva/oasis3-mct_Patch

#########################################################
main(){
echo "########################################"
echo "start patching oasis3-mct" 
echo "########################################"
echo " "

echo "cd in sva-root: ${SVAROOT}"
cd ${SVAROOT}

echo "remove exisiting oasis"
rm -rf oasis3

echo "load current oasis-mct 2.0 from svn"
svn checkout http://oasis3mct.cerfacs.fr/svn/branches/OASIS3-MCT_2.0_branch/oasis3-mct
#
echo "rename to oasis3"
mv oasis3-mct oasis3

echo "cd back to patch-root  ${patchroot}"
cd ${patchroot}

echo "run the script that renames the interface of mct in oasis..." 
echo "in order to run simultaniously with existing mct in clm"
./rename_mct_interface.ksh

echo "theres something wrong with the write_grids & write_corner interfaces..."
echo "replace it by a bugfix"
cp mod_oasis_grid.F90 ${SVAROOT}/oasis3/lib/psmile/src/

echo "cd in oasis make folder"
cd ${SVAROOT}/oasis3/util/make_dir

#add support for other machines here
if [[ ${host} = "juqueen" ]]; then  
	echo "copy juqueen makefile from patchfolder in the oasis make folder"
	cp ${patchroot}/make.xlf_juqueen_oa3 .
	echo "set the juqueen makefile active in make.inc"
 	sed -i "s@include.*@include ${SVAROOT}/oasis3/util/make_dir/make.xlf_juqueen_oa3@" make.inc
fi
if [[ ${host} = "gnu" ]]; then
        echo "copy gnu makefile from patchfolder in the oasis make folder"
        cp ${patchroot}/make.gnu_cluma_oa3 .
        echo "set the gnu makefile active in make.inc"
        sed -i "s@include.*@include ${SVAROOT}/oasis3/util/make_dir/make.gnu_cluma_oa3@" make.inc
fi
#
#
#


echo "cd back to patch-root  ${patchroot}"
cd ${patchroot}

echo "change all includes from mod_prism_proto in mod_prism"
sed -i "s/mod_prism_proto/mod_prism/" ${SVAROOT}/clm/src/oas3/oas_clm_vardef.F90
sed -i "s/mod_prism_proto/mod_prism/" ${SVAROOT}/parflow/pfsimulator/amps/oas3/oas_pfl_vardef.F90
sed -i "s/mod_prism_proto/mod_prism/" ${SVAROOT}/cosmo/src/oas3/oas_cos_vardef.F90

echo "remove all obsolete includes"
sed -i "s/USE mod_prism.*//" ${SVAROOT}/clm/src/oas3/oas_clm_rcv.F90
sed -i "s/USE mod_prism.*//" ${SVAROOT}/clm/src/oas3/oas_clm_snd.F90
sed -i "s/USE mod_prism.*//" ${SVAROOT}/cosmo/src/oas3/oas_cos_rcv.F90
sed -i "s/USE mod_prism.*//" ${SVAROOT}/cosmo/src/oas3/oas_cos_snd.F90
sed -i "s/USE mod_prism.*//" ${SVAROOT}/cosmo/src/oas3/oas_cos_define.F90
sed -i "s/USE mod_prism.*//" ${SVAROOT}/parflow/pfsimulator/amps/oas3/oas_pfl_define.F90
//sed -i "s/USE mod_prism.*//" ${SVAROOT}/bldsva/ibm_changes/nrw_data-set/oas_pfl_define.F90
//sed -i "s/USE mod_prism.*//" ${SVAROOT}/bldsva/ibm_changes/idealized-test-case_data-set/oas_pfl_define.F90
sed -i "s/USE mod_prism.*//" ${SVAROOT}/parflow/pfsimulator/amps/oas3/oas_pfl_snd.F90
sed -i "s/USE mod_prism.*//" ${SVAROOT}/parflow/pfsimulator/amps/oas3/oas_pfl_rcv.F90

echo "prism_get_freq is not supported at the moment (copy module with own implementation)"
#sed -i "s/CALL prism_get_freq.*/cplfreq = 900/" ${SVAROOT}/cosmo/src/oas3/oas_cos_define.F90
cp mod_oasis_method.F90 ${SVAROOT}/oasis3/lib/psmile/src/
sed -i "/oasis_get_debug/a   use mod_oasis_method ,only: prism_get_freq            => oasis_get_freq" ${SVAROOT}/oasis3/lib/psmile/src/mod_prism.F90  	# critical anchor

echo "more difficulties for clm"
cp oas_clm_define.F90 ${SVAROOT}/clm/src/oas3/

echo "adjust PSMILE libs in build_oas3"
sed -i "s/setenv LIBOASIS.*//" ${SVAROOT}/bldsva/build_oas3
sed -i 's@setenv PSMILE_INCDIR.*@setenv PSMILE_INCDIR "-I${LIBBUILD}/psmile.${CHAN}"@' ${SVAROOT}/bldsva/build_oas3
sed -i 's@setenv LIBPSMILE.*@setenv LIBPSMILE "${ARCHDIR}/lib/libpsmile.${CHAN}.a ${ARCHDIR}/lib/libmct.a ${ARCHDIR}/lib/libmpeu.a ${ARCHDIR}/lib/libscrip.a"@' ${SVAROOT}/bldsva/build_oas3

echo "adjust jobscript creation"
sed -i "s/set nproc_oas = 1/set nproc_oas = 0/" ${SVAROOT}/bldsva/build_oas3
sed -i "s/set bgs_oas.*/set bgs_oas = 0/" ${SVAROOT}/bldsva/build_oas3
sed -i 's@set exec_line = "$MPIRUN -p $nppn -n $mpitasks --mapping mapfile.txt $flags_runjob : ./oasis3.MPI1.x "@set exec_line = "$MPIRUN -p $nppn -n $mpitasks --mapping mapfile.txt $flags_runjob : ./clm "@' ${SVAROOT}/bldsva/build_oas3
sed -i 's@ -np $nproc_oas ./oasis3.MPI1.x :@@' ${SVAROOT}/bldsva/build_oas3

echo "rename BILINEA rmp-files to BILINEAR"
sed -i '/cp $forcingdir\/oasis3\/\* $rundir\/\./a foreach x (*BILINEA*)' ${SVAROOT}/bldsva/build_oas3         	# critical anchor
sed -i '/foreach x (\*BILINEA\*)/a mv $x \`echo $x | sed "s\/BILINEA\/BILINEAR\/"\` ; end' ${SVAROOT}/bldsva/build_oas3

echo "adjust mapfile creation on juqueen"
sed -i 's/cat mh0/cat/' ${SVAROOT}/bldsva/build_oas3
sed -i 's/cat m0/cat/' ${SVAROOT}/bldsva/build_oas3

echo "cd in ${SVAROOT}/bldsva/data_oas3/"
cd ${SVAROOT}/bldsva/data_oas3/

echo "changing cosmo resolution in namcoupl from X / Y to X-8 / Y-8 since ghostcells can't be counted"
sed -i '/endif    #$compile_option == 5/a sed -i "s@$ngcosx $ngcosy@`expr $ngcosx - 8` `expr $ngcosy - 8`@g" \*namcouple\*' ${SVAROOT}/bldsva/build_oas3      # critical anchor

# to build oasis3-mct standalone:

#export MPIDIR=/homea/slts/slts06/local/juqueen/mpi
#export NETCDF=/homea/slts/sltsP3cp/local_juqueen/NetCDF3.6.3_xlc_xlf_64
#export PREP_C=""
#export OPT_C=""
#make -f TopMakefileOasis3 realclean
#make -f TopMakefileOasis3 oasis3_psmile

echo " "
echo "########################################"
echo "finished patching oasis3-mct" 
echo "########################################"

}


logDate=$(date +"%Y%m%d%H%M")
main 2>&1 | tee -a patch${logDate}.log

