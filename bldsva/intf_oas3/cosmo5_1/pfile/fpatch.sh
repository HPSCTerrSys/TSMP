#!/bin/bash

curr_dir=`pwd`
TSMP_root=`echo "$curr_dir" | sed 's@/bldsva/intf_oas3/cosmo5_1/pfile@@'` #remove bldsva/intf_oas3/cosmo5_1/pfile to get rootpath

# ptsmp: the directory contain the changed source files for coupling
ptsmp=$curr_dir
# porig: the directory of the original cosmo5_1 source code (the fresh copy located in the TSMP root directory)
porig=$TSMP_root/cosmo5_1

echo "porig::" $porig
echo "ptsmp::" $ptsmp

ltsmp=$ptsmp/*90
#ltsmp1=$(echo $ltsmp| sed -e 's/\s\s*/\n/g')
#lorig=$(find $opth -name '*.F90')                  
#echo $lorig | sed -e 's/\s\s*/\n/g'
for il in $ltsmp
do
   ntsmp=`basename $il`
   norig=$(find $porig -name $ntsmp)
   echo $ntsmp
   echo $norig
   diff -u $norig  $il > patch_$ntsmp.diff
done
                                   
