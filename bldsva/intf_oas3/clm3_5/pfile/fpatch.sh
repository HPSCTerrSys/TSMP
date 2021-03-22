#!/bin/bash


ptsmp=/p/home/jusers/ghasemi1/jureca/project/RTDO1/terrsysmp/bldsva/intf_oas3/clm3_5/tsmp
porig=/p/home/jusers/ghasemi1/jureca/project/RTDO1/terrsysmp/clm3_5/src

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
                                   
