#!/bin/csh
echo "----------------------------------------"
echo "Usage: ./archive_smresponse.csh clm-cos vg03"
echo "----------------------------------------"
echo " "
#
set modID  = $1
set smID   = $2
set outdir = ${modID}_${smID}
set indir  = $WORK/run_${outdir}
#
if ( -d ./$outdir ) then
  echo Directory already exists! $dir
  exit 1
else
  mkdir $outdir
endif
#
cd $outdir
#
cp $indir/lmrun_uc .
#
mkdir cosout
cd cosout
mkdir ivr
mv $indir/cosout/lfff00000000c.nc ivr/
mv $indir/cosout/lfff* . 
cd ..
#
if (${modID} != "cos") then
  cp $indir/lnd.stdin .  
  mkdir clmout
  cd clmout   
  mv $indir/clmoas.clm2.h0.* . 
endif
#
cd 

exit 0
