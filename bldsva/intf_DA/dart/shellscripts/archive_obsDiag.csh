#!/bin/csh

# usage: ./archive_obsDiag.csh archiveDir numDays 
#
echo "-------------------------------------------------------------------"
echo "Obs diagnostic files will be archived in the "$WORK" directory"
echo "usage: ./archive_obsDiag.csh archiveDir numDays"
echo "-------------------------------------------------------------------"
echo " "

set archiveDir = $1
set numDays = $2
set dir = "obsDIAG"
cd $archiveDir

if ( -d ./$dir ) then
  echo Directory already exists! $dir
else
  mkdir $dir
endif

foreach iday (`seq 1 $numDays`)
  set rundir = `printf $WORK/rundart%02d $iday`
  set outputdir = $rundir:t
  cp $rundir/obs_seq.final $dir/"obs_seq_"$outputdir".final"
  echo "Processed diagnostic " $rundir/obs_seq.final
  echo $archiveDir/$dir/"obs_seq_"$outputdir".final" >> obs_file_list.txt
end
exit 0
