#!/bin/csh
foreach instance (`seq 1 90`)
  set model_dir = `printf $WORK/rundart%02d $instance`
 ./archive_hist.csh $model_dir 48
end
exit 0
