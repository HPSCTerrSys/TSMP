#!/bin/csh
#To use for renaming Parflow files with different counter
#What a mesh
set filepath = "/home/pshrestha/terrsysmp/run/CLUMA2_1.2.0MCT_clm-cos-pfl_idealRTD_dart02"
set fileprfx = "rurlaf.out.satur"
set filext   = ".dist"
cd $filepath
foreach inst (`seq 0 1 47`)
  cd $filepath/"tsmp_instance_"$inst
  set ctr = 8
  foreach filename (`ls $fileprfx*$filext`)
    set fname=$filename:r
    set nindx = `echo $ctr | awk '{printf("%05d",$1);}'`
    if ($filext == ".dist") then
      echo $fname$filext $fileprfx.$nindx.pfb$filext
      mv $fname$filext $fileprfx.$nindx.pfb$filext
    else
      echo $fname$filext $fileprfx.$nindx$filext
      mv $fname$filext $fileprfx.$nindx$filext
    endif
    @ ctr = $ctr + 1
  end
  cd ..
end


