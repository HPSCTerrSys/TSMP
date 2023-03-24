#! /bin/bash

rmpfile="rmp_gclm_to_gicon_DISTWGT.nc"
echo "Preparing: $rmpfile"

# CONSERV
#ncap2 -s 'defdim("ys",1);gicon_lat($dst_grid_size,$ys)=dst_grid_center_lat;gicon_lon($dst_grid_size,$ys)=dst_grid_center_lon;gicon_imask($dst_grid_size,$ys)=dst_grid_imask;gicon_clat($dst_grid_size,$ys)=dst_grid_corner_lat;gicon_clon($dst_grid_size,$ys)=dst_grid_corner_lon' rmp_gicon_to_gclm_DISTWGT_FRACAREA.nc d0000.nc
#ncap2 -s 'defdim("yd",1);gclm_lat($src_grid_size,$yd)=src_grid_center_lat;gclm_lon($src_grid_size,$yd)=src_grid_center_lon;gclm_imask($src_grid_size,$yd)=src_grid_imask;gclm_clat($src_grid_size,$ys)=src_grid_corner_lat;gclm_clon($src_grid_size,$ys)=src_grid_corner_lon' d0000.nc d0001.nc
#ncks -x -v dst_grid_center_lat,dst_grid_center_lon,dst_grid_imask,src_grid_center_lat,src_grid_center_lon,src_grid_imask d0001.nc d0002.nc # 2 more tp delete

# DISTWGT
ncap2 -s 'defdim("ys",1);gicon_lat($dst_grid_size,$ys)=dst_grid_center_lat;gicon_lon($dst_grid_size,$ys)=dst_grid_center_lon;gicon_imask($dst_grid_size,$ys)=dst_grid_imask' $rmpfile d0000.nc
ncap2 -s 'defdim("yd",1);gclm_lat($src_grid_size,$yd)=src_grid_center_lat;gclm_lon($src_grid_size,$yd)=src_grid_center_lon;gclm_imask($src_grid_size,$yd)=src_grid_imask' d0000.nc d0001.nc
ncks -x -v dst_grid_center_lat,dst_grid_center_lon,dst_grid_imask,src_grid_center_lat,src_grid_center_lon,src_grid_imask d0001.nc d0002.nc

# fix mask value for OASIS3-MCT
ncap2 -s 'where(gclm_imask==1) gclm_imask=2' d0002.nc m0000.nc
ncap2 -s 'where(gclm_imask==0) gclm_imask=1' m0000.nc m0001.nc
ncap2 -s 'where(gclm_imask==2) gclm_imask=0' m0001.nc m0002.nc
ncap2 -s 'where(gicon_imask==1) gicon_imask=2' m0002.nc m0003.nc
ncap2 -s 'where(gicon_imask==0) gicon_imask=1' m0003.nc m0004.nc
ncap2 -s 'where(gicon_imask==2) gicon_imask=0' m0004.nc m0005.nc

ncrename -v gicon_lat,gicon.lat m0005.nc r0000.nc
ncrename -v gclm_lat,gclm.lat r0000.nc r0001.nc
ncrename -v gicon_lon,gicon.lon r0001.nc r0002.nc
ncrename -v gclm_lon,gclm.lon r0002.nc r0003.nc
ncrename -v gicon_imask,gicon.msk r0003.nc r0004.nc
ncrename -v gclm_imask,gclm.msk r0004.nc r0005.nc

mv -f r0005.nc $rmpfile
rm -f d00??.nc m00??.nc r00??.nc
