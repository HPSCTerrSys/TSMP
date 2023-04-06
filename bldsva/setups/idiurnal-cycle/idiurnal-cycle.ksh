#! /bin/ksh

source /p/software/$SYSTEMNAME/lmod/lmod/init/ksh
module load GCC OpenMPI CDO NCO

static_files=tsmp_idiurnal-cycle/input

defaultFDCLM="$rootdir/$static_files/clm"
defaultFDCOS="$rootdir/$static_files/cosmo"
defaultFDOAS="$rootdir/$static_files/oasis3"
defaultFDPFL="$rootdir/$static_files/parflow"
defaultFDICON="$rootdir/$static_files/icon"

defaultNLICON="$rootdir/bldsva/setups/$refSetup/icon_master.namelist $rootdir/bldsva/setups/$refSetup/NAMELIST_icon"
defaultICONProc=316

StartDate="2008-07-08 00" 
InitDate="2008-07-08 00" 

Runhours=72
defaultEndDate="$(date -d "${InitDate[0]} $dRunhours hours" "+%Y-%m-%d %H")"
  
DumpCLM=1
DumpCOS=1
DumpPFL=1

gxCLM=90
gyCLM=90
dtCLM=90
resCLM="0090x0090"

gxCOS=9800
gyCOS=
dt_icon=10.0
nboundlinesCOS=

gxPFL=90
gyPFL=90
dtPFL=0.025 
runnamePFL="idiurnal"
basePFL=0.0025

freq1OAS=90
freq2OAS=90

deltaobs=1

indPFL=
indPFL2=				#Name change in cordex setup 
pfsolPFL=
inipress=
slope=

gicon=$(find $defaultFDICON -name torus_\*.nc)
gclm=$(find $defaultFDCLM -name griddata_\*.nc)

  if [ "$gicon"x == "x" ]; then
    echo "  ERROR: Please provide ICON grid torus_*.nc in $defaultFDICON."
    exit 1
  fi
  if [ "$gclm"x == "x" ]; then
    echo "  ERROR: Please provide CLM grid griddata_*.nc in $defaultFDCLM."
    exit 1
  fi

  rmp_folder="$rootdir/bldsva/setups/icon-ccs/create_remap"
  rm -rf $rmp_folder
  mkdir $rmp_folder

  # setup pfl-clm remapping
  gclm_lon="$(cdo infov ${gclm} 2>> $err_file | grep LONGXY)"
  gclm_lat="$(cdo infov ${gclm} 2>> $err_file | grep LATIXY)"
  gclm_lon=($gclm_lon)
  gclm_lat=($gclm_lat)
  gx_clm=$(echo "scale=0;sqrt(${gclm_lon[5]})"|bc 2>> $err_file)
  gy_clm=$gx_clm  # !!!
  lat1=${gclm_lat[8]}
  lat2=${gclm_lat[10]}
  lon1=${gclm_lon[8]}
  lon2=${gclm_lon[10]}
  cat > $rmp_folder/tgrid.txt <<EOF
  gridtype  = lonlat
  xsize     = ${gx_clm}
  ysize     = ${gy_clm}
  xfirst    = ${lon1}
  xinc      = $(echo "scale=20;(${lon2}-(${lon1}))/${gx_clm}"|bc 2>> $err_file)
  yfirst    = ${lat1}
  yinc      = $(echo "scale=20;(${lat2}-(${lat1}))/${gy_clm}"|bc 2>> $err_file)
EOF
  res_x=$(printf "%04d" "${gx_clm}")
  res_y=$(printf "%04d" "${gy_clm}")
  res="${res_x}x${res_y}"

  # setup icon-clm remapping
  gicon_lon="$(cdo infov ${gicon} 2>> $err_file | grep lon_cell_centre)"
  gicon_lon=($gicon_lon)
  gx_icon=${gicon_lon[5]}


  cd $rmp_folder
  cdo setgrid,tgrid.txt $gclm clmgrid1.nc >> $log_file 2>> $err_file
  ncks -O -x -v EDGEN,EDGES,EDGEW,EDGEE,NUMLON clmgrid1.nc clmgrid.nc >> $log_file 2>> $err_file


  # perform pfl-clm remapping
  ln -s clmgrid.nc pflgrid.nc >> $log_file 2>> $err_file
  comment "  Generating weight file: clm2pfl DISTWGT"
  cdo -P 8 gendis,pflgrid.nc pflgrid.nc rmp_gclm_to_gpfl_DISTWGT.nc >> $log_file 2>> $err_file
  check
  comment "  Generating weight file: pfl2clm DISTWGT"
  cdo -P 8 gendis,pflgrid.nc pflgrid.nc rmp_gpfl_to_gclm_DISTWGT.nc >> $log_file 2>> $err_file
  check



  # perform icon-clm remapping
  ncks -O -x -v edge_of_cell,cell_domain_id,adjacent_cell_of_edge,cell_index,cell_no_of_domains,cells_of_vertex,child_cell_id,child_cell_index,child_edge_id,child_edge_index,dual_area_p,edge_cell_distance,edge_index,edge_orientation,edge_parent_type,edge_vert_distance,edge_vertices,edges_of_vertex,end_idx_c,end_idx_e,end_idx_v,index_c_list,index_e_list,index_v_list,meridional_normal_dual_edge,meridional_normal_primal_edge,neighbor_cell_index,orientation_of_normal,parent_cell_index,parent_cell_type,parent_edge_index,parent_vertex_index,refin_c_ctrl,refin_e_ctrl,refin_v_ctrl,start_idx_c,start_idx_e,start_idx_v,vertex_index,vertex_of_cell,vertices_of_vertex,zonal_normal_dual_edge,zonal_normal_primal_edge,cartesian_x_vertices,cartesian_y_vertices,cartesian_z_vertices $gicon icogrid.nc >> $log_file 2>> $err_file
  comment "  Generating weight file: clm2icon DISTWGT"
  cdo -P 8 gendis,clmgrid.nc icogrid.nc rmp_gicon_to_gclm_DISTWGT.nc >> $log_file 2>> $err_file
  check
  comment "  Generating weight file: icon2clm DISTWGT"
  cdo -P 8 gendis,icogrid.nc clmgrid.nc rmp_gclm_to_gicon_DISTWGT.nc >> $log_file 2>> $err_file
  check


  cd ..
