#! /bin/ksh

initSetup(){
 # sh $rootdir/bldsva/machines/${platform}/loadenvs.Intel
  source /gpfs/software/juwels/lmod/lmod/init/ksh
#  module load GCC ParaStationMPI CDO NCO
 module load GCC OpenMPI CDO NCO

  defaultFDOAS="$rootdir/bldsva/setups/icon-dcycl-nwp/input"
  defaultFDCLM="$rootdir/bldsva/setups/icon-dcycl-nwp/input"
  defaultFDPFL="$rootdir/bldsva/setups/icon-dcycl-nwp/input"
  defaultFDICON="$rootdir/bldsva/setups/icon-dcycl-nwp/input"

  #defaultFDOAS="/work/jibg34/jibg3401/data/TestCases/icon-ccs"
  #defaultFDCLM="/work/jibg34/jibg3401/data/TestCases/icon-ccs"
  #defaultFDPFL="/work/jibg34/jibg3401/data/TestCases/icon-ccs"
  #defaultFDICON="/work/jibg34/jibg3401/data/TestCases/icon-ccs"

#  defaultNLICON="$rootdir/bldsva/setups/icon-ccs-snwp/icon_master.namelist $rootdir/bldsva/setups/icon-ccs-snwp/NAMELIST_icon"
  defaultNLICON="$rootdir/bldsva/setups/icon-diurnal-snwp/icon_master.namelist $rootdir/bldsva/setups/icon-diurnal-snwp/NAMELIST_icon"
  defaultNLCLM=$rootdir/bldsva/setups/icon-ccs-snwp/lnd.stdin 
  defaultNLPFL=$rootdir/bldsva/setups/icon-ccs-snwp/coup_oas.tcl 

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

  defaultNppn=48
  defaultICONProc=92 #1456
  defaultCLMProcX=2 #4
  defaultCLMProcY=2 #4
  defaultPFLProcX=8
  defaultPFLProcY=8

  defaultStartDate="2008-07-08 00" #"2008-05-08 06"
  defaultInitDate="2008-07-08 00" #"2008-05-08 06"
  defaultRunhours=72 #5 #12
  defaultEndDate="$(date -d "${defaultInitDate[0]} $defaultRunhours hours" "+%Y-%m-%d %H")"

  defaultDumpCLM=1
  defaultDumpPFL=0 #1

  # local folder for remapping
  rmp_folder="$rootdir/bldsva/setups/icon-ccs-snwp/create_remap"
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
  dt_clm=60 # 900 #200 #3000 #45 # in seconds
  
  # setup parflow
  gx_pfl=${gx_clm}
  gy_pfl=${gx_clm}
  dt_pfl=0.8333333 # 5min 0.0125 # 45s in hours
  pflrunname="ccs"
  base_pfl=0.0025

  # setup icon-clm remapping
  gicon_lon="$(cdo infov ${gicon} 2>> $err_file | grep lon_cell_centre)"
  gicon_lon=($gicon_lon)
  gx_icon=${gicon_lon[5]}
  dt_icon=10. #0.5

#  # load needed modules for preprocessing
#  module load Intel/2019.3.199-GCC-8.3.0 IntelMPI/2018.5.288 CDO/1.9.5 NCO/4.7.9

  cd $rmp_folder

  if [[ $withPFL == "true" && $withCLM == "true" ]]; then
  # perform pfl-clm remapping
  cdo setgrid,tgrid.txt $gclm clmgrid1.nc >> $log_file 2>> $err_file
  ncks -O -x -v EDGEN,EDGES,EDGEW,EDGEE,NUMLON clmgrid1.nc clmgrid.nc >> $log_file 2>> $err_file
  ln -s clmgrid.nc pflgrid.nc >> $log_file 2>> $err_file
  comment "  Generating weight file: clm2pfl DISTWGT"
  cdo -P 8 gendis,pflgrid.nc pflgrid.nc rmp_gclm_to_gpfl_DISTWGT.nc >> $log_file 2>> $err_file
  check
  comment "  Generating weight file: pfl2clm DISTWGT"
  cdo -P 8 gendis,pflgrid.nc pflgrid.nc rmp_gpfl_to_gclm_DISTWGT.nc >> $log_file 2>> $err_file
  check
  fi

  if [[ $withICON == "true" && $withCLM == "true" ]]; then
  # SPo perform clmgrid for icon-clm only
  cdo setgrid,tgrid.txt $gclm clmgrid1.nc >> $log_file 2>> $err_file
  ncks -O -x -v EDGEN,EDGES,EDGEW,EDGEE,NUMLON clmgrid1.nc clmgrid.nc >> $log_file 2>> $err_file
  # perform icon-clm remapping
  ncks -O -x -v edge_of_cell,cell_domain_id,adjacent_cell_of_edge,cell_index,cell_no_of_domains,cells_of_vertex,child_cell_id,child_cell_index,child_edge_id,child_edge_index,dual_area_p,edge_cell_distance,edge_index,edge_orientation,edge_parent_type,edge_vert_distance,edge_vertices,edges_of_vertex,end_idx_c,end_idx_e,end_idx_v,index_c_list,index_e_list,index_v_list,meridional_normal_dual_edge,meridional_normal_primal_edge,neighbor_cell_index,orientation_of_normal,parent_cell_index,parent_cell_type,parent_edge_index,parent_vertex_index,refin_c_ctrl,refin_e_ctrl,refin_v_ctrl,start_idx_c,start_idx_e,start_idx_v,vertex_index,vertex_of_cell,vertices_of_vertex,zonal_normal_dual_edge,zonal_normal_primal_edge,cartesian_x_vertices,cartesian_y_vertices,cartesian_z_vertices $gicon icogrid.nc >> $log_file 2>> $err_file
  comment "  Generating weight file: clm2icon DISTWGT"
  cdo -P 8 gendis,clmgrid.nc icogrid.nc rmp_gicon_to_gclm_DISTWGT.nc >> $log_file 2>> $err_file
  check
  comment "  Generating weight file: icon2clm DISTWGT"
  cdo -P 8 gendis,icogrid.nc clmgrid.nc rmp_gclm_to_gicon_DISTWGT.nc >> $log_file 2>> $err_file
  check
  fi

  cd ..

  cplfreq1=`echo "$dt_clm*10" | bc` #9000 #3000 #450  # icon <-> clm, in tenths of seconds
  cplfreq2=3000 #450  # clm <-> pfl, in tenths of seconds

  if [[ $withPFL == "false" && $withICON == "true" ]]; then
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_icon_clm
  fi
  if [[ $withPFL == "true" && $withICON == "false" ]]; then
    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_pfl_clm
  fi
  if [[ $withPFL == "true" && $withICON == "true" ]]; then
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_icon_clm_pfl
  fi

}

finalizeSetup(){
route "${cblue}>> finalizeSetup${cnormal}"
  if [[ $withOAS == "true" ]] then

    if [[ $withOASMCT == "true" ]] then
      comment "   rename oasis3 remapping files (SBr: skip for now!)"
      #for x in $rundir/*BILINEA* ;do 
      #  comment "   rename oasis3 remapping files" 
      #    mv $x $(echo $x | sed "s/BILINEA/BILINEAR/") >> $log_file 2>> $err_file
      #  check 
      #done
    fi  
  fi  

  comment "  copy icon grid"
    ln -s $gicon $rundir/ >> $log_file 2>> $err_file
  check

  if [[ $withOAS == "true" ]] then
    comment "  copy remap-files"
      cp $rmp_folder/rmp* $rundir/ >> $log_file 2>> $err_file
#SPo    rm -r $rmp_folder
    check
  fi

route "${cblue}<< finalizeSetup${cnormal}"
}
