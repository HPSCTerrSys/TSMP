#! /bin/bash

cat > tgrid.txt <<EOF
gridtype  = lonlat
xsize     = 300 
ysize     = 300 
xfirst    = -3.166136
xinc      = 0.02097122
yfirst    = -0.1739269
yinc      = 0.00115951266666667
EOF

cdo setgrid,tgrid.txt griddata_0300x0300.nc clmgrid1.nc
ncks -x -v EDGEN,EDGES,EDGEW,EDGEE,NUMLON clmgrid1.nc clmgrid.nc

ncks -x -v edge_of_cell,cell_domain_id,adjacent_cell_of_edge,cell_index,cell_no_of_domains,cells_of_vertex,child_cell_id,child_cell_index,child_edge_id,child_edge_index,dual_area_p,edge_cell_distance,edge_index,edge_orientation,edge_parent_type,edge_vert_distance,edge_vertices,edges_of_vertex,end_idx_c,end_idx_e,end_idx_v,index_c_list,index_e_list,index_v_list,meridional_normal_dual_edge,meridional_normal_primal_edge,neighbor_cell_index,orientation_of_normal,parent_cell_index,parent_cell_type,parent_edge_index,parent_vertex_index,refin_c_ctrl,refin_e_ctrl,refin_v_ctrl,start_idx_c,start_idx_e,start_idx_v,vertex_index,vertex_of_cell,vertices_of_vertex,zonal_normal_dual_edge,zonal_normal_primal_edge torus_triangles_192x192_50m.nc icogrid.nc

#echo "Generating clm2icon conservative weights"
#date
#cdo gencon,clmgrid.nc icon.nc clm2icon_weights.nc
#date
#echo "Generating icon2clm conservative weights"
#cdo gencon,ml.nc clmgrid1.nc icon2clm_weights.nc
#date

echo "Generating weight file: clm2icon DISTWGT"
date
cdo -P 8 gendis,clmgrid.nc icogrid.nc rmp_gclm_to_gicon_DISTWGT.nc
date

echo "Generating weight file: icon2clm DISTWGT"
cdo -P 8 gendis,icogrid.nc clmgrid.nc rmp_gicon_to_gclm_DISTWGT.nc
date

rm -f clmgrid.nc clmgrid1.nc
rm -f icogrid.nc
rm -f tgrid.txt
exit 0
