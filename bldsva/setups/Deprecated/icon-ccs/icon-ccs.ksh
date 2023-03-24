#! /bin/ksh

static_files=/../../../../../../../../../p/project/cslts/local/data/terrsysmp/icon-ccs


StartDate="2008-05-08 06"
InitDate="2008-05-08 06"
  
DumpCLM=1
DumpCOS=1
DumpPFL=1
  
Runhours=1

gxCLM=300
gyCLM=300
dtCLM=900
resCLM="300x300"

gxCOS=15472
gyCOS=
dtCOS=
nboundlinesCOS=

gxPFL=300
gyPFL=300
dtPFL=0.25
runnamePFL="ccs"
basePFL=0.0025

freq1OAS=450
freq2OAS=450

deltaobs=1

indPFL=parflow_436x424x15_cosmomask_indicator_FAOonly.sa
indPFL2=pfl_ind.sa 					#Name change in cordex setup 
pfsolPFL=geom_cordex0.11_436x424.pfsol
inipress=
slope=slope*.sa  

