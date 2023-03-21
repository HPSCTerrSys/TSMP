#! /bin/ksh

static_files=/../../../../../../../../../p/project/cslts/poll1/data/


StartDate="2013-04-24 00"
InitDate="2013-04-24 00"
  
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
runnamePFL="nrw"
basePFL=0.0025

freq1OAS=900
freq2OAS=900

deltaobs=1

indPFL=parflow_436x424x15_cosmomask_indicator_FAOonly.sa
indPFL2=pfl_ind.sa 					#Name change in cordex setup 
pfsolPFL=geom_cordex0.11_436x424.pfsol
inipress=
slope=slope*.sa  

