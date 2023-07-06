#! /bin/ksh

static_files=$rootdir/tsmp_eur11_eraint_eval_v2/input

StartDate="2021-06-24 12"
InitDate="2021-06-24 12"
  
DumpCLM=1
DumpCOS=1
DumpPFL=1
  
Runhours=3

gxCLM=436
gyCLM=424
dtCLM=3600
resCLM="436x424"

gxCOS=444
gyCOS=432
dtCOS=60
nboundlinesCOS=4

gxPFL=436
gyPFL=424
dtPFL=1.0
runnamePFL="cordex0.11"
basePFL=0.0025

freq1OAS=3600
freq2OAS=3600

deltaobs=1

indPFL=parflow_436x424x15_cosmomask_indicator_FAOonly.sa
indPFL2=pfl_ind.sa 					#Name change in cordex setup 
pfsolPFL=geom_cordex0.11_436x424.pfsol
inipress=
slope=slope*.sa  

