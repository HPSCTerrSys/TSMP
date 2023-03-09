#! /bin/ksh

static_files=tsmp_idealscal/input_600300

StartDate="2008-05-08 00"
InitDate="2008-05-08 00"
  
DumpCLM=1
DumpCOS=1
DumpPFL=1
  
Runhours=3

gxCLM=600
gyCLM=600
dtCLM=900
resCLM="0600x0600"

gxCOS=300
gyCOS=300
dtCOS=10
nboundlinesCOS=4

gxPFL=600
gyPFL=600
dtPFL=0.25
runnamePFL="ideal"
basePFL=0.0025

freq1OAS=900
freq2OAS=900

deltaobs=1

indPFL=
pfsolPFL=
inipress= 
slope= 
