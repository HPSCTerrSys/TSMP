#! /bin/ksh

static_files=tsmp_idealscal/input_24001200

StartDate="2008-05-08 00"
InitDate="2008-05-08 00"
  
DumpCLM=1
DumpCOS=1
DumpPFL=1
  
Runhours=3

gxCLM=2400
gyCLM=2400
dtCLM=900
resCLM="2400x2400"

gxCOS=1200
gyCOS=1200
dtCOS=10
nboundlinesCOS=4

gxPFL=2400
gyPFL=2400
dtPFL=0.25
runnamePFL="ideal"
basePFL=0.0025

freq1OAS=900
freq2OAS=900

deltaobs=1

