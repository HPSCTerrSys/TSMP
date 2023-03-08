#! /bin/ksh

static_files=tsmp_idealscal/input_1200600

StartDate="2008-05-08 00"
InitDate="2008-05-08 00"
  
DumpCLM=1
DumpCOS=1
DumpPFL=1
  
Runhours=3

gxCLM=1200
gyCLM=1200
dtCLM=900
resCLM="1200x1200"

gxCOS=600
gyCOS=600
dtCOS=10
nboundlinesCOS=4

gxPFL=1200
gyPFL=1200
dtPFL=0.25
runnamePFL="ideal"
basePFL=0.0025

freq1OAS=900
freq2OAS=900

deltaobs=1

