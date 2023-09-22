#! /bin/ksh

static_files=$rootdir/tsmp_idealscal/input_300150

StartDate="2008-05-08 00"
InitDate="2008-05-08 00"
  
DumpCLM=1
DumpCOS=1
DumpPFL=1
  
Runhours=3

gxCLM=300
gyCLM=300
dtCLM=900
resCLM="0300x0300"

gxCOS=150
gyCOS=150
dtCOS=10
nboundlinesCOS=4

gxPFL=300
gyPFL=300
dtPFL=0.25
runnamePFL="ideal"
basePFL=0.0025

freq1OAS=900
freq2OAS=900

deltaobs=1

indPFL=
indPFL2=$indPFL
pfsolPFL=
inipress= *.pfb
slope= 

rtimeFactor=1
