#! /bin/ksh

static_files=$rootdir/tsmp_idealrtd/input

StartDate="2015-08-07 00"
InitDate="2015-08-07 00"
  
DumpCLM=1
DumpCOS=1
DumpPFL=1
  
Runhours=24

gxCLM=54
gyCLM=24
dtCLM=18
resCLM="0024x0054"

gxCOS=60
gyCOS=30
dtCOS=18
nboundlinesCOS=3

gxPFL=54
gyPFL=24
dtPFL=0.005
runnamePFL="idealRTD"
basePFL=0.0025

freq1OAS=18
freq2OAS=18

deltaobs=1

indPFL=
indPFL2=$indPFL
pfsolPFL=
inipress=*.pfb
slope=

rtimeFactor=1   
