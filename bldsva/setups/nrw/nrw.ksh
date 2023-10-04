#! /bin/ksh

static_files=$rootdir/tsmp_nrw/input

StartDate="2008-05-08 00"
InitDate="2008-05-08 00"
  
DumpCLM=1
DumpCOS=1
DumpPFL=1
  
Runhours=24

gxCLM=300
gyCLM=300
dtCLM=90
resCLM="0300x0300"

gxCOS=150
gyCOS=150
dtCOS=10
nboundlinesCOS=3

gxPFL=300
gyPFL=300
dtPFL=0.025
runnamePFL="rurlaf"
basePFL=0.0025

freq1OAS=90
freq2OAS=90

deltaobs=1

indPFL=rurSoil.sa
indPFL2=$indPFL
pfsolPFL=geom_cordex0.11_436x424.pfsol 
inipress=*.ics.pfb*
slope=*slope.pfb*

rtimeFactor=10

