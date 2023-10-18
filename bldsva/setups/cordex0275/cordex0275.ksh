#! /bin/ksh

static_files=$rootdir/tsmp_eur0275/static

StartDate="2020-01-01 12"
InitDate="2020-01-01 12"
  
DumpCLM=1
DumpCOS=1
DumpPFL=1
  
Runhours=1

gxCLM=1592
gyCLM=1544
dtCLM=900
resCLM="1592x1544"

gxCOS=1600
gyCOS=1552
dtCOS=25
nboundlinesCOS=4

gxPFL=1592
gyPFL=1544
dtPFL=0.25
runnamePFL="cordex0.0275"
basePFL=0.0025

freq1OAS=900
freq2OAS=900

deltaobs=1

indPFL=EUR-0275_TSMP_FZJ-IBG3_CLMPFLDomain_1592x1544_INDICATOR_regridded_rescaled_SoilGrids250-v2017_BGR3_alv.sa
indPFL2=EUR-0275_TSMP_FZJ-IBG3_CLMPFLDomain_1592x1544_INDICATOR_regridded_rescaled_SoilGrids250-v2017_BGR3_alv.sa				#Name change in cordex setup 
pfsolPFL=PfbMask4SolidFile.pfsol
inipress=
slope=EUR-0275_TSMP_FZJ-IBG3_CLMPFLDomain_1592x1544_?SLOPE_TPS_HydroRIVER_sea_streams_corr.sa 

rtimeFactor=1

