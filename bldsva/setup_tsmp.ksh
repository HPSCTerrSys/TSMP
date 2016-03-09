#! /bin/ksh

cblue=$(tput setaf 4)
cnormal=$(tput setaf 9)
cred=$(tput setaf 1)
cgreen=$(tput setaf 2)

getDefaults(){
  def_platform="JURECA"                #This should be correct - change with caution
  def_version="MCT"                 #This should be correct - change with caution
  def_rootdir="/homea/slts/slts06/newtsmp/terrsysmp"    #This should be correct - change with caution
  def_combination="cosclmpfl"           #This should be correct - change with caution
  def_bindir="$def_rootdir/bin/${def_platform}_$date"  #This should be correct - change with caution
  def_rundir=""  			#Will be set to $rootdir/run/$platform_$date if empty
  def_oasdir=""                         #Will be set to $rootdir/XXX_$platform_$date if empty
  def_cosdir=""
  def_clmdir=""
  def_pfldir=""
  # parameters  will be set to tested platform defaults if empty
  nppn=""
  wtime=""
  queue=""
  px_clm=""
  py_clm=""
  px_cos=""
  py_cos=""
  px_pfl=""
  py_pfl=""

  refSetup=""
  
  #profiling ("yes" , "no") - will use machine standard
  def_profiling="no"
}

setDefaults(){
  platform=$def_platform
  version=$def_version
  rootdir=$def_rootdir
  bindir=$def_bindir
  rundir=$def_rundir
  profiling=$def_profiling

  oasdir=$def_oasdir
  clmdir=$def_clmdir
  cosdir=$def_cosdir
  pfldir=$def_pfldir

  log_file=$rootdir/bldsva/log_all.txt
  err_file=$rootdir/bldsva/err_all.txt
  rm -f $log_file $err_file
 
  combination=$def_combination

}

setMachine(){

  if [[ $nppn == "" ]] then ; nppn=$defaultNppn ; fi
  if [[ $wtime == "" ]] then ; wtime=$defaultwtime  ; fi
  if [[ $queue == "" ]] then ; queue=$defaultQ ; fi
  if [[ $px_clm == "" ]] then ; px_clm=$defaultCLMProcX ; fi
  if [[ $py_clm == "" ]] then ; py_clm=$defaultCLMProcY ; fi
  if [[ $px_cos == "" ]] then ; px_cos=$defaultCOSProcX ; fi
  if [[ $py_cos == "" ]] then ; py_cos=$defaultCOSProcY ; fi
  if [[ $px_pfl == "" ]] then ; px_pfl=$defaultPFLProcX ; fi
  if [[ $py_pfl == "" ]] then ; py_pfl=$defaultPFLProcY ; fi

  if [[ $refSetup == "" ]] then ; refSetup=$defaultRefSetup ; fi
}

check(){
 if [[ $? == 0  ]] then
    print "    ... ${cgreen}passed!${cnormal}" 
 else
    print "    ... ${cred}error!!! - aborting...${cnormal}"
    print "See $log_file and $err_file"
    exit 1
  fi
}



  yyyy=2008
  mm=05
  dd=01
  hh=00
  restart=0
  restDate=2015-11-08
  runhours=3



debug(){
print "##############################"
print "            debug:            "
print "##############################"
print "platform=$platform"
print "version=$version"
print "root=$rootdir"
print "mpi path: $mpiPath"
print "silo path: $siloPath"
print "hypre path: $hyprePath"
print "tcl path: $tclPath"
print "grib1 path: $grib1Path"
print "ncdf path: $ncdfPath"
print "profiling: $profiling"
}


#######################################
#               Main
#######################################

  date=`date +%d%m%Y`
  getDefaults
  setDefaults


  #GetOpts definition
  USAGE=$'[-?\n@(#)$Id: TerrSysMP setup script 1.0 - '
  USAGE+=$' date: 10.10.2015 $\n]'
  USAGE+="[-author?Fabian Gasper <f.gasper@fz-juelich.de>]"
  USAGE+="[+NAME?TerrSysMP setup script]"
  USAGE+="[+DESCRIPTION?sets up TSMP run by handling namelists and copying necessary files into a run directory]"
  USAGE+="[b:bash?Bash mode - set command line arguments will overwrite default values (no interactive mode) (This is the default with arguments).]"
  USAGE+="[i:interactive?Interactive mode - command line arguments are ignored and defaults will be overwritten during the interactive session (This is the default without arguments).]"
#  USAGE+="[a:avail?Prints a listing of every machine with available versions. The script will exit afterwards.]"
  USAGE+="[t:tutorial?Prints a tutorial/description on how to add new versions and platforms to this script. The script will exit afterwards.]"
  USAGE+="[r:rootdir?Absolut path to TerrSysMP root directory.]:[path:='$def_rootdir']"
  USAGE+="[B:bindir?Absolut path to bin directory for the builded executables.]:[path:='$def_bindir']"
  USAGE+="[R:rundir?Absolut path to run directory for TSMP.]:[path:='$def_rundir']"
  USAGE+="[v:version?Tagged TerrSysMP version. Note that not every version might be implemented on every machine. Run option -a, --avail to get a listing.]:[version:='$version']{"
  for v in ${!versions[@]} ; do
        str=$(printf "%-20s #%s" "$v" "${versions[$v]}")
        USAGE+="[?$str]"
  done
  USAGE+="}"
  USAGE+="[m:machine?Target Platform. Possible selection:]:[machine:='$platform']{"
  for p in ${!platforms[@]} ; do
        str=$(printf "%-20s #%s" "$p" "${platforms[$p]}")
        USAGE+="[?$str]"
  done
  USAGE+="}"

  USAGE+="[p:profiling?Makes necessary changes to compile with a profiling tool if available.]:[profiling:='$def_profiling']"
#  USAGE+="[C:combination? Combination of component models.]:[combination:='$def_combination']"
#  USAGE+="[W:optoas?Build option for Oasis.]:[optoas:='${def_options["oas"]}']{"
#  USAGE+=$(printf "[?%-12s #%s]" "build" "build from scratch")
#  USAGE+=$(printf "[?%-12s #%s]" "make" "only (resume) make and make install - no make clean and configure")
#  USAGE+=$(printf "[?%-12s #%s]" "configure" "only make clean and configure - no make")
#  USAGE+=$(printf "[?%-12s #%s]" "skip" "no build")
#  USAGE+="}"
#  USAGE+="[X:optcos?Build option for Cosmo.]:[optcos:='${def_options["cos"]}']"
#  USAGE+="[Y:optclm?Build option for CLM.]:[optclm:='${def_options["clm"]}']"
#  USAGE+="[Z:optpfl?Build option for Parflow.]:[optpfl:='${def_options["pfl"]}']"
#  USAGE+="[H:hyprepath?Include Path for Hypre.]:[hyprepath:='$hyprePath']"
#  USAGE+="[S:silopath?Include Path for Silo.]:[silopath:='$siloPath']"
#  USAGE+="[T:tclpath?Include Path for TCL.]:[tclpath:='$tclPath']"
#  USAGE+="[G:grib1path?Include Path for Grib1.]:[grib1path:='$grib1Path']"
#  USAGE+="[M:mpipath?Include Path for MPI.]:[mpipath:='$mpiPath']"
#  USAGE+="[N:ncdfpath?Include Path for NetCDF.]:[ncdfpath:='$ncdfPath']"
  USAGE+=$'\n\n\n\n'


  mode=0
  args=0
  # parsing the command line arguments
  while getopts "$USAGE" optchar ; do
    case $optchar in
    i)  mode=2 ;;
    b)  mode=1 ;;
    m)  #if [[ "${platforms[${OPTARG}]}" != "" ]] then 
                platform=$OPTARG
        #else
        #    print "The selected platform '${OPTARG}' is not available. run '.$call --man' for help"
        #    terminate    
        #fi
        args=1
        ;;
    p)  profling="${OPTARG}" ; args=1 ;;
    v)  if [[ "${!versions[${OPTARG}]}" != ""  ]] then
            version=$OPTARG
        else
            print "The selected version '${OPTARG}' is not available. run '.$call --man' for help"
            terminate
        fi
        args=1 ;;
    a)  listAvailabilities ;;
    t)  listTutorial ;;
    r)  new_rootdir="$OPTARG" ; args=1 ;;
    B)  bindir="$OPTARG" ; args=1 ;;
    R)  rundir="$OPTARG" ; args=1 ;;
    C)  combination="$OPTARG" ; args=1 ;;
#    W)  options+=(["oas"]="$OPTARG") ; args=1 ;;
#    X)  options+=(["cos"]="$OPTARG") ; args=1 ;;
#    Y)  options+=(["clm"]="$OPTARG") ; args=1 ;;
#    Z)  options+=(["pfl"]="$OPTARG") ; args=1 ;;
#    M)  mpiPath="$OPTARG" ; args=1 ;;
#    N)  ncdfPath="$OPTARG" ; args=1 ;;
#    G)  grib1Path="$OPTARG" ; args=1 ;;
#    T)  tclPath="$OPTARG" ; args=1 ;;
#    H)  hyprePath="$OPTARG" ; args=1 ;;
#    S)  siloPath="$OPTARG" ; args=1 ;;
    esac
  done

  if [[ $rundir == "" ]] then
     rundir=$rootdir/run/${platform}_${date}
  else
     rundir="$rundir/${platform}_${date}"
  fi
  mkdir -p $rundir
  rm -rf $rundir/*

  . ${rootdir}/bldsva/machines/${platform}/build_interface_${platform}.ksh

  getMachineDefaults
  setMachine

  withOAS="false"
  withCOS="false"
  withPFL="false"
  withCLM="false"
  withOASMCT="false"
  withPCLM="false"

  case "$combination" in *clm*) withCLM="true" ;; esac
  case "$combination" in *cos*) withCOS="true" ;; esac
  case "$combination" in *pfl*) withPFL="true" ;; esac
  if [[ $withCLM == "true" && ( $withCOS == "true" || $withPFL == "true" )  ]]; then
    withOAS="true"
    case "$version" in *MCT*) withOASMCT="true" ;; esac
  fi



  if [[ $oasdir == "" && withOAS  ]] then ;  oasdir=$rootdir/oasis3${platform}_${date} ; fi
  if [[ $oasdir == "" && withOASMCT  ]] then ;  oasdir=$rootdir/oasis3-mct${platform}_${date} ; fi
  if [[ $cosdir == "" ]] then ;  cosdir=$rootdir/cosmo${platform}_${date} ; fi
  if [[ $clmdir == "" ]] then ;  clmdir=$rootdir/clm${platform}_${date} ; fi
  if [[ $pfldir == "" ]] then ;  pfldir=$rootdir/parflow${platform}_${date} ; fi


  nproc_oas=0
  nproc_pfl=0
  nproc_clm=0
  nproc_cos=0
  if [[ $withPFL == "true" ]] ; then ; nproc_pfl=$((${px_pfl}*${py_pfl})) ;fi
  if [[ $withCOS == "true" ]] ; then ; nproc_cos=$((${px_cos}*${py_cos})) ;fi
  if [[ $withCLM == "true" ]] ; then ; nproc_clm=$((${px_clm}*${py_clm})) ;fi
  if [[ $withOAS == "true" ]] ; then ; nproc_oas=1 ;fi
  if [[ $withOASMCT == "true" ]] ; then ; nproc_oas=0 ;fi

#  start setup


  setup=$refSetup

  if [[ $setup != "" ]] ;then ; . ${rootdir}/bldsva/setups/$setup/${setup}_setup.ksh ;fi

  initSetup

  if [[ $withCLM == "true" ]] ; then ;. ${rootdir}/bldsva/intf_oas3/clm3_5/arch/${platform}/build_interface_clm3_5_${platform}.ksh ;  setup_clm ; cp $bindir/clm $rundir ; fi
  if [[ $withCOS == "true" ]] ; then ;. ${rootdir}/bldsva/intf_oas3/cosmo4_21/arch/${platform}/build_interface_cosmo4_21_${platform}.ksh ;  setup_cos ; cp $bindir/lmparbin_pur $rundir ; fi
  if [[ $withPFL == "true" ]] ; then ;. ${rootdir}/bldsva/intf_oas3/parflow/arch/${platform}/build_interface_parflow_${platform}.ksh ;  setup_pfl ; cp $bindir/parflow $rundir ; fi
  if [[ $withOAS == "true" ]] ; then ;. ${rootdir}/bldsva/intf_oas3/oasis3-mct/arch/${platform}/build_interface_oasis3-mct_${platform}.ksh ;  setup_oas ; fi

  finalizeSetup


  createRunscript

  cp $err_file $rundir
  cp $log_file $rundir

#  debug

