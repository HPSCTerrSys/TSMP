#! /bin/ksh

cblue=$(tput setaf 4)
cnormal=$(tput setaf 9)
cred=$(tput setaf 1)
cgreen=$(tput setaf 2)
#combinations_all=("clm" "cos" "pfl" "clmcos" "clmpfl" "clmcospfl")
typeset -A options
typeset -A def_options

getDefaults(){
  def_platform="JURECA" #This should be correct - change with caution
  def_version="1.1.0MCT" #This should be correct - change with caution
  def_rootdir="/homea/slts/slts06/newtsmp/terrsysmp" #This should be correct - change with caution
  combinations_limitation=("pflclmcos") #This should be correct - change with caution

  def_bindir=""				#Will be set to $rootdir/bin/$platform_$date if empty
  def_oasdir=""				#Will be set to $rootdir/XXX_$platform_$date if empty
  def_cosdir=""
  def_clmdir=""
  def_pfldir=""

  # pathes will be set to tested platform defaults if empty
  mpiPath=""
  ncdfPath=""
  grib1Path=""
  tclPath=""
  hyprePath=""
  siloPath=""


  #compiler optimization
  def_optComp=""   # will be set to platform defaults if empty

  #profiling ("yes" , "no") - will use machine standard
  def_profiling="no"

  # build=build from scratch
  # make=(resume) make - no configure
  # configure= only configure - no make
  # skip=do nothing
  def_options+=(["oas"]="fresh")
  def_options+=(["clm"]="fresh")
  def_options+=(["pfl"]="fresh")
  def_options+=(["cos"]="fresh")

  def_combination="${combinations_limitation[0]}"
}

setDefaults(){
  platform=$def_platform
  version=$def_version
  rootdir=$def_rootdir
  bindir=$def_bindir
  optComp=$def_optComp
  profiling=$def_profiling
  oasdir=$def_oasdir
  clmdir=$def_clmdir
  cosdir=$def_cosdir
  pfldir=$def_pfldir


  log_file=$rootdir/bldsva/log_all.txt
  err_file=$rootdir/bldsva/err_all.txt
  rm -f $log_file $err_file

  options+=(["oas"]=${def_options["oas"]})
  options+=(["cos"]=${def_options["cos"]})
  options+=(["clm"]=${def_options["clm"]})
  options+=(["pfl"]=${def_options["pfl"]})

  combination=$def_combination

}

setMachine(){

  if [[ $mpiPath == "" ]] then ; mpiPath=$defaultMpiPath ; fi
  if [[ $ncdfPath == "" ]] then ; ncdfPath=$defaultNcdfPath  ; fi
  if [[ $grib1Path == "" ]] then ; grib1Path=$defaultGrib1Path ; fi
  if [[ $tclPath == "" ]] then ; tclPath=$defaultTclPath ; fi
  if [[ $hyprePath == "" ]] then ; hyprePath=$defaultHyprePath ; fi
  if [[ $siloPath == "" ]] then ; siloPath=$defaultSiloPath ; fi

  #compiler optimization
  if [[ $mpiPath == "" ]] then ; def_optComp=$defaultOptC ; fi
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


compileClm(){
print "${cblue}> c_compileClm${cnormal}"
print -n "  source clm interface script"
  . ${rootdir}/bldsva/intf_oas3/clm3_5/arch/${platform}/build_interface_clm3_5_${platform}.ksh >> $log_file 2>> $err_file
check
always_clm
if [[ ${options["clm"]} == "skip" ]] ; then ; return  ;fi 
if [[ ${options["clm"]} == "fresh" ]] ; then 
print -n "  backup clm dir to: $clmdir"
    rm -rf $clmdir >> $log_file 2>> $err_file
check
    cp -rf ${rootdir}/clm $clmdir >> $log_file 2>> $err_file
check
    substitutions_clm  
fi
if [[ ${options["clm"]} == "configure" || ${options["clm"]} == "build" || ${options["clm"]} == "fresh" ]] ; then
    configure_clm    
fi

if [[ ${options["clm"]} == "make" || ${options["clm"]} == "build" || ${options["clm"]} == "fresh" ]] ; then
    make_clm
fi
print "${cblue}< c_compileClm${cnormal}"
}

compileCosmo(){
print "${cblue}> c_compileCosmo${cnormal}"
print -n "  source cos interface script"
  . ${rootdir}/bldsva/intf_oas3/cosmo4_21/arch/${platform}/build_interface_cosmo4_21_${platform}.ksh >> $log_file 2>> $err_file
check
always_cos
if [[ ${options["cos"]} == "skip" ]] ; then ; return  ;fi 
if [[ ${options["cos"]} == "fresh" ]] ; then 
print -n "  backup cos dir to: $cosdir"
    rm -rf $cosdir >> $log_file 2>> $err_file
check
    cp -rf ${rootdir}/cosmo $cosdir >> $log_file 2>> $err_file
check
    substitutions_cos  
fi
if [[ ${options["cos"]} == "configure" || ${options["cos"]} == "build" || ${options["cos"]} == "fresh" ]] ; then
    configure_cos    
fi

if [[ ${options["cos"]} == "make" || ${options["cos"]} == "build" || ${options["cos"]} == "fresh" ]] ; then
    make_cos
fi
print "${cblue}< c_compileCosmo${cnormal}"
}

compileOasis(){
print "${cblue}> c_compileOasis${cnormal}"
print -n "  source oas interface script"
  . ${rootdir}/bldsva/intf_oas3/oasis3-mct/arch/${platform}/build_interface_oasis3-mct_${platform}.ksh >> $log_file 2>> $err_file
check
always_oas
if [[ ${options["oas"]} == "skip" ]] ; then ; return  ;fi 
if [[ ${options["oas"]} == "fresh" ]] ; then 
print -n "  backup oas dir to: $oasdir"
    rm -rf $oasdir >> $log_file 2>> $err_file
check
    cp -rf ${rootdir}/oasis3-mct $oasdir >> $log_file 2>> $err_file
check
    substitutions_oas  
fi
if [[ ${options["oas"]} == "configure" || ${options["oas"]} == "build" || ${options["oas"]} == "fresh" ]] ; then
    configure_oas    
fi

if [[ ${options["oas"]} == "make" || ${options["oas"]} == "build" || ${options["oas"]} == "fresh" ]] ; then
    make_oas
fi
print "${cblue}< c_compileOasis${cnormal}"
}

compileParflow(){
print "${cblue}> c_compileParflow${cnormal}"
print -n "  source pfl interface script"
  . ${rootdir}/bldsva/intf_oas3/parflow/arch/${platform}/build_interface_parflow_${platform}.ksh >> $log_file 2>> $err_file
check
always_pfl
if [[ ${options["pfl"]} == "skip" ]] ; then ; return  ;fi 
if [[ ${options["pfl"]} == "fresh" ]] ; then 
print -n "  backup pfl dir to: $pfldir"
    rm -rf $pfldir >> $log_file 2>> $err_file
check
    cp -rf ${rootdir}/parflow $pfldir >> $log_file 2>> $err_file
check
    substitutions_pfl
fi
if [[ ${options["pfl"]} == "configure" || ${options["pfl"]} == "build" || ${options["pfl"]} == "fresh" ]] ; then
    configure_pfl    
fi

if [[ ${options["pfl"]} == "make" || ${options["pfl"]} == "build" || ${options["pfl"]} == "fresh" ]] ; then
    make_pfl
fi
print "${cblue}< c_compileParflow${cnormal}"
}


runCompilation(){
# run model compilation
  if [[ $withOAS == "true" ]] ; then ; compileOasis ; fi # must be called first bc of dependecies

  if [[ $withCLM == "true" ]] ; then ; compileClm ; fi
  if [[ $withCOS == "true" ]] ; then ; compileCosmo ; fi
  if [[ $withPFL == "true" ]] ; then ; compileParflow ; fi
}

debug(){
print "##############################"
print "            debug:            "
print "##############################"
print "platform=$platform"
print "version=$version"
print "root=$rootdir"
print "possible combinations=${combinations_limitation[@]}"
print "mpi path: $mpiPath"
print "silo path: $siloPath"
print "hypre path: $hyprePath"
print "tcl path: $tclPath"
print "grib1 path: $grib1Path"
print "ncdf path: $ncdfPath"
print "optComp: $optComp"
print "profiling: $profiling"
print "combination: $combination"
}


terminate(){
 print ""
 print "Terminating $call. No changes were made..."
 exit 0
}

warning(){
 print "Warning!!! - $wmessage"
 print "This could lead to errors or wrong behaviour. Would you like to continue anyway?"
 PS3="Your selection(1-2)?"
 select ret in "yes" "no, exit"
 do
   if [[ -n $ret ]]; then
      case $ret in
         "yes") break ;;
         "no, exit") terminate ;;
      esac
      break
   fi
 done
}



#######################################
#		Main
#######################################


  getDefaults
  setDefaults


  #GetOpts definition
  USAGE=$'[-?\n@(#)$Id: TerrSysMP build script 1.0 - '
  USAGE+=$' date: 07.09.2015 $\n]'
  USAGE+="[-author?Fabian Gasper <f.gasper@fz-juelich.de>]"
  USAGE+="[+NAME?TerrSysMP build script]"
  USAGE+="[+DESCRIPTION?builds TSMP based on decisions for included modules]"
  USAGE+="[b:bash?Bash mode - set command line arguments will overwrite default values (no interactive mode) (This is the default with arguments).]"
  USAGE+="[i:interactive?Interactive mode - command line arguments are ignored and defaults will be overwritten during the interactive session (This is the default without arguments).]"
  USAGE+="[a:avail?Prints a listing of every machine with available versions. The script will exit afterwards.]"
  USAGE+="[t:tutorial?Prints a tutorial/description on how to add new versions and platforms to this script. The script will exit afterwards.]"
  USAGE+="[r:rootdir?Absolut path to TerrSysMP root directory.]:[path:='$def_rootdir']"
  USAGE+="[B:bindir?Absolut path to bin directory for the builded executables.]:[path:='$def_bindir']"
   
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
  USAGE+="[o:optimization?Compiler optimisation flags.]:[optimization:='$optComp']"
  USAGE+="[C:combination? Combination of component models.]:[combination:='$def_combination']"
  USAGE+="[W:optoas?Build option for Oasis.]:[optoas:='${def_options["oas"]}']{"
  USAGE+=$(printf "[?%-12s #%s]" "build" "build from scratch")
  USAGE+=$(printf "[?%-12s #%s]" "make" "only (resume) make and make install - no make clean and configure")
  USAGE+=$(printf "[?%-12s #%s]" "configure" "only make clean and configure - no make")
  USAGE+=$(printf "[?%-12s #%s]" "skip" "no build")
  USAGE+="}"
  USAGE+="[X:optcos?Build option for Cosmo.]:[optcos:='${def_options["cos"]}']"
  USAGE+="[Y:optclm?Build option for CLM.]:[optclm:='${def_options["clm"]}']"
  USAGE+="[Z:optpfl?Build option for Parflow.]:[optpfl:='${def_options["pfl"]}']"
  USAGE+="[H:hyprepath?Include Path for Hypre.]:[hyprepath:='$hyprePath']"
  USAGE+="[S:silopath?Include Path for Silo.]:[silopath:='$siloPath']"
  USAGE+="[T:tclpath?Include Path for TCL.]:[tclpath:='$tclPath']"
  USAGE+="[G:grib1path?Include Path for Grib1.]:[grib1path:='$grib1Path']"
  USAGE+="[M:mpipath?Include Path for MPI.]:[mpipath:='$mpiPath']"
  USAGE+="[N:ncdfpath?Include Path for NetCDF.]:[ncdfpath:='$ncdfPath']"
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
    o)  optComp="${OPTARG}" ; args=1 ;; 
    v)  if [[ "${!versions[${OPTARG}]}" != ""  ]] then
            version=$OPTARG
        else
            print "The selected version '${OPTARG}' is not available. run '.$call --man' for help"
            terminate
        fi   
        args=1 ;;
    a)  listAvailabilities ;;
    t)  listTutorial ;;
    r)  rootdir="$OPTARG" ; args=1 ;;
    B)  bindir="$OPTARG" ; args=1 ;;
    C)  case "${combinations_limitation[@]}" in *$OPTARG*)  combination="$OPTARG" ; args=1 ; valid="true" ;; esac 
        if [[ $valid != "true" ]] then; wmessage="This combination is not supported in this version" ; warning ; combination="$OPTARG" ; args=1  ;fi
        ;;
    W)  options+=(["oas"]="$OPTARG") ; args=1 ;;
    X)  options+=(["cos"]="$OPTARG") ; args=1 ;;
    Y)  options+=(["clm"]="$OPTARG") ; args=1 ;;
    Z)  options+=(["pfl"]="$OPTARG") ; args=1 ;;

    w)  oasdir="$OPTARG" ; args=1 ;;
    x)  cosdir="$OPTARG"; args=1 ;;
    y)  clmdir="$OPTARG"; args=1 ;;
    z)  pfldir="$OPTARG"; args=1 ;;

    M)  mpiPath="$OPTARG" ; args=1 ;;
    N)  ncdfPath="$OPTARG" ; args=1 ;;
    G)  grib1Path="$OPTARG" ; args=1 ;;
    T)  tclPath="$OPTARG" ; args=1 ;;
    H)  hyprePath="$OPTARG" ; args=1 ;;
    S)  siloPath="$OPTARG" ; args=1 ;;    
    esac
  done


  date=`date +%d%m%Y`


  if [[ $bindir == "" ]] then
     bindir="$rootdir/bin/${platform}_${date}" 
  else
     bindir="$bindir/${platform}_${date}" 
  fi


print -n "  create bindir: $bindir"
  mkdir -p $bindir >> $log_file 2>> $err_file
check
  rm -f $bindir/* >> $log_file 2>> $err_file
check

print -n "  source common interface"
  . ${rootdir}/bldsva/intf_oas3/common_build_interface.ksh >> $log_file 2>> $err_file
check
print -n "  source machine build interface for $platform"
  . ${rootdir}/bldsva/machines/${platform}/build_interface_${platform}.ksh >> $log_file 2>> $err_file
check


  getMachineDefaults 
  setMachine

  withOAS="false"
  withCOS="false"
  withPFL="false"
  withCLM="false"
  withOASMCT="false"

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



 runCompilation 

  mv -f $err_file $bindir
  mv -f $log_file $bindir

 # debug
