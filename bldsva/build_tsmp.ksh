#! /bin/ksh

getDefaults(){
  def_platform="JURECA" 
  def_version="1.1.0MCT" 
  def_combination="clm-cos-pfl"
  def_rootdir="$estdir" #This should be correct - change with caution

  def_bindir=""				#Will be set to $rootdir/bin/$platform_$version_$combination if empty
  def_oasdir=""				#Will be set to $rootdir/XXX_$platform_$combination if empty
  def_cosdir=""
  def_clmdir=""
  def_pfldir=""

  # pathes will be set to tested platform defaults if empty
  def_mpiPath=""
  def_ncdfPath=""
  def_grib1Path=""
  def_tclPath=""
  def_hyprePath=""
  def_siloPath=""

  #compiler optimization
  def_optComp=""   # will be set to platform defaults if empty

  #profiling ("yes" , "no") - will use machine standard
  def_profiling="no"

  # fresh build clean in new folder
  # build=build clean
  # make=(resume) make - no configure
  # configure= only configure - no make
  # skip=do nothing
  def_options+=(["oas"]="fresh")
  def_options+=(["clm"]="fresh")
  def_options+=(["pfl"]="fresh")
  def_options+=(["cos"]="fresh")
}

setDefaults(){
  #load the default values
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
  mpiPath=$def_mpiPath
  ncdfPath=$def_ncdfPath
  grib1Path=$def_grib1Path
  tclPath=$def_tclPath
  hyprePath=$def_hyprePath
  siloPath=$def_siloPath
  combination=$def_combination

  log_file=$cpwd/log_all_${date}.txt
  err_file=$cpwd/err_all_${date}.txt
  rm -f $log_file $err_file

  options+=(["oas"]=${def_options["oas"]})
  options+=(["cos"]=${def_options["cos"]})
  options+=(["clm"]=${def_options["clm"]})
  options+=(["pfl"]=${def_options["pfl"]})

}

setSelection(){

  if [[ $mpiPath == "" ]] then ; mpiPath=$defaultMpiPath ; fi
  if [[ $ncdfPath == "" ]] then ; ncdfPath=$defaultNcdfPath  ; fi
  if [[ $grib1Path == "" ]] then ; grib1Path=$defaultGrib1Path ; fi
  if [[ $tclPath == "" ]] then ; tclPath=$defaultTclPath ; fi
  if [[ $hyprePath == "" ]] then ; hyprePath=$defaultHyprePath ; fi
  if [[ $siloPath == "" ]] then ; siloPath=$defaultSiloPath ; fi

  #compiler optimization
  if [[ $optComp == "" ]] then ; optComp=$defaultOptC ; fi


  set -A mList ${modelVersion[$version]}
  if [[ $oasdir == "" ]] then ;  oasdir=$rootdir/${mList[0]}_${platform}_${combination} ; fi
  if [[ $cosdir == "" ]] then ;  cosdir=$rootdir/${mList[2]}_${platform}_${combination} ; fi
  if [[ $clmdir == "" ]] then ;  clmdir=$rootdir/${mList[1]}_${platform}_${combination} ; fi
  if [[ $pfldir == "" ]] then ;  pfldir=$rootdir/${mList[3]}_${platform}_${combination} ; fi
  if [[ $bindir == "" ]] then
     bindir="$rootdir/bin/${platform}_${version}_${combination}"
  fi

print -n "  create bindir: $bindir"
  mkdir -p $bindir >> $log_file 2>> $err_file
check
#  rm -f $bindir/* >> $log_file 2>> $err_file
#check


  if [[ $combination == "" ]] ; then
    set -A array ${combinations[$version]}
    combination=${array[0]}
  fi
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
    . ${rootdir}/bldsva/intf_oas3/${mList[1]}/arch/${platform}/build_interface_${mList[1]}_${platform}.ksh >> $log_file 2>> $err_file
  check
    always_clm
    if [[ ${options["clm"]} == "skip" ]] ; then ; return  ;fi 
    if [[ ${options["clm"]} == "fresh" ]] ; then 
  print -n "  backup clm dir to: $clmdir"
      rm -rf $clmdir >> $log_file 2>> $err_file
  check
      cp -rf ${rootdir}/${mList[1]} $clmdir >> $log_file 2>> $err_file
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
    . ${rootdir}/bldsva/intf_oas3/${mList[2]}/arch/${platform}/build_interface_${mList[2]}_${platform}.ksh >> $log_file 2>> $err_file
  check
    always_cos
    if [[ ${options["cos"]} == "skip" ]] ; then ; return  ;fi 
    if [[ ${options["cos"]} == "fresh" ]] ; then 
  print -n "  backup cos dir to: $cosdir"
      rm -rf $cosdir >> $log_file 2>> $err_file
  check
      cp -rf ${rootdir}/${mList[2]} $cosdir >> $log_file 2>> $err_file
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
    . ${rootdir}/bldsva/intf_oas3/${mList[0]}/arch/${platform}/build_interface_${mList[0]}_${platform}.ksh >> $log_file 2>> $err_file
  check
    always_oas
    if [[ ${options["oas"]} == "skip" ]] ; then ; return  ;fi 
    if [[ ${options["oas"]} == "fresh" ]] ; then 
  print -n "  backup oas dir to: $oasdir"
      rm -rf $oasdir >> $log_file 2>> $err_file
  check
      cp -rf ${rootdir}/${mList[0]} $oasdir >> $log_file 2>> $err_file
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
    . ${rootdir}/bldsva/intf_oas3/${mList[3]}/arch/${platform}/build_interface_${mList[3]}_${platform}.ksh >> $log_file 2>> $err_file
  check
    always_pfl
    if [[ ${options["pfl"]} == "skip" ]] ; then ; return  ;fi 
    if [[ ${options["pfl"]} == "fresh" ]] ; then 
  print -n "  backup pfl dir to: $pfldir"
      rm -rf $pfldir >> $log_file 2>> $err_file
  check
      cp -rf ${rootdir}/${mList[3]} $pfldir >> $log_file 2>> $err_file
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

interactive(){
  clear
  print "${cblue}##############################################${cnormal}"
  print "${cblue}         Interactive installation...          ${cnormal}"
  print "${cblue}##############################################${cnormal}"
  print "The following variables are needed:"
  printState
  print "${cblue}##############################################${cnormal}"
  PS3="Your selection(1-3)?"
  select ret in "!!!start!!!" "edit" "exit"
  do  
    if [[ -n $ret ]]; then
       case $ret in
          "!!!start!!!") break ;;
          "edit") 
		while true
		do
		  print "Please select a ${cred}number${cnormal} you want to change!"
		  print "Select '0' to go back."
		  read numb
		  if [[ $numb == 0 ]] ; then ; break ; fi
		  if [[ $numb == 1 ]] ; then 
			print "The following options are available:"
			for a in "${!platforms[@]}" ; do
    				printf "%-20s #%s\n" "$a" "${platforms[$a]}"
			done
			print "Please type in your desired value..."
			read platform 
		  fi
		  if [[ $numb == 2 ]] ; then 
                        print "The following versions are available for $platform:"
                        for a in ${availability[$platform]} ; do
                                printf "%-20s #%s\n" "$a" "${versions[$a]} consisting of: "
				for b in ${modelVersion[$a]} ; do
				   printf "%-20s - %s\n" "" "$b"
				done
                        done		
			print "Please type in your desired value..."
			read version 
		  fi
		  if [[ $numb == 3 ]] ; then  
                        print "The following combinations are available for $version:"
                        for a in ${combinations[$version]} ; do
                                printf "%-20s\n" "$a"
                        done
                        print "Please type in your desired value..."
			read combination 
		  fi
		  if [[ $numb == 4 ]] ; then ; read val ;options+=(["oas"]="$val") ; fi
		  if [[ $numb == 5 ]] ; then ; read val ;options+=(["clm"]="$val") ; fi
		  if [[ $numb == 6 ]] ; then ; read val ;options+=(["cos"]="$val") ; fi
		  if [[ $numb == 7 ]] ; then ; read val ;options+=(["pfl"]="$val") ; fi
		  if [[ $numb == 8 ]] ; then ; read rootdir ; fi
		  if [[ $numb == 9 ]] ; then ; read bindir ; fi
		  if [[ $numb == 10 ]] ; then ; read oasdir ; fi
		  if [[ $numb == 11 ]] ; then ; read clmdir ; fi
		  if [[ $numb == 12 ]] ; then ; read cosdir ; fi
		  if [[ $numb == 13 ]] ; then ; read pfldir ; fi
		  if [[ $numb == 14 ]] ; then ; read mpiPath ; fi
		  if [[ $numb == 15 ]] ; then ; read siloPath ; fi
		  if [[ $numb == 16 ]] ; then ; read hyprePath ; fi
	 	  if [[ $numb == 17 ]] ; then ; read tclPath ; fi
		  if [[ $numb == 18 ]] ; then ; read grib1Path ; fi
		  if [[ $numb == 19 ]] ; then ; read ncdfPath ; fi
		  if [[ $numb == 20 ]] ; then ; read optComp ; fi
		  if [[ $numb == 21 ]] ; then ; read profiling ; fi
		done	
		interactive
	  ;;
          "exit") terminate ;;
       esac
       break
    fi  
  done

}
printState(){
  print ""
  print "${cred}(1)${cnormal} platform (default=$def_platform): ${cgreen}$platform${cnormal}"
  print "${cred}(2)${cnormal} version (default=$def_version): ${cgreen}$version${cnormal}"
  print "${cred}(3)${cnormal} combination (default=$def_combination): ${cgreen}$combination${cnormal}"
  print ""
  print "${cred}(4)${cnormal} oasis build option (default=${def_options["oas"]}): ${cgreen}${options["oas"]}${cnormal}" 
  print "${cred}(5)${cnormal} clm build option (default=${def_options["clm"]}): ${cgreen}${options["clm"]}${cnormal}"
  print "${cred}(6)${cnormal} cosmo build option (default=${def_options["cos"]}): ${cgreen}${options["cos"]}${cnormal}"
  print "${cred}(7)${cnormal} parflow build option (default=${def_options["pfl"]}): ${cgreen}${options["pfl"]}${cnormal}"
  print ""
  print "${cred}(8)${cnormal} root dir (default=$def_rootdir): ${cgreen}$rootdir${cnormal}"
  print "${cred}(9)${cnormal} bin dir (default=$def_rootdir/bin/${def_platform}_${version}_${combination}): ${cgreen}$bindir${cnormal}"
  print "${cred}(10)${cnormal} oasis dir (default=$def_rootdir/${mList[0]}_${def_platform}_$combination): ${cgreen}$oasdir${cnormal}"
  print "${cred}(11)${cnormal} clm dir (default=$def_rootdir/${mList[1]}_${def_platform}_$combination): ${cgreen}$clmdir${cnormal}"
  print "${cred}(12)${cnormal} cosmo dir (default=$def_rootdir/${mList[2]}_${def_platform}_$combination): ${cgreen}$cosdir${cnormal}"
  print "${cred}(13)${cnormal} parflow dir (default=$def_rootdir/${mList[3]}_${def_platform}_$combination): ${cgreen}$pfldir${cnormal}"
  print ""
  print "${cred}(14)${cnormal} mpi path (default=$defaultMpiPath): ${cgreen}$mpiPath${cnormal}"
  print "${cred}(15)${cnormal} silo path (default=$defaultSiloPath): ${cgreen}$siloPath${cnormal}"
  print "${cred}(16)${cnormal} hypre path (default=$defaultHyprePath): ${cgreen}$hyprePath${cnormal}"
  print "${cred}(17)${cnormal} tcl path (default=$defaultTclPath): ${cgreen}$tclPath${cnormal}"
  print "${cred}(18)${cnormal} grib1 path (default=$defaultGrib1Path): ${cgreen}$grib1Path${cnormal}"
  print "${cred}(19)${cnormal} ncdf path (default=$defaultNcdfPath): ${cgreen}$ncdfPath${cnormal}"
  print ""
  print "${cred}(20)${cnormal} optComp (default=$defaultOptComp): ${cgreen}$optComp${cnormal}"
  print "${cred}(21)${cnormal} profiling (default=$def_profiling): ${cgreen}$profiling${cnormal}"
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

sanityCheck(){

  if [[ "${versions[${version}]}" == ""  ]] then
      print "The selected version '${version}' is not available. run '.$call --man' for help"
      terminate
  fi

  if [[ "${platforms[${platform}]}" == ""  ]] then
      print "The selected platform '${platform}' is not available. run '.$call --man' for help"
      terminate
  fi


  valid="false"
  case "${availability[${platform}]}" in *" ${version} "*) valid="true" ;; esac
  if [[ $valid != "true" ]] then; wmessage="This version is not supported on this machine" ; warning  ;fi

  valid="false"
  cstr="invalid"
  if [[ $withCLM == "true" &&  $withCOS == "true" && $withPFL == "true"  ]]; then ;cstr=" clm-cos-pfl " ; fi
  if [[ $withCLM == "true" &&  $withCOS == "true" && $withPFL == "false"  ]]; then ;cstr=" clm-cos " ; fi
  if [[ $withCLM == "true" &&  $withCOS == "false" && $withPFL == "true"  ]]; then ;cstr=" clm-pfl " ; fi
  if [[ $withCLM == "true" &&  $withCOS == "false" && $withPFL == "false"  ]]; then ;cstr=" clm " ; fi
  if [[ $withCLM == "false" &&  $withCOS == "true" && $withPFL == "false"  ]]; then ;cstr=" cos " ; fi
  if [[ $withCLM == "false" &&  $withCOS == "false" && $withPFL == "true"  ]]; then ;cstr=" pfl " ; fi
  case "${combinations[${version}]}" in *${cstr}*) valid="true" ;; esac
  if [[ $valid != "true" ]] then; wmessage="This combination is not supported in this version" ; warning  ;fi
}

listAvailabilities(){
   
  print ${cblue}"A list of all supported versions for a special platform."${cnormal}
  print ""
  for p in "${!platforms[@]}" ; do
    printf "%-20s #%s\n" "$p" "${platforms[$p]}"
    for a in ${availability[$p]} ; do
        print $'\t'"$a"
    done
  done
  print ""
  print ${cblue}"A list of details for each version."${cnormal}	
  print ""
  for v in "${!versions[@]}" ; do
    printf "%-20s #%s\n" "$v" "${versions[$v]}"
    print ${cgreen}$'\t possible combinations:'${cnormal}
    print $'\t'${combinations[$v]}		
    print ${cgreen}$'\t componentmodel versions:'${cnormal}	
    for a in ${modelVersion[$v]} ; do
        print $'\t '"$a"
    done
  done

  exit 0
}

listTutorial(){
  print "1) The required component models need to be loaded in beforehand"
  print "2) The models must be named with their version name as specified in the supported_versions.ksh and in the intf-folder: for example cosmo4_21, parflow, clm3_5 or oasis3-mct in the root directory"
  print "2b) You can download oasis3-mct with: svn checkout http://oasis3mct.cerfacs.fr/svn/branches/OASIS3-MCT_2.0_branch/oasis3-mct"
  print "3) If not specified other, the component models will be copied to a working version with the name: MODEL_PLATFORM_COMBINATION"
  print "4) If a new version or platform is supported, edit the supported_versions.ksh and reflect all dependencies and constraints"
	
  exit 0
}

getRoot(){
  #automatically determine root dir
  call=`echo $0 | sed 's@^\.@@'`                    #clean call from leading dot
  cpwd=`pwd` 
  call=`echo $call | sed 's@^/@@'`                  #clean call from leading /
  call=`echo "/$call" | sed 's@^/\./@/\.\./@'`      #if script is called without leading ./ replace /./ by /../
  curr=`echo $cpwd | sed 's@^/@@'`                   #current directory without leading /    
  call=`echo $call | sed "s@$curr@@"`               #remove current directory from call if absolute path was called
  estdir=`echo "/$curr$call" | sed 's@/bldsva/build_tsmp.ksh@@'` #remove bldsva/configure machine to get rootpath
}

#######################################
#		Main
#######################################

  cblue=$(tput setaf 4)
  cnormal=$(tput setaf 9)
  cred=$(tput setaf 1)
  cgreen=$(tput setaf 2)
  typeset -A options
  typeset -A def_options

  date=`date +%d%m%y-%H%M%S`
  getRoot 
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
  USAGE+="[R:rootdir?Absolut path to TerrSysMP root directory.]:[path:='$def_rootdir']"
  USAGE+="[B:bindir?Absolut path to bin directory for the builded executables. bin/MACHINE_DATE will be taken if ''.]:[path:='$def_bindir']"
   
  USAGE+="[v:version?Tagged TerrSysMP version. Note that not every version might be implemented on every machine. Run option -a, --avail to get a listing.]:[version:='$version']"
  USAGE+="[m:machine?Target Platform. Run option -a, --avail to get a listing.]:[machine:='$def_platform']"

  USAGE+="[p:profiling?Makes necessary changes to compile with a profiling tool if available.]:[profiling:='$def_profiling']"
  USAGE+="[o:optimization?Compiler optimisation flags.]:[optimization:='$def_optComp']"
  USAGE+="[c:combination? Combination of component models.]:[combination:='$def_combination']"
  USAGE+="[W:optoas?Build option for Oasis.]:[optoas:='${def_options["oas"]}']{"
  USAGE+=$(printf "[?%-12s #%s]" "fresh" "build from scratch in a new folder")
  USAGE+=$(printf "[?%-12s #%s]" "build" "build clean")
  USAGE+=$(printf "[?%-12s #%s]" "make" "only (resume) make and make install - no make clean and configure")
  USAGE+=$(printf "[?%-12s #%s]" "configure" "only make clean and configure - no make")
  USAGE+=$(printf "[?%-12s #%s]" "skip" "no build")
  USAGE+="}"
  USAGE+="[X:optcos?Build option for Cosmo.]:[optcos:='${def_options["cos"]}']"
  USAGE+="[Y:optclm?Build option for CLM.]:[optclm:='${def_options["clm"]}']"
  USAGE+="[Z:optpfl?Build option for Parflow.]:[optpfl:='${def_options["pfl"]}']"
   
  USAGE+="[w:oasdir?Source directory for Oasis3. oasis3_MACHINE_DATE will be taken if ''.]:[oasdir:='${def_oasdir}']"
  USAGE+="[x:cosdir?Source directory for Cosmo. cosmo_MACHINE_DATE will be taken if ''.]:[cosdir:='${def_cosdir}']"
  USAGE+="[y:clmdir?Source directory for CLM. clm_MACHINE_DATE will be taken if ''.]:[clmdir:='${def_clmdir}']"
  USAGE+="[z:pfldir?Source directory for Parflow. parflow_MACHINE_DATE will be taken if ''.]:[pfldir:='${def_pfldir}']"

  USAGE+="[H:hyprepath?Include Path for Hypre. The machine default will be taken if ''.]:[hyprepath:='$hyprePath']"
  USAGE+="[S:silopath?Include Path for Silo. The machine default will be taken if ''.]:[silopath:='$siloPath']"
  USAGE+="[T:tclpath?Include Path for TCL. The machine default will be taken if ''.]:[tclpath:='$tclPath']"
  USAGE+="[G:grib1path?Include Path for Grib1. The machine default will be taken if ''.]:[grib1path:='$grib1Path']"
  USAGE+="[M:mpipath?Include Path for MPI. The machine default will be taken if ''.]:[mpipath:='$mpiPath']"
  USAGE+="[N:ncdfpath?Include Path for NetCDF. The machine default will be taken if ''.]:[ncdfpath:='$ncdfPath']"
  USAGE+=$'\n\n\n\n'



  mode=0
  args=0
  # parsing the command line arguments
  while getopts "$USAGE" optchar ; do
    case $optchar in
    i)  mode=2 ;;  
    b)  mode=1 ;;
    m)  platform="$OPTARG" ; args=1 ;;
    p)  profling="${OPTARG}" ; args=1 ;;
    o)  optComp="${OPTARG}" ; args=1 ;; 
    v)  version="$OPTARG"  ;  args=1 ;;
    a)  listA="true" ;;
    t)  listTutorial ;;
    R)  rootdir="$OPTARG" ; args=1 ;;
    B)  bindir="$OPTARG" ; args=1 ;;
    c)  combination="$OPTARG" ; args=1 ;;

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


  #choose combination
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


print -n "  source list with supported machines and configurations"
  . $rootdir/bldsva/supported_versions.ksh
check

  if [[ $listA == "true" ]] ; then ; listAvailabilities ; fi
  sanityCheck

print -n "  source common interface"
  . ${rootdir}/bldsva/intf_oas3/common_build_interface.ksh >> $log_file 2>> $err_file
check
print -n "  source machine build interface for $platform"
  . ${rootdir}/bldsva/machines/${platform}/build_interface_${platform}.ksh >> $log_file 2>> $err_file
check

  getMachineDefaults 
  setSelection

  # determine whether or not to run interactive session
  if [[ $mode == 0 ]] then
    if [[ $args == 0 ]] then
        mode=2
    else
        mode=1
    fi  
  fi
  if [[ $mode == 2 ]] then ; interactive ; fi


  runCompilation | tee test.txt 

  printState >> $log_file
  print "$call $*">> $log_file
  mv -f $err_file $bindir
  mv -f $log_file $bindir

  print ${cgreen}"build script finished sucessfully"${cnormal}
  print "Rootdir: ${rootdir}"
  print "Bindir: ${bindir}"


