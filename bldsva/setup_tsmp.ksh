#! /bin/ksh

######################################################
##################### Defaults #######################
######################################################

getDefaults(){
  def_platform=""               
  def_version=""                 
  def_rootdir="$estdir"    #This should be correct - change with caution
  def_combination=""           
  def_refSetup=""
  def_bindir=""				#Will be set to $rootdir/bin/$platform_$version_$combination if empty
  def_rundir=""  			#Will be set to $rootdir/run/$platform_${version}_$combination_$refSetup_$exp_id if empty
  def_exp_id="__DATE__"				
  def_numInst=""
  def_startInst=""
  # following parameters  will be set to tested platform defaults if empty
  def_nppn=""
  def_ngpn=""
  def_wtime=""
  ## following parameters  will be set to tested setup defaults if empty
  def_namelist_clm=""
  def_namelist_cos=""
  def_namelist_icon=""
  def_namelist_oas=""
  def_namelist_pfl=""
  def_restart_script=""
#DA
  def_namelist_da=""

  def_forcingdir_clm=""
  def_forcingdir_cos=""
  def_forcingdir_icon=""
  def_forcingdir_oas=""
  def_forcingdir_pfl=""
  def_restfile_pfl=""
  def_restfile_clm=""
  def_restfile_cos=""
  def_restfile_icon=""
  def_dump_clm=""
  def_dump_cos=""
  def_dump_icon=""
  def_processor=""
  def_dump_pfl=""	  
  def_queue=""
  def_px_clm=""
  def_py_clm=""
  def_px_cos=""
  def_py_cos=""
  def_p_icon=""
  def_px_pfl=""
  def_py_pfl=""
  def_startDate=""
  def_initDate=""
  def_runhours=""
  #profiling ("yes" , "no") - will use machine standard
  def_profiling="no"
  def_cplscheme=""
  def_mode=""
  def_compiler="Intel" # set default Intel
}

#####################################################
# USERS SHOULD NOT EDIT BELOW THIS LINE
#####################################################

setDefaults(){
  platform=$def_platform
  compiler=$def_compiler
  if [[ $platform == "" ]] then ; platform="JUWELS" ; fi #We need a hard default here
  version=$def_combination
  if [[ $version == "" ]] then ; version="" ; fi #We need a hard default here
  bindir=$def_bindir
  rundir=$def_rundir
  exp_id=$def_exp_id
  profiling=$def_profiling
  rootdir=$def_rootdir

  cplscheme=$def_cplscheme
  if [[ $cplscheme == "" ]] then ; cplscheme="true" ; fi #We need a hard default here
  mode=$def_mode
  if [[ $mode == "" ]] then ; mode="0" ; fi #We need a hard default here

  log_file=$cpwd/log_all_${date}.txt
  err_file=$cpwd/err_all_${date}.txt
  stdout_file=$cpwd/stdout_all_${date}.txt
  rm -f $log_file $err_file $stdout_file
 
  refSetup=$def_refSetup
  combination=$def_combination
  nppn=$def_nppn
  ngpn=$def_ngpn
  numInst=$def_numInst
  if [[ $numInst == "" ]] then ; numInst="1" ; fi #We need a hard default here
  startInst=$def_startInst
  if [[ $startInst == "" ]] then ; startInst="0" ; fi #We need a hard default here
  queue=$def_queue
  wtime=$def_wtime
  px_clm=$def_px_clm
  py_clm=$def_py_clm
  px_cos=$def_px_cos
  py_cos=$def_py_cos
  p_icon=$def_p_icon
  px_pfl=$def_px_pfl
  py_pfl=$def_py_pfl
  forcingdir_clm=$def_forcingdir_clm
  forcingdir_cos=$def_forcingdir_cos
  forcingdir_icon=$def_forcingdir_icon
  forcingdir_oas=$def_forcingdir_oas
  forcingdir_pfl=$def_forcingdir_pfl
  restfile_pfl=$def_restfile_pfl
  restfile_clm=$def_restfile_clm
  restfile_cos=$def_restfile_cos
  restfile_icon=$def_restfile_icon
  dump_clm=$def_dump_clm
  dump_cos=$def_dump_cos
  dump_icon=$def_dump_icon
  processor=$def_processor
  dump_pfl=$def_dump_pfl
  namelist_clm=$def_namelist_clm
  namelist_cos=$def_namelist_cos
  namelist_icon=$def_namelist_icon
  namelist_oas=$def_namelist_oas
  namelist_pfl=$def_namelist_pfl
  restart_script=$def_restart_script
#DA
  namelist_da=$def_namelist_da

  startDate=$def_startDate
  runhours=$def_runhours
  initDate=$def_initDate

  #profiling
  profComp=""
  profRun=""
  profVar=""

  withCESM="false"
}

clearMachineSelection(){
  wtime=""
  queue=""
  clearSetupSelection
}

clearSetupSelection(){
  nppn=""
  px_clm=""
  py_clm=""
  px_cos=""
  py_cos=""
  p_icon=""
  px_pfl=""
  py_pfl=""
  startDate=""
  initDate=""
  runhours=""

  forcingdir_clm=""
  forcingdir_cos=""
  forcingdir_icon=""
  forcingdir_pfl=""
  forcingdir_oas=""
  clearPathSelection
}

clearPathSelection(){
  rundir=""
  bindir=""
  namelist_clm=""
  namelist_cos=""
  namelist_icon=""
  namelist_pfl=""
  namelist_oas=""
  restart_script=""
#DA
  namelist_da=""
}

setSelection(){

  if [[ $nppn == "" ]] then ; nppn=$defaultNppn ; fi
  if [[ $ngpn == "" ]] then ; ngpn=$defaultNgpn ; fi
  if [[ $wtime == "" ]] then ; wtime=$defaultwtime  ; fi
  if [[ $queue == "" ]] then ; queue=$defaultQ ; fi
  if [[ $px_clm == "" ]] then ; px_clm=$defaultCLMProcX ; fi
  if [[ $py_clm == "" ]] then ; py_clm=$defaultCLMProcY ; fi
  if [[ $px_cos == "" ]] then ; px_cos=$defaultCOSProcX ; fi
  if [[ $py_cos == "" ]] then ; py_cos=$defaultCOSProcY ; fi
  if [[ $p_icon == "" ]] then ; p_icon=$defaultICONProc ; fi
  if [[ $px_pfl == "" ]] then ; px_pfl=$defaultPFLProcX ; fi
  if [[ $py_pfl == "" ]] then ; py_pfl=$defaultPFLProcY ; fi
  
  if [[ $forcingdir_clm == "" ]] then ; forcingdir_clm=$defaultFDCLM ; fi
  if [[ $forcingdir_cos == "" ]] then ; forcingdir_cos=$defaultFDCOS ; fi
  if [[ $forcingdir_icon == "" ]] then ; forcingdir_icon=$defaultFDICON ; fi
  if [[ $forcingdir_oas == "" ]] then ; forcingdir_oas=$defaultFDOAS ; fi
  if [[ $forcingdir_pfl == "" ]] then ; forcingdir_pfl=$defaultFDPFL ; fi

  if [[ $namelist_clm == "" ]] then ; namelist_clm=$defaultNLCLM ; fi
  if [[ $namelist_cos == "" ]] then ; namelist_cos=$defaultNLCOS ; fi
  if [[ $namelist_icon == "" ]] then ; namelist_icon=$defaultNLICON ; fi
  if [[ $namelist_oas == "" ]] then ; namelist_oas=$defaultNLOAS ; fi
  if [[ $namelist_pfl == "" ]] then ; namelist_pfl=$defaultNLPFL ; fi
  if [[ $restart_script == "" ]] then ; restart_script=$defaultRST ; fi
#DA
  if [[ $namelist_da == "" ]] then ; namelist_da=$defaultNLDA ; fi

  if [[ $startDate == "" ]] then ; startDate=$defaultStartDate ; fi
  if [[ $initDate == "" ]] then ; initDate=$defaultInitDate ; fi
  if [[ $runhours == "" ]] then ; runhours=$defaultRunhours ; fi

  if [[ $dump_clm == "" ]] then ; dump_clm=$defaultDumpCLM ; fi
  if [[ $dump_cos == "" ]] then ; dump_cos=$defaultDumpCOS ; fi
  if [[ $dump_icon == "" ]] then ; dump_icon=$defaultDumpICON ; fi
  if [[ $processor == "" ]] then ; processor=$defaultprocessor ; fi
  if [[ $dump_pfl == "" ]] then ; dump_pfl=$defaultDumpPFL ; fi

  if [[ $exp_id == "__DATE__" ]] then
     exp_id="_${date}"
  fi

  if [[ $rundir == "" ]] then
     rundir="$rootdir/run/${platform}_${combination}_${refSetup}"
  fi
  rundir="${rundir}${exp_id}"


  if [[ $bindir == "" ]] then
     bindir="$rootdir/bin/${platform}_${combination}"
  fi
  
   if echo "$combination" | grep -q 'pdaf'; then

     if echo "$combination" | grep -q 'clm5'; then

       mListgen="clm5-cos5-pfl-pdaf"

     else

       if echo "$combination" | grep -q 'cos4'; then
	 mListgen="clm3-cos4-pfl-pdaf"
       else
         mListgen="clm3-cos5-pfl-pdaf"
       fi

     fi

   elif echo "$combination" | grep -q 'clm4' && echo "$combination" | grep -q 'cos4'; then
	mListgen="clm4-cos4-pfl"
   elif echo "$combination" | grep -q 'clm4'; then
	mListgen="clm4-cos5-pfl"
   
   elif echo "$combination" | grep -q 'icon21'; then
	mListgen="clm3-icon21-pfl"
   elif echo "$combination" | grep -q 'icon26'; then
	mListgen="clm3-icon26-pfl"
	   
   elif echo "$combination" | grep -q 'eclm'; then
	mListgen="eclm"
   elif echo "$combination" | grep -q 'eclm-mct'; then
	mListgen="eclm-mct"

   else 
	if echo "$combination" | grep -q 'cos4'; then
		mListgen="clm3-cos4-pfl"
	else
		mListgen="clm3-cos5-pfl"
	fi
  fi
  version=$mListgen 
  set -A mList ${modelVersion[$mListgen]}

  #Fix for cos /clm namelist because they differ in newer version.
  if [[ ${mList[1]} == clm4_0  ]] ; then ; namelist_clm+=4_0 ; fi
  if [[ ${mList[1]} == clm5_0  ]] ; then ; namelist_clm+=5_0 ; fi
  if [[ ${mList[2]} == cosmo5_1  ]] ; then ; namelist_cos+=5_1 ; fi


}


setCombination(){
  withOAS="false"
  withCOS="false"
  withICON="false"
  withPFL="false"
  withCLM="false"
  withOASMCT="false"
  withPCLM="false"
#DA
  withDA="false"
  withPDAF="false"

  case "$combination" in *icon*) withICON="true" ;; esac
  case "$combination" in *clm*) withCLM="true" ;; esac
  case "$combination" in *cos*) withCOS="true" ;; esac
  case "$combination" in *pfl*) withPFL="true" ;; esac
  if [[ $withCLM == "true" && ( $withCOS == "true" || $withICON == "true" || $withPFL == "true" )  ]]; then
    withOAS="true"
    withOASMCT="true"
  fi  
#DA
  case "$combination" in *pdaf*) withDA="true" ; withPDAF="true" ;; esac
}

finalizeSelection(){

  mkdir -p $rundir
  if [[ -n $(ls -A $rundir) ]] ; then
     mv ${rundir} ${rundir}_backup_${date}
     mkdir $rundir
  fi

  nproc_oas=0
  nproc_pfl=0
  nproc_clm=0
  nproc_cos=0
  nproc_icon=0
  if [[ $withPFL == "true" ]] ; then ; nproc_pfl=$((${px_pfl}*${py_pfl})) ;fi 
  if [[ $withCOS == "true" ]] ; then ; nproc_cos=$((${px_cos}*${py_cos})) ;fi 
  if [[ $withICON == "true" ]] ; then ; nproc_icon=${p_icon} ;fi 
  if [[ $withCLM == "true" ]] ; then ; nproc_clm=$((${px_clm}*${py_clm})) ;fi 
  if [[ $withOAS == "true" ]] ; then ; nproc_oas=1 ;fi 
  if [[ $withOASMCT == "true" ]] ; then ; nproc_oas=0 ;fi 


  yyyy=$(date '+%Y' -d "$startDate")
  mm=$(date '+%m' -d "$startDate")
  dd=$(date '+%d' -d "$startDate")
  hh=$(date '+%H' -d "$startDate")

}

terminate(){
 print ""
 print "Terminating $call. No changes were made...${cnormal}"
 rm -f $err_file
 rm -f $log_file
 rm -f $stdout_file
 exit 0
}


check(){
 if [[ $? == 0  ]] then
    print "    ... ${cgreen}OK!${cnormal}"  | tee -a $stdout_file
 else
    print "    ... ${cred}error!!! - aborting...${cnormal}" | tee -a $stdout_file
    print "See $log_file and $err_file" | tee -a $stdout_file
    exit 1
  fi
}
check_pfl(){
 if [[ $? == 0  ]] then
    print "    ... ${cgreen}OK!${cnormal}"  | tee -a $stdout_file
 else
    print " No file for this setup case" | tee -a $stdout_file
  fi
}

comment(){
  print -n "$1" | tee -a $stdout_file
}

route(){
  print "$1" | tee -a $stdout_file
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


hardSanityCheck(){
  if [[ "${platforms[${platform}]}" == ""  ]] then
      print "The selected platform '${platform}' is not available. run '$call --man' for help"
      terminate
  fi

  valid="false"
  if [[ ${refSetup} == "" ]] ; then ; valid="true" ; fi
  case "${setupsAvail[${platform}]}" in *" ${refSetup} "*) valid="true" ;; esac
  if [[ $valid != "true" ]] then; print "This setup is not supported on this machine" ; terminate  ;fi

}

deprecatedVersion(){
  if [[ "${version}" != ""  ]] then
      print "The use of the internal version with -v is deprecated. Please provide your desired combination with -c including version numbers (clm3-cos5-pfl). "
      terminate
  fi
}

softSanityCheck(){

  #check multi instance functionality (only working with Oasis3-MCT)
  if [[ $numInst > 1 && $withOAS == "true" && $withOASMCT == "false" ]]; then ; wmessage="The -N option is only supported with Oasis3-MCT. it will be ignored if you continue." ; warning  ;fi

  valid="false"
  case "${availability[${platform}]}" in *" ${version} "*) valid="true" ;; esac
  if [[ $valid != "true" ]] then; wmessage="This version is not supported on this machine" ; warning  ;fi 

  valid="false"
  cstr="invalid"
  if [[ $withCLM == "true" &&  $withICON == "false" &&  $withCOS == "true" && $withPFL == "true"  ]]; then ;cstr=" clm-cos-pfl " ; fi
  if [[ $withCLM == "true" &&  $withICON == "true" &&  $withCOS == "false" && $withPFL == "true"  ]]; then ;cstr=" clm-icon-pfl " ; fi
  if [[ $withCLM == "true" &&  $withICON == "false" &&  $withCOS == "true" && $withPFL == "false"  ]]; then ;cstr=" clm-cos " ; fi
  if [[ $withCLM == "true" &&  $withICON == "true" &&  $withCOS == "false" && $withPFL == "false"  ]]; then ;cstr=" clm-icon " ; fi
  if [[ $withCLM == "true" &&  $withICON == "false" &&  $withCOS == "false" && $withPFL == "true"  ]]; then ;cstr=" clm-pfl " ; fi
  if [[ $withCLM == "true" &&  $withICON == "false" &&  $withCOS == "false" && $withPFL == "false"  ]]; then ;cstr=" clm " ; fi
  if [[ $withCLM == "false" &&  $withICON == "false" &&  $withCOS == "true" && $withPFL == "false"  ]]; then ;cstr=" cos " ; fi
  if [[ $withCLM == "false" &&  $withICON == "true" &&  $withCOS == "false" && $withPFL == "false"  ]]; then ;cstr=" icon " ; fi
  if [[ $withCLM == "false" &&  $withICON == "false" &&  $withCOS == "false" && $withPFL == "true"  ]]; then ;cstr=" pfl " ; fi
  case "${combinations[${version}]}" in *${cstr}*) valid="true" ;; esac
  if [[ $valid != "true" ]] then; wmessage="This combination is not supported in this version" ; warning  ;fi

  valid="false"
  case "${profilingImpl}" in *" ${profiling} "*) valid="true" ;; esac
  if [[ $valid != "true" ]] then; wmessage="This profiling tool is not supported on this machine" ; warning  ;fi

}

interactive(){
  clear
  print "${cyellow}##############################################${cnormal}"
  print "${cyellow}         Interactive installation...          ${cnormal}"
  print "${cyellow}##############################################${cnormal}"
  print "The following variables are needed:"
  printState
  print "${cyellow}##############################################${cnormal}"
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
                        comment "  source machine build interface for $platform"
  			  . ${rootdir}/bldsva/machines/config_${platform}.ksh >> $log_file 2>> $err_file
      			check
			clearMachineSelection
      			getMachineDefaults
			#reset features if not supported by new machine selection
    			case "${availability[$platform]}" in 
				*" $version "*) ;;
		     		*)
	  			set -A array ${availability[$platform]}
          			version=${array[0]} ;; 
    			esac
    			case "${combinations[$version]}" in 
        			*" $combination "*);;
		         	*)      
          			set -A array ${combinations[$version]}
          			combination=${array[0]} ;;     
    			esac

    			case "${setupsAvail[$platform]}" in 
        			*" $refSetup "*) ;;
		      		*)     
          			set -A array ${setupsAvail[$platform]}
          			refSetup=${array[0]} ;;     
    			esac
      			setCombination
      			comment "  source setup for $refSetup on $platform"
        		  . ${rootdir}/bldsva/setups/$refSetup/${refSetup}.ksh >> $log_file 2>> $err_file
        		  . ${rootdir}/bldsva/setups/common_setup.ksh >> $log_file 2>> $err_file
      			check
      			initSetup
      			setSelection

                  fi
                  if [[ $numb == 3 ]] ; then
                        print "The following combinations are available for $version:"
                        for a in ${combinations[$version]} ; do
                                printf "%-20s\n" "$a"
                        done
                        print "Please type in your desired value..."
                        read combination
                        clearSetupSelection
                        setCombination
                        comment "  source setup for $refSetup on $platform"
                          . ${rootdir}/bldsva/setups/$refSetup/${refSetup}.ksh >> $log_file 2>> $err_file
                          . ${rootdir}/bldsva/setups/common_setup.ksh >> $log_file 2>> $err_file
                        check
                        initSetup
                        setSelection

                  fi
                  if [[ $numb == 4 ]] ; then 
			print "The following setups are available for $platform:"
                        for a in ${setupsAvail[$platform]} ; do
                                printf "%-20s\n" "$a"
                        done
                        print "Please type in your desired value..."
		        read refSetup
                        clearSetupSelection
                        setCombination
                        comment "  source setup for $refSetup on $platform"
                          . ${rootdir}/bldsva/setups/$refSetup/${refSetup}.ksh >> $log_file 2>> $err_file
                          . ${rootdir}/bldsva/setups/common_setup.ksh >> $log_file 2>> $err_file
                        check
                        initSetup
                        setSelection

		  fi


                  if [[ $numb == 5 ]] ; then ; read queue ; fi
                  if [[ $numb == 6 ]] ; then ; read wtime ; fi
                  if [[ $numb == 7 ]] ; then ; read startDate  ; fi
                  if [[ $numb == 8 ]] ; then ; read runhours ; fi
                  if [[ $numb == 9 ]] ; then ; read nppn ; fi
                  if [[ $numb == 10 ]] ; then ; read px_clm ; fi
                  if [[ $numb == 11 ]] ; then ; read py_clm ; fi
                  if [[ $numb == 12 ]] ; then ; read px_cos ; fi
                  if [[ $numb == 13 ]] ; then ; read py_cos ; fi
                  if [[ $numb == 40 ]] ; then ; read p_icon ; fi
                  if [[ $numb == 14 ]] ; then ; read px_pfl ; fi
                  if [[ $numb == 15 ]] ; then ; read py_pfl ; fi
                  if [[ $numb == 16 ]] ; then ; read rootdir ;clearPathSelection; setSelection ; fi
                  if [[ $numb == 17 ]] ; then ; read bindir ; fi
                  if [[ $numb == 18 ]] ; then ; read rundir ; fi
                  if [[ $numb == 19 ]] ; then ; read namelist_clm ; fi
                  if [[ $numb == 20 ]] ; then ; read forcingdir_clm ; fi
                  if [[ $numb == 21 ]] ; then ; read namelist_cos ; fi
                  if [[ $numb == 41 ]] ; then ; read namelist_icon ; fi
                  if [[ $numb == 22 ]] ; then ; read forcingdir_cos ; fi
                  if [[ $numb == 42 ]] ; then ; read forcingdir_icon ; fi
                  if [[ $numb == 23 ]] ; then ; read namelist_pfl ; fi
                  if [[ $numb == 24 ]] ; then ; read forcingdir_pfl ; fi
                  if [[ $numb == 25 ]] ; then ; read namelist_oas ; fi
                  if [[ $numb == 26 ]] ; then ; read forcingdir_oas ; fi
		  if [[ $numb == 27 ]] ; then ; read namelist_da ; fi
                  if [[ $numb == 28 ]] ; then 
                         print "The following profiling tools are available for $platform:"
                         for a in ${profilingImpl} ; do
                                printf "%-20s\n" "$a"
                         done
                         print "Please type in your desired value..."
                         read profiling
                  fi
                  if [[ $numb == 29 ]] ; then ; read numInst ; fi
  		  if [[ $numb == 30 ]] ; then ; read startInst ; fi
                  if [[ $numb == 31 ]] ; then ; read initDate ; fi
		  if [[ $numb == 32 ]] ; then ; read cplscheme ; fi
   		  if [[ $numb == 33 ]] ; then ; read exp_id ; fi
		  if [[ $numb == 34 ]] ; then ; read restfile_pfl ; fi
		  if [[ $numb == 35 ]] ; then ; read restfile_clm ; fi
		  if [[ $numb == 36 ]] ; then ; read restfile_cos ; fi
		  if [[ $numb == 43 ]] ; then ; read restfile_icon ; fi
                  if [[ $numb == 37 ]] ; then ; read dump_pfl ; fi
                  if [[ $numb == 38 ]] ; then ; read dump_clm ; fi
                  if [[ $numb == 39 ]] ; then ; read dump_cos ; fi
                  if [[ $numb == 44 ]] ; then ; read dump_icon ; fi
                  if [[ $numb == 45 ]] ; then ; read compiler ; fi
                  if [[ $numb == 46 ]] ; then ; read processor ; fi
                  if [[ $numb == 47 ]] ; then ; read ngpn ; fi
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
  print "${cred}(45)${cnormal} compiler (default=$def_compiler): ${cgreen}$compiler ${cnormal}"
  print "${cred}(2)${cnormal} Deprecated. Please specify your desired combination with the -c option."
  print "${cred}(3)${cnormal} combination (default=$def_combination): ${cgreen}$combination${cnormal}"
  print "${cred}(4)${cnormal} refSetup (default=$def_refSetup): ${cgreen}$refSetup${cnormal}"
  print ""
  print "${cred}(5)${cnormal} scheduler queue (default=$def_queue): ${cgreen}$queue${cnormal}"
  print "${cred}(6)${cnormal} wallclock time (default=$def_wtime): ${cgreen}$wtime${cnormal}"
  print "${cred}(7)${cnormal} startDate (default=$def_startDate): ${cgreen}$startDate${cnormal}"
  print "${cred}(8)${cnormal} run hours (default=$def_runhours): ${cgreen}$runhours${cnormal}"
  print "${cred}(9)${cnormal} number of processors per node (default=$def_nppn): ${cgreen}$nppn${cnormal}"
  print "${cred}(47)${cnormal} number of GPUs per node (default=$def_ngpn): ${cgreen}$ngpn${cnormal}"
  print ""
  print "${cred}(10)${cnormal} processors in X for clm (default=$def_px_clm): ${cgreen}$px_clm${cnormal}"
  print "${cred}(11)${cnormal} processors in Y for clm (default=$def_py_clm): ${cgreen}$py_clm${cnormal}"
  print "${cred}(12)${cnormal} processors in X for cos (default=$def_px_cos): ${cgreen}$px_cos${cnormal}"
  print "${cred}(13)${cnormal} processors in Y for cos (default=$def_py_cos): ${cgreen}$py_cos${cnormal}"
  print "${cred}(40)${cnormal} processors for icon     (default=$def_p_icon): ${cgreen}$p_icon${cnormal}"
  print "${cred}(14)${cnormal} processors in X for pfl (default=$def_px_pfl): ${cgreen}$px_pfl${cnormal}"
  print "${cred}(15)${cnormal} processors in Y for pfl (default=$def_py_pfl): ${cgreen}$py_pfl${cnormal}"
  print ""
  print "${cred}(16)${cnormal} root dir (default=$def_rootdir): ${cgreen}$rootdir${cnormal}"
  print "${cred}(17)${cnormal} bin dir (default=$def_rootdir/bin/${def_platform}_${version}_${combination}): ${cgreen}$bindir${cnormal}"
  print "${cred}(18)${cnormal} run dir (default=$def_rootdir/run/${def_platform}_${version}_${combination}_${refSetup}_${exp_id}): ${cgreen}$rundir${cnormal}"
  print ""
  print "${cred}(19)${cnormal} namelist dir for clm (default='$def_namelist_clm'): ${cgreen}$namelist_clm${cnormal}"
  print "${cred}(20)${cnormal} forcing dir for clm (default='$def_forcingdir_clm'): ${cgreen}$forcingdir_clm${cnormal}"
  print "${cred}(21)${cnormal} namelist dir for cos (default='$def_namelist_cos'): ${cgreen}$namelist_cos${cnormal}"
  print "${cred}(41)${cnormal} namelist dir for icon (default='$def_namelist_icon'): ${cgreen}$namelist_icon${cnormal}"
  print "${cred}(22)${cnormal} forcing dir for cos (default='$def_forcingdir_cos'): ${cgreen}$forcingdir_cos${cnormal}"
  print "${cred}(42)${cnormal} forcing dir for icon (default='$def_forcingdir_icon'): ${cgreen}$forcingdir_icon${cnormal}"
  print "${cred}(23)${cnormal} namelist dir for pfl (default='$def_namelist_pfl'): ${cgreen}$namelist_pfl${cnormal}"
  print "${cred}(24)${cnormal} forcing dir for pfl (default='$def_forcingdir_pfl'): ${cgreen}$forcingdir_pfl${cnormal}"
  print "${cred}(25)${cnormal} namelist dir for oas (default='$def_namelist_oas'): ${cgreen}$namelist_oas${cnormal}"
  print "${cred}(26)${cnormal} forcing dir for oas (default='$def_forcingdir_oas'): ${cgreen}$forcingdir_oas${cnormal}"
#DA
  print "${cred}(27)${cnormal} namelist dir for data assimilation (default='$def_namelist_da'): ${cgreen}$namelist_da${cnormal}"
  print ""
  print "${cred}(28)${cnormal} profiling (default=$def_profiling): ${cgreen}$profiling${cnormal}"
  print "${cred}(29)${cnormal} number of TerrSysMP instances (default=$def_numInst): ${cgreen}$numInst${cnormal}"
  print "${cred}(30)${cnormal} counter to start TerrSysMP instances with (default=$def_startInst): ${cgreen}$startInst${cnormal}"
  print "${cred}(31)${cnormal} init Date (default=$def_initDate): ${cgreen}$initDate${cnormal}"
  print "${cred}(32)${cnormal} Couple-Scheme (default=$def_cplscheme): ${cgreen}$cplscheme ${cnormal}"
  print "${cred}(33)${cnormal} Experiment ID. DATE will be taken if ''. (default=$def_exp_id): ${cgreen}$exp_id ${cnormal}"
  print ""
  print "${cred}(34)${cnormal} Path to restart file for pfl. No restart if ''. (default=$def_restfile_pfl): ${cgreen}$restfile_pfl ${cnormal}"
  print "${cred}(35)${cnormal} Path to restart file for clm. No restart if ''. (default=$def_restfile_clm): ${cgreen}$restfile_clm ${cnormal}"
  print "${cred}(36)${cnormal} Path to restart file for cos. No restart if ''. (default=$def_restfile_cos): ${cgreen}$restfile_cos ${cnormal}"  
  print "${cred}(43)${cnormal} Path to restart file for icon. No restart if ''. (default=$def_restfile_icon): ${cgreen}$restfile_icon ${cnormal}"  
  print "${cred}(37)${cnormal} Dump interval for pfl.  (default=$def_dump_pfl): ${cgreen}$dump_pfl ${cnormal}"
  print "${cred}(38)${cnormal} Dump interval for clm.  (default=$def_dump_clm): ${cgreen}$dump_clm ${cnormal}"
  print "${cred}(39)${cnormal} Dump interval for cos.  (default=$def_dump_cos): ${cgreen}$dump_cos ${cnormal}"
  print "${cred}(43)${cnormal} Dump interval for icon.  (default=$def_dump_icon): ${cgreen}$dump_icon ${cnormal}"
  print "${cred}(44)${cnormal} Architecture type: CPU, GPU, MSA.  (default=$def_processor): ${cgreen}$processor ${cnormal}"
}


listAvailabilities(){

  print ${cyellow}"A list of details for each available platform."${cnormal}
  print ""
  for p in "${!platforms[@]}" ; do
    printf "%-20s #%s\n" "$p" "${platforms[$p]}"
    print ${cgreen}$'\t supported versions:'${cnormal}
    for a in ${availability[$p]} ; do
        print $'\t '"$a"
    done
    print ${cgreen}$'\t supported setups:'${cnormal}
    for a in ${setupsAvail[$p]} ; do
        print $'\t '"$a"
    done
  done


  print ""
  print ${cyellow}"A list of details for each version."${cnormal} 
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
  print ""
  print ${cyellow}"A list of details for each setup."${cnormal}
  print ""
  for v in "${!setups[@]}" ; do
    printf "%-20s #%s\n" "$v" "${setups[$v]}"
  done



  exit 0
}

listTutorial(){
  print "1) TODO"
 
  exit 0
}


getRoot(){
  #automatically determine root dir
  cpwd=`pwd`
  if [[ "$0" == '/'*  ]] ; then
    #absolut path
    estdir=`echo "$0" | sed 's@/bldsva/setup_tsmp.ksh@@'` #remove bldsva/configure machine to get rootpath
    call=$0
  else
    #relative path
    call=`echo $0 | sed 's@^\.@@'`                    #clean call from leading dot
    call=`echo $call | sed 's@^/@@'`                  #clean call from leading /
    call=`echo "/$call" | sed 's@^/\./@/\.\./@'`      #if script is called without leading ./ replace /./ by /../
    curr=`echo $cpwd | sed 's@^/@@'`                   #current directory without leading /   
    call=`echo $call | sed "s@$curr@@"`               #remove current directory from call if absolute path was called
    estdir=`echo "/$curr$call" | sed 's@/bldsva/setup_tsmp.ksh@@'` #remove bldsva/configure machine to get rootpath
    call="$estdir/bldsva$call"
  fi 
  
}


############################ 
# ICON interface methods
############################

c_setup_icon(){
route "${cyellow}>>> c_setup_icon${cnormal}"

comment "  cp add_run_routines to rundir"
  cp $rootdir/bldsva/setups/idiurnal-cycle/common/add_run_routines $rundir >> $log_file 2>> $err_file
check

comment "  cp namelist to rundir"
  cp ${namelist_icon} $rundir >> $log_file 2>> $err_file
check

comment "  sed dt to namelist"
  sed "s,__dt_icon_bldsva__,$dt_icon," -i $rundir/NAMELIST_icon >> $log_file 2>> $err_file
check

comment "  sed start time to namelist"
  dSD=($defaultStartDate)
  sed "s,__starttime_icon_bldsva__,${dSD[0]}T${dSD[1]}:00:00Z," -i $rundir/icon_master.namelist >> $log_file 2>> $err_file
  sed "s,__starttime_icon_bldsva__,${dSD[0]}T${dSD[1]}:00:00Z," -i $rundir/NAMELIST_icon >> $log_file 2>> $err_file
check

comment "  sed end time to namelist"
  dED=($defaultEndDate)
  sed "s,__endtime_icon_bldsva__,${dED[0]}T${dED[1]}:00:00Z," -i $rundir/icon_master.namelist >> $log_file 2>> $err_file
  sed "s,__endtime_icon_bldsva__,${dED[0]}T${dED[1]}:00:00Z," -i $rundir/NAMELIST_icon >> $log_file 2>> $err_file
check

route "${cyellow}<<< c_setup_icon${cnormal}"
}

setup_icon(){
route "${cyellow}>> setupIcon${cnormal}"

  c_setup_icon

route "${cyellow}<< setupIcon${cnormal}" 
}

############################ 
#Cosmo interface methods
############################

c_setup_cos(){
route "${cyellow}>>> c_setup_cos${cnormal}"

comment "  cp namelist to rundir"
  cp ${namelist_cos} $rundir/lmrun_uc >> $log_file 2>> $err_file
check

nstop_cos=$((  ($runhours*3600 + ($(date -u '+%s' -d "${startDate}") - $(date -u '+%s' -d "${initDate}")) )  /$dt_cos  ))
#if [[ $withCESM == "false" ]] ; then ; nstop_cos=$(($nstop_cos-($cplfreq1/$dt_cos))) ; fi

comment "  sed dt to namelist"
  sed "s,dt_cos_bldsva,$dt_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check

comment "  sed number of procs to namelist"
  sed "s,nprocx_cos_bldsva,$px_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s,nprocy_cos_bldsva,$py_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check

comment "  sed gridpoints to namelist"
  sed "s,ie_tot_bldsva,$gx_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s,je_tot_bldsva,$gy_cos," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check

comment "  sed gridpoints to namelist"
  sed "s,nbdl_cos_bldsva,$nbndlines," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check


comment "  create input dir for cosmo"
  mkdir -p $rundir/cosmo_in >> $log_file 2>> $err_file
check
comment "  fill cosmo input dir with softlinks from cosmo forcing dir"
  ln -sf $forcingdir_cos/* $rundir/cosmo_in >> $log_file 2>> $err_file
check

comment "  sed forcingdir to namelist"
  sed "s,__forcingdir__,$rundir/cosmo_in," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed rundir to namelist"
  sed "s,__rundir__,$rundir," -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed stop time to namelist"
  sed "s/nstop_cos_bldsva/$nstop_cos/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed date to namelist"
  sed "s/init_y_bldsva/$(date '+%Y' -d "$initDate")/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed init_m_bldsva to namelist"  
  sed "s/init_m_bldsva/$(date '+%m' -d "$initDate")/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_d_bldsva/$(date '+%d' -d "$initDate")/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
  sed "s/init_h_bldsva/$(date '+%H' -d "$initDate")/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check

cnt=$(( ($(date -u '+%s' -d "${startDate}") - $(date -u '+%s' -d "${initDate}"))/3600))
comment "  sed start hour to namelist"
sed "s/__hstart__/$cnt/" -i $rundir/lmrun_uc >> $log_file 2>> $err_file
check
comment "  sed restart interval to namelist"
#sed "s/__nhour_restart_start__/$(($cnt+$runhours))/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
sed "s/__nhour_restart_start__/$cnt/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check
sed "s/__nhour_restart_stop__/$(($cnt+$runhours))/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check
sed "s/__nhour_restart_incr__/1/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check

cnts=$(( ( $(date -u '+%s' -d "${startDate}") - $(date -u '+%s' -d "${initDate}")) / ${dt_cos} ))
comment "  sed output interval to namelist"
sed "s/__ncomb_start__/$cnts/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check
sed "s/__dump_cos_interval__/ $(python -c "print ($dump_cos*(3600/$dt_cos))")/" -i $rundir/lmrun_uc  >> $log_file 2>> $err_file
check

if [[ $restfile_cos != "" ]] then
comment "  softlink restart file to input dir"
ln -s $restfile_cos $rundir/cosmo_in  >> $log_file 2>> $err_file
check
fi

comment "  cd to rundir"
  cd $rundir >> $log_file 2>> $err_file
check
comment "  run lmrun_uc clean"
  $rundir/lmrun_uc cleancluma >> $log_file 2>> $err_file
check
comment "  run lmrun_uc exe"
  $rundir/lmrun_uc execluma >> $log_file 2>> $err_file
check

if [[ $withPDAF == "true" ]] ; then
  cp $rundir/INPUT_IO $rundir/INPUT_IO_$(printf "%05d" $(($instance-$startInst)))     
fi

route "${cyellow}<<< c_setup_cos${cnormal}"
}

setup_cos(){
route "${cyellow}>> setupCos${cnormal}"

  c_setup_cos

route "${cyellow}<< setupCos${cnormal}" 
}

############################ 
#OASIS interface methods
############################

c_setup_oas(){
route "${cyellow}>>> c_setup_oas${cnormal}"

  comment "   copy cf_name_table to rundir"
    cp $rootdir/bldsva/data_oas3/cf_name_table.txt $rundir >> $log_file 2>> $err_file
  check
  comment "   copy oas namelist to rundir"
    cp $namelist_oas $rundir/namcouple >> $log_file 2>> $err_file
  check
  comment "   sed procs, gridsize & coupling freq into namcouple"

  ncpl_exe1=$nproc_cos
  ncpl_exe2=$nproc_pfl
  ncpl_exe3=1
  if [[ $withCESM == "true" || $withOASMCT == "true" ]] ; then ; ncpl_exe3=$nproc_clm ; fi


  if [[ $withICON == "true" ]]; then
    sed "s/ngiconx/$gx_icon/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq1/$cplfreq1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

    sed "s/ngclmx/$(($gx_clm*$gy_clm))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  fi


  if [[ $withPFL == "true" && $withCOS == "true" ]] then

    sed "s/nproc_exe1/$nproc_cos/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe3/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe3/$ncpl_exe3/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq1/$cplfreq1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq2/$cplfreq2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    
    sed "s/ngpflx/$gx_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngpfly/$gy_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmx/$(($gx_clm*$gy_clm))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngcosx/$(($gx_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngcosy/$(($gy_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

  fi  
if [[ $withPFL == "true" && $withCOS == "false" ]] then

    sed "s/nproc_exe1/$nproc_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe3/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq2/$cplfreq2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

    sed "s/ngpflx/$gx_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngpfly/$gy_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmx/$gx_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/$gy_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  fi
  if [[ $withPFL == "false" && $withCOS == "true" ]] then
    sed "s/nproc_exe1/$nproc_cos/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe3/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq1/$cplfreq1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

    sed "s/ngcosx/$(($gx_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngcosy/$(($gy_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmx/$gx_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/$gy_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

  fi
  rtime=$((($runhours*3600 + $cplfreq1)*$rtimeFactor))
  if [[ $withCESM == "true" ]] ; then ; rtime=$(($rtime+$cplfreq1)) ; fi
  comment "   sed sim time into namcouple"
    sed "s/totalruntime/$rtime/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  comment "   sed startdate into namcouple"
    sed "s/yyyymmdd/${yyyy}${mm}${dd}/" -i $rundir/namcouple  >> $log_file 2>> $err_file
  check


route "${cyellow}<<< c_setup_oas${cnormal}"
}

setup_oas(){
route "${cyellow}>> setupOas${cnormal}"

  c_setup_oas

route "${cyellow}<< setupOas${cnormal}"
}

############################ 
#CLM interface methods
############################

c_setup_clm(){
route "${cyellow}>>> c_setup_clm${cnormal}"
if [[ ${mList[1]} == "clm5_0" ]]; then
 comment "  CLM5.0 setup"

 comment "  sed rundir to namelist"
   sed "s,__rundir__,$rundir," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
 check

 comment "  sed dt to namelist"
   sed "s,__dt_clm_bldsva__,$dt_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
 check

 comment "  sed forcingdir to namelist"
   sed "s,__forcingdir__,$forcingdir_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
 check
 comment "  sed dump interval namelist"
   sed "s,__dump_clm_interval__,$dump_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
 check
 comment "  sed runtime to namelist"
   if [[ $withICON == "true" ]]; then
     runstep_clm=$((($runhours*3600 + $cplfreq1/10)/$dt_clm))
   else
     runstep_clm=$((($runhours*3600 + $cplfreq1)/$dt_clm))
   fi
   sed "s,__runstep_clm_bldsva__,$runstep_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
 check
 comment "  add axe rights to clm namelist"
     chmod 755 $rundir/lnd.stdin >> $log_file 2>> $err_file
 check
 comment "  run clm namelist"
     $rundir/lnd.stdin >> $log_file 2>> $err_file
 check
# Below default clm 3.5 setup
else

comment "  sed rundir to namelist"
  sed "s,__rundir__,$rundir," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check

comment "  sed starttime to namelist"
  sed "s,__seconds_clm_bldsva__,$seconds_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed dt to namelist"
  sed "s,__dt_clm_bldsva__,$dt_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed forcingdir to namelist"
  sed "s,__forcingdir__,$forcingdir_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed gridsize to namelist"
  sed "s,__gridsize__,$res," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  create rpointer dummy file"
  touch $rpointer >> $log_file 2>> $err_file
check
comment "  sed rpointer path to namelist"
  sed "s,__rundir_rpointerdir__,$rpointer," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed date to namelist"
  sed "s,__yyyymmdd_bldsva__,${yyyy}${mm}${dd}," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed dump interval namelist"
  sed "s,__dump_clm_interval__,$dump_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed runtime to namelist"

runstep_clm=$((($runhours*3600 + $cplfreq1)/$dt_clm))

  sed "s,__runstep_clm_bldsva__,$runstep_clm," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  sed restart file path to namelist"
    sed "s,__finidat__,${restfile_clm}," -i $rundir/lnd.stdin >> $log_file 2>> $err_file
check


comment "  add axe rights to clm namelist"
    chmod 755 $rundir/lnd.stdin >> $log_file 2>> $err_file
check
comment "  run clm namelist"
    $rundir/lnd.stdin >> $log_file 2>> $err_file
check

if [[ $withPDAF == "true" ]] ; then
  cp $rundir/lnd.stdin $rundir/lnd.stdin_$(printf "%05d" $(($instance-$startInst)))     
fi

fi
route "${cyellow}<<< c_setup_clm${cnormal}"
}

setup_clm(){
route "${cyellow}>> setupClm${cnormal}"
  seconds_clm=$(($hh*3600))
#  runstep_clm=$(($runhours*3600/$dt_clm))
  rpointer=$rundir/lnd.clmoas.rpointer

comment "  cp namelist to rundir"
  cp $namelist_clm $rundir/lnd.stdin >> $log_file 2>> $err_file
check


  c_setup_clm
}  

############################ 
# eCLM interface methods
############################

c_setup_eclm(){
  route "${cyellow}>>> c_setup_eclm${cnormal}"
  route "${cyellow}<<< c_setup_eclm${cnormal}"
}

############################ 
#Parflow interface methods
############################

c_setup_pfl(){
route "${cyellow}>>> c_setup_pfl${cnormal}"

  if [ $numInst > 1 ] || [! -f "$rundir/coup_oas.tcl" ]; then
    comment "  $rundir/coup_oas.tcl does not exist (or for ensemble runs: numInst > 1), is copied, see c_setup_pfl()"
    cp $namelist_pfl $rundir/coup_oas.tcl >> $log_file 2>> $err_file
    check
  fi
 
  if [[ $withICON == "true" ]]; then
   comment "  sed ccs_ic_press file to pfl namelist."
      sed "s,__ccs_ic_press__,$defaultFDPFL/ccs_ic_press.pfb," -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
   check
  fi
 
  comment "   sed nproc x to pfl namelist."
    sed "s/__nprocx_pfl_bldsva__/$px_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed nproc y to pfl namelist."
    sed "s/__nprocy_pfl_bldsva__/$py_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed gridpoints x to pfl namelist."
    sed "s/__ngpflx_bldsva__/$gx_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed gridpoints y to pfl namelist."
    sed "s/__ngpfly_bldsva__/$gy_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed forcingdir to pfl namelist."
    sed "s,__forcingdir__,$rundir," -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed dt to pfl namelist."
    sed "s/__dt_pfl_bldsva__/$dt_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed end time to pfl namelist."
#    sed "s/__stop_pfl_bldsva__/$runstep_clm/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
    sed "s/__stop_pfl_bldsva__/$(python -c "print (${runhours} + ${base_pfl})")/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed dump interval to pfl namelist."
    sed "s/__dump_pfl_interval__/$dump_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check

  comment "   sed timing base to pfl namelist."
    sed "s/__base_pfl__/$base_pfl/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check

  comment "   sed start counter to pfl namelist."
      cnt=$(( ($(date -u '+%s' -d "${startDate}") - $(date -u '+%s' -d "${initDate}"))))
      cnt=$(python -c "print ($cnt/($dump_pfl*3600.))")
      sed "s/__start_cnt_pfl__/$cnt/" -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check


    if [[ $restfile_pfl == "" ]] then
  comment "   sed initial condition to pfl namelist."
      sed "s/__pfl_ICPpressureType__/HydroStaticPatch/" -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file      # HydrostaticPatch > PFBFile
  check
      sed "s/__pfl_ICPpressureValue__/-5.0/" -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file      # delete this during restart run
  check
  comment "   sed delete restart file name from pfl namelist."
      sed '/__pfl_ICPpressureFileName__/d' -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file
  check
    else
  comment "   sed delete IC_Pressure value  from pfl namelist." 
      sed '/__pfl_ICPpressureValue__/d' -i $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check
  comment "   sed initial condition to pfl namelist."
      sed "s/__pfl_ICPpressureType__/PFBFile/" -i $rundir/coup_oas.tcl   >> $log_file 2>> $err_file      # HydrostaticPatch > PFBFile
  check
      sed "s,__pfl_ICPpressureFileName__,$restfile_pfl," -i $rundir/coup_oas.tcl  >> $log_file 2>> $err_file
  check
    fi
  
  export PARFLOW_DIR=$bindir
  comment "   cd to rundir."
    cd $rundir >> $log_file 2>> $err_file
  check

  comment "   create parflow db with tclsh from namelist."
  check
  tclsh $rundir/coup_oas.tcl >> $log_file 2>> $err_file
  check


route "${cyellow}<<< c_setup_pfl${cnormal}"
}

setup_pfl(){
route "${cyellow}>> setup_pfl${cnormal}"
  c_setup_pfl

route "${cyellow}<< setup_pfl${cnormal}"
}

############################ 
#PDAF interface methods
############################


c_setup_pdaf(){
route "${cyellow}>>> c_setup_pdaf${cnormal}"
  comment "   copy pdaf namelist to rundir."
    cp $namelist_da $rundir/enkfpf.par >> $log_file 2>> $err_file
  check 
  comment "   sed num instances into pdaf namelist."
    sed "s/__ninst__/$(($numInst-$startInst))/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed pflname into pdaf namelist."
    sed "s/__pflname__/$pflrunname/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed pflproc into pdaf namelist."
    sed "s/__pflproc__/$nproc_pfl/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed dt into pdaf namelist."
    sed "s/__dt__/$dt_pfl/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed endtime into pdaf namelist."
    sed "s/__endtime__/$(python -c "print (${runhours} + ${base_pfl})")/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed clmproc into pdaf namelist."
    sed "s/__clmproc__/$nproc_clm/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check
  comment "   sed cosproc into pdaf namelist."
    sed "s/__cosproc__/$nproc_cos/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check 
  comment "   sed dtmult into pdaf namelist."
    sed "s/__dtmult__/$(python -c "print (${dt_pfl} * 3600 / ${dt_cos})")/" -i $rundir/enkfpf.par >> $log_file 2>> $err_file
  check 

route "${cyellow}<<< c_setup_pdaf${cnormal}"
}

c_setup_rst(){

 comment " copy $restart_script to $rundir"
   cp $restart_script $rundir >> $log_file 2>> $err_file
 check

 comment "   sed startDate into restart template."
    sed 's/__startDate_bldsva__/"'"$startDate"'"/' -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
  check

 comment "   sed initDate into restart template."
    sed 's/__initDate_bldsva__/"'"$initDate"'"/' -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
  check

 comment "   sed dt_clm into restart template."
    sed "s/__dt_clm_bldsva__/$dt_clm/" -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
  check

 comment "   sed dt_cosmo into restart template."
    sed "s/__dt_cos_bldsva__/$dt_cos/" -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
  check

 comment "   sed PARFLOW_DIR into restart template $bindir."
#    sed "/__PARFLOW_DIR__/ \$bindir" -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
    sed -i "s|__PARFLOW_DIR__|$bindir|" $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
#    sed "s/__PARFLOW_DIR__/$bindir/" -i $rundir/tsmp_restart.sh >> $log_file 2>> $err_file
  check

}

setup_da(){
route "${cyellow}>> setup_da${cnormal}"
  c_setup_pdaf
route "${cyellow}<< setup_da${cnormal}"
}

#######################################
#               Main
#######################################

  cyellow=$(tput setaf 3)
  cnormal=$(tput sgr0)
  cred=$(tput setaf 1)
  cgreen=$(tput setaf 2)

  date=`date +%d%m%y-%H%M%S`
  getRoot 
  getDefaults
  setDefaults

  #GetOpts definition
  USAGE=$'[-?\n@(#)$Id: TerrSysMP setup script 1.0 - '
  USAGE+=$' date: 10.10.2015 $\n]'
  USAGE+="[-author?Fabian Gasper]"
  USAGE+="[+NAME?TerrSysMP setup script]"
  USAGE+="[+DESCRIPTION?sets up TSMP run by handling namelists and copying necessary files into a run directory]"
  USAGE+="[b:bash?Bash mode - set command line arguments will overwrite default values (no interactive mode) (This is the default with arguments).]"
  USAGE+="[i:interactive?Interactive mode - command line arguments and defaults will be overwritten during the interactive session (This is the default without arguments).]"
  USAGE+="[a:avail?Prints a listing of every machine with available versions. The script will exit afterwards.]"
  USAGE+="[t:tutorial?Prints a tutorial/description on how to add new versions and platforms to this script. The script will exit afterwards.]"
  USAGE+="[R:rootdir?Absolute path to TerrSysMP root directory.]:[path:='$def_rootdir']"
  USAGE+="[B:bindir?Absolute path to bin directory for the builded executables.]:[path:='$def_bindir']"
  USAGE+="[r:rundir?Absolute path to run directory for TSMP.]:[path:='$def_rundir']"
  USAGE+="[I:expid?Experiment ID suffix for rundir. DATE will be taken if empty.]:[expid:='$def_exp_id']"
  USAGE+="[v:version?Tagged TerrSysMP version. Note that not every version might be implemented on every machine. Run option -a, --avail to get a listing.]:[version:='$version']"
  USAGE+="[V:refsetup?Reference setup. This is a setup that is supported and tested on a certain machine. No further inputs are needed for this setup. If you leave it '' the machine default will be taken.]:[refsetup:='$def_refSetup']"
  USAGE+="[m:machine?Target Platform. Run option -a, --avail to get a listing.]:[machine:='$platform']"
  USAGE+="[p:profiling?Makes necessary changes to compile with a profiling tool if available.]:[profiling:='$def_profiling']"
  USAGE+="[c:combination? Combination of component models.]:[combination:='$def_combination']"
  USAGE+="[C:cplscheme? Couple-Scheme for CLM/COS coupling.]:[cplscheme:='$def_cplscheme']"
  USAGE+="[O:compiler? Compiler used.]:[compiler:='$def_compiler']"
  USAGE+="[P:nppn? Number of processors per node. If you leave it '' the machine default will be taken.]:[nppn:='$def_nppn']"
  USAGE+="[H:ngpn? Number of GPUs per node. If you leave it '' the machine default will be taken.]:[ngpn:='$def_ngpn']"
  USAGE+="[N:numinst? Number of instances of TerrSysMP. Currently only works with Oasis3-MCT - ignored otherwise.]:[numinst:='$def_numInst']"
  USAGE+="[n:startinst? Instance counter to start with. Currently only works with Oasis3-MCT - ignored otherwise.]:[startinst:='$def_startInst']"
  USAGE+="[q:queue? Scheduler Queue name. If you leave it '' the machine default will be taken.]:[queue:='$def_queue']"
  USAGE+="[Q:wtime? Wallclocktime for your run. If you leave it '' the machine default will be taken.]:[wtime:='$def_wtime']"
  USAGE+="[s:startdate? (Restart-) Date for your model run. Must be given in the format: 'YYYY-MM-DD_HH']:[startdate:='$def_startDate']"
  USAGE+="[S:initdate? Initial date for your model (restart) runs. Must be given in the format: 'YYYY-MM-DD_HH']:[initdate:='$def_initDate']" 
  USAGE+="[T:runhours? Number of simulated hours.]:[runhours:='$def_runhours']" 
  USAGE+="[w:pxclm? Number of tasks for clm in X direction. If you leave it '' the machine default will be taken.]:[pxclm:='$def_px_clm']"
  USAGE+="[W:pyclm? Number of tasks for clm in Y direction. If you leave it '' the machine default will be taken.]:[pyclm:='$def_py_clm']"
  USAGE+="[x:pxcos? Number of tasks for cosmo in X direction. If you leave it '' the machine default will be taken.]:[pxcos:='$def_px_cos']"
  USAGE+="[X:pycos? Number of tasks for cosmo in Y direction. If you leave it '' the machine default will be taken.]:[pycos:='$def_py_cos']"
  USAGE+="[z:picon? Number of tasks for icon. If you leave it '' the machine default will be taken.]:[picon:='$def_p_icon']"
  USAGE+="[y:pxpfl? Number of tasks for ParFlow in X direction. If you leave it '' the machine default will be taken.]:[pxpfl:='$def_px_pfl']"
  USAGE+="[Y:pypfl? Number of tasks for ParFlow in Y direction. If you leave it '' the machine default will be taken.]:[pypfl:='$def_py_pfl']"
  USAGE+="[d:namoas? Namelist for Oasis. This script will always try to substitute the placeholders by the reference setup values. Make sure your namelist and placeholders are compatible with the reference setup. If you don't wont the substitution remove placeholders from your namelist. This flag will replace the default namelist from the reference setup ]:[namoas:='']"
  USAGE+="[D:forcediroas? Forcing directory for oasis. This will replace the default forcing dir from the reference setup.]:[forcediroas:='']"
  USAGE+="[e:namclm? Namelist for clm. This script will always try to substitute the placeholders by the reference setup values. Make sure your namelist and placeholders are compatible with the reference setup. If you don't wont the substitution remove placeholders from your namelist. This flag will replace the default namelist from the reference setup ]:[namclm:='']"
  USAGE+="[E:forcedirclm? Forcing directory for clm. This will replace the default forcing dir from the reference setup.]:[forcedirclm:='']"
  USAGE+="[f:namcos? Namelist for Cosmo. This script will always try to substitute the placeholders by the reference setup values. Make sure your namelist and placeholders are compatible with the reference setup. If you don't wont the substitution remove placeholders from your namelist. This flag will replace the default namelist from the reference setup ]:[namcos:='']"
  USAGE+="[Z:namicon? Namelist for icon. This script will always try to substitute the placeholders by the reference setup values. Make sure your namelist and placeholders are compatible with the reference setup. If you don't wont the substitution remove placeholders from your namelist. This flag will replace the default namelist from the reference setup ]:[namicon:='']"
  USAGE+="[U:forcediricon? Forcing directory for ICON. This will replace the default forcing dir from the reference setup.]:[forcediricon:='']"
  USAGE+="[F:forcedircosmo? Forcing directory for Cosmo. This will replace the default forcing dir from the reference setup.]:[forcedircosmo:='']"
  USAGE+="[g:nampfl? Namelist for ParFlow. This script will always try to substitute the placeholders by the reference setup values. Make sure your namelist and placeholders are compatible with the reference setup. If you don't wont the substitution remove placeholders from your namelist. This flag will replace the default namelist from the reference setup ]:[nampfl:='']"
#DA
  USAGE+="[h:namda? Namelist for data assimilation. This script will always try to substitute the placeholders by the reference setup values. Make sure your namelist and placeholders are compatible with the reference setup. If you don't wont the substitution remove placeholders from your namelist. This flag will replace the default namelist from the reference setup ]:[namda:='']"

  USAGE+="[G:forcedirpfl? Forcing directory for ParFlow. This will replace the default forcing dir from the reference setup.]:[forcedirpfl:='']"
  USAGE+="[j:restfileclm? Restart file for CLM. No restart for CLM if ''.]:[restfileclm:='$def_restfile_clm']"
  USAGE+="[k:restfilecos? Restart file for Cosmo. No restart for Cosmo if ''.]:[restfilecos:='$def_restfile_cos']"
  USAGE+="[e:restfileicon? Restart file for icon. No restart for icon if ''.]:[restfileicon:='$def_restfile_icon']"
  USAGE+="[l:restfilepfl? Restart file for ParFlow. No restart for ParFlow if ''.]:[restfilepfl:='$def_restfile_pfl']"
  USAGE+="[J:dumpclm? Dump interval for CLM (in hours).]:[dumpclm:='$def_dump_clm']"
  USAGE+="[K:dumpcos? Dump interval for Cosmo (in hours).]:[dumpcos:='$def_dump_cos']"
  USAGE+="[E:dumpicon? Dump interval for Cosmo (in hours).]:[dumpicon:='$def_dump_icon']"
  USAGE+="[A:processor? Processor architecture: CPU, GPU, MSA.]:[processor:='$def_processor']" 
  USAGE+="[L:dumppfl? Dump interval for ParFlow (in hours).]:[dumppfl:='$def_dump_pfl']"


  USAGE+=$'\n\n\n\n'


  args=0
  # parsing the command line arguments
  while getopts "$USAGE" optchar ; do
    case $optchar in
    i)  mode=2 ;;
    b)  mode=1 ;;
    m)  platform="$OPTARG" ; args=1 ;; 
    p)  profiling="${OPTARG}" ; args=1 ;;
    v)  version="$OPTARG"  ;  args=1 ;;
    V)  refSetup="$OPTARG"  ;  args=1 ;; 
    a)  listA="true" ;;
    t)  listTutorial ;;
    q)  queue=$OPTARG ; args=1 ;; 
    Q)  wtime=$OPTARG ; args=1 ;;   
    P)  nppn=$OPTARG ; args=1 ;;
    H)  ngpn=$OPTARG ; args=1 ;;
    n)  startInst=$OPTARG ; args=1 ;;
    N)  numInst=$OPTARG ; args=1 ;;
    s)  startDate=$(echo $OPTARG | sed -e "s/_/ /g"); args=1 ;;
    S)  initDate=$(echo $OPTARG | sed -e "s/_/ /g") ; args=1 ;;
    T)  runhours=$OPTARG ; args=1 ;;    
    
    d)  namelist_oas=$OPTARG ; args=1 ;;
    D)  forcingdir_oas=$OPTARG ; args=1 ;;
    e)  namelist_clm=$OPTARG ; args=1 ;;
    E)  forcingdir_clm=$OPTARG ; args=1 ;;
    f)  namelist_cos=$OPTARG ; args=1 ;;
    Z)  namelist_icon=$OPTARG ; args=1 ;;
    F)  forcingdir_cos=$OPTARG ; args=1 ;;
    U)  forcingdir_icon=$OPTARG ; args=1 ;;
    g)  namelist_pfl=$OPTARG ; args=1 ;;
    G)  forcingdir_pfl=$OPTARG ; args=1 ;;
#DA
    h)  namelist_da=$OPTARG ; args=1 ;;

    w)  px_clm=$OPTARG ; args=1 ;;
    W)  py_clm=$OPTARG ; args=1 ;;
    x)  px_cos=$OPTARG ; args=1 ;;
    X)  py_cos=$OPTARG ; args=1 ;;
    z)  p_icon=$OPTARG ; args=1 ;;
    y)  px_pfl=$OPTARG ; args=1 ;;
    Y)  py_pfl=$OPTARG ; args=1 ;;
 
    R)  rootdir="$OPTARG" ; args=1 ;;
    B)  bindir="$OPTARG" ; args=1 ;;
    r)  rundir="$OPTARG" ; args=1 ;;
    I)  exp_id="$OPTARG" ; args=1 ;;
    c)  combination="$OPTARG" ; args=1 ;;
    C)  cplscheme="$OPTARG" ; args=1 ;;
    O)  compiler="$OPTARG" ; args=1 ;;

    l)  restfile_pfl="$OPTARG"; args=1 ;;
    j)  restfile_clm="$OPTARG"; args=1 ;;
    k)  restfile_cos="$OPTARG"; args=1 ;;
    e)  restfile_icon="$OPTARG"; args=1 ;;

    L)  dump_pfl="$OPTARG"; args=1 ;;
    J)  dump_clm="$OPTARG"; args=1 ;;
    K)  dump_cos="$OPTARG"; args=1 ;;
    E)  dump_icon="$OPTARG"; args=1 ;;
    A)  processor="$OPTARG" ; args=1 ;;


    esac
  done


deprecatedVersion


comment "  source list with supported machines and configurations"
  . $rootdir/bldsva/supported_versions.ksh
check

  if [[ $listA == "true" ]] ; then ; listAvailabilities ; fi

  hardSanityCheck

  #if no combination is set, load first as default
  if [[ $combination == "" ]] ; then
    set -A array ${combinations[$version]}
    combination=${array[0]}
  fi
  if [[ $refSetup == "" ]] ; then
    set -A array ${setupsAvail[$platform]}
    refSetup=${array[0]}
  fi
  
  comment "  source machine build interface for $platform"
    . ${rootdir}/bldsva/machines/config_${platform}.ksh >> $log_file 2>> $err_file
  check

  comment "  source setup for $refSetup on $platform"
	. ${rootdir}/bldsva/setups/$refSetup/${refSetup}.ksh >> $log_file 2>> $err_file
    . ${rootdir}/bldsva/setups/common_setup.ksh >> $log_file 2>> $err_file
    
  check

  setCombination
  initSetup

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


  softSanityCheck
  finalizeMachine
  finalizeSelection 

#  start setup
  origrundir=$rundir
  orignamelist_cos=$namelist_cos
  orignamelist_icon=$namelist_icon
  orignamelist_clm=$namelist_clm
  orignamelist_pfl=$namelist_pfl
  if [[ $withCLM == "true" ]] ; then 
    if [[ ${mList[1]} == "eclm" ]] ; then
        comment "   Creating symlink to eclm.exe"
        ln -s $bindir/bin/eclm.exe $rundir/eclm.exe >> $log_file 2>> $err_file
        check
    elif [[ $withPDAF == "false" ]] ; then
      comment "  cp clm exe to $rundir" 
        cp $bindir/clm $rundir >> $log_file 2>> $err_file
      check
    fi
  fi
  if [[ $withCOS == "true" ]] ; then
    if [[ $withPDAF == "false" ]] ; then
      comment "  cp cos exe to $rundir"
        cp $bindir/lmparbin_pur $rundir  >> $log_file 2>> $err_file
      check
    fi
  fi
  if [[ $withICON == "true" ]] ; then
    if [[ $withPDAF == "false" ]] ; then
      comment "  cp icon exe to $rundir"
        cp $bindir/icon $rundir  >> $log_file 2>> $err_file
      check
    fi
  fi
  if [[ $withPFL == "true" ]] ; then
    if [[ $withPDAF == "false" ]] ; then
      comment "  cp pfl exe to $rundir"
        cp $bindir/parflow $rundir  >> $log_file 2>> $err_file
      check
    fi
  fi
  if [[ $withOAS == "true" ]] ; then
    if [[ $withOASMCT == "false" ]] ; then
      comment "  cp oas exe to $rundir"
         cp $bindir/oasis3.MPI1.x $rundir  >> $log_file 2>> $err_file
      check
    fi
  fi

#DA
  if [[ $withPDAF == "true" ]] ; then
    comment "  cp da exe to $rundir"
      cp $bindir/tsmp-pdaf $rundir  >> $log_file 2>> $err_file
    check
  fi



  for instance in {$startInst..$(($startInst+$numInst-1))}
  do
  route ${cyellow}"> creating instance: $instance"${cnormal}
    # Ensemble only creation
    if [[ $numInst > 1 && ( $withOASMCT == "true" || $withOAS == "false"   ) && $withPDAF == "false" ]] ; then 
      rundir=$origrundir/tsmp_instance_$instance
      comment "  mkdir sub-directory for instance"
        mkdir -p $rundir >> $log_file 2>> $err_file
      check
      if [[ -e "${orignamelist_icon}_${instance}" ]] ; then 
	namelist_icon="${orignamelist_icon}_${instance}"
      else
        namelist_icon="$orignamelist_icon"
      fi
      if [[ -e "${orignamelist_cos}_${instance}" ]] ; then 
	namelist_cos="${orignamelist_cos}_${instance}"
      else
        namelist_cos="$orignamelist_cos"	
      fi
      if [[ -e "${orignamelist_clm}_${instance}" ]] ; then
        namelist_clm="${orignamelist_clm}_${instance}"
      else
        namelist_clm="$orignamelist_clm"
      fi
      if [[ -e "${orignamelist_pfl}_${instance}" ]] ; then
        namelist_pfl="${orignamelist_pfl}_${instance}"
      else
        namelist_pfl="$orignamelist_pfl"
      fi
    fi

 #DA
    #PDAF creation
    if [[ $withPDAF == "true" ]] ; then
      comment "  mkdir sub-directory for instance"
        mkdir -p $origrundir/tsmp_instance_$(printf "%05d" $instance) >> $log_file 2>> $err_file
      check
	
	namelist_icon="${orignamelist_icon}_${instance}"
	namelist_cos="${orignamelist_cos}_${instance}"
        namelist_clm="${orignamelist_clm}_${instance}"
        namelist_pfl="${orignamelist_pfl}_${instance}"
    fi

    if [[ $withCLM == "true" ]] ; then ; setup_clm ;  fi
    if [[ $withCOS == "true" ]] ; then ; setup_cos ;  fi
    if [[ $withICON == "true" ]] ; then ; setup_icon ;  fi
    if [[ $withPFL == "true" ]] ; then ; setup_pfl ;  fi
    if [[ $withOAS == "true" ]] ; then ; setup_oas ;  fi

    if [[ $withPDAF == "true" ]] ; then
      cp ${pflrunname}_$(printf "%05d" $instance).pfidb $origrundir/tsmp_instance_$(printf "%05d" $instance)   		
    else	
      finalizeSetup
    fi
  route ${cyellow}"< creating instance: $instance"${cnormal}
  done

#DA
  if [[ $withPDAF == "true" ]] ; then
    setup_da
    finalizeSetup
  fi

  rundir=$origrundir
  namelist_cos="$orignamelist_cos"
  namelist_clm="$orignamelist_clm"
  namelist_pfl="$orignamelist_pfl"

  createRunscript

  echo "Git:" >> $log_file
  cd $rootdir
  git rev-parse --abbrev-ref HEAD >> $log_file
  git rev-parse HEAD >> $log_file
  echo "Selection:" >> $log_file
  printState >> $log_file
  echo "Call:" >> $log_file
  print "$call $*">> $log_file
  #remove special charecters for coloring from logfiles
  sed -i "s,.\[32m,,g" $log_file
  sed -i "s,.\[31m,,g" $log_file
  sed -i "s,.\[34m,,g" $log_file
  sed -i "s,.[(]B.\[m,,g" $log_file

  sed -i "s,.\[32m,,g" $stdout_file
  sed -i "s,.\[31m,,g" $stdout_file
  sed -i "s,.\[34m,,g" $stdout_file
  sed -i "s,.[(]B.\[m,,g" $stdout_file

  mv -f $err_file $rundir
  mv -f $log_file $rundir
  mv -f $stdout_file $rundir

  print ${cgreen}"install script finished sucessfully"${cnormal}
  print "Rootdir: ${rootdir}"
  print "Bindir: ${bindir}"
  print "Rundir: ${rundir}" 

