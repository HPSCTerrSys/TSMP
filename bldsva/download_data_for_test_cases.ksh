#! /bin/ksh

# This is a script to download the input files needed to run the tsmp test cases.
# Version: V1.0
# Version: V1.1: adding -l 8 in wget to download the whole directory tree 
# Version: V1.2: Adding the script for downloading the input data for nrw test case, 04.02.2021
# Version: V1.3: Adding the script for downloading the input data and pre-processing-tools for IdealRTD test case, 15.06.2021
# Version: V1.4: Update cordex test case - change of path, 22.09.2022
# Version: V1.5: Adding the script for download the input data for idealscal test case, 26.09.2022
# Version: V1.6: Adding the script for download the input data for idiurnal-cycle test case, 03.04.2023
#--------------------------------------
# Initial owner: Abouzar Ghasemi, email: a.ghasemi@fz-juelich.de
# Current owner: Carl Hartick, email: c.hartick@fz-juelich.de
#-------------------------------------

if [ -z "$1" ]
  then
    echo "!! No arguments supplied !!"
    echo "Please nrw for nrw test case"
    echo "Please cordex for cordex test case"
    echo "Please idealrtd for idealised RTD test case"
    echo "Please idealscal for idealised scaling test case"
    echo "Please idiurnal-cycle for icon test case"
    exit 1
fi


testcase=$1
echo "input data downloads for $testcase"

if [ $testcase = cordex ]
then
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_eur11_eraint_eval_v2/input/
elif [ $testcase = nrw ]
then
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_nrw/input/
elif [ $testcase = idealrtd ]
then
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_idealrtd/input/
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_idealrtd/pre-processing-tools/
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_idealrtd/external/
elif [ $testcase = idealscal ]
then
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_idealscal/
elif [ $testcase = idiurnal-cycle ]
then
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_idiurnal-cycle/
elif [ $testcase = cordex0275 ]
then
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_eur0275/static/
fi

# Checking if the downloaded input files are corrupted.
if [ $testcase = cordex ]
then
     for i in $( ls -d tsmp_eur11_eraint_eval_v2/input/* )
     do
       cd $i
       md5sum -c checksums.md5
       cd ../../..
     done

elif [ $testcase = nrw ]
then
     for i in $( ls -d tsmp_nrw/input/* )
     do
       cd $i
       md5sum -c checksums.md5
       cd ../../..
     done
elif [ $testcase = idealrtd ]
then
     cd tsmp_idealrtd/pre-processing-tools
     md5sum -c checksums.md5
     cd ../..
elif [ $testcase = idealscal ]
then
     for i in $( find tsmp_idealscal/ -mindepth 2 -maxdepth 2 -type d )
     do
       cd $i
       md5sum -c checksums.md5
       cd ../../..
     done
elif [ $testcase = idiurnal-cycle ]
then
     for i in $( find tsmp_idiurnal-cycle/ -mindepth 2 -maxdepth 2 -type d )
     do
       cd $i
       md5sum -c checksums.md5
       cd ../../..
     done
fi

if [ $? -eq 0 ]; then
    echo "---------------------------------------"
    echo "Input files are downloaded succesfully"
    echo "---------------------------------------"
else
    echo "------------------------------------------------"
    echo "Download failed, please see teh checksum output"
    echo "------------------------------------------------"
    exit
fi

# move data to ../terrsysmp directory
if [ $testcase = cordex ]
then
     mv tsmp_eur11_eraint_eval_v2 ../
elif [ $testcase = nrw ]
then
     mv tsmp_nrw ../
elif [ $testcase = idealrtd ]
then
     mv tsmp_idealrtd ../
elif [ $testcase = idealscal ]
then
     mv tsmp_idealscal ../
elif [ $testcase = idiurnal-cycle ]
then
     mv tsmp_idiurnal-cycle ../
elif [ $testcase = cordex0275 ]
then
     mv tsmp_eur0275 ../
fi

