#! /bin/ksh

# This is a script to download the input files needed to run the cordex.11 test case.
# Version: V1.0
# Version: V1.1: adding -l 8 in wget to download the whole directory tree 
# Version: V1.2: Adding the script for downloading the input data for nrw test case, 04.02.2021
# Version: V1.3: Adding the script for downloading the input data and pre-processing-tools for IdealRTD test case, 015.06.2021
#--------------------------------------
# Abouzar Ghasemi
# Jülich Supercomputing Centre (JSC)
# Jülich, Germany
# email: a.ghasemi@fz-juelich.de
#-------------------------------------

if [ -z "$1" ]
  then
    echo "!! No arguments supplied !!"
    echo "Please nrw for nrw test case"
    echo "Please cordex for cordex test case"
    exit 1
fi


testcase=$1
echo "input data downloads for $testcase"

if [ $testcase = cordex ]
then
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_eur11_eraint_eval/input/
elif [ $testcase = nrw ]
then
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_nrw/input/
elif [ $testcase = idealrtd ]
then
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_idealrtd/input/
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_idealrtd/pre-processing-tools/
      wget -x -l 8 -nH --cut-dirs=3 -e robots=off --recursive  --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_idealrtd/external/
fi

# Checking if the downloaded input files are corrupted.
if [ $testcase = cordex ]
then
     for i in $( ls -d tsmp_eur11_eraint_eval/input/* )
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
     mv tsmp_eur11_eraint_eval ../
elif [ $testcase = nrw ]
then
     mv tsmp_nrw ../
elif [ $testcase = idealrtd ]
then
     mv tsmp_idealrtd ../
fi

