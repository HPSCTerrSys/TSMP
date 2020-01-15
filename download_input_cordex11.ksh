#! /bin/ksh

# This is a script to download the input files needed to run the cordex.11 test case.
# Version: V1.0
# Written by Abouzar Ghasemi



wget -r -nH --cut-dirs=3 -e robots=off --no-parent --reject="index.html*" https://datapub.fz-juelich.de/slts/tsmp_testcases/data/tsmp_eur11_eraint_eval/input/


# Checking if the downloaded input files are corrupted.
for i in $( ls -d tsmp_eur11_eraint_eval/input/* )
do
     cd $i
     md5sum -c checksums.md5
     cd ../../..
done

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

