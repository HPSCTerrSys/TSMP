#!/bin/sh

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Before installing the third party libraries please read carefully the following points: "
echo " "
echo "- If you want to install the third party libraries on a linux machine with softwares installed as module, please run: "
echo " "
echo "module --force purge "
echo "in order to clean env variables."
echo "  "
echo "- If you want to install it on a local pc, please clean env variables such as PATH and LD_LIBRARY_PATH from the installed Netcdf, GribApi, HDF5, Openmpi (or MVAPICH), Silo, Hypre. "
echo "  "
echo "Note that if you have them already installed on your system and they are exported in .bashrc or .profile, the following installation might be messed up with your previous installation "
echo " "
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo " "
echo "would you like to check the env variables?"

read -p "would you like to check the env variables (y/n)? " answer
case $answer in
    y|Y )
         echo "***** check below the env variables ****** "
         sleep 10
         env
    ;;
    * )
        echo "No!"
    ;;
esac

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
read -p "if every thing is fine with env variables then type y to continue " answer
case $answer in
    y|Y )
       echo "env variables are clean"
    ;;
    * )
      echo "clean env variables and then continue "
      exit
    ;;
esac

bldsva=$(pwd)
[ ! -d "../lib" ] && mkdir ../lib
cd ../lib && dlib=$(pwd) 
#[ ! -d "$dlib/openmpi/lib" ] && 
rm -rf *
[ ! -d "$dlib/src" ] && mkdir $dlib/src
cd $dlib/src
dsrc=$(pwd)
mkdir -p $dlib/openmpi
mkdir -p $dlib/hypre
mkdir -p $dlib/silo
mkdir -p $dlib/gribapi
mkdir -p $dlib/netcdf
mkdir -p $dlib/hdf5
mkdir -p $dlib/zlib
mkdir -p $dlib/curl
mkdir -p $dlib/eccodes
mkdir -p $dlib/pnetcdf

log_file=$dlib/lib_install1.out
err_file=$dlib/lib_install1.err

#*****************************************
err=0
pist=`which gcc`
if [ $? -ne 0 ]; then
    echo "gcc is not installed"
    echo "install it via sudo apt or consult with your sys admin"
    err+=1
fi
pist=`which gfortran`
if [ $? -ne 0 ]; then
    echo "gfortran is not installed"
    echo "install it via sudo apt or consult with your sys admin"
    err+=1
fi
pist=`which g++`
if [ $? -ne 0 ]; then
    echo "g++ is not installed"
    echo "install it via sudo apt or consult with your sys admin"
    err+=1    
fi
pinst=`which ksh`
if [ $? -ne 0 ]; then
    echo "ksh is not installed"
    echo "install it via sudo apt or consult with your sys admin"
    err+=1
fi

pinst=`which curl`
if [ $? -ne 0 ]; then
    echo "curl is not installed"
    echo "install it via sudo apt or consult with your sys admin"
    echo "if via sudo apt: you need to install both curl and libcurl4-openssl-dev"
    err+=1
fi
pinst=`which m4`
if [ $? -ne 0 ]; then
    echo "m4 is not installed"
    echo "install it via sudo apt or consult with your sys admin"
    err+=1
fi
pinst=`which python`
if [ $? -ne 0 ]; then
    echo "python is not installed"
    echo "install (recommended Python 2.7.xx) it via sudo apt or consult with your sys admin"
    err+=1
fi
pinst=`which make`
if [ $? -ne 0 ]; then
    echo "make is not installed"
    echo "install it via sudo apt or consult with your sys admin"
    echo "if via sudo apt: after installing make: "
    echo  " make this symbolic link: ln -s /usr/bin/make /usr/bin/gmake"
    err+=1 
    echo "err = $err "
    echo "please install the missing library befor proceeding"
    [$err -ne 0 ] && exit
fi

#********* download the source files of the Libs ******************

pwget=`wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.0.tar.gz`
if [ $? -ne 0 ]; then
    echo "openmpi can not be downloaded!!!check wget command"
    exit
fi
tar -xvf openmpi-4.1.0.tar.gz >> $log_file 2>> $err_file
#****
pwget=`wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.6/src/hdf5-1.10.6.tar.gz`
if [ $? -ne 0 ]; then
    echo "HDF5 can not be downloaded!!!check wget command"
    exit
fi
tar -xvf hdf5-1.10.6.tar.gz >> $log_file 2>> $err_file
#*****
pwget=`wget https://github.com/Unidata/netcdf-c/archive/v4.6.1.zip`
if [ $? -ne 0 ]; then
    echo "Netcdf-c can not be downloaded!!!check wget command"
    exit
fi
unzip v4.6.1.zip >> $log_file 2>> $err_file
#*****
pwget=`wget https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.5.3.zip`
if [ $? -ne 0 ]; then
    echo "Netcdf-fortran can not be downloaded!!!check wget command"
    exit
fi
unzip v4.5.3.zip >> $log_file 2>> $err_file
#*****
pwget=`wget https://datapub.fz-juelich.de/slts/tsmp_testcases/TSMP_lib/grib_api-1.25.0-Source.tar.gz%3fapi=v2`
if [ $? -ne 0 ]; then
    echo "Grib-api can not be downloaded!!!check wget command"
    exit
fi
tar -xvf grib_api-1.25.0-Source.tar.gz?api=v2 >> $log_file 2>> $err_file
#*****
pwget=`wget https://github.com/LLNL/Silo/archive/refs/tags/4.10.2.tar.gz`
if [ $? -ne 0 ]; then
    echo "Silo can not be downloaded!!!check wget command"
    exit
fi
tar -xvf 4.10.2.tar.gz >> $log_file 2>> $err_file
#*****
pwget=`wget https://github.com/hypre-space/hypre/archive/refs/tags/v2.20.0.zip`
if [ $? -ne 0 ]; then
    echo "Hypre can not be downloaded!!!check wget command"
    exit
fi
unzip v2.20.0.zip >> $log_file 2>> $err_file
#****
pwget=`wget https://prdownloads.sourceforge.net/tcl/tcl8.6.10-src.tar.gz`
if [ $? -ne 0 ]; then
    echo "TCL can not be downloaded!!!check wget command"
    exit
fi
tar -xvf tcl8.6.10-src.tar.gz >> $log_file 2>> $err_file
#****
pwget=`wget https://github.com/madler/zlib/archive/refs/tags/v1.2.11.tar.gz`
if [ $? -ne 0 ]; then
    echo "zlib can not be downloaded!!!check wget command"
    exit
fi
tar -xvf v1.2.11.tar.gz >> $log_file 2>> $err_file
#****
pwget=`wget https://curl.haxx.se/download/curl-7.71.1.tar.gz`
if [ $? -ne 0 ]; then
    echo "curl can not be downloaded!!!check wget command"
    exit
fi
tar -xvf curl-7.71.1.tar.gz >> $log_file 2>> $err_file
echo "****** LIBs are downloaded and extracted in $dsrc *****"

#****
pwget=`wget https://github.com/ecmwf/eccodes/archive/refs/tags/2.18.0.zip`
if [ $? -ne 0 ]; then
    echo "ecCodes can not be downloaded!!!check wget command"
    exit
fi
unzip 2.18.0.zip >> $log_file 2>> $err_file
echo "****** LIBs are downloaded and extracted in $dsrc *****"

#****
pwget=`wget https://github.com/Parallel-NetCDF/PnetCDF/archive/refs/tags/checkpoint.1.12.1.zip`
if [ $? -ne 0 ]; then
    echo "PnetCDF can not be downloaded!!!check wget command"
    exit
fi
unzip checkpoint.1.12.1.zip >> $log_file 2>> $err_file
echo "****** LIBs are downloaded and extracted in $dsrc *****"

#*************** Zlib *************************************
cd $dsrc/zlib* 
./configure --prefix=$dlib/zlib >> $log_file 2>> $err_file
make >> $log_file 2>> $err_file
make install >> $log_file 2>> $err_file

#
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$dlib/zlib/lib"
export PATH="$dlib/zlib/bin:$PATH"

#******************* curl ***************************************
echo "***********************************************************"
echo "installing curl in $dlib/curl ......"
#
cd $dsrc/curl*
./configure --prefix="$dlib/curl" --with-zlib=$dlib/zlib  >> $log_file 2>> $err_file
make  >> $log_file 2>> $err_file
make install  >> $log_file 2>> $err_file
#
export PATH="$dlib/curl/bin:$PATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$dlib/curl/lib"

#*************** openmpi *************************************
echo "installing openmpi ........ "

/usr/bin/env >> $log_file 2>> $err_file
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib:/usr/lib/x86_64-linux-gnu"
env > $dlib/env.log
#
cd $dsrc/openmpi-*
./configure  --prefix="$dlib/openmpi" --with-pmix=internal >> $log_file 2>> $err_file
make all >> $log_file 2>> $err_file
make install >> $log_file 2>> $err_file

# export opnempi env vairable for installtion of other Apps

export PATH="$dlib/openmpi/bin:$PATH"
export CC=$dlib/openmpi/bin/mpicc
export FC=$dlib/openmpi/bin/mpif90
export F77=$dlib/openmpi/bin/mpif77

#***************** HDF5 *****************************************
echo "***********************************************************"
echo "installing hdf5 in $dlib/hdf5 ......"
#
cd $dsrc/hdf5-*
./configure --prefix=$dlib/hdf5 --enable-production --enable-hl --enable-fortran  >> $log_file 2>> $err_file
make  >> $log_file 2>> $err_file
make install   >> $log_file 2>> $err_file
#
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$dlib/hdf5/lib"

#******************* netcdf-c *****************************************
echo "***********************************************************"
echo "installing netcdf-c in $dlib/netcdf ..."
#
export CFLAGS="-I$dlib/hdf5/include -I$dlib/curl/include"
export LDFLAGS="-L$dlib/hdf5/lib -L$dlib/curl/lib -L$dlib/zlib/lib"
export LIBS="-lhdf5 -lhdf5_hl -lcurl -lz"
#
cd $dsrc/netcdf-c*
./configure --prefix=$dlib/netcdf  >> $log_file 2>> $err_file
make >> $log_file 2>> $err_file
make install >> $log_file 2>> $err_file
#
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$dlib/netcdf/lib"
export CFLAGS="-I$dlib/netcdf/include -I$dlib/hdf5/include -I$dlib/curl/include"
export FCFLAGS="-I$dlib/hdf5/include -I$dlib/curl/include -I$dlib/netcdf/include"
export CPPFLAGS="-I$dlib/netcdf/include -I$dlib/hdf5/include -I$dlib/curl/include"
export LDFLAGS="-L$dlib/hdf5/lib -L$dlib/netcdf/lib -L$dlib/curl/lib -L$dlib/zlib/lib"
export LIBS="-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl"
#********************** netcdf-fortran *****************************
echo "***********************************************************"
echo "installing netcdf-fortran in $dlib/netcdf ..."
#
cd $dsrc/netcdf-f*
./configure --prefix=$dlib/netcdf  >> $log_file 2>> $err_file
make >> $log_file 2>> $err_file
make install >> $log_file 2>> $err_file
#****************************** grib-api **************************
echo "*************************************************************"
echo "installing gribapi in $dlib/gribapi ..."
#
cd $dsrc/grib*
./configure --prefix=$dlib/gribapi --disable-jpeg >> $log_file 2>> $err_file
make >> $log_file 2>> $err_file
make install >> $log_file 2>> $err_file

#******************************* Silo ***************************
echo "***********************************************************"
echo "installing silo in $dlib/silo ..."
export FCFLAGS="-g -O2 -I$dlib/netcdf/include -I$dlib/hdf5/include -L$dlib/hdf5/lib -lhdf5_fortran -lhdf5"
#
cd $dsrc/Silo-*
./configure --prefix=$dlib/silo --with-hdf5=$dlib/hdf5/include,$dlib/hdf5/lib --enable-fortran --enable-shared  >> $log_file 2>> $err_file
make  >> $log_file 2>> $err_file
make install  >> $log_file 2>> $err_file

# *************************** Hypre ****************************
echo "***********************************************************"
echo "installing hypre in $dlib/hypre ..."
#
cd $dsrc/hypre-*/src
./configure --prefix=$dlib/hypre --with-MPI-lib-dirs=$dlib/openmpi  --enable-fortran --enable-shared >> $log_file 2>> $err_file
make  >> $log_file 2>> $err_file
make install  >> $log_file 2>> $err_file

#******************************* TCL *******************************
echo "**************************************************************"
echo "installing tcl in $dlib/tcl ..."
export LDFLAGS=""
export LIBS=""
export CFLAGS=""
export FCFLAGS=""
export CPPFLAGS=""
#
cd $dsrc/tcl*/unix
./configure --prefix=$dlib/tcl --enable-shared >> $log_file 2>> $err_file
make  >> $log_file 2>> $err_file
make install  >> $log_file 2>> $err_file
ln -s $dlib/tcl/bin/tclsh8.5 $dlib/tcl/bin/tclsh  >> $log_file 2>> $err_file
#

#******************************* ecCodes *******************************
echo "**************************************************************"
echo "installing ecCódes in $dlib/eccodes ..."
cd $dsrc/eccodes*
cmake  src/
make
ctest
make install

#******************************* PnetCDF *******************************
echo "**************************************************************"
echo "installing ecCódes in $dlib/eccodes ..."
cd $dsrc/PnetCDF-checkpoint.1.12.1/
make
make install

cat > $bldsva/machines/loadenv_x86 << end_loadenv
export PATH="$dlib/openmpi/bin:\$PATH"
export PATH="$dlib/tcl/bin:\$PATH"
export LD_LIBRARY_PATH="\$LD_LIBRARY_PATH:$dlib/hdf5/lib"
export LD_LIBRARY_PATH="\$LD_LIBRARY_PATH:$dlib/netcdf/lib"
export LD_LIBRARY_PATH="\$LD_LIBRARY_PATH:$dlib/gribapi/lib"
export LD_LIBRARY_PATH="\$LD_LIBRARY_PATH:$dlib/hypre/lib"
export LD_LIBRARY_PATH="\$LD_LIBRARY_PATH:$dlib/silo/lib"
export LD_LIBRARY_PATH="\$LD_LIBRARY_PATH:$dlib/tcl/lib"

export GRIB_DEFINITION_PATH=$dlib/gribapi/share/grib_api/definitions/
export GRIB_SAMPLES_PATH=$dlib/gribapi/share/grib_api/samples/

export CC=$dlib/openmpi/bin/mpicc
export FC=$dlib/openmpi/bin/mpif90
export F77=$dlib/openmpi/bin/mpif77


end_loadenv

echo " installation of third party lib finished!!!"
