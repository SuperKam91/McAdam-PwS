#!/bin/bash

if [ "${1,,}" == "debug" ]; then
  echo "DEBUGGING MODE ENABLED"
  debug=true
fi

echo -n "Which system are you on (appci/CAMGRID/COSMOS/HPC/IOA/Universe): "
read which_sys
which_sys=`echo $which_sys | awk '{ print tolower($0) }'`

cp Makefile.template Makefile
cp McCompile.template McCompile

# Make sure these directories exist, otherwise it messes up the file system
mkdir -p bin lib 

if [ $debug ]; then
  # No optimization, turn on debugging
  IPO=""
  DEBUGFLAGS='-O0 -g -traceback -debug all -fp-speculation=off -ftrapuv'
  FDEBUGFLAGS='$(DEBUGFLAGS) -warn all -check noarg_temp_created -implicitnone -gen-interfaces -warn interfaces -fpe0 -fpe-all=0'
  CDEBUGFLAGS='$(DEBUGFLAGS) -check-uninit'
  AR="ar r"
else
  IPO="-ipo -xHost -O3 -w"
  DEBUGFLAGS=""
  FDEBUGFLAGS=""
  CDEBUGFLAGS=""
  AR="xiar r"
fi

# Some extra setup stuff for some systems
rm -f extra.setup
touch extra.setup
echo "Do 'source extra.setup' to load any modules required"
if [ $which_sys == "appci" ]; then
  echo module load intel/fce/15.0.0 >> extra.setup
  echo module load mpich/intel/2.4.1p1 >> extra.setup
  echo module load intel/mkle/15.0.0 >> extra.setup
elif [ $which_sys == "hpc" ]; then
  #echo module load cfitsio/intel/3.300 >> extra.setup
  #echo module unload intel/fce/12.1.10.319 >> extra.setup
  #echo module unload intel/cce/12.1.10.319 >> extra.setup
  #echo module load intel/fce/15.0.1.133 >> extra.setup
  #echo module load intel/cce/15.0.1.133 >> extra.setup
  # New for CSD3
  echo module load cfitsio-3.410-intel-17.0.4-4qrgkot >> extra.setup
  echo module unload intel/impi/2017.4/intel >> extra.setup
  echo module load intel/impi/2018.1/intel >> extra.setup
fi

# Parameters for headers for Makefile, customised for each system
if [ $which_sys == "appci" ]; then
  FC="mpif90 -fpp"
  CC="mpicc"
  CXX="mpiCC"
  LAPACKLIB="-mkl=parallel -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread"
  FFTWLIB="-lfftw3"
  EXTRAFLAGS="-shared-intel -openmp"
elif [ $which_sys == "camgrid" ]; then
  echo "Warning: Makefiles haven't been tested on CAMGRID"
  FC="mpif90 -fpp"
  CC="mpicc"
  CXX="mpicxx"
  LAPACKLIB="-mkl=parallel -L/usr/lib64 -llapack -lblas -lgfortran"
  FFTWLIB="-L/soft/star/lib -lfftw3"
  EXTRAFLAGS="-openmp -Bstatic"
elif [[ $which_sys == "cosmos" || $which_sys == "universe" ]]; then
  # Not sure these flags are all necessary but they're recommended on the COSMOS user-guide page so go with it
  FC="ifort -lmpi -g -align -ansi-alias -mcmodel=medium -traceback"
  CC="icc -lmpi -g -align -ansi-alias -mcmodel=medium -restrict"
  CXX="icpc -lmpi++ lmpi"
  LAPACKLIB="-mkl=parallel"
  FFTWLIB="-lfftw3"
  EXTRAFLAGS="-openmp"
  IPO=`echo $IPO | sed 's/-xHost/-axAVX,SSE4.2 -msse4.2/g'`
elif [ $which_sys == "hpc" ]; then
  FC="mpif90 -fpp"
  CC="mpicc"
  CXX="mpiicpc"
  LAPACKLIB="-mkl=parallel"
  FFTWLIB="-lfftw3"
  EXTRAFLAGS="-qopenmp"
elif [ $which_sys == "ioa" ]; then
  echo "Warning: Makefiles haven't been tested on IoA systems"
  FC="mpif90 -fpp"
  CC="icc -lmpi"
  CXX="icpc -lmpi"
  LAPACKLIB="-mkl=parallel"
  FFTWLIB="-lfftw3"
  EXTRAFLAGS="-openmp"
else
  echo "Don't recognise system name $which_sys"
  exit 1
fi

sed -i "s/THIS_FC/$FC/" Makefile
sed -i "s/THIS_CC/$CC/" Makefile
sed -i "s/THIS_CXX/$CXX/" Makefile
sed -i "s/THIS_AR/$AR/" Makefile
sed -i "s/THIS_LAPACKLIB/$LAPACKLIB/" Makefile
sed -i "s/THIS_FFTWLIB/$FFTWLIB/" Makefile
sed -i "s/THIS_EXTRAFLAGS/$EXTRAFLAGS/" Makefile
sed -i "s/THIS_DEBUGFLAGS/$DEBUGFLAGS/" Makefile
sed -i "s/THIS_FDEBUGFLAGS/$FDEBUGFLAGS/" Makefile
sed -i "s/THIS_CDEBUGFLAGS/$CDEBUGFLAGS/" Makefile
sed -i "s/THIS_IPO/$IPO/" Makefile

# Don't need to do any customising in the McCompile currently

exit 0

