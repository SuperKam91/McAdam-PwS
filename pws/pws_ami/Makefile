

CXX= gcc
CXXL= gcc
AR = ar rs

LDFLAGS= -L../../Healpix_3.11/src/cxx/generic_gcc/lib -L/usr/lib64 -O3 -fopenmp -s

#LDFLAGS=  -O2 -s

#CFLAGS=-DLINUX -DNORANDOMBACKGROUND -DMYCONTOURSASS -Wno-write-strings -Wformat -I/usr/local/healpix/Healpix_3.11/src/cxx/cxxsupport/ -I/usr/local/healpix/Healpix_3.11/src/cxx/Healpix_cxx/ -I/usr/local/healpix/Healpix_3.11/src/cxx/libsharp/ -I./include -O2 -s -fomit-frame-pointer

#CFLAGSMAPCUT= -DLINUX -DNORANDOMBACKGROUND -DMYCONTOURSASS  -DFORCEDETECTIONINFORMATIVEPRIORS -Wno-write-strings -Wformat -I../../Healpix_3.11/src/cxx/generic_gcc/include -I./include -O3 -s -fomit-frame-pointer -fopenmp

CFLAGS= -DAMI -DLINUX -DNORANDOMBACKGROUND -DMYCONTOURSASS  -DFORCEDETECTIONINFORMATIVEPRIORS -Wno-write-strings -Wformat -I../../Healpix_3.11/src/cxx/generic_gcc/include -I./include -O3 -s -fomit-frame-pointer -fopenmp

#CFLAGS=-DLINUX -DNORANDOMBACKGROUND -DMYCONTOURSASS -DFORCEDETECTIONINFORMATIVEPRIORS -Wno-write-strings -Wformat  -I./include -O2 -s -fomit-frame-pointer

OUTPUT_FILE_PWS=PwSCore
OUTPUT_FILE_MAPC=MapCutter
OUTPUT_FILE_AMILIB=libpwsamilike.a

CPPFLAGS= $(CFLAGS)

OBJFILESPWS= bandmat.o cholesky.o evalue.o fft.o hholder.o jacobi.o mersenne.o mother.o newfft.o newmat1.o newmat2.o newmat3.o newmat4.o newmat5.o newmat6.o\
newmat7.o newmat8.o newmat9.o newmatex.o newmatnl.o newmatrm.o sfmt.o solution.o sort.o submat.o svd.o FFTW_factory.o FFTW_Storage.o FFTW_Traits.o myexcept.o\
ZEUS_Debug.o ZEUS_FourierMachine.o ZEUS_GaussianRandomGen.o ZEUS_General.o ZEUS_GlobalValuesStore.o ZEUS_Histogram.o ZEUS_InOut.o ZEUS_InOutBinFile.o\
ZEUS_InOutHealpixFits.o ZEUS_InOutPipeline.o ZEUS_InOutTxtFile.o ZEUS_Object.o ZEUS_ObjGauss.o ZEUS_ObjSZ.o ZEUS_ObjSZNgai.o ZEUS_PhysicsMath.o ZEUS_Priors.o\
ZEUS_StorageBaseManip.o ZEUS_WorkSpace.o\
PWS_Antenna.o PWS_AntGauss.o PWS_BackgroundProcessor.o PWS_BrentMinimizer.o PWS_CatalogueComp.o PWS_GlobalInfoStore.o PWS_Globals.o PWS_Main.o PWS_NestedSampler.o\
PWS_NestSamplCases.o PWS_NestSamplEllipsBound.o PWS_OddsEval.o PWS_PatchProcessor.o PWS_PowellSnakes.o PWS_Zone.o qa_assess_contours.o

OBJFILESMAPC=  bandmat.o cholesky.o evalue.o fft.o hholder.o jacobi.o mersenne.o mother.o newfft.o newmat1.o newmat2.o newmat3.o newmat4.o newmat5.o newmat6.o\
newmat7.o newmat8.o newmat9.o newmatex.o newmatnl.o newmatrm.o sfmt.o solution.o sort.o submat.o svd.o FFTW_factory.o FFTW_Storage.o FFTW_Traits.o myexcept.o\
ZEUS_Debug.o ZEUS_FourierMachine.o ZEUS_GaussianRandomGen.o ZEUS_General.o ZEUS_GlobalValuesStore.o ZEUS_Histogram.o ZEUS_InOut.o ZEUS_InOutBinFile.o\
ZEUS_InOutHealpixFits.o ZEUS_InOutPipeline.o ZEUS_InOutTxtFile.o ZEUS_Object.o ZEUS_ObjGauss.o ZEUS_ObjSZ.o ZEUS_ObjSZNgai.o ZEUS_PhysicsMath.o ZEUS_Priors.o\
ZEUS_StorageBaseManip.o ZEUS_WorkSpace.o\
MC_Antenna.o MC_AntennaGaussian.o MC_HealpixCutter.o MC_HPixMapPatchCutting.o MC_ObjFactory.o MC_PatchFinder.o MC_PixExtractProc.o

all:$(OUTPUT_FILE_AMILIB) $(OUTPUT_FILE_MAPC)

mapcut:$(OUTPUT_FILE_MAPC)

ami:$(OUTPUT_FILE_AMILIB)

$(OUTPUT_FILE_AMILIB): $(OBJFILESPWS)
	$(AR) $(OUTPUT_FILE_AMILIB) $(OBJFILESPWS) 

$(OUTPUT_FILE_MAPC): $(OBJFILESMAPC)
	$(CXX) $(OBJFILESMAPC) $(LDFLAGS) -lstdc++ -lm -lfftw3 -lhealpix_cxx -lcxxsupport -lcfitsio -lsharp -lc_utils -lfftpack -o $(OUTPUT_FILE_MAPC)

#$(OUTPUT_FILE_PWS): $(OBJFILESPWS)
#	$(CXX) $(OBJFILESPWS) $(LDFLAGS) -lstdc++ -lm -lfftw3 -lhealpix_cxx -lcxxsupport -lcfitsio -lsharp -lc_utils -lfftpack -o $(OUTPUT_FILE_PWS)

clean:
	rm -f *.o $(OUTPUT_FILE_AMILIB) $(OUTPUT_FILE_MAPC)

bandmat.o: bandmat.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatrc.h include/controlw.h
cholesky.o: cholesky.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h
evalue.o: evalue.cpp include/include.h include/newmatap.h \
 include/newmat.h include/include.h include/boolean.h include/myexcept.h \
 include/newmatrm.h include/precisio.h
fft.o: fft.cpp include/include.h include/newmatap.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h
FFTW_factory.o: FFTW_factory.cpp include/FFTW_Factory.h \
 include/FFTW_Storage.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/TypeManip.h include/static_check.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_Exceptions.h \
 include/ZEUS_Strings.h include/ZEUS_General.h include/FFTW_Traits.h \
 /usr/include/fftw3.h \
 include/ZEUS_SingletonSerialize.h include/FFTW_Globals.h \
 include/FFTW_Strings.h
FFTW_Storage.o: FFTW_Storage.cpp include/FFTW_Storage.h \
 include/ZEUS_StorageBaseManip.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/FFTW_Traits.h \
 /usr/include/fftw3.h \
 include/ZEUS_SingletonSerialize.h include/FFTW_Globals.h \
 include/FFTW_Strings.h
FFTW_Traits.o: FFTW_Traits.cpp include/ZEUS_SingletonSerialize.h \
 include/Threads.h include/SmartPtr.h include/SmallObj.h \
 include/Singleton.h include/TypeManip.h include/static_check.h \
 include/ZEUS_Strings.h include/ZEUS_Exceptions.h include/ZEUS_General.h \
 include/FFTW_Traits.h \
 /usr/include/fftw3.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_SingletonSerialize.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Device.h \
 include/FFTW_Storage.h include/FFTW_Traits.h
hholder.o: hholder.cpp include/include.h include/newmatap.h \
 include/newmat.h include/include.h include/boolean.h include/myexcept.h
jacobi.o: jacobi.cpp include/include.h include/newmatap.h \
 include/newmat.h include/include.h include/boolean.h include/myexcept.h \
 include/precisio.h include/newmatrm.h
MC_Antenna.o: MC_Antenna.cpp include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/MC_Antenna.h include/ZEUS_WorkSpace.h
MC_AntennaGaussian.o: MC_AntennaGaussian.cpp include/MC_AntennaGaussian.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Strings.h \
 include/ZEUS_Exceptions.h include/ZEUS_General.h include/ZEUS_Debug.h \
 include/MC_Antenna.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h
MC_HealpixCutter.o: MC_HealpixCutter.cpp \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Strings.h \
 include/ZEUS_Exceptions.h include/ZEUS_General.h include/ZEUS_Debug.h \
 include/ZEUS_InOutTxtFile.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_StorageBaseManip.h \
 include/ZEUS_LokiVariant.h include/Typelist.h include/NullType.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/ZEUS_InOutBinFile.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alm_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/xcomplex.h \
 include/MC_Defaults.h include/MC_HealpixCutter.h \
 include/MC_HPixMapPatchCutting.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 include/ZEUS_GlobalValuesStore.h include/MC_HPixMapPatchCutting.h
MC_HPixMapPatchCutting.o: MC_HPixMapPatchCutting.cpp \
 include/ZEUS_Exceptions.h include/ZEUS_Strings.h include/ZEUS_General.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Exceptions.h include/ZEUS_Debug.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/FFTW_Traits.h \
 include/ZEUS_SingletonSerialize.h include/FFTW_Globals.h \
 include/FFTW_Strings.h include/ZEUS_Strings.h \
 include/MC_HPixMapPatchCutting.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_WorkSpace.h include/ZEUS_GlobalValuesStore.h \
 include/MC_HealpixCutter.h include/MC_HPixMapPatchCutting.h \
 include/MC_PatchFinder.h include/MC_PixExtractProc.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Storage.h include/FFTW_Device.h \
 include/ZEUS_GaussianRandomGen.h include/randomc.h include/MC_Defaults.h \
 include/MC_Antenna.h include/MC_ObjFactory.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/ZEUS_ObjGauss.h \
 include/ZEUS_ObjSZ.h include/ZEUS_ObjSZNgai.h
MC_ObjFactory.o: MC_ObjFactory.cpp include/MC_ObjFactory.h \
 include/ZEUS_Object.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_Factory.h include/Factory.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_ObjGauss.h include/ZEUS_ObjSZ.h include/ZEUS_ObjSZNgai.h
MC_PatchFinder.o: MC_PatchFinder.cpp include/ZEUS_PhysicsMath.h \
 include/PWS_Strings.h include/ZEUS_Strings.h include/ZEUS_Exceptions.h \
 include/ZEUS_General.h include/ZEUS_Debug.h include/MC_PatchFinder.h \
 include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_StorageBaseManip.h \
 include/ZEUS_LokiVariant.h include/Typelist.h include/NullType.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/MC_HPixMapPatchCutting.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 include/ZEUS_GlobalValuesStore.h include/MC_HPixMapPatchCutting.h
MC_PixExtractProc.o: MC_PixExtractProc.cpp \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alm.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alm_powspec_tools.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/xcomplex.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alm_healpix_tools.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 include/MC_HPixMapPatchCutting.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_General.h include/ZEUS_Strings.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_Exceptions.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h \
 include/ZEUS_GlobalValuesStore.h include/MC_HealpixCutter.h \
 include/MC_HPixMapPatchCutting.h include/MC_PixExtractProc.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_GaussianRandomGen.h include/randomc.h \
 include/MC_Defaults.h include/MC_Antenna.h include/MC_ObjFactory.h \
 include/ZEUS_Object.h include/ZEUS_Factory.h include/Factory.h \
 include/ZEUS_ObjGauss.h include/ZEUS_ObjSZ.h include/ZEUS_ObjSZNgai.h \
 include/MC_AntennaGaussian.h include/MC_AntennaPixWindow.h \
 include/ZEUS_Debug.h
mersenne.o: mersenne.cpp include/randomc.h
mother.o: mother.cpp include/randomc.h
myexcept.o: myexcept.cpp include/include.h include/myexcept.h \
 include/include.h
newfft.o: newfft.cpp include/newmatap.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h
newmat1.o: newmat1.cpp include/newmat.h include/include.h \
 include/boolean.h include/myexcept.h
newmat2.o: newmat2.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatrc.h include/controlw.h
newmat3.o: newmat3.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatrc.h include/controlw.h
newmat4.o: newmat4.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatrc.h include/controlw.h
newmat5.o: newmat5.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatrc.h include/controlw.h
newmat6.o: newmat6.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatrc.h include/controlw.h
newmat7.o: newmat7.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatrc.h include/controlw.h
newmat8.o: newmat8.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatrc.h include/controlw.h include/precisio.h
newmat9.o: newmat9.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatio.h include/newmat.h include/newmatrc.h \
 include/controlw.h
newmatex.o: newmatex.cpp include/include.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h
newmatnl.o: newmatnl.cpp include/newmatap.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatnl.h
newmatrm.o: newmatrm.cpp include/newmat.h include/include.h \
 include/boolean.h include/myexcept.h include/newmatrm.h
PWS_Antenna.o: PWS_Antenna.cpp include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/PWS_Antenna.h include/ZEUS_WorkSpace.h \
 include/ZEUS_Factory.h include/Factory.h \
 include/ZEUS_SingletonSerialize.h include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_GlobalValuesStore.h include/ZEUS_FourierMachine.h \
 include/FFTW_Traits.h include/FFTW_Globals.h include/FFTW_Strings.h \
 include/FFTW_Storage.h include/FFTW_Device.h include/ZEUS_Bins.h \
 include/ZEUS_Object.h include/PWS_Defaults.h include/ZEUS_Debug.h
PWS_AntGauss.o: PWS_AntGauss.cpp include/PWS_AntGauss.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Strings.h \
 include/ZEUS_Exceptions.h include/ZEUS_General.h include/ZEUS_Debug.h \
 include/PWS_Antenna.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_Factory.h include/Factory.h \
 include/ZEUS_SingletonSerialize.h include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_GlobalValuesStore.h include/ZEUS_FourierMachine.h \
 include/FFTW_Traits.h include/FFTW_Globals.h include/FFTW_Strings.h \
 include/FFTW_Storage.h include/FFTW_Device.h include/ZEUS_Bins.h \
 include/ZEUS_Object.h include/PWS_Defaults.h
PWS_BackgroundProcessor.o: PWS_BackgroundProcessor.cpp \
 include/PWS_BackgroundProcessor.h include/ZEUS_GaussianRandomGen.h \
 include/Threads.h include/SmartPtr.h include/SmallObj.h \
 include/Singleton.h include/TypeManip.h include/static_check.h \
 include/randomc.h include/ZEUS_Strings.h include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_Exceptions.h include/ZEUS_General.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_PhysicsMath.h \
 include/PWS_Strings.h include/ZEUS_Debug.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_GlobalValuesStore.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_Bins.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/PWS_Defaults.h \
 include/PWS_GlobalInfoStore.h include/PWS_Zone.h include/PWS_Antenna.h \
 include/ZEUS_InOut.h include/ZEUS_WorkSpace.h include/ZEUS_Debug.h \
 include/PWS_PatchProcessor.h include/ZEUS_Priors.h \
 include/PWS_BackgroundProcessor.h
PWS_BrentMinimizer.o: PWS_BrentMinimizer.cpp include/PWS_BrentMinimizer.h \
 include/PWS_OddsEval.h include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_GlobalValuesStore.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_Bins.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/PWS_Defaults.h \
 include/PWS_GlobalInfoStore.h include/PWS_Zone.h include/PWS_Antenna.h \
 include/PWS_PatchProcessor.h include/ZEUS_Priors.h \
 include/PWS_BackgroundProcessor.h include/ZEUS_GaussianRandomGen.h \
 include/randomc.h
PWS_CatalogueComp.o: PWS_CatalogueComp.cpp include/PWS_CatalogueComp.h \
 include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_GlobalValuesStore.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_Bins.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/PWS_Defaults.h \
 include/PWS_GlobalInfoStore.h include/PWS_Zone.h include/PWS_Antenna.h \
 include/ZEUS_InOut.h
PWS_GlobalInfoStore.o: PWS_GlobalInfoStore.cpp include/ZEUS_ObjSZNgai.h \
 include/ZEUS_Object.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_Factory.h include/Factory.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_InOutPipeline.h include/PWS_CatalogueComp.h \
 include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 include/ZEUS_GlobalValuesStore.h include/ZEUS_FourierMachine.h \
 include/FFTW_Traits.h include/FFTW_Globals.h include/FFTW_Strings.h \
 include/FFTW_Storage.h include/FFTW_Device.h include/ZEUS_Bins.h \
 include/PWS_Defaults.h include/PWS_GlobalInfoStore.h include/PWS_Zone.h \
 include/PWS_Antenna.h include/PWS_GlobalInfoStore.h
PWS_Globals.o: PWS_Globals.cpp include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_GlobalValuesStore.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_Bins.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/PWS_Defaults.h
PWS_Main.o: PWS_Main.cpp include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_GlobalValuesStore.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_Bins.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/PWS_Defaults.h \
 include/PWS_BackgroundProcessor.h include/ZEUS_GaussianRandomGen.h \
 include/randomc.h include/PWS_Globals.h include/PWS_GlobalInfoStore.h \
 include/PWS_Zone.h include/PWS_Antenna.h include/PWS_PatchProcessor.h \
 include/ZEUS_Priors.h include/PWS_BackgroundProcessor.h \
 include/ZEUS_Priors.h include/ZEUS_Debug.h
PWS_NestedSampler.o: PWS_NestedSampler.cpp include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_General.h include/ZEUS_Strings.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_Exceptions.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h \
 include/ZEUS_Debug.h include/PWS_NestedSampler.h \
 include/ZEUS_GaussianRandomGen.h include/randomc.h
PWS_NestSamplCases.o: PWS_NestSamplCases.cpp include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_GlobalValuesStore.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_Bins.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/PWS_Defaults.h \
 include/PWS_NestSamplCases.h include/ZEUS_Priors.h \
 include/PWS_OddsEval.h include/PWS_Globals.h \
 include/PWS_GlobalInfoStore.h include/PWS_Zone.h include/PWS_Antenna.h \
 include/PWS_PatchProcessor.h include/PWS_BackgroundProcessor.h \
 include/ZEUS_GaussianRandomGen.h include/randomc.h \
 include/PWS_NestedSampler.h include/PWS_NestSamplEllipsBound.h \
 include/newmatap.h include/newmat.h include/include.h include/boolean.h \
 include/myexcept.h
PWS_NestSamplEllipsBound.o: PWS_NestSamplEllipsBound.cpp \
 include/PWS_NestSamplEllipsBound.h include/newmatap.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/PWS_NestedSampler.h include/ZEUS_PhysicsMath.h \
 include/PWS_Strings.h include/ZEUS_Strings.h include/ZEUS_Exceptions.h \
 include/ZEUS_General.h include/ZEUS_Debug.h \
 include/ZEUS_GaussianRandomGen.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/randomc.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h
PWS_OddsEval.o: PWS_OddsEval.cpp include/PWS_PatchProcessor.h \
 include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_GlobalValuesStore.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_Bins.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/PWS_Defaults.h \
 include/ZEUS_Priors.h include/PWS_GlobalInfoStore.h include/PWS_Zone.h \
 include/PWS_Antenna.h include/PWS_BackgroundProcessor.h \
 include/ZEUS_GaussianRandomGen.h include/randomc.h \
 include/PWS_OddsEval.h include/PWS_PatchProcessor.h
PWS_PatchProcessor.o: PWS_PatchProcessor.cpp include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_General.h include/ZEUS_Strings.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_Exceptions.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h \
 include/PWS_BackgroundProcessor.h include/ZEUS_GaussianRandomGen.h \
 include/randomc.h include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 include/ZEUS_InOut.h include/ZEUS_GlobalValuesStore.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_Bins.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/PWS_Defaults.h \
 include/PWS_GlobalInfoStore.h include/PWS_Zone.h include/PWS_Antenna.h \
 include/PWS_PatchProcessor.h include/ZEUS_Priors.h \
 include/PWS_BackgroundProcessor.h include/PWS_OddsEval.h \
 include/PWS_PatchProcessor.h include/PWS_BrentMinimizer.h \
 include/PWS_OddsEval.h include/PWS_PowellSnakes.h \
 include/PWS_BrentMinimizer.h include/PWS_NestSamplCases.h \
 include/PWS_NestedSampler.h include/PWS_NestSamplEllipsBound.h \
 include/newmatap.h include/newmat.h include/include.h include/boolean.h \
 include/myexcept.h include/ZEUS_Debug.h include/qa_assess_contours.h
PWS_PowellSnakes.o: PWS_PowellSnakes.cpp include/PWS_PowellSnakes.h \
 include/PWS_OddsEval.h include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_GlobalValuesStore.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_Bins.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/PWS_Defaults.h \
 include/PWS_GlobalInfoStore.h include/PWS_Zone.h include/PWS_Antenna.h \
 include/PWS_PatchProcessor.h include/ZEUS_Priors.h \
 include/PWS_BackgroundProcessor.h include/ZEUS_GaussianRandomGen.h \
 include/randomc.h include/PWS_BrentMinimizer.h
PWS_Zone.o: PWS_Zone.cpp include/PWS_GlobalInfoStore.h \
 include/PWS_Globals.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_GlobalValuesStore.h \
 include/ZEUS_FourierMachine.h include/FFTW_Traits.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h include/ZEUS_Bins.h include/ZEUS_Object.h \
 include/ZEUS_Factory.h include/Factory.h include/PWS_Defaults.h \
 include/PWS_Zone.h include/PWS_Antenna.h include/PWS_AntGauss.h \
 include/ZEUS_ObjGauss.h include/ZEUS_ObjSZNgai.h include/ZEUS_ObjSZ.h \
 include/ZEUS_Debug.h include/PWS_Zone.h
qa_assess_contours.o: qa_assess_contours.cpp include/ZEUS_PhysicsMath.h \
 include/PWS_Strings.h include/ZEUS_Strings.h include/ZEUS_Exceptions.h \
 include/ZEUS_General.h include/ZEUS_Debug.h include/qa_assess_contours.h
sfmt.o: sfmt.cpp include/sfmt.h include/randomc.h
solution.o: solution.cpp include/include.h include/boolean.h \
 include/myexcept.h include/include.h include/solution.h \
 include/boolean.h include/myexcept.h
sort.o: sort.cpp include/include.h include/newmatap.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h
submat.o: submat.cpp include/include.h include/newmat.h include/include.h \
 include/boolean.h include/myexcept.h include/newmatrc.h \
 include/controlw.h
svd.o: svd.cpp include/include.h include/newmatap.h include/newmat.h \
 include/include.h include/boolean.h include/myexcept.h \
 include/newmatrm.h include/precisio.h
ZEUS_Debug.o: ZEUS_Debug.cpp \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 include/ZEUS_Debug.h
ZEUS_FourierMachine.o: ZEUS_FourierMachine.cpp \
 include/ZEUS_FourierMachine.h include/ZEUS_StorageBaseManip.h \
 include/Threads.h include/SmartPtr.h include/SmallObj.h \
 include/Singleton.h include/TypeManip.h include/static_check.h \
 include/ZEUS_LokiVariant.h include/Typelist.h include/NullType.h \
 include/ZEUS_Exceptions.h include/ZEUS_Strings.h include/ZEUS_General.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h \
 include/FFTW_Traits.h include/ZEUS_SingletonSerialize.h \
 include/FFTW_Globals.h include/FFTW_Strings.h include/FFTW_Storage.h \
 include/FFTW_Device.h
ZEUS_GaussianRandomGen.o: ZEUS_GaussianRandomGen.cpp \
 include/ZEUS_GaussianRandomGen.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/randomc.h include/ZEUS_Strings.h
ZEUS_General.o: ZEUS_General.cpp include/ZEUS_General.h \
 include/ZEUS_Strings.h
ZEUS_GlobalValuesStore.o: ZEUS_GlobalValuesStore.cpp \
 include/ZEUS_Exceptions.h include/ZEUS_Strings.h include/ZEUS_General.h \
 include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_Exceptions.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h \
 include/ZEUS_InOutTxtFile.h include/ZEUS_InOut.h \
 include/ZEUS_GlobalValuesStore.h
ZEUS_Histogram.o: ZEUS_Histogram.cpp include/ZEUS_Histogram.h \
 include/ZEUS_GaussianRandomGen.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/randomc.h include/ZEUS_Strings.h \
 include/ZEUS_General.h
ZEUS_InOutBinFile.o: ZEUS_InOutBinFile.cpp include/ZEUS_InOutBinFile.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alm_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/xcomplex.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 include/ZEUS_General.h include/ZEUS_Strings.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_Exceptions.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h
ZEUS_InOut.o: ZEUS_InOut.cpp include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_General.h include/ZEUS_Strings.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_Exceptions.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h \
 include/ZEUS_InOutBinFile.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alm_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/xcomplex.h \
 include/ZEUS_InOut.h include/ZEUS_InOutHealpixFits.h \
 include/ZEUS_InOutPipeline.h include/ZEUS_InOutTxtFile.h
ZEUS_InOutHealpixFits.o: ZEUS_InOutHealpixFits.cpp \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 include/ZEUS_InOutHealpixFits.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 include/ZEUS_General.h include/ZEUS_Strings.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_Exceptions.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h
ZEUS_InOutPipeline.o: ZEUS_InOutPipeline.cpp \
 include/ZEUS_GaussianRandomGen.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/randomc.h include/ZEUS_Strings.h \
 include/ZEUS_InOutPipeline.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_General.h include/ZEUS_SingletonSerialize.h \
 include/ZEUS_Exceptions.h include/ZEUS_StorageBaseManip.h \
 include/ZEUS_LokiVariant.h include/Typelist.h include/NullType.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h
ZEUS_InOutTxtFile.o: ZEUS_InOutTxtFile.cpp include/ZEUS_InOutTxtFile.h \
 include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_General.h include/ZEUS_Strings.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_Exceptions.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h
ZEUS_Object.o: ZEUS_Object.cpp include/ZEUS_Object.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_Factory.h include/Factory.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h
ZEUS_ObjGauss.o: ZEUS_ObjGauss.cpp include/ZEUS_ObjGauss.h \
 include/ZEUS_Object.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_Factory.h include/Factory.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h
ZEUS_ObjSZ.o: ZEUS_ObjSZ.cpp include/ZEUS_ObjSZ.h include/ZEUS_Object.h \
 include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_Factory.h include/Factory.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h
ZEUS_ObjSZNgai.o: ZEUS_ObjSZNgai.cpp include/ZEUS_ObjSZNgai.h \
 include/ZEUS_Object.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h include/ZEUS_Factory.h include/Factory.h \
 include/ZEUS_SingletonSerialize.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h
ZEUS_PhysicsMath.o: ZEUS_PhysicsMath.cpp include/ZEUS_Exceptions.h \
 include/ZEUS_Strings.h include/ZEUS_General.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_Exceptions.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h \
 include/ZEUS_PhysicsMath.h include/ZEUS_WorkSpace.h include/ZEUS_Debug.h
ZEUS_Priors.o: ZEUS_Priors.cpp include/ZEUS_Priors.h include/ZEUS_InOut.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/trafos.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/vec3.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/datatypes.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/error_handling.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/pointing.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rotmatrix.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/arr.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/alloc_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/math_utils.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_base.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_tables.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/rangeset.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/safe_cast.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/healpix_map_fitsio.h \
 ../../Healpix_3.11/src/cxx/generic_gcc/include/fitshandle.h \
 include/ZEUS_General.h include/ZEUS_Strings.h \
 include/ZEUS_SingletonSerialize.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_Exceptions.h \
 include/ZEUS_StorageBaseManip.h include/ZEUS_LokiVariant.h \
 include/Typelist.h include/NullType.h include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/ZEUS_PhysicsMath.h include/PWS_Strings.h include/ZEUS_Debug.h \
 ZEUS_Priors.inc
ZEUS_StorageBaseManip.o: ZEUS_StorageBaseManip.cpp \
 include/ZEUS_StorageBaseManip.h include/Threads.h include/SmartPtr.h \
 include/SmallObj.h include/Singleton.h include/TypeManip.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h
ZEUS_WorkSpace.o: ZEUS_WorkSpace.cpp include/ZEUS_WorkSpace.h \
 /usr/include/fftw3.h \
 include/TypeManip.h include/ZEUS_Exceptions.h include/ZEUS_Strings.h \
 include/ZEUS_General.h include/ZEUS_StorageBaseManip.h include/Threads.h \
 include/SmartPtr.h include/SmallObj.h include/Singleton.h \
 include/static_check.h include/ZEUS_LokiVariant.h include/Typelist.h \
 include/NullType.h include/ZEUS_PhysicsMath.h include/PWS_Strings.h \
 include/ZEUS_Debug.h
