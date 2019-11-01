project (solarmain)
cmake_minimum_required(VERSION 3.12 )
add_subdirectory (libsolar-fortran)
add_subdirectory (libsolar-c)
add_subdirectory (libsolar-cc)
add_executable(solarmain solarmain.c)
set (C_COMPILER gcc-7)
if(EXISTS "$ENV{C_COMPILER}") 
	set (C_COMPILER $ENV{C_COMPILER})
endif(EXISTS "$ENV{C_COMPILER}") 
	
set(C_FLAGS -std=c99 -m64 -O2 -I./include -DUSE_SAFELIB -fexceptions)
target_compile_options (solarmain PRIVATE ${C_FLAGS})
SET(CMAKE_C_LINK_EXECUTABLE ${C_COMPILER})
SET(CMAKE_C_COMPILER ${C_COMPILER})
set (MKL $ENV{MKLROOT})
set(LINK_OPTIONS -static  -std=c++11 -m64 -L${MKL}/lib/intel64 -fexceptions -pthread -fopenmp)
link_directories(BEFORE ./lib ${MKL}/lib/intel64 )


add_library(mkl_intel_lp64 STATIC IMPORTED)
set_property(TARGET mkl_intel_lp64 PROPERTY IMPORTED_LOCATION "${MKL}/lib/intel64/libmkl_intel_lp64.a")

add_library(mkl_gnu_thread STATIC IMPORTED)
set_property(TARGET mkl_gnu_thread PROPERTY IMPORTED_LOCATION "${MKL}/lib/intel64/libmkl_gnu_thread.a")

add_library(mkl_core STATIC IMPORTED)
set_property(TARGET mkl_core PROPERTY IMPORTED_LOCATION "${MKL}/lib/intel64/libmkl_core.a")

target_link_libraries(solarmain -Wl,--start-group solar-fortran solar-c solar-cc -Wl,--end-group  tcl8.4  z rt   -Wl,--start-group c gfortran gcc stdc++  m -Wl,--end-group dl  RicVolume  plinkio mkl_intel_lp64 pthread   mkl_gnu_thread mkl_core gomp) 

