project(solar-c)
set (C_COMPILER gcc-7)
if(DEFINED ENV{C_COMPILER}) 
	set (C_COMPILER $ENV{C_COMPILER})
endif() 

SET(CMAKE_C_COMPILER ${C_COMPILER})
include_directories(./include)
set(C_FLAGS -std=c99 -m64 -O2 -I../include -DUSE_SAFELIB -fexceptions)
set(C_FILES ccuserid.c       cplink_converter.c  fdate.c        mehd.c  ncdf.c  nrutil.c    plinkio.c   safelib.c    tclgr.c     w_mehd.c
cpedifromsnps.c  eigstruc.c          getusername.c  mt.c    npdf.c  pipeback.c  plotpipe.c    testuidt.c)
add_library (solar-c STATIC ${C_FILES})
target_compile_options (solar-c PRIVATE ${C_FLAGS})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../lib)
