project(solar-fortran)
enable_language(Fortran)
set(FORTRAN_FLAGS -m64 -O2 -fno-second-underscore -fexceptions)
set(FORTRAN_FILES alnorm.f  calc.f    clear.f    dcopyped.f  dfun.f    dppfa.f   fdist.f      fun_mehd.f  hpsort.f  linpack.f  mvncdf.f  pedtst.f   preopt.f  resid.f   smpoutput.f  tdist.f     yesno.f
astop.f   cdfchi.f  copyped.f  dcovar.f    digam.f   edftst.f  fit2dp.f     gamlog.f    inital.f  logo.f     newlik.f  phidens.f  preped.f  scor.f    spmpar.f     trigam.f
asycov.f  cdfnor.f  covar.f    ddfun.f     direct.f  evdlik.f  fortfiles.f  gaucdf.f    ipmpar.f  loop.f     optima.f  pinput.f   qdprog.f  search.f  sweep.f      upcase.f
beta.f    chisq.f   dasycov.f  delta7.f    dmean.f   evdout.f  fun.f        hessin.f    kin.f     mblank.f   output.f  ppnd.f     random.f  settab.f  symeig.f     writestr.f
)
add_library(solar-fortran STATIC ${FORTRAN_FILES})
target_compile_options(solar-fortran PRIVATE ${FORTRAN_FLAGS})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../lib)
