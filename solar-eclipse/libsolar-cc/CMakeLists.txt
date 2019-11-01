project(solar-cc)
set (CC_COMPILER g++-7)
if(EXISTS "$ENV{CC_COMPILER}") 
	set (CC_COMPILER $ENV{CC_COMPILER})
endif(EXISTS "$ENV{CC_COMPILER}") 

SET(CMAKE_CXX_COMPILER ${CC_COMPILER})
include_directories(../include)
set(CC_FLAGS -O2 -std=c++11 -fopenmp  -DGPU_TOOLS  -m64  -I../include -DUSE_SAFELIB -fexceptions -DTR1 -DLINUX64 -DUSE_SAFELIB  -I${MKLROOT}/include/ -DEIGEN_USE_MKL_ALL)
set(CC_FILES alnormc.cc               define.cc      fisherpedigree.cc  ibd.cc          loglike.cc     model.cc           parameter.cc          scale.cc        solarfilecmd.cc        token.cc
annotate_gwas.cc         drand.cc       fphi.cc            ibdoption.cc    map.cc         mu.cc              pedifromsnps.cc       scan.cc         solar_mle_setup.cc     trait.cc
ccsearch.cc              eqvar.cc       freq.cc            ibs.cc          marker.cc      nifti_assemble.cc  pedigree.cc           scratchfile.cc  solar-trait-reader.cc  verbosity.cc
chi.cc                   evd.cc         function.cc        inorm_nifti.cc  mathmatrix.cc  nifti_to_csv.cc    phenotypes.cc         simqtl.cc       split_pheno_file.cc    voxel.cc
constraint.cc               gwas.cc            key.cc          matrix.cc      normal.cc          plink_converter.cc    snp.cc          sporadic.cc            zscore.cc
covariate.cc             expression.cc  help.cc            listmaker.cc    maximize.cc    omega.cc           power.cc              solar.cc        tablefile.cc
create_fake_pedigree.cc  field.cc       howclose.cc        loadsave.cc     mibd.cc        option.cc          reorder_phenotype.cc  solarfile.cc    tablefilecmd.cc)
add_library (solar-cc STATIC ${CC_FILES})
target_compile_options (solar-c PRIVATE ${CC_FLAGS})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../lib)

