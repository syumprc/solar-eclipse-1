//
//  pedifromsnps.cpp
//  
//
//  Created by Brian Donohue on 8/25/18.
//

#include <cstdlib>
#include <iostream>
#include "plinkio.h"
#include <fstream>
#include "solar.h"
#include <string>
using namespace std;
extern "C" void create_empirical_pedigree_allele_sharing(struct pio_file_t * input_file, const char * output_filename, const int per_chromosomet);
extern "C" void calculate_correlation_empirical_pedigree(struct pio_file_t * plink_file, const char * output_filename, const double alpha, const int per_chromosome);
extern "C" void write_epedigree_to_file(char * name_buffer, double * epedigree, pio_file_t * plink_file, const unsigned int n_subjects){
    ofstream output_stream(name_buffer);
    
    pio_sample_t * sample_i;
    pio_sample_t * sample_j;
    output_stream << "IDA,IDB,KIN\n";
    for(size_t j = 0; j < n_subjects ;j++){
        sample_j = fam_get_sample(&plink_file->fam_file, j);
        output_stream << sample_j->iid << "," << sample_j->iid <<  "," <<  epedigree[j*n_subjects + j] << endl;
        for(size_t i = j + 1; i < n_subjects ;i++){
            sample_i = fam_get_sample(&plink_file->fam_file, i);
            output_stream << sample_j->iid << "," << sample_i->iid <<  "," <<  epedigree[j*n_subjects + i] << endl;
            output_stream << sample_i->iid << "," << sample_j->iid <<  "," <<  epedigree[j*n_subjects + i] << endl;
        }
    }
    
    output_stream.close();
}
static void print_help(Tcl_Interp * interp){
    Solar_Eval(interp, "help pedifromsnps");
    
}
extern "C" int pedfromsnpsCmd(ClientData clientData, Tcl_Interp *interp,
                              int argc,const char *argv[]){
    
    int use_allele_sharing = 1;
    int per_chromosome = 0;
    const char * plink_filename = 0;
    const char * output_filename = 0;
    double alpha;
    
    for(unsigned arg = 1; arg < argc; arg++){
        if((!StringCmp(argv[arg], "--i", case_ins) || \
            !StringCmp(argv[arg], "--input", case_ins) \
            || !StringCmp(argv[arg], "-i", case_ins)\
            || !StringCmp(argv[arg], "-input", case_ins)) && arg + 1 < argc){
            plink_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--o", case_ins) || \
                  !StringCmp(argv[arg], "--output", case_ins) \
                  || !StringCmp(argv[arg], "-o", case_ins)\
                  || !StringCmp(argv[arg], "-output", case_ins)) && arg + 1 < argc){
            output_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "-corr", case_ins) || \
                  !StringCmp(argv[arg], "--corr", case_ins)) && arg + 1 < argc){
            alpha = atof(argv[++arg]);
            use_allele_sharing = 0;
        }else if(!StringCmp(argv[arg], "-king", case_ins) || \
                  !StringCmp(argv[arg], "--king", case_ins) ){
    
            use_allele_sharing = 1;
        }else if(!StringCmp(argv[arg], "-per-chromo", case_ins) || \
                 !StringCmp(argv[arg], "--per-chromo", case_ins)){
            per_chromosome = 1;
        }else if(!StringCmp(argv[arg], "-help", case_ins) || \
                 !StringCmp(argv[arg], "--help", case_ins) || \
                 !StringCmp(argv[arg], "help", case_ins) ){
            print_help(interp);
            return TCL_OK;
        }else{
            string error_message = "Invalid argument was entered: " + string(argv[arg]);
            RESULT_BUF(error_message.c_str());
            return TCL_ERROR;
        }
    }
    
    if(!plink_filename){
        RESULT_LIT("No plink file set specified with --i");
        return TCL_ERROR;
    }
    
    if(!output_filename){
        RESULT_LIT("No output filename specified with --o");
        return TCL_ERROR;
    }

    pio_file_t * plink_file = new pio_file_t;
    pio_status_t status;
    status = pio_open(plink_file, plink_filename);
    if(status != PIO_OK){
        if(status == P_FAM_IO_ERROR){
            RESULT_LIT("Error in loading .fam file");
            return TCL_ERROR;
        }else if (status == P_BIM_IO_ERROR){
            RESULT_LIT("Error in loading .bim file");
            return TCL_ERROR;
        }else if (status == P_BED_IO_ERROR){
            RESULT_LIT("Error in loading .bed file");
            return TCL_ERROR;
        }else{
            RESULT_LIT("Error loading plink file");
            return TCL_ERROR;
        }
    }
    int version = plink_file->bed_file.header.version;
    if (plink_file->bed_file.header.snp_order == BED_UNKNOWN_ORDER){
        pio_close(plink_file);
        RESULT_LIT("Error in the .bed snp order. Retry creation of file using a different version of plink");
        
        return TCL_ERROR;
    }/*else if (plink_file->bed_file.header.snp_order == BED_ONE_SAMPLE_PER_ROW){
        pio_close(plink_file);
        printf("In order to read efficiently the transpose of specified plink file must be performed\n");
        string transpose_filename = string(plink_filename) + string(".trans");
        string message = "Filename of transposed plink file is " + transpose_filename + "\n";
        printf(message.c_str());
        status = pio_transpose(plink_filename, transpose_filename.c_str());
        if(status != PIO_OK){
            RESULT_LIT("Error in creating transpose");
            return TCL_ERROR;
        }
        
        status = pio_open(plink_file, transpose_filename.c_str());
        
        if(status != PIO_OK){
            printf("Error in opening transposed plink file\n");
            return TCL_ERROR;
        }
    }*/
                   
    

    if(use_allele_sharing){
        create_empirical_pedigree_allele_sharing( plink_file,  output_filename, per_chromosome);
    }else{
        string message = "Creating correlation empirical pedigree with alpha=" + to_string(alpha);
        RESULT_BUF(message.c_str());
        try{
            calculate_correlation_empirical_pedigree(plink_file, output_filename, alpha, per_chromosome);
        }catch(...){
            RESULT_LIT("An error happened during empirical pedigree creation");
            pio_close(plink_file);
            return TCL_ERROR;
        }
        
    }
    pio_close(plink_file);
    delete plink_file;
    RESULT_LIT("Empirical pedigree creatipn is complete");
    return TCL_OK;
}
