//  
//
//  Created by Brian Donohue on 8/25/18.
//

#include"plinkio.h"
#include<stdio.h>
#include "safelib.h"
#include <math.h>
#include <stdio.h>
void write_epedigree_to_file(char * name_buffer, double * epedigree, struct pio_file_t * plink_file, const unsigned int n_subjects);
static void calculate_empirical_pedigree_correlation(char * name_buffer, struct pio_file_t * input_file,  unsigned int * missing_snp_count,\
                                                                   double * empirical_pedigree, const double alpha, const size_t start_index,\
                                                                   const size_t n_subjects, const size_t n_snps){
    
    snp_t * buffer = (snp_t*)malloc(sizeof(snp_t)*n_subjects);
    for(size_t col = 0; col < n_subjects; col++){
        for(size_t row = col ; row < n_subjects; row++){
            empirical_pedigree[col*n_subjects + row] = 0.0;
           missing_snp_count[col*n_subjects + row] = 0;

        }
    }
    
    for(size_t snp_index = 0; snp_index < n_snps; snp_index++){
        pio_next_row(input_file,buffer);
        unsigned int missing_snps[n_subjects];
        for(unsigned i = 0; i < n_subjects; i++) missing_snps[i] = 0;
        unsigned int  n_subjects_missing = 0;
        unsigned int  hetero_count = 0;
        unsigned int homo_count = 0;
        for(size_t col = 0; col < n_subjects; col++){
            switch(buffer[col]){
                case 1:
                    hetero_count++;
                    break;
		 case 2:
                    homo_count++;
                    break;	

	          case 3:
		     missing_snps[col] = 1;
                    n_subjects_missing++;
                    break;
		  default:
			break;
           }

        }
        for(size_t col = 0; col < n_subjects; col++){
	
            for(size_t row = col; row < n_subjects ; row++){
                if(missing_snps[col] == 1 || missing_snps[row] == 1) missing_snp_count[col*n_subjects + row]++;
            }
        }
        const  double freq = (double) (2*homo_count + hetero_count)/(2.0*(double)(n_subjects - n_subjects_missing));
	 if(freq == 1.0 || freq == 0.0) continue;
        const double variance = 2.0*freq*(1.0 - freq);
        const double factor =  pow(variance, alpha);
	printf("%f \n", factor);
	  const double mean = 2.0*freq;
        for(size_t col = 0; col < n_subjects; col++){
            const snp_t  snp_i = buffer[col];
            if(snp_i == 3) continue;
            const double component_i = (snp_i - mean)*factor;
            empirical_pedigree[col*n_subjects + col] += (snp_i - mean)*component_i;
            for(size_t row = col + 1; row < n_subjects; row++) if(buffer[row] != 3) empirical_pedigree[col*n_subjects + row] += component_i*(buffer[row] - mean);
            
        }
    }
    
    for(size_t col = 0; col < n_subjects; col++){
        for(size_t row = col; row < n_subjects; row++){
          printf("%f \n" , empirical_pedigree[col*n_subjects + row] );
	 empirical_pedigree[col*n_subjects + row] /= (n_snps - missing_snp_count[col*n_subjects + row]);
	}
    }
    double norms[n_subjects];
    for(size_t col = 0 ; col < n_subjects ; col++){
       norms[col] = sqrt(fabs(empirical_pedigree[col*n_subjects + col]));
    }
    for(size_t col = 0; col < n_subjects; col++){
        for(size_t row = col; row < n_subjects ; row++){
	 empirical_pedigree[col*n_subjects + row] /= norms[row]*norms[col];
   	 }
    }
write_epedigree_to_file(name_buffer, empirical_pedigree, input_file, n_subjects);    
    
    free(buffer);
}
void calculate_correlation_empirical_pedigree(struct pio_file_t * plink_file, const char * output_filename, const double alpha, const int per_chromosome){
    const size_t n_subjects = plink_file->bed_file.header.num_samples;
    double * epedigree = (double*)malloc(sizeof(double)*n_subjects*n_subjects);
    unsigned int * snp_count = (unsigned int *)malloc(sizeof(unsigned int)*n_subjects*n_subjects);
    if(per_chromosome){
        struct pio_locus_t * locus;
        locus = bim_get_locus(&plink_file->bim_file, 0);
        unsigned char chromosome = locus->chromosome;
        size_t start = 0;
        size_t end = 0;
        while(start < plink_file->bed_file.header.num_loci){
            size_t snp = start;
            while(locus->chromosome == chromosome && snp < plink_file->bed_file.header.num_loci) locus = bim_get_locus(&plink_file->bim_file, snp++);
            end = snp;
            size_t snp_batch_size = end - start;
            char name_buffer[200];
            sprintf(name_buffer, "%s.chr%u.csv", output_filename, (unsigned int )chromosome);
            calculate_empirical_pedigree_correlation(name_buffer,plink_file,  snp_count,\
                                                      epedigree,  alpha, start,\
                                                     plink_file->bed_file.header.num_samples, snp_batch_size);
            
            if(end != plink_file->bed_file.header.num_loci){
                chromosome = locus->chromosome;
            }
            start = end;
        }
    }else{
       char name_buffer[200];
        sprintf(name_buffer, "%s", output_filename);
        calculate_empirical_pedigree_correlation(name_buffer, plink_file,  snp_count,\
                                                 epedigree,  alpha, 0,\
                                                 plink_file->bed_file.header.num_samples, plink_file->bed_file.header.num_loci);
    }
    
    free(epedigree);
    free(snp_count);
    
}
static void calculate_empirical_pedigree_allele_sharing_quantities(struct pio_file_t * input_file,unsigned int * homo_share, unsigned int  * hetero_share,\
                                                                   unsigned int * hetero_count_i, unsigned int * hetero_count_j,\
                                                                     const size_t start_index,\
                                                                   const size_t n_subjects, const size_t n_snps){
    
    snp_t * buffer = (snp_t*)malloc(sizeof(snp_t)*n_subjects);
    for(size_t col = 0; col < n_subjects; col++){
        for(size_t row = col ; row < n_subjects; row++){
            homo_share[col*n_subjects + row] = 0;
            hetero_share[col*n_subjects + row] = 0;
            hetero_count_i[col*n_subjects + row] = 0;
            hetero_count_j[col*n_subjects + row] = 0;
           // freq_square_sum[col*n_subjects + row] = 0.0;
        }
    }
    
    for(size_t snp_index = 0; snp_index < n_snps; snp_index++){
        pio_next_row(input_file,buffer);
        

        
        //size_t n_subjects_missing = 0;
        for(size_t col = 0; col < n_subjects; col++){
            int snp_j = buffer[col];
            if(snp_j == 3) continue;
            for(size_t row = col ; row < n_subjects; row++){
                int snp_i = buffer[row];
                if(snp_i == 3) continue;
                if((snp_j  == 0 && snp_i  == 2) ||  (snp_i== 0 &&  snp_j == 2) ) homo_share[col*n_subjects + row]++;
                if(snp_j  == 1 && snp_i == 1) hetero_share[col*n_subjects + row]++;
                if(snp_i== 1)  hetero_count_i[col*n_subjects + row]++;
                if(snp_j == 1)  hetero_count_j[col*n_subjects + row]++;
            }
        }
        
        

    }
    
}


static inline double calculate_robust_kinship(const unsigned int hetero_share,const  unsigned int  homo_share,const unsigned int  hetero_count_i,const  unsigned int  hetero_count_j){
    unsigned int denom = hetero_count_i;
    if(hetero_count_i > hetero_count_j) denom = hetero_count_j;
    
    return ((hetero_share-2.0*homo_share)/(denom) + 1.0 - ((hetero_count_i + hetero_count_j)/(2.0*denom)) );
}

void create_empirical_pedigree_allele_sharing(struct pio_file_t * input_file, const char * output_filename, const int per_chromosome){
    struct pio_locus_t * locus;
    pio_status_t status;
    const unsigned int n_subjects = input_file->bed_file.header.num_samples;
    const unsigned int n_snps = input_file->bed_file.header.num_loci;
    unsigned int * hetero_count_i = (unsigned int * )malloc(sizeof(unsigned int)*n_subjects*n_subjects);
    unsigned int * hetero_count_j = (unsigned int * )malloc(sizeof(unsigned int)*n_subjects*n_subjects);
    unsigned int * hetero_share= (unsigned int * )malloc(sizeof(unsigned int)*n_subjects*n_subjects);
    unsigned int * homo_share = (unsigned int * )malloc(sizeof(unsigned int)*n_subjects*n_subjects);


   // double * squared_variance = (double * )malloc(sizeof(double)*n_subjects*n_subjects);
    double * empirical_pedigree = (double * )malloc(sizeof(double)*n_subjects*n_subjects);
    if(per_chromosome){
        struct pio_locus_t * locus;
        locus = bim_get_locus(&input_file->bim_file, 0);
        unsigned char chromosome = locus->chromosome;
        size_t start = 0;
        size_t end = 0;
        while(start < input_file->bed_file.header.num_loci){
            size_t snp = start;
            while(locus->chromosome == chromosome && snp < input_file->bed_file.header.num_loci) locus = bim_get_locus(&input_file->bim_file, snp++);
            end = snp;
            size_t snp_batch_size = end - start;
            calculate_empirical_pedigree_allele_sharing_quantities(input_file,homo_share,  hetero_share,\
              hetero_count_i,  hetero_count_j,  start, n_subjects, snp_batch_size);
            
            for(size_t col =0 ; col < n_subjects ; col++){
               // empirical_pedigree[col*n_subjects + col] = 1.0;
                for(size_t row = col ; row < n_subjects; row++) {
                    empirical_pedigree[col*n_subjects + row] = calculate_robust_kinship(hetero_share[col*n_subjects + row], homo_share[col*n_subjects + row],hetero_count_i[col*n_subjects + row],hetero_count_j[col*n_subjects + row]);
                }
            }
            
            char name_buffer[200];
            sprintf(name_buffer, "%s.chr%u.csv", output_filename, (unsigned int )chromosome);
            write_epedigree_to_file(name_buffer, empirical_pedigree, input_file, n_subjects);
            if(end != input_file->bed_file.header.num_loci){
                chromosome = locus->chromosome;
            }
            start = end;
        }
    }else{
        calculate_empirical_pedigree_allele_sharing_quantities(input_file,homo_share,  hetero_share,\
                                                               hetero_count_i,  hetero_count_j, 0, n_subjects, n_snps);

        
        for(size_t col =0 ; col < n_subjects ; col++){
            // empirical_pedigree[col*n_subjects + col] = 1.0;
            for(size_t row = col ; row < n_subjects; row++) {
                empirical_pedigree[col*n_subjects + row] = calculate_robust_kinship(hetero_share[col*n_subjects + row], homo_share[col*n_subjects + row],hetero_count_i[col*n_subjects + row],hetero_count_j[col*n_subjects + row]);
            }
        }
        
        char name_buffer[200];
        sprintf(name_buffer, "%s", output_filename);
        write_epedigree_to_file(name_buffer, empirical_pedigree, input_file, n_subjects);
    }



    free(hetero_count_i);
    free(hetero_count_j);
    free(hetero_share);
    free(homo_share);

    free(empirical_pedigree);
}

    
