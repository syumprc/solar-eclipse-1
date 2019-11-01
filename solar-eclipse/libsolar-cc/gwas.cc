#include "solar.h"
#include <Eigen/Dense>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <array>
#include <iterator>
#include "plinkio.h"
#include <unordered_map>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include "solar-trait-reader.h"
using namespace std;

void load_phi2_matrix(Tcl_Interp * interp);
static const unsigned GWAS_BATCH_SIZE = 6000;
static const int PERMUTATION_BATCH_SIZE = 1000;
extern "C" void cdfchi_ (int*, double*, double*, double*, double*,
			 int*, double*);
			 
			 
 double chicdf(double chi, double df){
	double p, q, bound;
	int status = 0;
	int which = 1;
	
	
	cdfchi_ (&which, &p, &q, &chi, &df, &status, &bound);
	
	return q;
}
double gwas_chicdf(double chi, double df){
return chicdf(chi, df);
}
extern "C" void symeig_ (int*, double*, double*, double*, double*, int*);
void calculate_eigenvectors (double * phi2, double * eigenvalues, double * eigenvectors ,int n)
{
    double* e =  new double[n];
    memset(e, 0, sizeof(double)*n);
    int * info = new int;
    *info  = 0;
    symeig_(&n, phi2, eigenvalues, e, eigenvectors, info);
    delete [] e;
    delete [] info;
}


static int create_gwas_matrices(Tcl_Interp * interp, pio_file_t * plink_file, int * & plink_index_map, Eigen::VectorXd & trait_vector, Eigen::MatrixXd & phi2,  int * & snp_values, const char * trait_name, const char * phenotype_filename, const char * plink_filename, vector<string> & snp_names){
    const char * errmsg = 0;
    cout << "Reading phenotype file...\n";
    SolarFile * file = SolarFile::open("gwas", phenotype_filename, &errmsg);
    
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    int field_count;
    const char ** fields = file->names(&field_count, &errmsg);
    
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    file->start_setup(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    file->setup("id", &errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    char ** file_data;
    vector<string> phenotype_ids;
    
    while (0 != (file_data = file->get (&errmsg))){
        
        phenotype_ids.push_back(string(file_data[0]));
        
        
    }
    file->rewind(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    file->start_setup(&errmsg);
    
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    file->setup(trait_name, &errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    double * trait_data = new double[phenotype_ids.size()];
    double row_value;
    size_t row_index = 0;

    

    vector<string> missing_trait_ids;
    while (0 != (file_data = file->get (&errmsg))){
        string str_data = file_data[0];
        if(str_data.length() != 0){
            trait_data[row_index] = strtod((const char*)file_data[0], NULL);
        }else{
            missing_trait_ids.push_back(phenotype_ids[row_index]);
            trait_data[row_index] = nan("");
        }
        row_index++;
    }
    int * raw_plink_data;
    
    size_t n_snps;
   // size_t n_snp_rows = phenotype_ids.size();
    vector<string> plink_ids;
    if(plink_file){
      //  pio_file_t * plink_file = new pio_file_t;
       // pio_status_t status = pio_open(plink_file, plink_filename);
	pio_status_t status;
       /* 
        if(status != PIO_OK){
            pio_close(plink_file);
            if(status == P_FAM_IO_ERROR){
                RESULT_LIT("Error in loading .fam file");
                return TCL_ERROR;
            }else if (status == P_BIM_IO_ERROR){
                RESULT_LIT("Error in loading .bim file");
                return TCL_ERROR;
            }else if (status == P_BED_IO_ERROR){
                RESULT_LIT("Error in loading .bed file");
                return TCL_ERROR;
            }
        }
        */
        n_snps  = plink_file->bed_file.header.num_loci;
        size_t n_snp_rows = plink_file->bed_file.header.num_samples;
        pio_sample_t * sample;
        pio_locus_t * locus;
       // int * snp_index_map = new int[n_snp_rows];
        row_index = 0;
        for(size_t i = 0; i < n_snp_rows; i++){
            sample = fam_get_sample(&plink_file->fam_file, i);
            string id  = string(sample->iid);
            vector<string>::iterator find_iter = find(phenotype_ids.begin(), phenotype_ids.end(), id);
            if(find_iter != phenotype_ids.end()){
                vector<string>::iterator trait_iter = find(missing_trait_ids.begin(),missing_trait_ids.end(), id);
                if(trait_iter == missing_trait_ids.end()){
                    plink_ids.push_back(id);
                   // snp_index_map[i] = row_index++;
                }//else{
                 //   snp_index_map[i] = -1;
              //  }
            }//else{
                //snp_index_map[i] = -1;
           // }
        }

        for(size_t i = 0; i < n_snps; i++){
            locus = bim_get_locus(&plink_file->bim_file, i);
            snp_names.push_back(string(locus->name));
        }
      //  raw_plink_data = new int[plink_ids.size()*n_snps];
      //  snp_t * snp_buffer = new snp_t[n_snp_rows];
       // for(size_t snp = 0 ; snp < n_snps; snp++){
       //     pio_next_row(plink_file,snp_buffer);
       //     for(size_t row =0 ;row < n_snp_rows; row++){
          //      int index = snp_index_map[row];
         //       if(index != -1) raw_plink_data[snp*plink_ids.size() +index] = snp_buffer[row];
         //   }
       // }
	//delete [] snp_index_map;
     //   delete [] snp_buffer;
      //  pio_close(plink_file);
     //   delete plink_file;
    }else{
        file->rewind(&errmsg);
        if(errmsg){
            RESULT_LIT(errmsg);
            return TCL_ERROR;
        }
        file->start_setup(&errmsg);
        
        if(errmsg){
            RESULT_LIT(errmsg);
            return TCL_ERROR;
        }
        for(size_t field = 0; field < field_count; field++){
            if(strstr(fields[field], "snp_")){
                string snp_name = string(fields[field]);
                snp_name = snp_name.substr(4, snp_name.length() - 4);
                snp_names.push_back(snp_name);
                file->setup(fields[field], &errmsg);
                if(errmsg){
                    RESULT_LIT(errmsg);
                    return TCL_ERROR;
                }
            }
        }
        n_snps = snp_names.size();
        size_t n_snp_rows = phenotype_ids.size();
        int * temp_plink_data = new int[n_snps*phenotype_ids.size()];
        size_t row_index = 0;
        vector<string> missing_plink_ids;
        while (0 != (file_data = file->get (&errmsg))){
            size_t missing_count = 0;
            for(size_t snp_index = 0; snp_index < n_snps; snp_index++){
                string str_data = file_data[snp_index];
		if(str_data.length() != 0){
                    int snp_value = stoi(string(file_data[snp_index]));
                    temp_plink_data[snp_index*phenotype_ids.size() + row_index] = snp_value;
                    if(snp_value == 3) missing_count++;
                }else{
                    temp_plink_data[snp_index*phenotype_ids.size() + row_index] = 3;
                    missing_count++;
                }
            }
            if(missing_count == n_snps){
                missing_plink_ids.push_back(phenotype_ids[row_index]);
            }
            row_index++;
        }
        int * plink_index_map = new int[phenotype_ids.size()];
        int map_index = 0;
        for(size_t row = 0; row < phenotype_ids.size(); row++){
            string id = phenotype_ids[row];
            vector<string>::iterator find_iter = find(missing_trait_ids.begin(), missing_trait_ids.end(), id);
            if(find_iter == missing_trait_ids.end()){
                find_iter = find(missing_plink_ids.begin(), missing_plink_ids.end(), id);
                if(find_iter == missing_plink_ids.end()){
                    plink_index_map[row] = map_index++;
                    plink_ids.push_back(id);
                }else{
                    plink_index_map[row] = -1;
                }
            }else{
                plink_index_map[row] = -1;
            }
        }
        raw_plink_data = new int[n_snps*plink_ids.size()];
        for(size_t row = 0 ;row < phenotype_ids.size(); row++){
            int index = plink_index_map[row];
            if(index != -1){
                for(size_t snp = 0; snp < n_snps ; snp++){
                    raw_plink_data[index + snp*plink_ids.size()] = temp_plink_data[row + snp*phenotype_ids.size()];
                }
            }
        }
        delete [] temp_plink_data;

    }

    vector<string> ids;
    vector<int> ibdids;
    SolarFile * pedindex_file = SolarFile::open("gwas", "pedindex.out", &errmsg);
    
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }

    pedindex_file->start_setup(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    pedindex_file->setup("id", &errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
	
    vector<string> shared_id_set = plink_ids;
  //  char ** file_data;
    int ibdid = 1;
    vector<string> ids_out;
    while (0 != (file_data = pedindex_file->get (&errmsg))){
        
        string id = string(file_data[0]);
        vector<string>::iterator find_iter = find(shared_id_set.begin(), shared_id_set.end(), id);
        if(find_iter != shared_id_set.end()){
            ids.push_back(id);
         //   shared_id_set.erase(find_iter);
            ibdids.push_back(ibdid);
        }
        
        ibdid++;
    }
    trait_vector.resize(ids.size());
    size_t id_index = 0;
    size_t row = 0;
    for(size_t row_index = 0; row_index < phenotype_ids.size(); row_index++){
        const double trait_value = trait_data[row_index];
        if(trait_value == trait_value){
            string trait_id = phenotype_ids[row_index];
            vector<string>::iterator find_iter = find(ids.begin(), ids.end(), trait_id);
            if(find_iter != ids.end()){
                trait_vector(distance(ids.begin(), find_iter)) = trait_value;
            }
        }
    }
    if(!plink_file){
    int * index_map = new int[plink_ids.size()];
    for(size_t i = 0; i < plink_ids.size(); i++){
        string plink_id = plink_ids[i];
        vector<string>::iterator find_iter = find(ids.begin(), ids.end(), plink_id);
        if(find_iter != ids.end()){
            index_map[i] = distance(ids.begin(), find_iter);
        }else{
            index_map[i] = -1;
        }
    }
    snp_values = new int[ids.size()*n_snps];
    for(size_t row = 0; row < plink_ids.size(); row++){
        int index = index_map[row];
        if(index != -1){
            for(size_t col = 0; col < n_snps; col++){
                snp_values[col*ids.size() + index] = raw_plink_data[col*plink_ids.size() + row];
            }
        }
    }
    delete [] raw_plink_data;
    }else{

	plink_index_map = new int[ids.size()];
	pio_sample_t * sample;
	vector<string> plink_id_list;
	for(unsigned i = 0; i < plink_file->bed_file.header.num_samples; i++){
		sample = fam_get_sample(&plink_file->fam_file, i);
		plink_id_list.push_back(string(sample->iid));
	}
	for(unsigned row = 0; row < ids.size(); row++){
		plink_index_map[row] =  distance(plink_id_list.begin(), find(plink_id_list.begin(), plink_id_list.end(), ids[row]));
	}
    }
 
    
    phi2.resize(ids.size(), ids.size());
    
    Matrix* pm2;
    
    pm2 = Matrix::find("phi2");
    if (!pm2) {
        Solar_Eval(interp, "matrix load phi2.gz phi2");
        pm2 = Matrix::find("phi2");
        if(!pm2){
            RESULT_LIT("phi2 matrix could not be loaded");
            return TCL_ERROR;
        }
    }

    for(int col_idx = 0; col_idx < ids.size(); col_idx++){
        for(int row_idx = col_idx; row_idx < ids.size(); row_idx++) {
            const double phi2_value = pm2->get(ibdids[row_idx], ibdids[col_idx]);
            phi2(row_idx, col_idx) = phi2_value;
            phi2(col_idx, row_idx) = phi2_value;
            
        }
        
    }
   
    return TCL_OK;
}

static int * get_permutation_indices(int * permutation_indices, int n_rows){
	iota (permutation_indices, permutation_indices + n_rows, 0);
	random_shuffle(permutation_indices, permutation_indices + n_rows);
//	return permutation_indices;

}
typedef struct gwas_data{
	double beta;
	double pvalue;
	double SE;
	double chi;
	double SD;
    double h2r;
    double loglik;
	gwas_data & operator = (const gwas_data & var) {
		beta = var.beta;
		pvalue = var.pvalue;
		SE = var.SE;
		chi = var.chi;
		SD= var.SD;
        h2r = var.h2r;
        loglik = var.loglik;
	}
	
}gwas_data;
/*
static inline double calculate_constraint(const double x){
	const double x_squared = x*x;
	return x_squared/(1 + x_squared);
	//return exp(x)/(1.0+exp(x));
}
static inline double calculate_dconstraint(const double x){
	return 2.0*x*pow((1.0 + x*x), -2);
	//return exp(x)*pow(1 + exp(x), -2);
}
static inline double calculate_ddconstraint(const double x){
	const double x_squared = x*x;
	return -2.0*(3.0*x_squared - 1.0)*pow((x_squared + 1.0), -3);
	//return -exp(x)*(exp(x) - 1.0)*pow(exp(x) + 1.0, -3);
}
static inline double calculate_dloglik(Eigen::VectorXd lambda_minus_one,Eigen::VectorXd residual_squared, Eigen::VectorXd sigma,const double variance){
	const double part_one = variance*lambda_minus_one.dot(sigma);
	const double part_two = variance*lambda_minus_one.dot(residual_squared.cwiseProduct(sigma.cwiseAbs2()));
	return -0.5*(part_one-part_two);
}
static inline double calculate_ddloglik(Eigen::VectorXd lambda_minus_one,Eigen::VectorXd residual_squared, Eigen::VectorXd sigma, const double variance){
	Eigen::VectorXd lambda_minus_one_squared = lambda_minus_one.cwiseAbs2();
	Eigen::VectorXd sigma_squared = sigma.cwiseAbs2();
	const double part_one = variance*variance*lambda_minus_one_squared.dot(sigma_squared);
	const double part_two = 2.0*variance*variance*lambda_minus_one_squared.dot(residual_squared.cwiseProduct(sigma.cwiseProduct(sigma_squared)));

	return -0.5*(-part_one + part_two);
}
static inline double calculate_ddloglik_with_constraint(const double t,const double dloglik, const double ddloglik){

	return pow(calculate_dconstraint(t), 2)*ddloglik + calculate_ddconstraint(t)*dloglik;
}

static inline double calculate_dloglik_with_constraint(const double t, const double dloglik){

	return calculate_dconstraint(t)*dloglik;
}

static inline Eigen::VectorXd calculate_sigma(const double t, const double variance, Eigen::MatrixXd aux){
	Eigen::VectorXd theta(2);
	const double h2 = calculate_constraint(t);
	theta(0) = 1.0 - h2;
	theta(1) = h2;
	theta = theta*variance;
	return aux*theta;
}
static inline double calculate_quick_loglik(Eigen::VectorXd Sigma, const unsigned N){
	return -0.5*(-log(Sigma.array()).sum() + N);
}

static void gwas_maximize_newton_raphson_method(gwas_data * result, Eigen::VectorXd Y, Eigen::MatrixXd X, Eigen::MatrixXd U, gwas_data null_result, const int precision){
    	Eigen::VectorXd raw_A = X.col(0).cwiseAbs2();
    	Eigen::VectorXd raw_B = X.col(1).cwiseProduct(X.col(0));
   	Eigen::VectorXd raw_C = X.col(1).cwiseAbs2();
   	Eigen::VectorXd raw_D = Y.cwiseProduct(X.col(0));
    	Eigen::VectorXd raw_E = Y.cwiseProduct(X.col(1));

	double t = 1;
	double h2 = 0.5;
	Eigen::VectorXd theta(2);
	theta(0) = 0.5;
	theta(1) = 0.5;
	double A,B,C,D,E;
	Eigen::VectorXd omega = (U*theta).cwiseInverse();
	A = raw_A.dot(omega);
	B = raw_B.dot(omega);
	C = raw_C.dot(omega);
	D = raw_D.dot(omega);
	E = raw_E.dot(omega);
	Eigen::VectorXd beta(2);
	double beta_denom  = A*C - B*B;
    	beta(0) = C*D - B*E;
    	beta(1) = A*E - B*D;
   	beta /= beta_denom;	
	Eigen::VectorXd residual = Y - X*beta;
	Eigen::VectorXd residual_squared = residual.cwiseAbs2();
	double variance = residual_squared.dot(omega)/Y.rows();
	Eigen::VectorXd sigma = omega*pow(variance, -1);
	double loglik = calculate_quick_loglik(sigma, Y.rows());
	Eigen::VectorXd lambda_minus_one = (U.col(1).array() - 1.0).matrix();
	double dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
	double ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
	double score = calculate_dloglik_with_constraint( t,dloglik);
	double hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
	double delta = -score/hessian;
	double new_h2 = 0;
	if (delta == delta){
		t += delta;
		new_h2 = calculate_constraint(t);
	}
	const double end = pow(10, -precision);

	while( delta == delta && fabs(new_h2 - h2) >= end){
		h2 = new_h2;

		theta(0) = 1.0 - h2;
		theta(1) = h2;
		omega = (U*theta).cwiseInverse();
		A = raw_A.dot(omega);
		B = raw_B.dot(omega);
		C = raw_C.dot(omega);
		D = raw_D.dot(omega);
		E = raw_E.dot(omega);
		beta_denom  = A*C - B*B;
    		beta(0) = C*D - B*E;
    		beta(1) = A*E - B*D;
   		beta /= beta_denom;	
		residual = Y - X*beta;
		residual_squared = residual.cwiseAbs2();
		variance = residual_squared.dot(omega)/Y.rows();
		sigma = omega*pow(variance, -1);
	 	dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
		ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
		score = calculate_dloglik_with_constraint( t,dloglik);
		hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
		delta = -score/hessian;
		if(delta == delta){
			t += delta;
			new_h2 = calculate_constraint(t);
		}

	}
	if(delta == delta){
		h2 = new_h2;
		theta(0) = 1.0 - h2;
		theta(1) = h2;
		omega = (U*theta).cwiseInverse();
		A = raw_A.dot(omega);
		B = raw_B.dot(omega);
		C = raw_C.dot(omega);
		D = raw_D.dot(omega);
		E = raw_E.dot(omega);
		beta_denom  = A*C - B*B;
    		beta(0) = C*D - B*E;
    		beta(1) = A*E - B*D;
   		beta /= beta_denom;	
		residual = Y - X*beta;
		residual_squared = residual.cwiseAbs2();
		variance = residual_squared.dot(omega)/Y.rows();
		sigma = omega/variance;
	}
	//result_loglik = loglik;
	//result_variance = variance;
    loglik = calculate_quick_loglik(sigma, Y.rows());
    double chi = 2.0*(loglik - null_result.loglik);
  //  Eigen::MatrixXd beta_covariance_matrix = X.transpose()*omega.asDiagonal()*X;
   // double denom = beta_covariance_matrix(0,0)*beta_covariance_matrix(1, 1) - beta_covariance_matrix(1,0)*beta_covariance_matrix(0,1);
    
    //  gwas_data result;
	//cout << "chi : " << chi << " delta: " << delta << " h2: " << h2 << endl;
    if(chi >= 0.0 && chi == chi && delta == delta){
        
        //        Eigen::MatrixXd hessian = calculate_REML_hessian(residual, SNP, Sigma, eigenvalues, variance);
        //      Eigen::VectorXd standard_errors = hessian.inverse().diagonal().cwiseSqrt();
        
       // double beta_se = beta_covariance(1,1);
        result->SE = sqrt(A/beta_denom);
        result->beta = beta(1);
        result->chi = chi;
        result->SD = sqrt(variance);
//        result->SE = 0.0;//beta_se;//standard_errors(0);
        result->pvalue = chicdf(result->chi , 1);
        result->h2r = h2;
        result->loglik = loglik;
    }else{
        result->beta = 0.0;
        result->chi = 0.0;
        result->pvalue = 1.0;
        result->h2r = null_result.h2r;
        result->loglik = null_result.loglik;
        result->SD = null_result.SD;
        result->SE = 0.0;

    }
	

}
double calculate_GWAS_loglik(Eigen::VectorXd SY2, Eigen::VectorXd lambda, Eigen::VectorXd Sigma){
    // double var_sum  = 0.0;
    // double log_sum = 0.0;
    // Eigen::VectorXd theta(2);
    // theta(0) = 1.0-h2r;
    // theta(1) = h2r;
    // Eigen::VectorXd Sigma  = (lambda*h2r).array() + e2;
    // const double variance_sum  = SY2.cwiseQuotient(variance);
    
    return -0.5*(log(Sigma.array()).sum() + SY2.rows());
    // for(int i = 0 ; i < SY.rows(); i++)
    
}

static gwas_data gwas_maximize_newton_raphson_method_null_model(Eigen::VectorXd Y, Eigen::VectorXd mean_column, Eigen::MatrixXd U, const int precision){
	gwas_data result;
	Eigen::VectorXd raw_A = Y.cwiseProduct(mean_column);
	Eigen::VectorXd raw_B = mean_column.cwiseAbs2();

	double t = 1;
	double h2 = 0.5;
	Eigen::VectorXd theta(2);
	theta(0) = 0.5;
	theta(1) = 0.5;
	double A,B;
	Eigen::VectorXd omega = (U*theta).cwiseInverse();
	A = raw_A.dot(omega);
	B = raw_B.dot(omega);
	Eigen::VectorXd beta(1);

    	beta(0) = A/B;	
	Eigen::VectorXd residual = Y - mean_column*beta;
	Eigen::VectorXd residual_squared = residual.cwiseAbs2();
	double variance = residual_squared.dot(omega)/Y.rows();
	Eigen::VectorXd sigma = omega*pow(variance, -1);
	double loglik = calculate_quick_loglik(sigma, Y.rows());
	Eigen::VectorXd lambda_minus_one = (U.col(1).array() - 1.0).matrix();
	double dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
	double ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
	double score = calculate_dloglik_with_constraint( t,dloglik);
	double hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
	double delta = -score/hessian;
	double new_h2 = 0;
	if (delta == delta){
		t += delta;
		new_h2 = calculate_constraint(t);
	}
	const double end = pow(10, -precision);

	while( delta == delta && fabs(new_h2 - h2) >= end){
		h2 = new_h2;

		theta(0) = 1.0 - h2;
		theta(1) = h2;
		omega = (U*theta).cwiseInverse();
		A = raw_A.dot(omega);
		B = raw_B.dot(omega);
		//Eigen::VectorXd beta(2);
    		beta(0) = A/B;
		residual = Y - mean_column*beta;
		residual_squared = residual.cwiseAbs2();
		variance = residual_squared.dot(omega)/Y.rows();
		sigma = omega*pow(variance, -1);
	 	dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
		ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
		score = calculate_dloglik_with_constraint( t,dloglik);
		hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
		delta = -score/hessian;
		if(delta == delta){
			t += delta;
			new_h2 = calculate_constraint(t);
		}

	}
	if(delta == delta){
		h2 = new_h2;
		theta(0) = 1.0 - h2;
		theta(1) = h2;
		omega = (U*theta).cwiseInverse();
		A = raw_A.dot(omega);
		B = raw_B.dot(omega);
    		beta(0) = A/B;
		residual = Y - mean_column*beta;
		residual_squared = residual.cwiseAbs2();
		variance = residual_squared.dot(omega)/Y.rows();
		sigma = omega/variance;
	}

    loglik = calculate_quick_loglik(sigma, Y.rows());

    if(delta == delta){
        result.beta = 0.0;
        result.pvalue = 1.0;
        result.SE = 0.0;
        result.loglik = loglik;
        result.chi = 0.0;
        result.h2r = h2;
        result.SD = sqrt(variance);

    }else{
        result.beta = 0.0;
        result.chi = 0.0;
        result.pvalue = 1.0;
        result.h2r = 0.0;
        result.loglik = 0.0;
        result.SD = 0.0;
        result.SE = 0.0;

    }
    return result;
}*/
double calculate_GWAS_loglik(Eigen::VectorXd SY2, Eigen::VectorXd lambda, Eigen::VectorXd Sigma){
    // double var_sum  = 0.0;
    // double log_sum = 0.0;
    // Eigen::VectorXd theta(2);
    // theta(0) = 1.0-h2r;
    // theta(1) = h2r;
    // Eigen::VectorXd Sigma  = (lambda*h2r).array() + e2;
    // const double variance_sum  = SY2.cwiseQuotient(variance);
    
    return -0.5*(log(Sigma.array()).sum() + SY2.rows());
    // for(int i = 0 ; i < SY.rows(); i++)
    
}

gwas_data  compute_null_model_MLE(Eigen::VectorXd Y, Eigen::VectorXd mean_column, Eigen::MatrixXd U, const unsigned precision){
    gwas_data result;
    double SD = 0.0;
    double loglik = 0.0;
    Eigen::VectorXd theta(2);
    theta(0) = 1.0;
    theta(1) = 0.0;
    Eigen::VectorXd Sigma = U*theta;
    Eigen::VectorXd Omega = Sigma.cwiseInverse();
    Eigen::VectorXd Y_mean_column = Y.cwiseProduct(mean_column);
    Eigen::VectorXd mean_column_squared = mean_column.cwiseAbs2();
    double mean = Omega.dot(Y_mean_column)/Omega.dot(mean_column_squared);
    
    Eigen::VectorXd residual = Y  - mean_column*mean;
    Eigen::VectorXd residual2 = residual.cwiseAbs2();
    double variance = residual2.dot(Omega)/residual2.rows();
    theta = theta*variance;
    Sigma = Sigma*variance;
    Eigen::VectorXd eigenvalues = U.col(1);
    double max_loglik = calculate_GWAS_loglik(residual2, eigenvalues, Sigma);
    double h2r = 0.0;
    for(double decimal = 0.1; decimal <= 1.0; decimal += 0.1){
        Eigen::VectorXd test_theta(2);
        test_theta(0) = 1.0 - decimal;
        test_theta(1) = decimal;
        Eigen::VectorXd test_sigma = U*test_theta;
        Eigen::VectorXd test_omega = test_sigma.cwiseInverse();
        double test_mean = test_omega.dot(Y_mean_column)/test_omega.dot(mean_column_squared);
        
        Eigen::VectorXd test_residual = Y - mean_column*mean;
        Eigen::VectorXd test_residual2 = test_residual.cwiseAbs2();
        double test_variance = test_residual2.dot(test_omega)/test_residual2.rows();
        test_sigma = test_sigma*test_variance;
        test_theta = test_theta*test_variance;
        double test_loglik = calculate_GWAS_loglik(test_residual2, eigenvalues, test_sigma);
        if(test_loglik > max_loglik && test_loglik == test_loglik){
            max_loglik = test_loglik;
            theta = test_theta;
            Sigma = test_sigma;
            variance = test_variance;
            residual = test_residual;
            h2r = decimal;
            mean = test_mean;
        }
    }
    
    for(int decimal_place = 2; decimal_place <= precision; decimal_place++){
        double scale = pow(10, -decimal_place);
        Eigen::VectorXd next_theta = theta;
        Eigen::VectorXd next_sigma = Sigma;
        Eigen::VectorXd next_residual = residual;
        double next_mean = mean;
        double next_variance = variance;
        double next_h2r= h2r;
        double next_loglik = max_loglik;
        for(int i = -6; i <= 6; i++){
            double value = i*scale;
            double test_h2r = h2r + value;
            if(test_h2r > 1.0 || test_h2r < 0.0 || test_h2r == h2r){
                continue;
            }
            Eigen::VectorXd test_theta(2);
            test_theta(0) = 1.0 - test_h2r;
            test_theta(1) = test_h2r;
            Eigen::VectorXd test_sigma = U*test_theta;
            Eigen::VectorXd test_omega = test_sigma.cwiseInverse();
            double test_mean = test_omega.dot(Y_mean_column)/test_omega.dot(mean_column_squared);
            
            Eigen::VectorXd test_residual = Y - mean_column*test_mean;
            Eigen::VectorXd test_residual2 = test_residual.cwiseAbs2();
            double test_variance = test_residual2.dot(test_omega)/test_residual2.rows();
            test_sigma = test_sigma*test_variance;
            test_theta = test_theta*test_variance;
            double test_loglik = calculate_GWAS_loglik(test_residual2, eigenvalues, test_sigma);
            if(test_loglik > next_loglik && test_loglik == test_loglik){
                next_loglik = test_loglik;
                next_theta = test_theta;
                next_sigma = test_sigma;
                next_variance = test_variance;
                next_residual = test_residual;
                next_mean = test_mean;
                next_h2r = test_h2r;
            }
            
            
        }
        
        max_loglik = next_loglik;
        theta = next_theta;
        Sigma = next_sigma;
        variance = next_variance;
        residual = next_residual;
        mean = next_mean;
        
        h2r = next_h2r;
        
    }
    //    cout << "null h2r " << h2r << " null loglik  " << max_loglik << endl;
    SD = 0.0;
    if(max_loglik == max_loglik){
        loglik = max_loglik;
        SD = sqrt(variance);
        result.beta = 0.0;
        result.pvalue = 1.0;
        result.SE = 0.0;
        result.loglik = loglik;
        result.chi = 0.0;
        result.h2r = h2r;
        result.SD = sqrt(variance);

    }else{
        result.beta = 0.0;
        result.pvalue = 0.0;
        result.chi = 0.0;
        result.h2r = 0.0;
        result.SD = 0.0;
        result.SE = 0.0;
        result.loglik = 0.0;
    }
    return result;
    
}

//double calculate_GWAS_loglik(Eigen::VectorXd SY2, Eigen::VectorXd lambda, Eigen::VectorXd Sigma);
static void  MLE_GWAS(gwas_data * result, Eigen::VectorXd Y, Eigen::MatrixXd X, Eigen::MatrixXd U, gwas_data null_result, const unsigned precision){
    
    Eigen::VectorXd raw_A = X.col(0).cwiseAbs2();
    Eigen::VectorXd raw_B = X.col(1).cwiseProduct(X.col(0));
    Eigen::VectorXd raw_C = X.col(1).cwiseAbs2();
    Eigen::VectorXd raw_D = Y.cwiseProduct(X.col(0));
    Eigen::VectorXd raw_E = Y.cwiseProduct(X.col(1));

    Eigen::VectorXd theta(2);
    theta(0) = 1.0;
    theta(1) = 0.0;
    
    Eigen::VectorXd Sigma = U*theta;
    Eigen::VectorXd Omega = Sigma.cwiseInverse();
    double A  = raw_A.dot(Omega);
    double B = raw_B.dot(Omega);
    double C = raw_C.dot(Omega);
    double D = raw_D.dot(Omega);
    double E = raw_E.dot(Omega);
    Eigen::VectorXd beta(2);
    double beta_denom  = A*C - B*B;
    beta(0) = C*D - B*E;
    beta(1) = A*E - B*D;
    beta /= beta_denom;
    Eigen::VectorXd residual = Y - X*beta;
    Eigen::VectorXd residual2 = residual.cwiseAbs2();
    double variance = residual2.dot(Omega)/residual2.rows();
    theta = theta*variance;
    Sigma = Sigma*variance;
    Eigen::VectorXd eigenvalues = U.col(1);
    double max_loglik = calculate_GWAS_loglik(residual2, eigenvalues, Sigma);
    double h2r = 0.0;
    for(double decimal = 0.1; decimal <= 1.0; decimal += 0.1){
        Eigen::VectorXd test_theta(2);
        test_theta(0) = 1.0 - decimal;
        test_theta(1) = decimal;
        Eigen::VectorXd test_sigma = U*test_theta;
        Eigen::VectorXd test_omega = test_sigma.cwiseInverse();
        double test_A = raw_A.dot(test_omega);
        double test_B = raw_B.dot(test_omega);
        double test_C = raw_C.dot(test_omega);
        double test_D = raw_D.dot(test_omega);
        double test_E = raw_E.dot(test_omega);
        double test_beta_denom = test_A*test_C - test_B*test_B;
        Eigen::VectorXd test_beta(2);
        test_beta(0) = (test_C*test_D - test_B*test_E)/test_beta_denom;
        test_beta(1) = (test_A*test_E - test_B*test_D)/test_beta_denom;
        Eigen::VectorXd test_residual = Y - X*test_beta;
        Eigen::VectorXd test_residual2 = test_residual.cwiseAbs2();
        double test_variance = test_residual2.dot(test_omega)/test_residual2.rows();
        test_sigma = test_sigma*test_variance;
        test_theta = test_theta*test_variance;
        double test_loglik = calculate_GWAS_loglik(test_residual2, eigenvalues, test_sigma);
        if(test_loglik > max_loglik && test_loglik == test_loglik){
            max_loglik = test_loglik;
            theta = test_theta;
            Sigma = test_sigma;
            variance = test_variance;
            residual = test_residual;
            h2r = decimal;
            beta = test_beta;
        }
    }
    
    for(int decimal_place = 2; decimal_place <= precision; decimal_place++){
        double scale = pow(10, -decimal_place);
        Eigen::VectorXd next_theta = theta;
        Eigen::VectorXd next_sigma = Sigma;
        Eigen::VectorXd next_residual = residual;
        Eigen::VectorXd next_beta = beta;
        double next_variance = variance;
        double next_h2r= h2r;
        double next_loglik = max_loglik;
        for(int i = -5; i <= 6; i++){
            double value = i*scale;
            double test_h2r = h2r + value;
            if(test_h2r > 1.0 || test_h2r < 0.0 || test_h2r == h2r){
                continue;
            }
            Eigen::VectorXd test_theta(2);
            test_theta(0) = 1.0 - test_h2r;
            test_theta(1) = test_h2r;
            Eigen::VectorXd test_sigma = U*test_theta;
            Eigen::VectorXd test_omega = test_sigma.cwiseInverse();
            double test_A = raw_A.dot(test_omega);
            double test_B = raw_B.dot(test_omega);
            double test_C = raw_C.dot(test_omega);
            double test_D = raw_D.dot(test_omega);
            double test_E = raw_E.dot(test_omega);
            double test_beta_denom = test_A*test_C - test_B*test_B;
            Eigen::VectorXd test_beta(2);
            test_beta(0) = (test_C*test_D - test_B*test_E)/test_beta_denom;
            test_beta(1) = (test_A*test_E - test_B*test_D)/test_beta_denom;
            Eigen::VectorXd test_residual = Y - X*test_beta;
            Eigen::VectorXd test_residual2 = test_residual.cwiseAbs2();
            double test_variance = test_residual2.dot(test_omega)/test_residual2.rows();
            test_sigma = test_sigma*test_variance;
            test_theta = test_theta*test_variance;
            double test_loglik = calculate_GWAS_loglik(test_residual2, eigenvalues, test_sigma);
            if(test_loglik > next_loglik && test_loglik == test_loglik){
                next_loglik = test_loglik;
                next_theta = test_theta;
                next_sigma = test_sigma;
                next_variance = test_variance;
                next_residual = test_residual;
                next_h2r = test_h2r;
                next_beta = test_beta;
            }
            
            
        }
        
        max_loglik = next_loglik;
        theta = next_theta;
        Sigma = next_sigma;
        variance = next_variance;
        residual = next_residual;
        h2r = next_h2r;
        beta = next_beta;
        
    }
    double chi = 2.0*(max_loglik - null_result.loglik);
    Eigen::MatrixXd beta_covariance_matrix = X.transpose()*Sigma.cwiseInverse().asDiagonal()*X;
    double denom = beta_covariance_matrix(0,0)*beta_covariance_matrix(1, 1) - beta_covariance_matrix(1,0)*beta_covariance_matrix(0,1);
    
    //  gwas_data result;
    if(chi >= 0.0 && chi == chi && denom > 0.0 && denom == denom){
        
        //        Eigen::MatrixXd hessian = calculate_REML_hessian(residual, SNP, Sigma, eigenvalues, variance);
        //      Eigen::VectorXd standard_errors = hessian.inverse().diagonal().cwiseSqrt();
        
       // double beta_se = beta_covariance(1,1);
        result->SE = sqrt(beta_covariance_matrix(0, 0)/denom);
        result->beta = beta(1);
        result->chi = chi;
        result->SD = sqrt(variance);
//        result->SE = 0.0;//beta_se;//standard_errors(0);
        result->pvalue = chicdf(result->chi , 1);
        result->h2r = h2r;
        result->loglik = max_loglik;
    }else{
        result->beta = 0.0;
        result->chi = 0.0;
        result->pvalue = 1.0;
        result->h2r = null_result.h2r;
        result->loglik = null_result.loglik;
        result->SD = null_result.SD;
        result->SE = 0.0;

    }
    
    
    // return result;
    
    
    
}

static void calculate_pvalue(vector<double>::iterator & pvalue_iterator, Eigen::MatrixXd syP_Sigma_P, Eigen::MatrixXd sigmaP, Eigen::MatrixXd snps_matrix){
	
	Eigen::ArrayXXd Ts_numerator = (snps_matrix*syP_Sigma_P).array();	
	
	Eigen::ArrayXXd Ts_denominator = ((snps_matrix.array()*snps_matrix.array()).matrix()*sigmaP).array();
	
    
	Eigen::ArrayXXd Ts = Ts_numerator*Ts_numerator/Ts_denominator;
	//double compare_value;
	/*double value;
	for(int row = 0 ; row < n_snps;row++){
	pvalue = 0.0;
	compare_value = Ts(0, 0);
	pvalue = (Ts >= compare_value).count()/(double)sigmaP.cols();	
	*/
//#pragma omp parallel for
	for(int row = 0; row < Ts.rows(); row++){
		*pvalue_iterator++ = (Ts.row(row) >= Ts(row, 0)).count()/double(Ts.cols());
	}
}	
static void print_gwas_help(Tcl_Interp * interp){
	Solar_Eval(interp, "help gwas");
}
static void calculate_eigen_data_null(Eigen::VectorXd & Y,  Eigen::MatrixXd & eigenvectors_transposed,  Eigen::MatrixXd & U,  Eigen::MatrixXd  phi2){
    const int n_subjects = Y.rows();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(phi2);
    
    Eigen::MatrixXd temp_eigenvectors = es.eigenvectors();
    Eigen::VectorXd temp_eigenvalues = es.eigenvalues();
    const int n_subjects_reml = n_subjects - (temp_eigenvalues.array() == 0.0).count();
    gwas_data result;

        eigenvectors_transposed = temp_eigenvectors.transpose();
        U.resize(n_subjects, 2);
        U.col(0) = Eigen::ArrayXd::Ones(n_subjects).matrix();
        U.col(1) = temp_eigenvalues;
        Y = eigenvectors_transposed * Y;//(Y.array() - Y.mean()).matrix();

    
}


static vector<gwas_data> GWAS_MLE_fix_missing_run(gwas_data null_result, Eigen::VectorXd default_Y, Eigen::VectorXd default_mean,\
						Eigen::MatrixXd default_eigenvectors_transposed,Eigen::MatrixXd default_U, \
                                              int * SNP_data, const int n_snps, const int n_subjects, const unsigned precision){
    // const int n_subjects = Y.rows();

    vector<gwas_data> results(n_snps);
#pragma omp parallel for schedule(dynamic)
    for(int iteration = 0; iteration < n_snps; iteration++){
        Eigen::MatrixXd local_X(n_subjects , 2);
        local_X.col(0) = default_mean;
        int local_n_subjects = n_subjects;
        double SNP_mean = 0.0;
        for(int row = 0 ; row < n_subjects; row++){
            local_X(row, 1) = SNP_data[iteration*n_subjects + row];
            if(local_X(row, 1) != 3.0){
                SNP_mean += local_X(row, 1);
            }else{
                local_n_subjects--;
            }
        }
        SNP_mean /= local_n_subjects;
        for(int row = 0; row < n_subjects; row++){
            if(local_X(row, 1) != 3.0){
                local_X(row, 1) -= SNP_mean;
            }else{
                local_X(row, 1) = 0.0;
            }
        }
        local_X.col(1) = default_eigenvectors_transposed*local_X.col(1);
        gwas_data result;

        MLE_GWAS(&result,default_Y, local_X, default_U, null_result, precision);
       // if(result.chi == 0.0) result.SD = default_SD;
        results[iteration] = result;
        
    }
    
    return results;
}


static vector<gwas_data> GWAS_MLE_run(gwas_data default_null_result,Eigen::VectorXd Y, Eigen::VectorXd default_Y, Eigen::VectorXd default_mean,\
				Eigen::MatrixXd default_U, Eigen::MatrixXd default_eigenvectors_transposed, Eigen::MatrixXd default_phi2,\
				int * SNP_data, const int n_snps, const int n_subjects, const unsigned precision){
    //const int n_subjects = Y.rows();

    vector<gwas_data> results(n_snps);
#pragma omp parallel for schedule(dynamic)
    for(int iteration = 0; iteration < n_snps; iteration++){
        Eigen::VectorXd local_Y;
        Eigen::MatrixXd local_U;
        Eigen::MatrixXd local_X;
        Eigen::MatrixXd local_eigenvectors_transposed;
        gwas_data local_null_result = default_null_result;
        double local_SD;
        Eigen::MatrixXd local_phi2;
        Eigen::VectorXd local_SNP(n_subjects);
        int local_n_subjects = n_subjects;
        for(int row = 0 ; row < n_subjects; row++){
            local_SNP(row) = SNP_data[iteration*n_subjects + row];
            if(local_SNP(row) == 3.0) local_n_subjects--;
        }
        if(local_n_subjects == 0){
            gwas_data result;
            result.beta = 0.0;
            result.SD = 0.0;
            result.SE = 0.0;
            result.chi = 0.0;
            result.pvalue = 0.0;
            results[iteration] = result;
        }else{
        if(local_n_subjects == n_subjects){
            local_U = default_U;
            local_eigenvectors_transposed = default_eigenvectors_transposed;
            local_Y = default_Y;
            local_SNP = local_eigenvectors_transposed*local_SNP;
            Eigen::MatrixXd temp_X(n_subjects, 2);
            temp_X.col(0) = default_mean;
            temp_X.col(1) = local_SNP;
            local_X = temp_X;
            local_SD = default_null_result.SD;
	    local_null_result = default_null_result;
        }else{
            Eigen::VectorXd temp_local_Y(local_n_subjects);
            Eigen::MatrixXd local_phi2(local_n_subjects,local_n_subjects);
            local_X = Eigen::ArrayXXd::Ones(local_n_subjects, 2).matrix();
            
            int local_row = 0;
            for(int row = 0; row < n_subjects; row++){
                if(local_SNP(row) != 3.0){
                    temp_local_Y(local_row) = Y(row);
                    local_X(local_row, 1) = local_SNP(row);
                    int local_col = local_row;
                    for(int col = row; col < n_subjects; col++){
                        if( local_SNP(col) != 3.0){
                            local_phi2(local_row, local_col) = default_phi2(row,col);
                            local_phi2(local_col++, local_row) = default_phi2(col, row);
                        }
                    }
                    local_row++;
                }
            }
            local_Y = temp_local_Y;
            calculate_eigen_data_null(local_Y,  local_eigenvectors_transposed, local_U,   local_phi2);
            
            local_X = local_eigenvectors_transposed*local_X;
            Eigen::VectorXd local_mean = local_X.col(0);
            local_null_result = compute_null_model_MLE(local_Y, local_mean, local_U, precision);
            
            
            
        }
        gwas_data result;
        if(local_null_result.loglik == local_null_result.loglik){
            MLE_GWAS(&result, local_Y, local_X, local_U, local_null_result, precision);
        }else{
            result.chi  = 0.0;
            result.beta = 0.0;
            result.loglik = 0.0;
            result.SD = 0.0;
            result.SE = 0.0;
            result.h2r = 0.0;
            result.pvalue = 0.0;
        }
        results[iteration] = result;
        }
    }
    
    return results;
    
}
static Eigen::VectorXd calculate_theta(Eigen::VectorXd Y, Eigen::MatrixXd aux){
    Eigen::VectorXd theta(2);
    
    theta(0) = 0.0;
    theta(1) = 0.0;
    
    Eigen::VectorXd F = Y.cwiseAbs2();
    
    const double F_mean = F.mean();
    
    const double score = aux.col(1).dot(((F/F_mean).array() - 1.0).matrix())/F_mean;
    
    if(score <= 0.0) return theta;
    
    theta = aux.colPivHouseholderQr().solve(F);
    
    if(theta(0) < 0.0) theta(0) = 0.0;
    
    if(theta(1) < 0.0) theta(1) = 0.0;
    
    if(theta(0) == 0.0 && theta(1) == 0.0) return theta;
    
    Eigen::VectorXd Sigma = aux*theta;
    
    Eigen::MatrixXd Omega = Sigma.cwiseAbs2().cwiseInverse().asDiagonal();
    
    theta = (Omega*aux).colPivHouseholderQr().solve(Omega*F);
    
    if(theta(0) < 0.0) theta(0) = 0.0;
    
    if(theta(1) < 0.0) theta(1) = 0.0;
    
  
    
    if(theta(0) != theta(0) || theta(1) != theta(1)){
        theta(0) = 0.0;
        theta(1) = 0.0;
    }
    
    
    return theta;
    
    
}
static void calculate_pvalues_permutation_method(vector<double>::iterator & pvalue_iterator,Eigen::MatrixXd Sigma_Y_permutated, Eigen::MatrixXd Sigma_permutated, Eigen::MatrixXd SNP_data){
    Eigen::ArrayXXd Ts_numerator = (SNP_data*Sigma_Y_permutated).cwiseAbs2().array();
    
    Eigen::ArrayXXd Ts_denominator = (SNP_data.cwiseAbs2()*Sigma_permutated).array();
    Eigen::ArrayXXd Ts = Ts_numerator/Ts_denominator;
   
    for(int row = 0; row < Ts.rows(); row++){
        *pvalue_iterator++ = (Ts.row(row) >= Ts(row, 0)).count()/double(Ts.cols());
    }
    
   
    
}


static void compute_pvalue_batch(vector<double>::iterator & pvalue_iterator, Eigen::MatrixXd Sigma_Y_permutated, Eigen::MatrixXd Sigma_permutated, Eigen::MatrixXd eigenvectors_transposed,\
                                 int * SNP_data, const int batch_size,const int n_subjects, const int n_permutations){
    Eigen::MatrixXd SNP_matrix(n_subjects, batch_size);
    
    for(int col = 0 ; col < batch_size; col++){
        double SNP_mean = 0.0;
        int subject_count = 0;
        for(int subject = 0; subject < n_subjects; subject++){
            int value = SNP_data[col*n_subjects + subject];
            if(value != 3){
                SNP_mean+= value;
                subject_count++;
            }
            SNP_matrix(subject, col) = value;
        }
        SNP_mean /= subject_count;
        for(int subject = 0; subject < n_subjects; subject++){
            if(SNP_matrix(subject,col) != 3.0){
                SNP_matrix(subject,col) -= SNP_mean;
            }else{
                SNP_matrix(subject, col) = 0.0;
            }
        }
    }
    SNP_matrix = eigenvectors_transposed*SNP_matrix;
    SNP_matrix.transposeInPlace();
    calculate_pvalues_permutation_method(pvalue_iterator, Sigma_Y_permutated, Sigma_permutated, SNP_matrix);
    
    
}

static void compute_permutation_data(Eigen::VectorXd Sigma_Y, Eigen::VectorXd Sigma, Eigen::MatrixXd & Sigma_Y_permutated, \
                                     Eigen::MatrixXd & Sigma_permutated, const int n_subjects, const int n_permutations){
    vector<int> indices(n_subjects);
    iota(indices.begin(), indices.end(), 0);
    
    Sigma_Y_permutated.col(0) = Sigma_Y;
    Sigma_permutated.col(0) = Sigma;
    
    
    for(int col = 1 ; col < n_permutations + 1; col++){
        
        random_shuffle(indices.begin(), indices.end());
        
        for(int row = 0; row < n_subjects; row++){
            Sigma_Y_permutated(row, col) = Sigma_Y(indices[row]);
            Sigma_permutated(row, col) = Sigma(indices[row]);
        }
        
    }
    
}

static vector<double> compute_pvalues_permutation(Eigen::VectorXd Y, Eigen::VectorXd Sigma, Eigen::MatrixXd eigenvectors_transposed, int * snp_data, \
                                                    const int n_snps, const int n_subjects, const  int n_permutations){
    
    Eigen::VectorXd Sigma_Y = Y.cwiseProduct(Sigma);
    
    Eigen::MatrixXd Sigma_Y_permutated(n_subjects, n_permutations + 1);
    Eigen::MatrixXd Sigma_permutated(n_subjects, n_permutations + 1);
    
    compute_permutation_data(Sigma_Y, Sigma, Sigma_Y_permutated, \
                             Sigma_permutated,  n_subjects,  n_permutations);
    vector<double> pvalues(n_snps);
    if(n_snps <= PERMUTATION_BATCH_SIZE){
        vector<double>::iterator pvalue_iterator = pvalues.begin();
        compute_pvalue_batch(pvalue_iterator, Sigma_Y_permutated, Sigma_permutated, eigenvectors_transposed, \
                             snp_data, n_snps, n_subjects, n_permutations);
    }else{
        const int batch_size = PERMUTATION_BATCH_SIZE;
        const int n_batches = ceil(double(n_snps)/batch_size);
        int remainder_snps = n_snps % batch_size;
        if(remainder_snps == 0)
            remainder_snps = batch_size;

#pragma omp parallel for
        for(int batch = 0; batch < n_batches; batch++){
            int  thread_batch_size = batch_size;
            if(batch == n_batches - 1)
                thread_batch_size = remainder_snps;
             vector<double>::iterator pvalue_iterator = pvalues.begin() + batch_size*batch;
            compute_pvalue_batch(pvalue_iterator, Sigma_Y_permutated, Sigma_permutated, eigenvectors_transposed, \
                                 snp_data + n_subjects*batch*batch_size, thread_batch_size, n_subjects, n_permutations);
        }
    }
    

    
    return pvalues;
    
}


vector<string> read_trait_list(const char * list_filename){
	ifstream input_stream(list_filename);
	vector<string> output;
	if(!input_stream.is_open()) return output;
	string trait;
	while(input_stream >> trait){
		output.push_back(trait);
	}

	return output;

}
static const char * run_gwas_list(const char * phenotype_filename, const char * list_filename, const char * plink_filename, const bool fix_missing, const unsigned precision){
	vector<string> trait_list;// = read_trait_list(list_filename);
	if(list_filename) {
		trait_list = read_trait_list(list_filename);
	}else{
		string trait_name = string(Trait::Name(0));
		trait_list.push_back(trait_name);
	}

	
	if(trait_list.size() == 0){
		return "No traits read from list file";
	}
	
	pio_file_t * plink_file = new pio_file_t;
	if (pio_open(plink_file, plink_filename) != PIO_OK){
		return "Error opening plink file";
	}
	const unsigned num_plink_samples = pio_num_samples(plink_file);
	const unsigned n_snps = pio_num_loci(plink_file);
	unsigned batch_size = GWAS_BATCH_SIZE;
	if(batch_size >= n_snps){
		batch_size = n_snps;
	}
	unsigned iterations = ceil(n_snps/batch_size);
	unsigned final_batch_size = n_snps % batch_size;
	if(final_batch_size == 0) final_batch_size = batch_size;
	if(n_snps >= GWAS_BATCH_SIZE){
		batch_size = GWAS_BATCH_SIZE;
		iterations = ceil(n_snps/GWAS_BATCH_SIZE);
		final_batch_size = n_snps % batch_size;
		if(final_batch_size == 0 ) final_batch_size = batch_size;
	}
	pio_sample_t * sample;
	vector<string> plink_ids;
	for(unsigned i = 0; i < num_plink_samples; i++){
		sample = pio_get_sample(plink_file, i);
		plink_ids.push_back(string(sample->iid));
	}
	vector<string> snp_names;
	for(unsigned snp = 0; snp < n_snps; snp++){
		pio_locus_t * locus = pio_get_locus(plink_file, snp);
		snp_names.push_back(string(locus->name));
	}
	snp_t * snp_buffer = new snp_t[num_plink_samples];
	Solar_Trait_Reader * trait_reader = new Solar_Trait_Reader(phenotype_filename, trait_list, plink_ids);
	for(unsigned set = 0; set < trait_reader->get_n_sets(); set++){
		Eigen_Data * eigen_data = trait_reader->get_eigen_data_set(set);
		vector<string> ids = eigen_data->get_ids();
		unsigned plink_index_map[ids.size()];
		for(unsigned index = 0; index < ids.size(); index++){
			string id = ids[index];
			plink_index_map[index] = distance(plink_ids.begin(), find(plink_ids.begin(), plink_ids.end(), id));
		}
		Eigen::VectorXd eigenvalues = Eigen::Map<Eigen::VectorXd>(eigen_data->get_eigenvalues(), ids.size());
		Eigen::MatrixXd eigenvectors_transposed = Eigen::Map<Eigen::MatrixXd>(eigen_data->get_eigenvectors_transposed(), ids.size(), ids.size());
		Eigen::MatrixXd phi2 = eigenvectors_transposed.transpose()*eigenvalues.asDiagonal()*eigenvectors_transposed;
		int * snp_data = new int[batch_size*ids.size()];
		for(unsigned trait = 0; trait < eigen_data->get_n_phenotypes(); trait++){
			string trait_name = eigen_data->get_trait_name(trait);
			const char * output_filename = string(trait_name + "-gwas.out").c_str();
			ofstream output_stream(output_filename);
			output_stream << "SNP,h2r,loglik,SD,beta_snp,beta_snp_se,chi2,p-value\n";

			Eigen::VectorXd trait_vector = Eigen::Map<Eigen::VectorXd>(eigen_data->get_phenotype_column(trait), ids.size());
			vector<gwas_data> results;
			pio_reset_row(plink_file);
			unsigned current_batch_size = batch_size;
			
    		
			Eigen::VectorXd default_Y = trait_vector;
    			Eigen::VectorXd mean = Eigen::ArrayXd::Ones(default_Y.rows()).matrix();
    			Eigen::MatrixXd default_U;
   		        Eigen::MatrixXd default_eigenvectors_transposed;
   			Eigen::MatrixXd default_phi2 = phi2;
      			calculate_eigen_data_null(default_Y, default_eigenvectors_transposed, default_U, default_phi2);
   			Eigen::VectorXd default_mean = default_eigenvectors_transposed*mean;
    		
    			gwas_data default_null_result = compute_null_model_MLE(default_Y, default_mean, default_U, precision);
			for(unsigned iteration = 0; iteration < iterations; iteration++){
				if(iteration == iterations - 1) current_batch_size = final_batch_size;
				for(unsigned snp = 0; snp < current_batch_size; snp++){
					pio_next_row(plink_file, snp_buffer);
					for(unsigned id = 0; id < ids.size(); id++){
						snp_data[ids.size()*snp + id] = snp_buffer[plink_index_map[id]];
					}
				}

				if(fix_missing){
					results = GWAS_MLE_fix_missing_run(default_null_result, default_Y, default_mean,\
						default_eigenvectors_transposed, default_U, \
                                              snp_data, current_batch_size, ids.size(),precision);
				}else{
					results = GWAS_MLE_run(default_null_result,trait_vector,  default_Y, default_mean,\
				 default_U,  default_eigenvectors_transposed, default_phi2,\
				snp_data, current_batch_size, ids.size(), precision);
				}

				for(unsigned snp = 0; snp < current_batch_size; snp++){
					output_stream << snp_names[iteration*batch_size + snp] << "," << results[snp].h2r << "," << \
					results[snp].loglik << "," << results[snp].SD << "," << results[snp].beta \
					 << "," << results[snp].SE << "," << results[snp].chi << "," << results[snp].pvalue << "\n";

			
				}		
					

			}
			output_stream.close();
		}
		delete [] snp_data;

	}
	delete trait_reader;

	return 0;
}

extern "C" int gwaCmd(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[]){
//	int n_permutations = 0;
	const char *  plink_filename = 0;
	
    bool correct_missing = false;
	const char * list_filename = 0;

    unsigned precision = 6;
    for(unsigned arg = 1; arg < argc; arg++){
	if(!StringCmp(argv[arg], "-help", case_ins) || !StringCmp(argv[arg], "--help", case_ins) || !StringCmp(argv[arg], "help", case_ins)){
		print_gwas_help(interp);
		return TCL_OK;
	}else if((!StringCmp(argv[arg], "-plink", case_ins) || !StringCmp(argv[arg], "--plink", case_ins)) && arg + 1 < argc){
		plink_filename = argv[++arg];
	}else if (!StringCmp(argv[arg], "-fix", case_ins) || !StringCmp(argv[arg], "--fix", case_ins) || !StringCmp(argv[arg], "-f", case_ins)){
		correct_missing = true;
	}else if((!StringCmp(argv[arg], "-precision", case_ins) || !StringCmp(argv[arg], "--precision", case_ins)) && arg + 1 < argc){
		precision = atoi(argv[++arg]);
	}else if((!StringCmp(argv[arg], "-list", case_ins) || !StringCmp(argv[arg], "--list", case_ins)) && arg + 1 < argc){
		list_filename = argv[++arg];
	}else{
		RESULT_LIT("Invalid argument entered");
		return TCL_ERROR;
	}
   }


    if((precision < 1 || precision > 6)){
        RESULT_LIT("Precision must be between 1 and 6");
        return TCL_ERROR;
    }
	Eigen::VectorXd trait_vector;

	int success;
	string phenotype_filename = Phenotypes::filenames(); 
	if(phenotype_filename.length() == 0){
		RESULT_LIT("No phenotype file is currently loaded");
		return TCL_ERROR;
	}
	if(!plink_filename){
		RESULT_LIT("No plink file was specified");
		return TCL_ERROR;
	}				
	
	if(!list_filename && Trait::Number_Of() == 0){
		RESULT_LIT("No trait selected with trait command or specified with list file");
		return TCL_ERROR;
	}
	try{
		load_phi2_matrix(interp);
	}catch(...){
		RESULT_LIT("phi2 matrix could not be loaded.  Check to see if pedigree has been properly loaded.");
		return TCL_ERROR;
	}
	const char * error = 0;
	error = run_gwas_list(phenotype_filename.c_str(), list_filename,\
				 plink_filename, correct_missing, precision);	
	if(!error){
		RESULT_LIT(error);
		return TCL_ERROR;
	}
										 
											
    return TCL_OK;
    
}
