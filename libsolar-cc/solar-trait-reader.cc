//
//  Solar_Trait_Reader.cc
//  
//
//  Created by Brian Donohue on 8/19/18.
//

#include "solar-trait-reader.h"
#include "solar.h"
#include <vector>
#include <string>
#include <iterator>
#include <algorithm>
#include <omp.h>
#include <cmath>

using namespace std;

extern "C" void symeig_ (int*, double*, double*, double*, double*, int*);
static Matrix * static_phi2 = 0;
void load_phi2_matrix(Tcl_Interp * interp){
    static_phi2 = Matrix::find("phi2");
    if (!static_phi2) {
        Solar_Eval(interp, "matrix load phi2.gz phi2");
        static_phi2 = Matrix::find("phi2");
        if(!static_phi2){
            throw Solar_Trait_Reader_Exception(std::string("Phi2 matrix could not be loaded"));
        }
    }
    
}
Eigen_Data::~Eigen_Data(){
    ibdids.clear();
    ids.clear();
    trait_names.clear();
    index_map.clear();
    if(eigenvalues){
        delete [] eigenvalues;
        eigenvalues = 0;
    }
    if(eigenvectors_transposed){
        delete [] eigenvectors_transposed;
        eigenvectors_transposed = 0;
    }
    if(phenotype_buffer){
        delete [] phenotype_buffer;
        phenotype_buffer = 0;
    }
}
Eigen_Data::Eigen_Data(vector<string> all_ids, vector<string> skip_ids, const size_t _n_sets){
    n_phenotypes = _n_sets;
    const char * errmsg = 0;
    SolarFile * pedindex_in = SolarFile::open("Eigen_Data", "pedindex.out", &errmsg);
    if(errmsg){
        throw Solar_Trait_Reader_Exception(string(errmsg));
    }
    pedindex_in->start_setup(&errmsg);
    if(errmsg){
        throw Solar_Trait_Reader_Exception(string(errmsg));
    }
    pedindex_in->setup("id", &errmsg);
    char ** file_data;
    size_t ibdid = 1;
    
    while (0 != (file_data = pedindex_in->get (&errmsg))){
        string id = string(file_data[0]);
        vector<string>::iterator  find_iter = find(skip_ids.begin(), skip_ids.end(), id);
        if(find_iter == skip_ids.end()){
            find_iter = find(all_ids.begin(), all_ids.end(), id);
            if(find_iter != all_ids.end()){
                ids.push_back(id);
                ibdids.push_back(ibdid);
            }
        }
        ibdid++;
    }
    index_map.resize(all_ids.size());
    for(size_t row = 0; row < all_ids.size(); row++){
        string id = all_ids[row];
        vector<string>::iterator find_iter = std::find(ids.begin(), ids.end(), id);
        if(find_iter != ids.end()){
            index_map[row] = distance(ids.begin(), find_iter);
        }else{
            index_map[row] = -1;
        }
    }
    if(static_phi2 == 0){
        throw Solar_Trait_Reader_Exception(string("Load phi2 matrix was not called prior to running Eigen_Data constructor"));
    }
    n_subjects = ids.size();
    phenotype_buffer = new double[n_subjects*n_phenotypes];
    trait_names.resize(n_phenotypes);
    double * eigenvectors = new double[n_subjects*n_subjects];
    eigenvectors_transposed = new double[n_subjects*n_subjects];
    eigenvalues = new double[n_subjects];
    double * phi2 = new double[n_subjects*n_subjects];
    double phi2_value;
    for(int col = 0; col < n_subjects; col++){
	try{
        phi2_value = static_phi2->get(ibdids[col], ibdids[col]);
	}catch(...){
	phi2_value = 0;
	}
        phi2[col*n_subjects + col] = phi2_value;
        for(int row = col+1; row < n_subjects; row++){
	    try{
            phi2_value = static_phi2->get(ibdids[row], ibdids[col]);
	    }catch(...){
		phi2_value = 0;
	    }
            phi2[col*n_subjects + row] = phi2_value;
            phi2[row*n_subjects + col] = phi2_value;
        }
    }
    

    calculate_eigenvectors_and_eigenvalues(phi2,eigenvectors, n_subjects);
    for(size_t col = 0; col < n_subjects; col++){
        for(size_t row = 0; row < n_subjects ; row++){
            eigenvectors_transposed[col*n_subjects + row] = eigenvectors[row*n_subjects + col];
        }
    }
    delete [] eigenvectors;
    delete [] phi2;

}

void Eigen_Data::calculate_eigenvectors_and_eigenvalues (double * phi2, double * eigenvectors, int n)
{
    double* e =  new double[n];
    memset(e, 0, sizeof(double)*n);
    int * info = new int;
    *info  = 0;
    symeig_(&n, phi2, eigenvalues, e, eigenvectors, info);
    delete [] e;
    delete [] info;
}


void Eigen_Data::set_phenotype_column(const size_t column_index, string name,  double *  input_buffer){
    
    for(size_t row = 0 ; row < index_map.size(); row++){
        int index = index_map[row];
        if(index != -1){
            phenotype_buffer[column_index*n_subjects + index] = input_buffer[row];
        }
    }
    trait_names[column_index] = name;
}

void Solar_Trait_Reader::Sort_ID_Sets(const double * const phenotype_buffer, vector<string> all_ids, vector<string> skip_ids, \
                                    vector<size_t> & n_sets_vector, vector<size_t> & eigen_data_indices, vector< vector <string> > & skip_id_sets,\
                                    const size_t start, const size_t end){

    for(size_t data_index = start; data_index < end; data_index++){
        vector<string> local_skip_ids;
        for(size_t row_index = 0; row_index < all_ids.size() ; row_index++){
            const double value = phenotype_buffer[data_index*all_ids.size() + row_index];
            if(value != value){
                string current_id = all_ids[row_index];
                local_skip_ids.push_back(current_id);
            }
        }
        vector<string> all_skip_ids = skip_ids;
        for(vector<string>::iterator id_iter = local_skip_ids.begin(); id_iter != local_skip_ids.end(); id_iter++){
            vector<string>::iterator find_iter = find(skip_ids.begin(), skip_ids.end(), *id_iter);
            if(find_iter == skip_ids.end()){
                all_skip_ids.push_back(*id_iter);
            }
        }
#pragma omp critical
        {
            if(skip_id_sets.size() == 0){
                skip_id_sets.push_back(all_skip_ids);
                n_sets_vector.push_back(1);
                eigen_data_indices[data_index] = 0;
            }else{
                int eigen_index = -1;
                for(size_t index  = 0 ; index < skip_id_sets.size(); index++){
                    if(skip_id_sets[index].size() != all_skip_ids.size()) continue;
		    bool is_equal = true;
                    vector<string> current_id_set = skip_id_sets[index];
		    
                    for(size_t  id_index  = 0; id_index <  current_id_set.size();id_index++){
			if(current_id_set[id_index]  != all_skip_ids[id_index]){
			    is_equal = false;
			    break;
                         }
                    }
		    if(is_equal){
			eigen_index = index;
			break;
	            }
                }
                if(eigen_index != -1){
                    eigen_data_indices[data_index] = eigen_index;
                    n_sets_vector[eigen_index]++;
                }else{
                    skip_id_sets.push_back(all_skip_ids);
                    n_sets_vector.push_back(1);
                    eigen_data_indices[data_index] =  skip_id_sets.size() - 1;
                }
            }
        }
    }
}
void  Solar_Trait_Reader::Read_Phenotype_Data_Thread_Launch(SolarFile * const file, double * const phenotype_buffer, vector<string> field_names,const size_t n_rows, const size_t start, const size_t end){
    const char * errmsg = 0;
    file->start_setup(&errmsg);
    if(errmsg){
        throw Solar_Trait_Reader_Exception(string(errmsg));
    }
    for(size_t trait_index = start; trait_index < end; trait_index++){
        file->setup(field_names[trait_index].c_str(), &errmsg);
        if(errmsg){
            throw Solar_Trait_Reader_Exception(string(errmsg));
        }
    }
    char ** file_data;
    const size_t n_fields = end - start;
    size_t row_index = 0;
    while (0 != (file_data = file->get (&errmsg))){
        
        for(size_t index = 0; index < n_fields; index++){
            if(StringCmp(file_data[index], 0, case_ins)){
                phenotype_buffer[row_index + index*n_rows] = atof (file_data[index]);
            }else{
                phenotype_buffer[row_index + index*n_rows] = nan("");
            }
        }
        
        row_index++;
        
    }
    
}
Solar_Trait_Reader::~Solar_Trait_Reader(){
    if(eigen_data){
        for(size_t index = 0; index < n_sets; index++){
            delete eigen_data[index];
        }
        delete [] eigen_data;
        eigen_data = 0;
    }
}
Solar_Trait_Reader::Solar_Trait_Reader(const char * phenotype_filename, vector<string> trait_names, vector<string> external_id_set){
    vector<string> all_ids;
    n_phenotypes = trait_names.size();
    const char * initial_errmsg = 0;
    SolarFile * initial_reader =  SolarFile::open("Multi-Trait-Read", phenotype_filename, &initial_errmsg);
    if(initial_errmsg){
        throw Solar_Trait_Reader_Exception(string(initial_errmsg));
    }
    
    initial_reader->start_setup(&initial_errmsg);
    if(initial_errmsg){
        throw Solar_Trait_Reader_Exception(string(initial_errmsg));
    }
    initial_reader->setup("id", &initial_errmsg);
    if(initial_errmsg){
        throw Solar_Trait_Reader_Exception(string(initial_errmsg));
    }
    char ** file_data;
    while (0 != (file_data = initial_reader->get (&initial_errmsg))){
        
        all_ids.push_back(string(file_data[0]));
        
        
    }
    delete initial_reader;
    double * phenotype_buffer = new double[all_ids.size()*trait_names.size()];
    
    vector<string> skip_ids;
    if(external_id_set.size() != 0){
        for(vector<string>::iterator id_iter = all_ids.begin(); id_iter != all_ids.end(); id_iter++){
            vector<string>::iterator find_iter = find(external_id_set.begin(), external_id_set.end(), *id_iter);
            if(find_iter != external_id_set.end()){
                external_id_set.erase(find_iter);
            }else{
                skip_ids.push_back(*id_iter);
            }
        }
    }
    
    if(n_phenotypes < 100){
        const char * errmsg = 0;
        SolarFile * file = SolarFile::open("Multi-Trait-Read", phenotype_filename, &errmsg);
        if(errmsg){
            throw Solar_Trait_Reader_Exception(string(errmsg));
        }
        Read_Phenotype_Data_Thread_Launch(file, phenotype_buffer, trait_names, all_ids.size(), 0, trait_names.size());
        delete file;
    }else{
#pragma omp parallel
        {
            int n_threads = omp_get_num_threads();
            int thread_index = omp_get_thread_num();
            const char * errmsg = 0;
            SolarFile * file = SolarFile::open("Multi-Trait-Read", phenotype_filename, &errmsg);
            if(errmsg){
                throw Solar_Trait_Reader_Exception(string(errmsg));
            }
            int batch_size = ceil(double(n_phenotypes)/n_threads);
            int end = (1 +thread_index)*batch_size;
            if(thread_index + 1 == n_threads)
                end = n_phenotypes;
            
            
            Read_Phenotype_Data_Thread_Launch(file, phenotype_buffer + thread_index*batch_size*all_ids.size(), trait_names, all_ids.size(), thread_index*batch_size, end);
            
            delete file;
        }
    }
    

    vector<size_t> eigen_data_indices(n_phenotypes);
    vector< vector<string> > skip_id_sets;
    vector<size_t> n_sets_vector;
    if(n_phenotypes < 100){
        Sort_ID_Sets(phenotype_buffer, all_ids, skip_ids, \
                    n_sets_vector, eigen_data_indices, skip_id_sets,\
                     0, n_phenotypes);
    }else{
#pragma omp parallel
        {
            int n_threads = omp_get_num_threads();
            int thread_index = omp_get_thread_num();
            int batch_size = ceil(double(n_phenotypes)/n_threads);
            int end = (1 +thread_index)*batch_size;
            if(thread_index + 1 == n_threads)
                end = n_phenotypes;
            Sort_ID_Sets(phenotype_buffer, all_ids, skip_ids, \
                         n_sets_vector, eigen_data_indices, skip_id_sets,\
                         thread_index*batch_size, end);
        }
    }
    n_sets = skip_id_sets.size();
    std::cout << n_sets << std::endl;
    eigen_data = new Eigen_Data*[skip_id_sets.size()];
//#pragma omp parallel for
    for(size_t data_index = 0; data_index < n_sets; data_index++){
	eigen_data[data_index] = new Eigen_Data(all_ids,skip_id_sets[data_index], n_sets_vector[data_index]);
    }
  /*  if(skip_id_sets.size() == 1){
        eigen_data[0] = new Eigen_Data(all_ids, skip_id_sets[0], n_sets_vector[0]);
    }else{
#pragma omp parallel
        {
            int n_threads = omp_get_num_threads();
            int thread_index = omp_get_thread_num();
            int batch_size = ceil(double(skip_id_sets.size())/n_threads);
            if(batch_size == 0) batch_size = 1;
            int end = (1 +thread_index)*batch_size;
            if(thread_index*batch_size < skip_id_sets.size() && (1 + thread_index)*batch_size >= skip_id_sets.size()){
                end = skip_id_sets.size();
            }
            for(size_t index = thread_index*batch_size; index < end ; index++){
                eigen_data[index] = new Eigen_Data(all_ids, skip_id_sets[index], n_sets_vector[index]);
            }
        
        }
        
    }*/
    size_t column_indices[skip_id_sets.size()];
    for(unsigned i = 0; i < skip_id_sets.size(); i++) column_indices[i] = 0;
  //  if(n_phenotypes < 100){
  #pragma omp parallel for
        for(int index = 0; index < n_phenotypes; index++){
          size_t column_index;   
	 #pragma omp critical
	    {
	    column_index = column_indices[eigen_data_indices[index]]++;
           }
 	     eigen_data[eigen_data_indices[index]]->set_phenotype_column(column_index, trait_names[index],phenotype_buffer + index*all_ids.size());
        }
   // }else{
/*#pragma omp parallel
        {
            int n_threads = omp_get_num_threads();
            int thread_index = omp_get_thread_num();
            int batch_size = floor(double(n_phenotypes)/n_threads);
            int end = (1 +thread_index)*batch_size;
            if(thread_index + 1 == n_threads)
                end = n_phenotypes;
            for(int index = thread_index*batch_size; index < end; index++){
                size_t column_index;
                #pragma omp critical
                {
                    column_index = column_indices[eigen_data_indices[index]]++;
                }
                eigen_data[eigen_data_indices[index]]->set_phenotype_column(column_index, trait_names[index],phenotype_buffer + index*all_ids.size());
            }
        }*/
  //  }
    delete [] phenotype_buffer;
}
