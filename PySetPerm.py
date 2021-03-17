import pandas as pd
import pyranges as pr
import numpy as np
import pickle
from itertools import repeat
import random
from scipy.stats import rankdata
import inspect


#--- classes
class FeaturesClass:
    #constructor
    def __init__(self,feature_file,range_modification):
        self.feature_file = feature_file
        self.range_modification = range_modification
        self.features = self._read_feature_file()
        self.features['idx'] = np.arange(len(self.features))
        self.features_user_def = self._feature_definition()
    
    def _read_feature_file(self):
        feature_table=pd.read_table(
            self.feature_file,
            header=0,
            names=['Chromosome',"Start","End","feature"],
            dtype={"Chromosome":str,"Start":int,"End":int,"feature":str}
            )
        return(feature_table)
    
    def _feature_definition(self):
        ftable=self.features.copy(deep=True)
        ftable['Start']=ftable['Start']-self.range_modification
        ftable['End']=ftable['End']+self.range_modification
        return(ftable)

class AnnotationSetClass:
    # constructor
    def __init__(self,annotation_file,features,min_size):
        # initializing instance variable
        self.annotation_file = annotation_file
        self.min_size = min_size
        self.annotation_sets = self._get_annotation_sets(features)
        self.annotation_array, self.annotation_array_ids = self._annotation_sets_to_array(features)
    
    def _get_annotation_sets(self,features):
        annotation_sets=pd.read_table(
            self.annotation_file,
            dtype={"id":str,"feature":str,"name":str}
            )
        #annotation_sets =annotation_sets.set_index('feature')
        return(annotation_sets)
    
    def _annotation_sets_to_array(self,features):
        sets = self.annotation_sets.join(features.set_index('feature'),on='feature').groupby('id')['idx'].apply(list)
        set_array =  [s for s in sets if len(s) >= self.min_size]
        set_array = np.sort(listnp_to_padded_nparray(set_array))
        set_names =  [i for i,s in enumerate(sets) if len(s) >= self.min_size]
        set_names = sets.index[set_names]
        return set_array, set_names



class InputClass:
    def __init__(self,candidate_file,background_file,features,annotation_array):
        # initializing instance variable
        self.candidate_file = candidate_file
        self.background_file = background_file
        self.candidates = _read_variant_file(candidate_file)
        self.background = _read_variant_file(background_file)
        self.background_features = _intersect_variants_features(self.background,features)
        self.background_features = _feature_list(self.background_features)
        self.candidate_features = _intersect_variants_features(self.candidates,features)
        self.candidate_array = _candidate_array(self.candidate_features)
        self.n_candidate_features = permutation_fset_intersect( (self.candidate_array,annotation_array) )
      
    def _read_variant_file(input_file):
        variant_table=pd.read_table(
            input_file,
            header=0,
            names=['Chromosome',"Start","End"],
            dtype={"Chromosome":str,"Start":int,"End":int}
            )
        return(variant_table)
    
    def _intersect_variants_features(variants, features):
        vtable=pr.PyRanges(variants).join(pr.PyRanges(features)).df
        return(vtable)
    
    def _feature_list(ftable):
        ftable['id']=ftable.Chromosome.astype(str).str.cat(ftable.Start.astype(str), sep='_')
        ftable=eastern_bg_genes.groupby('id')['idx'].apply(list).tolist()
        return(ftable)
    
    def _candidate_array(candidate_feature_mapped):
        features=np.asarray(np.unique(candidate_feature_mapped['idx']))
        out=np.ndarray((1,np.size(features)),dtype='uint16')
        out[0]=features
        return(out.astype('uint16'))



class PermutationClass:
	# constructor
	def __init__(self,background_features,n_candidate_features,n_permutations,n_cores)
	



	

#--- complete functions
    
def permutation_fset_intersect(args):
    start = time.perf_counter()
    permutation_array = args[0]
    function_array = args[1]
    max_z = max(permutation_array.max(), function_array.max()) + 1

    def csr_sparse(A, z):
    	m, n = A.shape
    	indptr = np.arange(0, m*n+1, n)
    	data = np.ones(m*n, dtype=np.uint16)
    	return csr_matrix((data, A.ravel(), indptr), shape=(m,z))

    intersection = csr_sparse(permutation_array,max_z ) * csr_sparse(function_array, max_z).T
    intersection = intersection.todense()
    finish = time.perf_counter()
    #print(f'This process took {round(finish-start, 2)} second(s)')
    return np.squeeze(np.asarray(intersection))



def sample_from_feature_list(feature_list,n_total):
    out = pd.unique([item for sublist in sample(feature_list,n_total) for item in sublist])
    while len(out) < n_total:
        out=np.append(out,pd.unique([item for sublist in sample(feature_list,n_total) for item in sublist]))
        out=pd.unique(out)
    out = out[:n_total]
    out = np.sort(out)
    return(out.astype('uint16'))

def array_of_resamples(feature_list,n_total,n_reps):
    out=np.ndarray((n_reps,n_total),dtype='uint16')
    for i in range(n_reps):
        out[i] = sample_from_feature_list(feature_list, n_total)
    return(out)

def array_of_resamples_tup(args):
    feature_list,n_total,n_reps = args[0],args[1],args[2]
    out=np.ndarray((n_reps,n_total),dtype='uint16')
    for i in range(n_reps):
        out[i] = sample_from_feature_list(feature_list, n_total)
    return(out)

def n_jobs_core_list(n_reps,n_cores):
	quotient, remainder = divmod(n_reps, n_cores)
	n_per_core = [quotient]*n_cores
	for i in range(remainder):
		n_per_core[i] = n_per_core[i]+1
	return(n_per_core)

def multicore_resample(feature_list,n_features,n_reps,n_cores):
	n_per_core = n_jobs_core_list(n_reps, n_cores)
	with concurrent.futures.ProcessPoolExecutor(max_workers=n_cores) as executor:
		results = executor.map(array_of_resamples_tup,zip(repeat(feature_list),repeat(n_features),n_per_core))
	results = list(results)
	return(np.concatenate(results))

def listnp_to_padded_nparray(listnp):
    max_width=np.max([np.size(l) for l in listnp])
    padded_array=np.asarray([np.pad(l,(0,max_width-np.size(l)),mode='constant',constant_values=(0,0) ) for l in listnp])
    return(padded_array.astype('uint16'))

def annotation_sets_to_array(annotation,features,min_size=3):
	sets = annotation.join(features.set_index('feature'),on='feature').groupby('id')['idx'].apply(list)
	set_array =  [s for s in sets if len(s) >= min_size]
	set_array = np.sort(listnp_to_padded_nparray(set_array))
	set_names =  [i for i,s in enumerate(sets) if len(s) >= min_size]
	set_names = sets.index[set_names]
	return set_array, set_names

def random_check_intersection(n_per_set,perms,sets,check_n):
	check_idxs = []
	n_perms=np.shape(perms)[0]
	n_sets=np.shape(sets)[0]
	for i in range(check_n):
		j=sample(range(0,n_perms-1),1)[0]
		k=sample(range(0,n_sets-1),1)[0]
		check_idxs.append( len(set(perms[j]).intersection(set(sets[k])))==n_per_set[j][k] )
	return(check_idxs)



def calculate_p_values(c_set_n, p_set_n):
	p_e = []
	p_d = []
	n_perm =  p_set_n.shape[0]
	if len(n_per_set.shape)==1:
		p_e.append((np.size(np.where(p_set_n>=c_set_n))+1)/(n_perm+1))
		p_d.append((np.size(np.where(p_set_n<=c_set_n))+1)/(n_perm+1))
	else: 
		for i in range(p_set_n.shape[1]):	
			p_e.append((np.size(np.where(p_set_n[:,i]>=c_set_n[i]))+1)/(n_perm+1))
			p_d.append((np.size(np.where(p_set_n[:,i]<=c_set_n[i]))+1)/(n_perm+1))
	return([p_e,p_d])

def make_results_table(cand_features,annotation,filtered_annotation_ids,mean_per_set,p_lists,n_per_set):
	cand_set=set(cand_features['feature'].values)
	cand_genes_in_sets = annotation.groupby('id')['feature'].apply(lambda x: np.unique(list(set(x).intersection(cand_set))))
	cand_genes_in_sets = pd.DataFrame(cand_genes_in_sets[pd.Index(filtered_annotation_ids)])
	cand_genes_in_sets = cand_genes_in_sets.reset_index(level=['id'])
	cand_genes_in_sets.columns = ['id','candidate_features']
	n_genes_in_sets = pd.DataFrame(annotation.groupby(['id','name'])['feature'].nunique()[filtered_annotation_ids])
	n_genes_in_sets = n_genes_in_sets.reset_index(level=['id', 'name'])
	out=n_genes_in_sets.join(cand_genes_in_sets.set_index('id'),on='id')
	out['n_candidates_in_set'] = [len(l) for l in out['candidate_features']]
	out['mean_permutation_n'] = mean_per_set
	out_col_names = out.columns.values
	out_col_names[2] = 'feature_set_n'
	out.columns = out_col_names
	out['emp_p_e'] = p_lists[0]
	out['emp_p_d'] = p_lists[1]
	out['fdr_e'] = fdr_from_p_matrix(n_per_set,out['emp_p_e'],method='enrichment')
	out['fdr_d'] = fdr_from_p_matrix(n_per_set,out['emp_p_d'],method='depletion')
	out['BH_fdr_e'] = p_adjust_bh(out['emp_p_e'])
	out['BH_fdr_d'] = p_adjust_bh(out['emp_p_d'])
	out = out.sort_values('emp_p_e')
	return(out)


def perm_p_matrix(perm_n_per_set,method='enrichment'):
	n_perms, n_sets = perm_n_per_set.shape
	out=np.ndarray((n_perms,n_sets),dtype='float64')
	method_int = 1
	if method=='enrichment':
		method_int = -1
	for i in range(n_sets):
		out[:,i]=rankdata(method_int*n_per_set[:,i],method='max')/n_perms
	return(out)

def fdr_from_p_matrix(perm_n_per_set,obs_p,method='enrichment'):
	p_matrix=perm_p_matrix(perm_n_per_set,method)
	obs_p_arr = np.asarray(obs_p)
	n_perm = p_matrix.shape[0]
	fdr_p = np.empty(len(obs_p),dtype='float64')
	obs_order = np.argsort(obs_p_arr)
	p_val, p_counts = np.unique(p_matrix,return_counts=True)
	current_max_fdr = 0
	for i, p_idx in enumerate(obs_order):
		if current_max_fdr==1:
			fdr_p[p_idx]=1
		else:
			obs=np.size(np.where(obs_p_arr<=obs_p_arr[p_idx]))
			exp=np.sum(p_counts[np.where(p_val <= obs_p_arr[p_idx])])/n_perm
			i_fdr = exp/obs
			if i_fdr > current_max_fdr and i_fdr < 1:
				fdr_p[p_idx] = i_fdr
				current_max_fdr = i_fdr
			elif i_fdr < current_max_fdr and i_fdr < 1:
				fdr_p[p_idx] = current_max_fdr
			else:
				fdr_p[p_idx] = 1
				current_max_fdr = 1
	return(fdr_p)


