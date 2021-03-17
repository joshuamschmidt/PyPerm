import pandas as pd
import pyranges as pr
import numpy as np
import concurrent.futures as cf
from itertools import repeat
from scipy.stats import rankdata
from scipy.sparse import csr_matrix
# import inspect
from random import sample
# import pickle


# --- global functions
def permutation_fset_intersect(args):
    permutation_array = args[0]
    function_array = args[1]
    max_z = max(permutation_array.max(), function_array.max()) + 1

    def csr_sparse(a, z):
        m, n = a.shape
        indptr = np.arange(0, m * n + 1, n)
        data = np.ones(m * n, dtype=np.uint16)
        return csr_matrix((data, a.ravel(), indptr), shape=(m, z))

    intersection = csr_sparse(permutation_array, max_z) * csr_sparse(function_array, max_z).T
    intersection = intersection.todense()
    return np.squeeze(np.asarray(intersection))


def listnp_to_padded_nparray(listnp):
    max_width = np.max([np.size(sublist) for sublist in listnp])
    padded_array = np.asarray(
        [np.pad(sublist, (0, max_width - np.size(sublist)), mode='constant', constant_values=(0, 0))
         for sublist
         in listnp])
    return padded_array.astype('uint16')


def annotation_sets_to_array(annotation, features, min_size=3):
    sets = annotation.join(features.set_index('feature'), on='feature').groupby('id')['idx'].apply(list)
    set_array = [s for s in sets if len(s) >= min_size]
    set_array = np.sort(listnp_to_padded_nparray(set_array))
    set_names = [i for i, s in enumerate(sets) if len(s) >= min_size]
    set_names = sets.index[set_names]
    return set_array, set_names


def sample_from_feature_list(feature_list, n_total):
    out = pd.unique([item for sublist in sample(feature_list, n_total) for item in sublist])
    while len(out) < n_total:
        out = np.append(out, pd.unique([item for sublist in sample(feature_list, n_total) for item in sublist]))
        out = pd.unique(out)
    out = out[:n_total]
    out = np.sort(out)
    return out.astype('uint16')


def array_of_resamples_tup(args):
    feature_list, n_total, n_reps = args[0], args[1], args[2]
    out = np.ndarray((n_reps, n_total), dtype='uint16')
    for i in range(n_reps):
        out[i] = sample_from_feature_list(feature_list, n_total)
    return out


def n_jobs_core_list(n_reps, n_cores):
    quotient, remainder = divmod(n_reps, n_cores)
    n_per_core = [quotient] * n_cores
    for i in range(remainder):
        n_per_core[i] = n_per_core[i] + 1
    return n_per_core


def multicore_resample(n_features, n_reps, n_cores, feature_list):
    n_per_core = n_jobs_core_list(n_reps, n_cores)
    with cf.ProcessPoolExecutor(max_workers=n_cores) as executor:
        results = executor.map(array_of_resamples_tup, zip(repeat(feature_list), repeat(n_features), n_per_core))
    results = list(results)
    return np.concatenate(results)


def multicore_intersect(permutation_array, functionalset_array, n_cores):
    split_permutation_array = np.array_split(permutation_array, n_cores)
    with cf.ProcessPoolExecutor(max_workers=n_cores) as executor:
        results = executor.map(permutation_fset_intersect, zip(split_permutation_array, repeat(functionalset_array)))
    results = list(results)
    return np.concatenate(results)


def calculate_p_values(c_set_n, p_set_n):
    p_e = []
    p_d = []
    n_perm = p_set_n.shape[0]
    if n_perm == 1:
        p_e.append((np.size(np.where(p_set_n >= c_set_n)) + 1) / (n_perm + 1))
        p_d.append((np.size(np.where(p_set_n <= c_set_n)) + 1) / (n_perm + 1))
    else:
        for i in range(p_set_n.shape[1]):
            p_e.append((np.size(np.where(p_set_n[:, i] >= c_set_n[i])) + 1) / (n_perm + 1))
            p_d.append((np.size(np.where(p_set_n[:, i] <= c_set_n[i])) + 1) / (n_perm + 1))
    return p_e, p_d


def make_results_table(cand_features, annotations, set_per_perm):
    cand_set = set(cand_features['feature'].values)
    cand_genes_in_sets = annotations.annotation_sets.groupby('id')['feature'].apply(
        lambda x: np.unique(list(set(x).intersection(cand_set))))
    cand_genes_in_sets = pd.DataFrame(cand_genes_in_sets[pd.Index(annotations.annotation_array_ids)])
    cand_genes_in_sets = cand_genes_in_sets.reset_index(level=['id'])
    cand_genes_in_sets.columns = ['id', 'candidate_features']
    n_genes_in_sets = pd.DataFrame(
        annotations.annotation_sets.groupby(['id', 'name'])['feature'].nunique()[annotations.annotation_array_ids])
    n_genes_in_sets = n_genes_in_sets.reset_index(level=['id', 'name'])
    out = n_genes_in_sets.join(cand_genes_in_sets.set_index('id'), on='id')
    out['n_candidates_in_set'] = [len(features) for features in out['candidate_features']]
    out['mean_permutation_n'] = set_per_perm.mean_per_set
    out_col_names = out.columns.values
    out_col_names[2] = 'feature_set_n'
    out.columns = out_col_names
    out['emp_p_e'] = set_per_perm.p_enrichment
    out['emp_p_d'] = set_per_perm.p_depletion
    out['fdr_e'] = fdr_from_p_matrix(set_per_perm.setNperPerm, out['emp_p_e'], method='enrichment')
    out['fdr_d'] = fdr_from_p_matrix(set_per_perm.setNperPerm, out['emp_p_d'], method='depletion')
    out['BH_fdr_e'] = p_adjust_bh(out['emp_p_e'])
    out['BH_fdr_d'] = p_adjust_bh(out['emp_p_d'])
    out = out.sort_values('emp_p_e')
    return out


def fdr_from_p_matrix(perm_n_per_set, obs_p, method='enrichment'):
    p_matrix = perm_p_matrix(perm_n_per_set, method)
    obs_p_arr = np.asarray(obs_p)
    n_perm = p_matrix.shape[0]
    fdr_p = np.empty(len(obs_p), dtype='float64')
    obs_order = np.argsort(obs_p_arr)
    p_val, p_counts = np.unique(p_matrix, return_counts=True)
    current_max_fdr = 0
    for i, p_idx in enumerate(obs_order):
        if current_max_fdr == 1:
            fdr_p[p_idx] = 1
        else:
            obs = np.size(np.where(obs_p_arr <= obs_p_arr[p_idx]))
            exp = np.sum(p_counts[np.where(p_val <= obs_p_arr[p_idx])]) / n_perm
            i_fdr = exp / obs
            if current_max_fdr < i_fdr < 1:
                fdr_p[p_idx] = i_fdr
                current_max_fdr = i_fdr
            elif current_max_fdr > i_fdr and i_fdr < 1:
                fdr_p[p_idx] = current_max_fdr
            else:
                fdr_p[p_idx] = 1
                current_max_fdr = 1
    return fdr_p


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


# --- classes
class Features:
    # constructor
    def __init__(self, feature_file, range_modification):
        self.feature_file = feature_file
        self.range_modification = range_modification

        def _read_feature_file(file):
            feature_table = pd.read_table(
                file,
                header=0,
                names=['Chromosome', "Start", "End", "feature"],
                dtype={"Chromosome": str, "Start": int, "End": int, "feature": str}
            )
            return feature_table

        self.features = _read_feature_file(self.feature_file)
        self.features['idx'] = np.arange(len(self.features))

        def _feature_definition(features, range_mod):
            ftable = features.copy(deep=True)
            ftable['Start'] = ftable['Start'] - range_mod
            ftable['End'] = ftable['End'] + range_mod
            return ftable

        self.features_user_def = _feature_definition(self.features, self.range_modification)


class AnnotationSet:
    # constructor
    def __init__(self, annotation_file, features, min_size):
        # initializing instance variable
        self.annotation_file = annotation_file
        self.min_size = min_size

        def _get_annotation_sets(a_file):
            annotation_sets = pd.read_table(
                a_file,
                dtype={"id": str, "feature": str, "name": str}
            )
            return annotation_sets

        def _annotation_sets_to_array(annotation_sets, feature, min_s):
            sets = annotation_sets.join(feature.set_index('feature'), on='feature').groupby('id')['idx'].apply(list)
            set_array = [s for s in sets if len(s) >= min_s]
            set_array = np.sort(listnp_to_padded_nparray(set_array))
            set_names = [i for i, s in enumerate(sets) if len(s) >= min_s]
            set_names = sets.index[set_names]
            return set_array, set_names

        self.annotation_sets = _get_annotation_sets(self.annotation_file)
        self.annotation_array, self.annotation_array_ids = _annotation_sets_to_array(self.annotation_sets, features,
                                                                                     self.min_size)


class Input:
    # constructor
    def __init__(self, candidate_file, background_file, features, annotation_array):
        # initializing instance variable
        self.candidate_file = candidate_file
        self.background_file = background_file

        def _read_variant_file(input_file):
            variant_table = pd.read_table(
                input_file,
                header=0,
                names=['Chromosome', "Start", "End"],
                dtype={"Chromosome": str, "Start": int, "End": int}
            )
            return variant_table

        def _intersect_variants_features(variants, features_table):
            vtable = pr.PyRanges(variants).join(pr.PyRanges(features_table)).df
            return vtable

        def _feature_list(ftable):
            ftable['id'] = ftable.Chromosome.astype(str).str.cat(ftable.Start.astype(str), sep='_')
            ftable = ftable.groupby('id')['idx'].apply(list).tolist()
            return ftable

        def _candidate_array(candidate_feature_mapped):
            feature_array = np.asarray(np.unique(candidate_feature_mapped['idx']))
            out = np.ndarray((1, np.size(feature_array)), dtype='uint16')
            out[0] = feature_array
            return out.astype('uint16')

        self.candidates = _read_variant_file(self.candidate_file)
        self.background = _read_variant_file(self.background_file)
        self.background_features = _intersect_variants_features(self.background, features)
        self.background_features = _feature_list(self.background_features)
        self.candidate_features = _intersect_variants_features(self.candidates, features)
        self.candidate_array = _candidate_array(self.candidate_features)
        self.n_candidates = np.size(self.candidate_array)
        self.n_candidate_per_function = permutation_fset_intersect((self.candidate_array, annotation_array))


class Permutation:
    # constructor
    def __init__(self, background_features, n_candidate_features, n_permutations, n_cores):
        self.n_candidate_features = n_candidate_features
        self.n_permutations = n_permutations
        self.n_cores = n_cores
        self.permutations = multicore_resample(self.n_candidate_features, self.n_permutations, self.n_cores,
                                               background_features)


class SetPerPerm:
    # constructor
    def __init__(self, permutations, annotation_sets, n_cores, n_candidates_per_set):
        self.setNperPerm = multicore_intersect(permutations, annotation_sets, n_cores)
        self.mean_per_set = np.array(np.mean(self.setNperPerm, axis=0))
        self.p_enrichment, self.p_depletion = calculate_p_values(n_candidates_per_set, self.setNperPerm)


# --- redundant and/or not used anymore

def perm_p_matrix(perm_n_per_set, method='enrichment'):
    n_perms, n_sets = perm_n_per_set.shape
    out = np.ndarray((n_perms, n_sets), dtype='float64')
    method_int = 1
    if method == 'enrichment':
        method_int = -1
    for i in range(n_sets):
        out[:, i] = rankdata(method_int * perm_n_per_set[:, i], method='max') / n_perms
    return out


def array_of_resamples(feature_list, n_total, n_reps):
    out = np.ndarray((n_reps, n_total), dtype='uint16')
    for i in range(n_reps):
        out[i] = sample_from_feature_list(feature_list, n_total)
    return out


def random_check_intersection(n_per_set, perms, sets, check_n):
    check_idxs = []
    n_perms = np.shape(perms)[0]
    n_sets = np.shape(sets)[0]
    for i in range(check_n):
        j = sample(range(0, n_perms - 1), 1)[0]
        k = sample(range(0, n_sets - 1), 1)[0]
        check_idxs.append(len(set(perms[j]).intersection(set(sets[k]))) == n_per_set[j][k])
    return check_idxs
