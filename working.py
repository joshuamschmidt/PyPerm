#----
species = ['eastern','central','internal']
start=time.perf_counter()
eastern_background_file='/Users/joshuaschmidt/Projects/pypermr/eastern_background.txt'
eastern_background=pypermr.read_variant_file(eastern_background_file)
eastern_bg_genes=pypermr.intersect_variants_features(eastern_background,genes_2kb)
eastern_bg_genes=feature_list(eastern_bg_genes)
eastern_candidate_genes=pypermr.intersect_variants_features(eastern_candidates,genes_2kb)
n_eastern_candidates=pypermr.n_features(eastern_candidate_genes)

n_perms=100000
n_cores=4




#--- annotation sets
annotation_file ='/Users/joshuaschmidt/Projects/pypermr/go.txt'
annotation_sets = get_annotation_sets(annotation_file,genes)
annotation_array, annotation_array_ids = annotation_sets_to_array(annotation_sets,genes,min_size=5)

#--- permutations
eastern_bg_perms=multicore_resample(eastern_bg_genes,n_eastern_candidates,n_perms,n_cores)

n_per_set = permutation_fset_Mintersect((eastern_bg_perms,annotation_array))
mean_per_set = np.array(np.mean(n_per_set,axis=0))
# check against slower set.intersection()
#np.all(random_check_intersection(n_per_set,eastern_bg_perms,annotation_array,1000))
candidate_feature_array = candidate_array(eastern_candidate_genes)
candidate_set_n = permutation_fset_Mintersect( (candidate_feature_array,annotation_array) )
p_lists = calculate_p_values(candidate_set_n, n_per_set)
results_table = make_results_table(eastern_candidate_genes,annotation_sets,annotation_array_ids,mean_per_set,p_lists,n_per_set)
end=time.perf_counter()
print(f'total time taken was {(end-start):.2f} seconds')

#--- prototype functions
feature_file = '/Users/joshuaschmidt/Projects/pypermr/genes.txt'
features = FeaturesClass(feature_file, 2000)
annotation_file ='/Users/joshuaschmidt/Projects/pypermr/go.txt'
annotation_set = AnnotationSetClass(annotation_file,genes_2kb,5)


