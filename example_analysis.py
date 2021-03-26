import pysetperm as psp
 import numpy as np
 import pandas as pd
# used for all sub analyses
n_perms = 100000
cores = 10
annotations = psp.AnnotationSet(annotation_file='data/genes.txt', range_modification=2000)
function_sets = psp.FunctionSets(function_set_file='data/kegg.txt', min_set_size=10, annotation_obj=annotations)
# specific inputs
e_candidates = psp.Variants(variant_file='data/eastern_candidates.txt')
e_candidates.annotate_variants(annotation_obj=annotations)
e_background = psp.Variants(variant_file='data/eastern_background.txt.gz')
e_background.annotate_variants(annotation_obj=annotations)

#c_candidates = psp.Variants(variant_file='data/central_candidates.txt')
#c_background = psp.Variants(variant_file='data/central_background.txt.gz')

e_test_obj = psp.TestObject(e_candidates,
                            e_background,
                            function_sets,
                            n_cores=cores)

# permutations
e_permutations = psp.Permutation(e_test_obj, n_perms, cores)

# distributions across permutations
e_per_set = psp.SetPerPerm(e_permutations,
                           function_sets,
                           e_test_obj,
                           cores)

# results tables
e_results = psp.make_results_table(e_input, annotations, e_per_set)


