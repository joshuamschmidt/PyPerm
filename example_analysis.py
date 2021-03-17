import pysetperm as psp
import numpy as np
import pandas as pd
# used for all sub analyses
features = psp.Features('data/genes.txt', 2000)
annotations = psp.AnnotationSet('data/kegg.txt', features.features_user_def, 5)
n_perms = 50000
cores = 4
# specific inputs
e_input = psp.Input('data/eastern_candidates.txt',
                    'data/eastern_background.txt.gz',
                    features.features_user_def,
                    annotations.annotation_array)

# permutations
e_permutations = psp.Permutation(e_input.background_features, e_input.n_candidates, n_perms, cores)

e_per_set = psp.SetPerPerm(e_permutations.permutations,
                           annotations.annotation_array,
                           cores,
                           e_input.n_candidate_per_function)
