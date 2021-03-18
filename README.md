PySetPerm
================
Joshua Schmidt

-   [Introduction.](#introduction)
-   [A worked example](#a-worked-example)

## Introduction.

PySetPerm is a python package to facilitate gene set enrichment tests
i.e. GO set enrichment. It is heavily inspired by the excellent GOWINDA
tool, developed by Robert Kofler.

While achieving similar results as GOWINDA, namely feature (gene) lengh
bias correction, PySetPerm adds additional functionality.

First, PySetPerm enables testing of depletion in addition to enrichment.

Second, PySetPerm enables testing of the enrichment/depletion of
combined candidate sets (joint distribution). The idea is best
illustrated by an hypothetical example.

Suppose a method to detect loci under positive selection is used in popA
and popB, identifying several candidate loci in each population. We
could test if candidates in either of the populations are enriched for
particular functions. This would be a common GO test for popA and popB
individually. Perhaps a more interesting question, biologically, is to
ask if the joint distribution of popA and popB candidate loci is
enriched for a particular function - this could suggest that a common
pathway has been under selection.

The aim for PySetPerm is to develop a simple, extensible framework for
developing these tests. At the moment simple testing of joint candidate
set distributions is implemented. Future plans include equivalent
explicit union and intersect tests, varied block permutation options and
more!

## A worked example

``` python
import pysetperm as psp
features = psp.Features('data/genes.txt', 2000)
annotations = psp.AnnotationSet('data/kegg.txt', features.features_user_def, 5)
n_perms = 100000
cores = 10
```

``` python
e_input = psp.Input('data/eastern_candidates.txt',
                    'data/eastern_background.txt.gz',
                    features,
                    annotations)

c_input = psp.Input('data/central_candidates.txt',
                    'data/central_background.txt.gz',
                    features,
                    annotations)

i_input = psp.Input('data/internal_candidates.txt',
                    'data/internal_background.txt.gz',
                    features,
                    annotations)
```

``` python
e_permutations = psp.Permutation(e_input, n_perms, cores)
c_permutations = psp.Permutation(c_input, n_perms, cores)
i_permutations = psp.Permutation(i_input, n_perms, cores)

e_per_set = psp.SetPerPerm(e_permutations,
                           annotations,
                           e_input,
                           cores)
c_per_set = psp.SetPerPerm(c_permutations,
                           annotations,
                           c_input,
                           cores)
i_per_set = psp.SetPerPerm(i_permutations,
                           annotations,
                           i_input,
                           cores)
```

.join\_objects() functions for both Input and SetPerPerm objects enable
fast creation of union sets.

``` python
# combine sims
ec_input = psp.Input.add_objects(e_input, c_input)
ec_per_set = psp.SetPerPerm.add_objects(e_per_set, c_per_set)
ei_input = psp.Input.add_objects(e_input, i_input)
ei_per_set = psp.SetPerPerm.add_objects(e_per_set, i_per_set)
ci_input = psp.Input.add_objects(c_input, i_input)
ci_per_set = psp.SetPerPerm.add_objects(c_per_set, i_per_set)
eci_input = psp.Input.add_objects(ec_input, i_input)
eci_per_set = psp.SetPerPerm.add_objects(ec_per_set, i_per_set)
```

``` python
# results
e_results = psp.make_results_table(e_input, annotations, e_per_set)
c_results = psp.make_results_table(c_input, annotations, c_per_set)
i_results = psp.make_results_table(i_input, annotations, i_per_set)
ec_results = psp.make_results_table(ec_input, annotations, ec_per_set)
ei_results = psp.make_results_table(ei_input, annotations, ei_per_set)
ci_results = psp.make_results_table(ci_input, annotations, ci_per_set)
eci_results = psp.make_results_table(eci_input, annotations, eci_per_set)
```

``` r
library(data.table)
as.data.table(py$eci_results)
```
