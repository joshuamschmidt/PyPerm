set_perm
================
Joshua Schmidt

-   [Introduction.](#introduction)
-   [A worked example](#a-worked-example)
-   [Basic usage](#basic-usage)
## Introduction.

set_perm is a python package to facilitate gene set enrichment tests
i.e. GO set enrichment. It is heavily inspired by the excellent GOWINDA
tool, developed by Robert Kofler.

While achieving similar results as GOWINDA, namely feature (gene) lengh
bias correction, PySetPerm adds additional functionality.

First, set_perm enables testing of depletion in addition to enrichment.

Second, set_perm enables testing of the enrichment/depletion of
combined candidate sets (joint distribution). The idea is best
illustrated by an hypothetical example.

Suppose a method to detect loci under positive selection is used in popA
and popB, identifying several candidate loci in each population. We
could test if candidates in either of the populations are enriched for
particular functions. This would be a run of the mill GO test for popA and popB
individually, where we might be lucky to find a significant enrichment in the
same gene set for both taxa. But this is likely underpowered.

Perhaps a more interesting question, biologically and statistically, is to
ask if the joint distribution of popA and popB candidate loci is
enriched for a particular function - this could suggest that a common
pathway has been under selection.

The aim for set_perm is to develop a simple, extensible framework for
developing these tests. At the moment simple testing of joint candidate
set distributions is implemented. Future plans include equivalent
explicit union and intersect tests, varied block permutation options and
more!

## A worked example
Go to the notebook here for an example workflow: [chimp example](/test_anlaysis.ipynb)


## Basic usage
`set_perm \
--candidates "central,data/central-0.000192-candidate.snps.bed.gz" "eastern,data/eastern-0.000228-candidate.snps.bed.gz" "ancestral,data/ancestral-0.005-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/kegg.txt \
--min_set_size 10 \
--n_perms 10000 \
--prefix "test" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 10`

### Units of analysis

### Create objects

### join objects
.join\_objects() functions for both Input and SetPerPerm objects enable
fast creation of union sets.


