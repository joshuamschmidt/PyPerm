PySetPerm
================
Joshua Schmidt

-   [Introduction.](#introduction)
-   [A worked example](#a-worked-example)
-   [Basic usage](#basic-usage)
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
particular functions. This would be a run of the mill GO test for popA and popB
individually, where we might be lucky to find a significant enrichment in the
same gene set for both taxa. But this is likely underpowered.

Perhaps a more interesting question, biologically and statistically, is to
ask if the joint distribution of popA and popB candidate loci is
enriched for a particular function - this could suggest that a common
pathway has been under selection.

The aim for PySetPerm is to develop a simple, extensible framework for
developing these tests. At the moment simple testing of joint candidate
set distributions is implemented. Future plans include equivalent
explicit union and intersect tests, varied block permutation options and
more!

## A worked example
Go to the notebook here for an example workflow: [chimp example](/test_anlaysis.ipynb)


## Basic usage

### Units of analysis

### Create objects

### join objects
.join\_objects() functions for both Input and SetPerPerm objects enable
fast creation of union sets.


