PySetPerm
================
Joshua Schmidt
18/03/2021

-   [Introduction.](#introduction.)
-   [A worked example](#a-worked-example)

## Introduction.

PySetPerm is a python package to facilitate gene set enrichment tests
i.e. GO set enrichment. It is heavily inspired by the excellent GOWINDA
tool, developed by Robert Kofler.

While achieving similar results as GOWINDA, namely feature (gene) lengh
bias correction, PySetPerm adds additional functionality.

First, PySetPerm enables testing of depletion in addition to enrichment.

Second, PySetPerm enables testing of the enrichment/depletion of
combined candidate sets. The idea is best illustrated by an hypothetical
example.

Suppose a method to detect loci under positive selection is used in popA
and popB, identifying several candidate loci in each population. We
could test if candidates in either of the populations are enriched for
particular functions. This would be a common GO test for popA and popB
individually. Perhaps a more interesting question, biologically, is to
ask if the union of popA and popB candidate loci is enriched for a
particular function - this could suggest that a common pathway has been
under selection.

The aim for PySetPerm is to develop a simple, extensible framework for
developing these tests. At the moment simple testing of candidate set
unions is implemented. Future plans include equivalent intersect tests,
block permutations and more!

## A worked example
