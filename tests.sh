#--- kegg

bin/set_perm \
--candidates "central,data/central-0.000192-candidate.snps.bed.gz" "eastern,data/eastern-0.000228-candidate.snps.bed.gz" "ancestral,data/ancestral-0.005-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/kegg.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "match-ancestral_0.5%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

bin/set_perm \
--candidates "central,data/central-4e-05-candidate.snps.bed.gz" "eastern,data/eastern-3.5e-05-candidate.snps.bed.gz" "ancestral,data/ancestral-0.001-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/kegg.txt \
--min_set_size 10 \
--n_perms 100000 \
--prefix "match-ancestral_0.1%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

bin/set_perm \
--candidates "central,data/central-2e-05-candidate.snps.bed.gz" "eastern,data/eastern-1.5e-05-candidate.snps.bed.gz" "ancestral,data/ancestral-5e-04-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/kegg.txt \
--min_set_size 10 \
--n_perms 100000 \
--prefix "match-ancestral_0.05%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

#--- vip
bin/set_perm \
--candidates "central,data/central-0.000192-candidate.snps.bed.gz" "eastern,data/eastern-0.000228-candidate.snps.bed.gz" "ancestral,data/ancestral-0.005-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/all_vips.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "match-ancestral_0.5%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

bin/set_perm \
--candidates "central,data/central-4e-05-candidate.snps.bed.gz" "eastern,data/eastern-3.5e-05-candidate.snps.bed.gz" "ancestral,data/ancestral-0.001-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/vip.txt \
--min_set_size 10 \
--n_perms 100000 \
--prefix "match-ancestral_0.1%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

bin/set_perm \
--candidates "central,data/central-2e-05-candidate.snps.bed.gz" "eastern,data/eastern-1.5e-05-candidate.snps.bed.gz" "ancestral,data/ancestral-5e-04-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/vip.txt \
--min_set_size 10 \
--n_perms 100000 \
--prefix "match-ancestral_0.05%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

#--- SIV

bin/set_perm \
--candidates "central,data/central-0.000192-candidate.snps.bed.gz" "eastern,data/eastern-0.000228-candidate.snps.bed.gz" "ancestral,data/ancestral-0.005-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/siv.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "match-ancestral_0.5%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

bin/set_perm \
--candidates "central,data/central-4e-05-candidate.snps.bed.gz" "eastern,data/eastern-3.5e-05-candidate.snps.bed.gz" "ancestral,data/ancestral-0.001-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/siv.txt \
--min_set_size 10 \
--n_perms 100000 \
--prefix "match-ancestral_0.1%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

bin/set_perm \
--candidates "central,data/central-2e-05-candidate.snps.bed.gz" "eastern,data/eastern-1.5e-05-candidate.snps.bed.gz" "ancestral,data/ancestral-5e-04-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/siv.txt \
--min_set_size 10 \
--n_perms 100000 \
--prefix "match-ancestral_0.05%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

#--- SIV modules
bin/set_perm \
--candidates "central,data/central-0.000192-candidate.snps.bed.gz" "eastern,data/eastern-0.000228-candidate.snps.bed.gz" "ancestral,data/ancestral-0.005-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/siv_modules.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "match-ancestral_0.5%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

bin/set_perm \
--candidates "central,data/central-4e-05-candidate.snps.bed.gz" "eastern,data/eastern-3.5e-05-candidate.snps.bed.gz" "ancestral,data/ancestral-0.001-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/siv_modules.txt \
--min_set_size 10 \
--n_perms 100000 \
--prefix "match-ancestral_0.1%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

bin/set_perm \
--candidates "central,data/central-2e-05-candidate.snps.bed.gz" "eastern,data/eastern-1.5e-05-candidate.snps.bed.gz" "ancestral,data/ancestral-5e-04-candidate.snps.bed.gz" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.bed.gz" \
--feature_def data/genes.txt \
--function_def data/siv_modules.txt \
--min_set_size 10 \
--n_perms 100000 \
--prefix "match-ancestral_0.05%" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6


## only bottom 90 windows

bin/set_perm \
--candidates "central,data/central-0.000192-candidate.snps.bed.gz" "eastern,data/eastern-0.000228-candidate.snps.bed.gz" "ancestral,data/ancestral-0.005-candidate.snps-b90.bed" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.snps.bed" \
--feature_def data/genes.txt \
--function_def data/kegg.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "match-ancestral_0.5%-b90" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6



bin/set_perm \
--candidates "central,data/central-0.000192-candidate.snps.bed.gz" "eastern,data/eastern-0.000228-candidate.snps.bed.gz" "ancestral,data/ancestral-0.005-candidate.snps-b90.bed" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.snps.bed" \
--feature_def data/genes.txt \
--function_def data/all_vips.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "match-ancestral_0.5%-b90" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

bin/set_perm \
--candidates "central,data/central-0.000192-candidate.snps.bed.gz" "eastern,data/eastern-0.000228-candidate.snps.bed.gz" "ancestral,data/ancestral-0.005-candidate.snps-b90.bed" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.snps.bed" \
--feature_def data/genes.txt \
--function_def data/siv.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "match-ancestral_0.5%-b90" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6

bin/set_perm \
--candidates "central,data/central-0.000192-candidate.snps.bed.gz" "eastern,data/eastern-0.000228-candidate.snps.bed.gz" "ancestral,data/ancestral-0.005-candidate.snps-b90.bed" \
--background "central,data/pbsnj-bg.snps.bed.gz" "eastern,data/pbsnj-bg.snps.bed.gz" "ancestral,data/ancestral-bg.snps.bed" \
--feature_def data/genes.txt \
--function_def data/siv_modules.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "match-ancestral_0.5%-b90" \
--joint "central,eastern" "central,ancestral" "eastern,ancestral" "central,eastern,ancestral" \
--gene_def 2000 \
--threads 6




# nested ancestral

bin/set_perm_nested \
--candidates "'0.05%',data/ancestral-5e-04-candidate.snps.bed.gz" "'0.5%',data/ancestral-0.005-candidate.snps.bed.gz" "'0.1%',data/ancestral-0.001-candidate.snps.bed.gz" \
--background "data/ancestral-bg.snps.bed" \
--feature_def data/genes.txt \
--function_def data/kegg.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "ancestral_nested" \
--gene_def 2000 \
--threads 4



bin/set_perm_nested \
--candidates "0.05%,data/ancestral-5e-04-candidate.snps.bed.gz" "0.5%,data/ancestral-0.005-candidate.snps.bed.gz" "0.1%,data/ancestral-0.001-candidate.snps.bed.gz" \
--background "data/ancestral-bg.snps.bed" \
--feature_def data/genes.txt \
--function_def data/all_vips.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "ancestral_nested" \
--gene_def 2000 \
--threads 4


bin/set_perm_nested \
--candidates "0.05%,data/ancestral-5e-04-candidate.snps.bed.gz" "0.5%,data/ancestral-0.005-candidate.snps.bed.gz" "0.1%,data/ancestral-0.001-candidate.snps.bed.gz" \
--background "data/ancestral-bg.snps.bed" \
--feature_def data/genes.txt \
--function_def data/siv.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "ancestral_nested" \
--gene_def 2000 \
--threads 4


bin/set_perm_nested \
--candidates "0.05%,data/ancestral-5e-04-candidate.snps.bed.gz" "0.5%,data/ancestral-0.005-candidate.snps.bed.gz" "0.1%,data/ancestral-0.001-candidate.snps.bed.gz" \
--background "data/ancestral-bg.snps.bed" \
--feature_def data/genes.txt \
--function_def data/siv_modules.txt \
--min_set_size 3 \
--n_perms 100000 \
--prefix "ancestral_nested" \
--gene_def 2000 \
--threads 4
