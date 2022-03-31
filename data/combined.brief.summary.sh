# combined branches - brief summary

#########################################################################################################################
## Rationale:
## JS:
## Some analysis suggests that there is significant gene length bias in calculating min-p (but not mean-p) per gene. 
# However, i think min-p is the most relevant approach, as it is more similar to running a GO on outlier SNPs. 
# So I have produced some length adjusted, input files for gowinda. I think I replicated the right SNPs in the files. 
# I did not produce the siv-responsive set. ## whether to generate this also following josh's procedure?
# I suggest running GOWINDA only on the new length corrected min-p candidate files. 
# I did a self coded test in R, on a few combinations of statistics (I don't need to calculate a gene min-p-value) and see that Th cell pathway is fdr significant - suggesting that Harvi's results do have signal. So lets try these, and see if we can enrich for signal in these results

#Analysis here: [genelength_stats.html]
(https://liveuclac.sharepoint.com/sites/ChimpsInternalBranchAdaptation/Shared Documents/General/genelength_stats.html)

#This file also has locations of the new candidate snp files in dropbox, but the root is
 Dropbox/Gene_lists/GOWINDA_input/length.corrected

"threepclr.focal.snps.genes.minimum.p.i.txt"

 due to the process of assigning the min-p per gene: the longer the gene is, the more draws it takes of single SNP p-values, 
 the more likely it samples a very low p-value.

Lets address this by binning genes by log10(length), and then z-transform the -log10(p-values) in each bin of length.


#########################################################################################################################

# example gowinda run -
#-----------------------------------------------------------------------------------------------------------------------
## A + C - merged lt + ht vips
## directory holding Gowinda-1.12.jar
gowinda="/Volumes/Ultra USB 3.0/UCL/anc.ce.selection/data/gene_enrichment_analyses/chimp_go-2" 
## directory with the GO input files
genesets="/Volumes/Ultra USB 3.0/UCL/anc.ce.selection/data/gene_enrichment_analyses/chimp_go-2" 
## directory with your background and candidate SNP list
input="/Volumes/Ultra USB 3.0/UCL/anc.ce.selection/data/gene_enrichment_analyses/KEGG.analysis/GOWINDA_input/length.corrected/minimum_p/candidate_genes/merged"
## directory to output results.
output="/Volumes/Ultra USB 3.0/UCL/anc.ce.selection/data/gene_enrichment_analyses/further.analysis"
## snp file dir
inputsnp="/Volumes/Ultra USB 3.0/UCL/anc.ce.selection/data/gene_enrichment_analyses/KEGG.analysis/GOWINDA_input/total_SNP_files/all_genes/merged"

java -Xmx6g -jar ${gowinda}/Gowinda-1.12.jar \
--snp-file ${inputsnp}/internal.daughter.all.genes.midpoint.uniq.txt \
--candidate-snp-file ${input}/internal.central.min.p.tail.midpoint.length.corrected.uniq.txt \
--gene-set-file ${genesets}/mergedlthtvip.gene.set_and.all.txt \
--annotation-file ${genesets}/chimp.genes.human.names.homologs.and.amb.homology_mar.17.gtf \
--simulations 50000 \
--min-significance 1 \
--gene-definition updownstream2000 \
--output-file ${output}/mergedlthtvip.all.candidates.min.p.internal.central.length.corrected.txt \
--mode gene \
--min-genes 3 \
--threads 4
#-----------------------------------------------------------------------------------------------------------------------
