# Follows validation steps used in the GOWINDA package: https://sourceforge.net/p/gowinda/wiki/Validation/

conda create -n python2 python=2.7 bedtools
conda activate python2
gowinda_v='/mnt/c/Users/jschmi06/Downloads/validation_files/validation_files'
work_dir='/mnt/c/Users/jschmi06/Documents/Projects/set_perm/validation'
# gowinda_v='/Users/joshuaschmidt/Downloads/validation_files/
# work_dir="/Users/joshuaschmidt/Projects/set_perm/validation"
python "$gowinda_v"/scripts/get_nonoverlapping_geneids.py --gtf "$gowinda_v"/Flybase.gtf \
> "$work_dir"/nonoverlapping_genes.txt

python "$gowinda_v"/scripts/get_geneids_havinggocategory.py \
--go "$gowinda_v"/association_gominer.txt \
--genelist "$work_dir"/nonoverlapping_genes.txt \
> "$work_dir"/geneids_withgoterm_nooverlap.txt

python "$gowinda_v"/scripts/create_snps_for_genes.py \
--genelist "$work_dir"/geneids_withgoterm_nooverlap.txt \
--gtf "$gowinda_v"/Flybase.gtf \
> "$work_dir"/snps_5pgene.txt

cat  "$work_dir"/snps_5pgene.txt|perl -ne 'print if rand()<0.03'|head -1000 >  "$work_dir"/rand_snps_1k.txt


# GOWINDA GO
java -Xmx4g -jar /mnt/c/Users/jschmi06/Downloads/Gowinda-1.12.jar \
--snp-file "$work_dir"/snps_5pgene.txt  \
--candidate-snp-file "$work_dir"/rand_snps_1k.txt \
--gene-set-file "$gowinda_v"/association_gominer.txt \
--annotation-file "$work_dir"/flybase_genes.gtf \
--simulations 100000 \
--min-significance 1 \
--gene-definition exon \
--threads 6 \
--output-file "$work_dir"/gowinda_res.txt \
--mode gene \
--min-genes 5

# CANDIATE SNPS to genes

cat "$work_dir"/snps_5pgene.txt| awk '{print $3}'|sort |uniq > "$work_dir"/gominer_total.txt
cat "$work_dir"/rand_snps_1k.txt| awk '{print $3}'|sort |uniq > "$work_dir"/gominer_totest.txt

# USED gprofiler webtool to convert to FBgn ids
# RESULTS IN 
# PANTHER_gominer_GO_process.txt

#---- Gene length biased results
python "$gowinda_v"/scripts/extract_geneids_fromgoaassociation.py \
--go "$gowinda_v"/association_gominer.txt \
> "$work_dir"/total_genes.txt

python "$gowinda_v"/scripts/create_length_biased_genes.py \
--gtf "$gowinda_v"/Flybase.gtf \
--genelist "$work_dir"/total_genes.txt \
>  "$work_dir"/snps_lengthbiased.txt

cat "$work_dir"/snps_lengthbiased.txt|perl -ne 'print if rand()<0.004'|head -1000> "$work_dir"/rand_1k_snps.txt

java -Xmx4g -jar /mnt/c/Users/jschmi06/Downloads/Gowinda-1.12.jar \
--snp-file "$work_dir"/snps_lengthbiased.txt  \
--candidate-snp-file "$work_dir"/rand_1k_snps.txt \
--gene-set-file "$gowinda_v"/association_gominer.txt \
--annotation-file "$gowinda_v"/Flybase.gtf \
--simulations 100000 \
--min-significance 1 \
--gene-definition exon \
--threads 6 \
--output-file "$work_dir"/gowinda_res_biased.txt \
--mode gene \
--min-genes 5

cat "$work_dir"/rand_1k_snps.txt | awk '{print $3}'|sort |uniq > "$work_dir"/gominer_totest_biased.txt

# PANTHER saved as PANTHER_gominer_GO_process_biased.txt


# MAKE set-perm compatible GO definitions
echo -e "id\tfeature\tname" > "$work_dir"/association_gominer_setperm.txt

while read p; do
	id=$(echo "$p" | cut -f1 );
	name=$(echo "$p" | cut -f2);
	genes=$(echo "$p" | cut -f3);
	gene_arr=(`echo ${genes}`);
	for gene in "${gene_arr[@]}"; do
		echo -e "$id\t$gene\t$name" >> "$work_dir"/association_gominer_setperm.txt;
  	done;
done < "$gowinda_v"/association_gominer.txt;

# set-perm gene def formats
awk 'BEGIN{OFS="\t"; print "chr", "start", "end", "gene";} {if(NR >1 && $2!="-") print $2"\t"$3"\t"$4"\t"$1}' "$work_dir"/FlyBase_Fields_download.txt \
> "$work_dir"/flybase_genes_setperm.txt

#---- I think these gene defs are too different to the GOWINDA rovided flybase gtf. Plus the mists of time and
# flybase upgrades....
# easier to use the gowinda GTF and convert it to set-perm format.
# take start of first exon to end of last exon as gene definiton.

# SNPs chr start end
awk 'BEGIN{OFS="\t"; print "chr", "start", "end";} {print $1"\t"$2-1"\t"$2}' "$work_dir"/rand_1k_snps.txt \
> "$work_dir"/rand_1k_snps.bed

awk 'BEGIN{OFS="\t"; print "chr", "start", "end";} {print $1"\t"$2-1"\t"$2}' "$work_dir"/snps_lengthbiased.txt \
> "$work_dir"/snps_lengthbiased.bed

awk 'BEGIN{OFS="\t"; print "chr", "start", "end";} {print $1"\t"$2-1"\t"$2}' "$work_dir"/rand_snps_1k.txt \
> "$work_dir"/rand_snps_1k.bed

awk 'BEGIN{OFS="\t"; print "chr", "start", "end";} {print $1"\t"$2-1"\t"$2}' "$work_dir"/snps_5pgene.txt \
> "$work_dir"/snps_5pgene.bed

conda deactivate
#conda activate pysetperm
cd $work_dir
../bin/set_perm \
--candidates "unbiased,rand_snps_1k.bed" \
--background "unbiased,snps_5pgene.bed" \
--feature_def "$work_dir"/flybase_genes_setperm.txt \
--function_def "$work_dir"/association_gominer_setperm.txt \
--min_set_size 5 \
--n_perms 100000 \
--prefix "unbiased" \
--gene_def 0 \
--threads 6

../bin/set_perm \
--candidates "biased,rand_1k_snps.bed" \
--background "biased,snps_lengthbiased.bed" \
--feature_def flybase_genes_setperm.txt \
--function_def association_gominer_setperm.txt \
--min_set_size 5 \
--n_perms 100000 \
--prefix "biased" \
--gene_def 0 \
--threads 6
