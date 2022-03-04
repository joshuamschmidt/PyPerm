setwd('C:/Users/jschmi06/Documents/Projects/ancestral_chimp_paper/')
library(data.table)

pbsnj <- readRDS('pbsnj.subspecies.pvalues.rds')
anc.windows <- readRDS('threepclr.anc.pvalues.rds')
genes <- fread('../PySetPerm/data/genes.txt')
genes[,start:=start-5000]
genes[,end:=end+5000]
setkey(genes, chr, start, end)
pbsnj[,end:=start]
anc.windows[,end:=start]
setkey(pbsnj, chr, start, end)
setkey(anc.windows, chr, start, end)

# ANC windows cuttofs
# 4090, 818, 409 windows
foverlaps(anc.windows[anc.p <= 0.005], genes, nomatch = 0)[,uniqueN(gene)]
# 817
foverlaps(anc.windows[anc.p <= 0.001], genes, nomatch = 0)[,uniqueN(gene)]
# 200
foverlaps(anc.windows[anc.p <= 0.0005], genes, nomatch = 0)[,uniqueN(gene)]
# 104
foverlaps(pbsnj[e.p <= 0.000228], genes, nomatch = 0)[,uniqueN(gene)]
# 806
foverlaps(pbsnj[e.p <= 0.000035], genes, nomatch = 0)[,uniqueN(gene)]
# 199
foverlaps(pbsnj[e.p <= 0.000015], genes, nomatch = 0)[,uniqueN(gene)]
# 104
foverlaps(pbsnj[c.p <= 0.000192], genes, nomatch = 0)[,uniqueN(gene)]
# 816
foverlaps(pbsnj[c.p <= 0.00004], genes, nomatch = 0)[,uniqueN(gene)]
# 200
foverlaps(pbsnj[c.p <= 0.00002], genes, nomatch = 0)[,uniqueN(gene)]
# 106

cutoffs <- c(0.005,0.001,0.0005, 0.000228, 0.000035, 0.000015, 0.000192, 0.00004, 0.00002)
names(cutoffs) <- c(rep.int("ancestral",times=3), rep.int("eastern",times=3), rep.int("central",times=3) )

cols <- c("anc.p","e.p","c.p")
names(cols) <- c("ancestral","eastern","central")


candidates <- list()

for (i in 1:length(cutoffs)) {
  c <- cutoffs[i]
  lineage <- names(cutoffs)[i]
  col <- cols[lineage]
  
  if(lineage=="ancestral"){
    tmp <- unique(foverlaps(anc.windows[get(col) <= c], genes, nomatch = 0)[,.(chr, start=i.start, end=i.end)])
    candidates[[paste(lineage,c,"candidate.snps.txt",sep="-")]] <- tmp
  } else {
    tmp <- unique(foverlaps(pbsnj[get(col) <= c], genes, nomatch = 0)[,.(chr, start=i.start, end=i.end)])
    candidates[[paste(lineage,c,"candidate.snps.txt",sep="-")]] <- tmp
  }
}
lapply(1:length(candidates), function(x) fwrite(x=candidates[[x]],file = names(candidates)[x], quote = F, sep = "\t",col.names = T, row.names = F))

# backgrounds.
pbs.bg <- unique(foverlaps(pbsnj, genes, nomatch = 0)[,.(chr,start=i.start,end=i.end)])
pclr.bg <- unique(foverlaps(anc.windows, genes, nomatch = 0)[,.(chr,start=i.start,end=i.end)])

fwrite(x=unique(pbs.bg[,.(chr,start,end)]),file = "pbsnj-bg.snps.txt", quote = F, sep = "\t",col.names = T, row.names = F)
fwrite(x=unique(pbs.bg[,.(chr,start,end)]),file = "ancestral-bg.snps.txt", quote = F, sep = "\t",col.names = T, row.names = F)




