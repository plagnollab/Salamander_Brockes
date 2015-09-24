library(DESeq)
library("RColorBrewer")
library("gplots")
count_table = read.csv("Genecounts_final_table.txt",header=TRUE,row.names=1)
cds = newCountDataSet(count_table,c("D0","D150","D15","D0","D150","D15","D0","D150","D15","D0","D150","D15","D0","D150","D15"))
rs = rowSums ( counts ( cds))
use = (rs > quantile(rs, probs=0.4))
cds = cds[ use, ]
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds )
pdf("Axo2Isoforms_plots.pdf")
plotDispEsts( cds )


vsdFull = varianceStabilizingTransformation( cds)
print(plotPCA(vsdFull))
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:100]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6),main="Expression Number genes:100")
heatmap.2(counts(cds)[select,], col = hmcol, trace="none", margin=c(10,6),main="Raw counts Number genes:100")
vsdFull = varianceStabilizingTransformation( cds)

select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:250]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(250)
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6),main="Expression Number genes:250")
heatmap.2(counts(cds)[select,], col = hmcol, trace="none", margin=c(10,6),main="Raw counts Number genes:250")
res = nbinomTest( cds, "D0", "D150" )
plotMA(res,main="D0_D150")
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="D0_D150")
resSig = res[ res$padj < 0.1, ]
write.csv( na.omit(res),file="D0_D150_fullresults.csv")
write.csv( na.omit(resSig),file="D0_D150_padj_results.csv")
res = nbinomTest( cds, "D0", "D15" )
plotMA(res,main="D0_D15")
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="D0_D15")
resSig = res[ res$padj < 0.1, ]
write.csv( na.omit(res),file="D0_D15_fullresults.csv")
write.csv( na.omit(resSig),file="D0_D15_padj_results.csv")
res = nbinomTest( cds, "D150", "D15" )
plotMA(res,main="D150_D15")
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="D150_D15")
resSig = res[ res$padj < 0.1, ]
write.csv( na.omit(res),file="D150_D15_fullresults.csv")
write.csv( na.omit(resSig),file="D150_D15_padj_results.csv")
dev.off()




===============

  count_tableDesign = data.frame(  row.names = colnames( count_table ), condition = c("D0","D150","D15","D0","D150","D15","D0","D150","D15","D0","D150","D15","D0","D150","D15"), libType = c('s_310', 's_310', 's_310', 's_312', 's_312', 's_312', 's_32', 's_32', 's_32', 's_33', 's_33', 's_33', 's_34', 's_34', 's_34'))
  
  
                                 
                                 
                                 
#cds1 = newCountDataSet(count_table,colnames(count_table))

cds1 = newCountDataSet(count_table, cddd$libType)
cds1 = estimateSizeFactors( cds1 )
cds1 = estimateDispersions( cds1,method="blind")
plotDispEsts( cds1 )


vsdFull1 = varianceStabilizingTransformation( cds1)
print(plotPCA(vsdFull1))



range(h[which(h[,14]>quantile(probs=0.9,h[,14])),14])


comp1338396_c0,comp16090_c0,comp69666_c0,comp587243_c0,comp31104_c0,comp2147_c0,comp66837_c0



write.csv(norm_counts[which((rownames(norm_counts) %in% c('comp1338396_c0', 'comp16090_c0', 'comp69666_c0', 'comp587243_c0', 'comp31104_c0', 'comp2147_c0', 'comp66837_c0'))),],"Contig_counts.csv")

norm_count_ordered= norm_counts[,order(c("D0","D150","D15","D0","D150","D15","D0","D150","D15","D0","D150","D15","D0","D150","D15"))]

write.csv(norm_count_ordered[which((rownames(norm_count_ordered) %in% c('comp1338396_c0', 'comp16090_c0', 'comp69666_c0', 'comp587243_c0', 'comp31104_c0', 'comp2147_c0', 'comp66837_c0'))),],"Contig_counts.csv")

  