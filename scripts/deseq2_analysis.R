library(DESeq2)

library("RColorBrewer")
library("ggplot2")
count_table = read.csv("data/Genecounts_final_table.txt.gz",header=TRUE,row.names=1)
count_table <- aggregate(count_table, FUN = sum, by = list(gsub(row.names(count_table), pattern = "_.*", replacement = "")))
row.names(count_table) <- paste0(count_table$Group.1)
count_table <- count_table[, -1]

time.points <- gsub(pattern = ".*_", replacement = "", names(count_table))
support <- data.frame(time.point = time.points, stringsAsFactors = TRUE)


support$D15.D150 <- factor(ifelse(support$time.point %in% c("D15", "D150"), support$time.point, NA))
support$D0.D15 <- factor(ifelse(support$time.point %in% c("D0", "D15"), support$time.point, NA))

dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData = support,
                              design = as.formula("~ time.point"))

good.rows <- apply(counts(dds), MAR = 1, FUN = max) > 20 ### here I require at least one sample to have at least a read count of 20


message("Number of genes being considered: ", nrow(dds))

dds <- estimateSizeFactors(dds)
sizeFactors <- sizeFactors(dds)
norm_counts <- counts(dds, normalized=TRUE)

################# now run the differential expression analysis
formula0 <- as.formula("~ 1")


good.data <- !is.na(support$D15.D150)
dds.D15.D150 <- DESeqDataSetFromMatrix(countData = count_table[good.rows, good.data], colData = support[ good.data,], design = as.formula("~ D15.D150"))
dds.D15.D150 <- DESeq(dds.D15.D150, test = "LRT", reduced = formula0, minReplicatesForReplace = 7 )
dds.D15.D150.results <- results(dds.D15.D150)
dds.D15.D150.results <- dds.D15.D150.results [ order(dds.D15.D150.results$pvalue),]
dds.D15.D150.results$contig <- row.names(dds.D15.D150.results) 
write.csv(x = dds.D15.D150.results, file = "data/diff_expression/D15_D150.csv", row.names = FALSE)


good.data <- !is.na(support$D0.D15)
dds.D0.D15 <- DESeqDataSetFromMatrix(countData = count_table[good.rows, good.data], colData = support[ good.data,], design = as.formula("~ D0.D15"))
dds.D0.D15 <- DESeq(dds.D0.D15, test = "LRT", reduced = formula0, minReplicatesForReplace = 7 )
dds.D0.D15.results <- results(dds.D0.D15)
dds.D0.D15.results <- dds.D0.D15.results [ order(dds.D0.D15.results$pvalue),]
dds.D0.D15.results$contig <- row.names(dds.D0.D15.results)
write.csv(x = dds.D0.D15.results, file = "data/diff_expression/D0_D15.csv", row.names = FALSE)





stop()



##################### plots happening below
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




##############################

count_tableDesign = data.frame(row.names = colnames( count_table ),
  condition = c("D0","D150","D15","D0","D150","D15","D0","D150","D15","D0","D150","D15","D0","D150","D15"),
  libType = c('s_310', 's_310', 's_310', 's_312', 's_312', 's_312', 's_32', 's_32', 's_32', 's_33', 's_33', 's_33', 's_34', 's_34', 's_34'))
                                 
                                 
#cds1 = newCountDataSet(count_table,colnames(count_table))

cds1 = newCountDataSet(count_table, cddd$libType)
cds1 = estimateSizeFactors( cds1 )
cds1 = estimateDispersions( cds1,method="blind")
plotDispEsts( cds1 )


vsdFull1 = varianceStabilizingTransformation( cds1)
print(plotPCA(vsdFull1))



range(h[which(h[,14]>quantile(probs=0.9,h[,14])),14])


##comp1338396_c0,comp16090_c0,comp69666_c0,comp587243_c0,comp31104_c0,comp2147_c0,comp66837_c0



write.csv(x = norm_counts[which((rownames(norm_counts) %in% c('comp1338396_c0', 'comp16090_c0', 'comp69666_c0', 'comp587243_c0', 'comp31104_c0', 'comp2147_c0', 'comp66837_c0'))),],
          file = "Contig_counts.csv")

norm_count_ordered= norm_counts[,order(c("D0","D150","D15","D0","D150","D15","D0","D150","D15","D0","D150","D15","D0","D150","D15"))]

write.csv(norm_count_ordered[which((rownames(norm_count_ordered) %in% c('comp1338396_c0', 'comp16090_c0', 'comp69666_c0', 'comp587243_c0', 'comp31104_c0', 'comp2147_c0', 'comp66837_c0'))),],"Contig_counts.csv")

  
