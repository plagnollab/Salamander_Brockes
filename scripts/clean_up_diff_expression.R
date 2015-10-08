D0.D15 <- read.csv("data/diff_expression/D0_D15.csv_annotated.csv")

D15.D150 <- read.csv("data/diff_expression/D15_D150.csv_annotated.csv")


my.annotations <- c("Associated.Gene.Name", "Description.Refseq", "RefSeq.Protein.ID", "Percent.match", "Evalue",
                    "Xenopus.Associated.Gene.Name", "Xenopus.Description", "Xenopus.Percent.match", "Xenopus.Evalue",
                    "Zebrafish.Associated.Gene.Name", "Zebrafish.Description", "Zebrafish.Percent.match", "Zebrafish.Evalue")
                    


combined.table <- D0.D15[, c("contig", "pvalue", "log2FoldChange", "padj", my.annotations)]
names(combined.table)[ 2 ] <- "pvalue.D0.D15"
names(combined.table)[ 3 ] <- "log2FoldChange.D0.D15"
names(combined.table)[ 4 ] <- "FDR.D0.D15"



combined.table$pvalue.D15.D150 <- D15.D150$pvalue[ match(combined.table$contig, table = D15.D150$contig) ]
combined.table$log2FoldChange.D15.D150 <- D15.D150$log2FoldChange[ match(combined.table$contig, table = D15.D150$contig) ]
combined.table$FDR.D15.D150 <- D15.D150$padj[ match(combined.table$contig, table = D15.D150$contig) ]


combined.table.final <- combined.table[, c("contig", "pvalue.D0.D15", "log2FoldChange.D0.D15", "FDR.D0.D15", "pvalue.D15.D150", "log2FoldChange.D15.D150", "FDR.D15.D150", my.annotations)]
save(list = "combined.table.final", file = "data/combined.table.final.RData")

combined.table.final.small <- combined.table.final[ (!is.na(combined.table.final$FDR.D0.D15) & combined.table.final$FDR.D0.D15 < 0.05) | (!is.na(combined.table.final$FDR.D15.D150) & combined.table.final$FDR.D15.D150 < 0.05),]
write.csv(x = combined.table.final.small, row.names = FALSE, file = "data/combined_D0_D15_D150_with_annotations_5pcFDR.csv")
