library(ggplot2)


D0.D15 <- read.csv("data/diff_expression/D0_D15.csv_annotated.csv")
D15.D150 <- read.csv("data/diff_expression/D15_D150.csv_annotated.csv")
#D0.D15 <- read.csv("D0_D15.csv", stringsAsFactors = FALSE)
#D15.D150 <- read.csv("D15_D150.csv", stringsAsFactors = FALSE)


contig.names <- intersect(D0.D15$contig, D15.D150$contig)



table.logFC <- data.frame(contig = contig.names,
                          description = D0.D15$Description.Refseq[ match(contig.names, table = D0.D15$contig) ],
                          gene.name = D0.D15$Associated.Gene.Name[ match(contig.names, table = D0.D15$contig) ],
                          D0.D15 =  D0.D15$log2FoldChange[ match(contig.names, table = D0.D15$contig) ],
                          D15.D150 =  D15.D150$log2FoldChange[ match(contig.names, table = D15.D150$contig) ],
                          D0.D15.pvalue =  D0.D15$pvalue[ match(contig.names, table = D0.D15$contig) ],
                          D15.D150.pvalue =  D15.D150$pvalue[ match(contig.names, table = D15.D150$contig) ])

table.logFC$comb.pvalue <- pmin( table.logFC$D0.D15.pvalue, table.logFC$D15.D150.pvalue)
table.logFC$comb.pvalue.discrete <- table.logFC$comb.pvalue < 10^-5



  
table.logFC$D0.D15 <- pmin(7, table.logFC$D0.D15)
table.logFC$D15.D150 <- pmin(7, table.logFC$D15.D150)
table.logFC$D0.D15 <- pmax(-7, table.logFC$D0.D15)
table.logFC$D15.D150 <- pmax(-7, table.logFC$D15.D150)


table.logFC$gene.name.clean <- ifelse(table.logFC$D15.D150 > -7 & (table.logFC$D15.D150 < -4 | table.logFC$D15.D150 > 4), as.character(table.logFC$gene.name), NA)
print(table(!is.na(table.logFC$gene.name.clean)))
table.logFC$signif <- ifelse(!is.na(table.logFC$comb.pvalue.discrete) & table.logFC$comb.pvalue.discrete, "red", "black")
                                    
g <- ggplot2::ggplot(data = table.logFC,
                     ggplot2::aes(x = table.logFC$D0.D15,
                                  y = table.logFC$D15.D150))
g <- g + ggplot2::geom_point(colour = table.logFC$signif)
g <- g + ggplot2::xlab("log FC, Day 0 to Day 15") + ggplot2::ylab("log FC, Day 15 to Day 150") + ggplot2::geom_text(ggplot2::aes(x = table.logFC$D0.D15,
                                                                                                                                 y = table.logFC$D15.D150,
                                                                                                                                 label = gene.name.clean),
                                                                                                                                 size = 3, vjust = -1.1)
ggplot2::ggsave(g, file = "figs/clustering_salamander.pdf", width = 6, height = 5)


write.csv(x = table.logFC, file = "results/logFC_table.csv", row.names = FALSE)
