library(readxl)
library(clusterProfiler)
library(enrichplot)
library(gprofiler2)

diff_table <- read_excel("dPSI_withScarpa_all.xlsx", col_names = T)

# filterings
diff_table$Islets_Tumor_p.adj <- as.numeric(diff_table$Islets_Tumor_p.adj)
res_significant <- subset(diff_table, diff_table$Islets_Tumor_p.adj <=.05) #only pvalue adjusted cutoff

res_significant <- res_significant[grep("HsaEX", res_significant$EVENT),]
res_significant <- res_significant[res_significant$COMPLEX %in% c("MIC","S","ANN","C1","C2","C3") & res_significant$LENGTH<50,]

gene <- unique(res_significant$GENE)

library(organism, character.only = TRUE)

ego <- enrichGO(gene          = gene,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP", # CC, MF
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

dotplot(ego, showCategory=10, font.size = 17, label_format = 50) +  
         scale_size_area(max_size =10) + #max size 15 for BP
         theme(
         axis.ticks=element_line(size=1.5),
         axis.title.x = element_text(size = 13), #15 for BP
         #axis.line = element_line(colour = 'black', size = 1),
         legend.key.size = unit(0.8, 'cm'),
         legend.title = element_text(size=13), #15 for BP
         legend.text = element_text(size=13), #15 for BP
         panel.background = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank()) +
         scale_colour_continuous(name  ="Adjusted\np-value", low = "red", high =  
         "blue") +
         scale_shape_discrete(name  ="Count")

# cnetplot
cnt_enrichment <- cnetplot(ego,showCategory=5, color_category='steelblue', 
                  cex_category =5,cex_gene = 3,layout = "kk", color_gene='grey',    
                  node_label = "all", cex_label_category = 3.5,cex_label_gene =2.5)
# extract the data object from the ggplot object
dat <- cnt_enrichment$data
# change the colors
dat$color[dat$name %in% c("CADPS2", "DCTN1", "ERGIC3", "CACNA1D", "SNX2", "PTK2")] <- "orange"
dat$size[dat$name %in% c("CADPS2", "DCTN1", "ERGIC3", "CACNA1D", "SNX2", "PTK2")] <- 10

# put it back
cnt_enrichment$data <- dat
# check it out
options(ggrepel.max.overlaps = Inf)
cnt_enrichment + theme(legend.text = element_text(size = 32), 
        legend.title = element_text(size = 40, ),
        legend.position = c(0.92, 0.8), # FOR BP: c(0.92,0.5)
        )

# Revigo
library(rrvgo)

simMatrix <- calculateSimMatrix(ego$ID,                       #calculating semantic similarity
                                orgdb="org.Hs.eg.db",
                                ont="CC",
                                method="Rel")

scores <- setNames(-log10(ego$qvalue), ego$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

# pick what to visualize in BP revigo
to_keep <- c("vesicle", "endocytosis","exocytosis", "GTP", "synap")
reducedTerms <- reducedTerms[grep(paste(to_keep,collapse="|"), reducedTerms$term),]

ggplot(reducedTerms,
       aes(x = term, y = -log10(size), size = -log10(size), fill = -log10(size))) +
       expand_limits(y = 1) +
       geom_point(shape = 21) +
       scale_size(range = c(2.5,12.5)) +
       scale_fill_continuous(low = 'royalblue', high = 'red4') +
       xlab('') + ylab('Enrichment score') +
       labs(
       title = 'GO Analysis')+
       geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),
             colour = c("black", "black", "black"),
             size = c(0.5, 1.5, 3)) +
       theme_bw(base_size = 24) +
       theme(
       legend.position = 'right',
       legend.background = element_rect(),
       plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
       plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
       plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
       axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
       axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
       axis.title = element_text(size = 12, face = 'bold'),
       axis.title.x = element_text(size = 12, face = 'bold'),
       axis.title.y = element_text(size = 12, face = 'bold'),
       axis.line = element_line(colour = 'black'),
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) + 
    coord_flip()
