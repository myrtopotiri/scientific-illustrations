seqdata <- read.delim(".../counts.txt",dec = ",",header = TRUE,sep = "\t", 
           stringsAsFactors = FALSE)
sampleinfo <- read.delim("sample_info.txt")

# rename the table
df <- seqdata

# keep the gene of interest
#SRRM3 - ENSG00000177679
gene <- df[rownames(df2) == "ENSG00000177679",]
gene <- t(gene)
gene=merge(gene,sampleinfo,by = "row.names",all  = F)

#plot the gene 
gene$Tumor.Grade[grep("1", gene$Tumor.Grade)] = "G1"
gene$Tumor.Grade[grep("2", gene$Tumor.Grade)] = "G2"
gene$Tumor.Grade[grep("3", gene$Tumor.Grade)] = "G3"
 
#keep only tumor samples and filter
gene1 <- gene[gene$Group == "tumor",]
gene1 <- gene1[-which(gene1$pStage == ""), ]
gene1 <- gene1[,c(1,2,3,14)]
library(reshape2)
gene1=melt(gene1)

# final plot
library(viridis)
library(ggplot2)
library(ggpubr)
ggboxplot(gene1[gene1$variable %in% "ENSG00000177679",], x = "pStage", y = "value", 
          fill ="pStage",add = "jitter", #width = 0.9,add.params = list(size = 4), 
          size = 1, order = c("IA", "IIA", "IB", "IIB", "III", "IV"))+
          scale_fill_manual(values=c(  
          "#acc9ab","#71b08c","#63b3ba","#6396ba","#ba6380","#ba6363"))+
          theme_pubr()+# Add horizontal line at base  + 
          ylim(0,8000)+
          labs(title="SRRM3 based on Tumor Stage")+ scale_x_discrete(name ="")+  
          ylab(expression(CPM~values))+#ylab(expression(Log[2]*FoldChange~values))+
          theme(legend.position="none",
                axis.title=element_text(size=18),
                plot.title = element_text(size = 18),
                axis.text.x = element_text(size = 18))+
          stat_compare_means(ref.group = "IB",label = "p.signif", size =6, label.y = 
          c(7500))
