# read the tables
total <- read.delim(".../dPSI.tab",dec = ".",header = TRUE,sep = "\t",stringsAsFactors = FALSE)
impact <- read.delim(".../PROT_IMPACT-hg38-v2.3.tab",dec = ".",header = TRUE,sep = "\t",stringsAsFactors = FALSE)
colnames(impact)[1] <- "EVENT"

# merge them in one table
joined <- merge(total, impact, all.x = T, all.y = F) 

# filter as wished
joined$dPSI_T_vs_Ilets <- as.numeric(joined$dPSI_T_vs_Ilets)
joined$Islets_Tumor_p.adj <- as.numeric(joined$Islets_Tumor_p.adj)
joined=joined[abs(joined$dPSI_T_vs_Ilets)>15 & joined$Islets_Tumor_p.adj<0.05,] #3930
table(joined$COMPLEX)

# keep only EX events
joined=as.data.frame(joined[joined$COMPLEX %in% c("S","C1","C2","C3","ANN", "MIC"),])
AS_genes <- unique(joined$GENE)

impact=data.frame(joined[!is.na(joined$ONTO),])
table(factor(impact$ONTO))

# change the names of the protein impact
impact$ONTO[grep("NonCoding", impact$ONTO)] <- "Non coding"
impact$ONTO[grep("ORF disruption", impact$ONTO)] <- "ORF disruption"
#impact$ONTO[grep("ORF disruption upon sequence exclusion", impact$ONTO)] <- "ORF disruption upon excl."
#impact$ONTO[grep("ORF disruption when splice site is used", impact$ONTO)] <- "ORF disruption upon incl."
impact$ONTO[grep("Protein isoform when splice site is used", impact$ONTO)] <- "Alternative protein isoforms"
impact$ONTO[grep("Alternative protein isoforms", impact$ONTO)] <- "Alternative protein isoforms"
impact$ONTO[grep("In the CDS, with uncertain impact", impact$ONTO)] <- "Unknown impact"
table(impact$ONTO)

# filter again
impact_data <- data.frame(table(factor(impact$ONTO)))
impact=as.data.frame(impact[impact$COMPLEX %in% c("S","C1","C2","C3","ANN", "MIC"),])
impact$dPSI_T_vs_Ilets <- as.numeric(impact$dPSI_T_vs_Ilets)
impact$Islets_Tumor_p.adj <- as.numeric(impact$Islets_Tumor_p.adj)
impact=impact[abs(impact$dPSI_T_vs_Ilets)>15 & impact$Islets_Tumor_p.adj<0.05,]
table(impact$COMPLEX)

# find the numbers in every category
impact$value <- impact %>% group_by(ONTO) %>%  
  mutate(pos = sum(dPSI_T_vs_Ilets>0),
         neg = sum(dPSI_T_vs_Ilets<0))
unique(impact$value$pos)
unique(impact$value$neg)
table(impact$value$neg)

# prepare the final table
# 3UTR, 5UTR, Alternative isoforms, Non coding, ORF disruption, Uknown impact
impact_data$pos <- c(2, 30, 285, 13, 112, 26)
impact_data$neg <- c(0, 40, 87, 7, 90, 2)

colnames(impact_data) <- c("Impact", "Total", "Positive dPSI", "Negative dPSI")

impact_data <- impact_data[,-2]

mdfr <- melt(impact_data, id.vars = "Impact")

# make the plot
ggplot(mdfr, aes(x =reorder(Impact, value), y = value, fill = variable)) +
  geom_col(aes(fill = variable), width = 0.85)+
  scale_fill_manual(values = c("black", "grey")) +
  scale_y_continuous(expand = c(0, 0))+
  #ggtitle("Predicted Impact on protein")+
  ylab("Exonic Events") +
  xlab("") +
  coord_flip()+
  theme(text = element_text(size = 50), #axis.text.y = element_text(face = "bold"), 
        axis.line.x = element_line(color="black", size = 1),  
        axis.ticks=element_line(size=1),
        axis.line.y = element_line(color="black", size = 1), 
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black"),
        legend.title = element_blank(),
        legend.position = c(.7,.1),
        panel.background = element_rect(fill='white', colour='white'), 
        plot.title = element_text(size=15))
