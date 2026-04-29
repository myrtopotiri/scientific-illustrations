library(Hmisc)       # rcorr()
library(reshape2)    # melt()
library(ggplot2)

# 1) Input tables (placeholders)
# AS_Chan_changing_exons: table with significant (|∆PSI|>15 & Padj<0.05) exonic AS events (rows = events, cols = samples PSI)
# GE_SF: gene expression table restricted to splicing factors (rows = SF genes, cols = samples expression)

# --- Ensure row identifiers ---
rownames(AS_Chan_changing_exons) <- AS_Chan_changing_exons$EVENT
rownames(GE_SF)                 <- GE_SF$Symbol

# 2) Build matrices with the SAME sample columns and numeric mode

df1 <- AS_Chan_changing_exons[, 8:136]   # PSI matrix: events x samples
df2 <- GE_SF[, 18:146]                   # Expression matrix: SF x samples

# Optional sanity: align columns (samples) explicitly

df2 <- df2[, colnames(df1)]              # enforce identical sample order
stopifnot(all(colnames(df1) == colnames(df2)))

# donut pie chart
library(webr)
library(ggplot2)

table.m$Var2 = gsub('Core SF', 'Core SF \n (13)', table.m$Var2)
table.m$Var2 = gsub('RBP', 'RNA binding proteins \n (9)', table.m$Var2)
table.m$Var2 = gsub('other', 'Other \n (4)', table.m$Var2)
table.m$Var2 = gsub('Minor', 'Minor (4)', table.m$Var2)
table.m$Var1 = gsub('Minor','Minor Spliceosome', table.m$Var1)

table.m$Var1 = gsub('PRP19 complex', 'PRP19 \n complex\u2025', table.m$Var1)
table.m$Var1 = gsub('tri-snRNP', '\u2025tri-\u2025 \n snRNP', table.m$Var1)
table.m$Var2 = gsub('RNA binding proteins \n (9)', 'RNA binding \n proteins (9)', table.m$Var2)

PieDonut(table.m, aes(Var2,Var1,count=value),addPieLabel = T,title = "",r0 =  0,  
         showRatioThreshold = 0.001, 
         explodeDonut=T,start=3*pi/2,ratioByGroup = T, showRatioDonut=F,
         showRatioPie=F,
         showPieName=F,pieLabelSize = 15, donutLabelSize = 14,  
         labelpositionThreshold=0.5, maxx=2,
         labelposition = 0.7)

PieDonut(table.m, aes(Var2,Var1,count=value), start=3*pi/2)


# 3) Stack events + SF into one matrix, transpose so that samples are rows
df  <- rbind(df1, df2)                   # (events + SF) x samples
M   <- t(df)                             # samples x (events + SF)
mode(M) <- "numeric"

# 4) Pearson correlations with rcorr (returns r, P, n)
rc <- rcorr(M)                           # rc$r, rc$P, rc$n

# 5) Extract event x SF correlation sub-matrices
# rows = events (df1), cols = splicing factors (df2)

cor_r <- as.data.frame(rc$r)
cor_r <- cor_r[, rownames(df2)]
cor_r <- cor_r[rownames(df1), ]
cor_r$EVENT <- rownames(cor_r)
cor_r_m <- melt(cor_r)                   # long format: EVENT, variable(SF), value(r)

cor_p <- as.data.frame(rc$P)
cor_p <- cor_p[, rownames(df2)]
cor_p <- cor_p[rownames(df1), ]
cor_p$EVENT <- rownames(cor_p)
cor_p_m <- melt(cor_p)                   # long format: EVENT, variable(SF), value(p)

cor_n <- as.data.frame(rc$n)
cor_n <- cor_n[, rownames(df2)]
cor_n <- cor_n[rownames(df1), ]
cor_n$EVENT <- rownames(cor_n)
cor_n_m <- melt(cor_n)                   # long format: EVENT, variable(SF), value(n)

# 6) Combine r + p + FDR + n into one table
cor_m <- cor_r_m
cor_m$pvalue <- cor_p_m$value
cor_m$fdr    <- p.adjust(cor_m$pvalue, method = "fdr")  # global FDR across all SF×event tests
cor_m$n      <- cor_n_m$value

# 7) Define significant correlations 
cor_m$DIR <- "Positive"
cor_m$DIR[cor_m$value < 0] <- "Negative"
cor_m$DIR[cor_m$fdr > 0.05] <- "ns"
cor_m$DIR[abs(cor_m$value) < 0.6] <- "ns"

# 8) Keep strong, significant edges (for downstream counts / networks)
cor_sig_pos <- cor_m[cor_m$fdr < 0.05 & cor_m$value >  0.6, ]
cor_sig_neg <- cor_m[cor_m$fdr < 0.05 & cor_m$value < -0.6, ]
cor_sig     <- rbind(cor_sig_pos, cor_sig_neg)

# Example summaries:
# - number of correlated events per splicing factor
sf_degree <- as.data.frame(table(cor_sig$variable))     # variable = splicing factor
sf_degree <- sf_degree[order(sf_degree$Freq, decreasing=TRUE), ]

# - distribution plot (optional)
ggplot(cor_m, aes(x=value, fill=DIR)) +
  geom_density(alpha = 0.8) +
  theme_classic() +
  xlab("Pearson r (PSI_event vs Expr_SF)") +
  ylab("Density")

# final plot with correlation number per SF 
 y=as.data.frame(table(corM_SF_events_sig_neg$variable))
 x=as.data.frame(table(corM_SF_events_sig_pos$variable))
 x$type="Positive"
 y$type="Negative"
 x=rbind(x,y)
 x <- x[order(x$Freq,decreasing = T),]
 x$Class="Core SF"
 x$Class[x$Var1 %in% SF_changing[SF_changing$Class.simple %in% "RBP","NAME"]]="RBP"
 x$Class[x$Var1 %in% SF_changing[SF_changing$Class.simple %in%"other","NAME"]]="other"
 x$Class[x$Var1 %in% SF_changing[SF_changing$Class.simple %in%   
 "Minor","NAME"]]="Minor"

 as.character(z[z$Freq>100,"Var1"])

 ggplot(data=x[x$Var1 %in% as.character(z[z$Freq>100,"Var1"]),] , aes(x= reorder(Var1,   

       -Freq), y=Freq,color = Class,fill=Class)) +
       geom_col(alpha=0.9)+ggtitle("Number of Correlations")+
       xlab("")+facet_grid(type ~ .)+
       theme(legend.position="topright", panel.background = element_rect(fill =    
       "white",  colour = "grey", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ #size =15),
        #axis.title.y = element_text(size = 15),
        #axis.text.y = element_text(size = 15))+
        theme(legend.position="none",panel.background = element_rect(fill = "white",     
        colour = "white",   size = 2, linetype = "solid"))+
        xlab("") + ylab("Total significat pairs")
