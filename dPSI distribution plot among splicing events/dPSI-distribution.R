# read the table
total <- read_excel("./dPSI.xlsx", col_names = T)

# filter 
total$POS_NEG="NS"
total$POS_NEG[total$Islets_Tumor_p.adj<0.05 & total$dPSI_T_vs_Ilets>15]="pos_sig_changing"
total$POS_NEG[total$Islets_Tumor_p.adj<0.05 & total$dPSI_T_vs_Ilets<(-15)]  
     ="neg_sig_changing"
table(total$POS_NEG)

# calculate the number of the events
total1 <- total[,c("GROUP", "POS_NEG")]
mdfr <- melt(total1, id.vars = "GROUP")
mdfr <- mdfr[!(mdfr$value=="NS"),]
table(mdfr$GROUP, mdfr$value)

mdfr1 <- data.frame(neg = c(47,37,200,43, 31), pos = c(32,64,325,2439, 148))
row.names(mdfr1) <- c("Alt3", "Alt5", "exons", "introns", "MIC")

mdfr1$GROUP <- row.names(mdfr1)

# percent value
mdfr1$neg_percent = (mdfr1$neg / (mdfr1$neg + mdfr1$pos))
mdfr1$pos_percent = (mdfr1$pos / (mdfr1$neg + mdfr1$pos)) 

# round the numbers
mdfr1$neg_percent <- round(mdfr1$neg_percent, digits = 2)
mdfr1$pos_percent <- round(mdfr1$pos_percent, digits = 2)

# make the final plot
library(ggplot2)
library(ggstance)

# Provided dataset
data <- data.frame(Neg = c(19, 36), Pos = c(81, 64))
row.names(data) <- c("MIC", "Exons")

# Reshape the data for ggplot2
library(tidyr)
data_long <- gather(data, key = "response", value = "value")

# Custom horizontal stacked bar plot centered at 0
ggplot(data_long, aes(x = reorder(row.names(data_long), -value), y = value, fill = 
  response)) +
  geom_bar(stat = "identity", position = "identity", color = "white") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "Custom Stacked Bar Plot", x = "Categories", y = "Percentage") +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) 
