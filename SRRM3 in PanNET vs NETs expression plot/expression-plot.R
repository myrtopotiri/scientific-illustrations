library(readxl)
library(dplyr)
library(tidyr)   

# File and small constant
file <- "cRPKM SRRM3.xlsx"
eps  <- 0.01

# 1. Read Sheet1 and rename the columns we care about
dat <- read_xlsx(file, sheet = "Sheet1") %>%
  rename(
    group = source_name,           # cancer/normal group
    SRRM3 = `SRRM3 cRPKM`          # numeric cRPKM values
  )

# Quick sanity check
unique(dat$group)
# [1] "I-NET" "lung-NET" "MCC" "NEPC" "ilets" "p-NET" "Digestive Tract" "Lung" "Prostate" "Skin"

# 2. Map each NET cancer to its matching normal tissue
normal_map <- c(
  "I-NET"    = "Digestive Tract",
  "p-NET"    = "ilets",
  "lung-NET" = "Lung",
  "NEPC"     = "Prostate",
  "MCC"      = "Skin"
)

net_groups    <- names(normal_map)          # NET groups
normal_groups <- unname(normal_map)         # normal tissue groups

# 3. Compute mean SRRM3 in each normal tissue
normal_means <- dat %>%
  filter(group %in% normal_groups) %>%      # keep only normal tissues
  group_by(group) %>%
  summarise(
    normal_mean = mean(SRRM3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(normal_group = group)

# 4. Build a table that says, for each NET group, which normal mean to use
map_tbl <- tibble(
  group        = net_groups,                # NET group name
  normal_group = normal_groups              # matching normal group
) %>%
  left_join(normal_means, by = "normal_group")

# 5. Add log2FC per sample (only for NETs)
dat_log2 <- dat %>%
  left_join(map_tbl, by = "group") %>%      # adds normal_mean for NET groups
  mutate(
    is_NET = group %in% net_groups,
    log2FC = ifelse(
      is_NET,
      log2((SRRM3 + eps) / (normal_mean + eps)),
      NA_real_
    )
  )

# 'dat_log2' now contains a 'log2FC' column for each NET sample
# You can inspect just those:
net_samples <- dat_log2 %>% filter(is_NET)

# 6. Optional: average log2FC per cancer type
per_cancer_log2FC <- net_samples %>%
  group_by(group) %>%
  summarise(
    n_samples   = n(),
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    sd_log2FC   = sd(log2FC, na.rm = TRUE),
    .groups = "drop"
  )

######### ====================================================================

FoldChange_data2$Cancer_Type <- gsub("p-NET", "PanNET", FoldChange_data2$Cancer_Type)

# define colors so PanNET is one group, everything else another
FoldChange_data2$Color <- ifelse(
  FoldChange_data2$Cancer_Type == "PanNET",
  "PanNET",
  "Other"
)

ggdotchart(FoldChange_data2, x = "Cancer_Type", y = "logFC",
           color = "Color",                                
           palette = c("Other" = "#00AFBB",   # blue for all others
                       "PanNET" = "#E7B800"),# yellow for PanNET
           sorting = "descending",
           add = "segments",
           add.params = list(color = "lightgray", size = 3),
           dot.size = 15,
           label = round(FoldChange_data2$logFC,1),
           font.label = list(color = "white", size = 17, vjust = 0.5),
           ggtheme = theme_pubr(),
           rotate = FALSE,
           xlab = FALSE
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  theme(panel.grid.minor = element_blank(),
        strip.placement = "outside",
        axis.title = element_text(size=20),
        axis.text = element_text(size=22),
        legend.position = "none",
        plot.title = element_text(size = 25)) +
  ggtitle("SRRM3 Expression in NET") +
  ylim(-2,6) +
  labs(y = expression(Log[2]*" Fold Change"))
