library(Gviz)
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(ggseqlogo)
library(grid)

# ===================================================
#gtf subset gene

library(rtracklayer)

# Load full GTF (this may take time)
gtf <- import("Homo_sapiens.GRCh38.115.gtf.gz")

# Extract all unique gene names
unique_genes <- unique(gtf$gene_name)
head(gtf$type)
head(unique_genes)

# Search for a specific gene (e.g., contains "ACOX")
grep("ACAD", unique_genes, value = TRUE)

# Filter for gene ACOX2 and exons only
acadvl_exons <- gtf[gtf$type == "exon" & !is.na(gtf$gene_name) & gtf$gene_name == "ACADVL"]
table(acadvl_exons$exon_number)
unique(acadvl_exons$transcript_id)
table(acadvl_exons$transcript_id)

# Save to new GTF
export(acadvl_exons, "ACADVL_exons.gtf")

#================== best transcript with most exons

library(dplyr)

# Group by transcript and count how many exons each one has
gtf[gtf$gene_name == "ACOX2" & !is.na(gtf$gene_name) & gtf$type == "exon",] %>%
  as.data.frame() %>%
  group_by(transcript_id) %>%
  summarize(n_exons = n()) %>%
  arrange(desc(n_exons))
best_tx <- "ENST00000900712"  # for ACOX2
best_tx <- "ENST00000543245"  # for ACADVL
acadvl_best_exons <- gtf[gtf$type == "exon" &
                           gtf$transcript_id == best_tx]
export(acadvl_best_exons, "acox2_best_exons.gtf")

##################


trim_exons_to_fasta_range <- function(fasta_path, exon_file) {
  # Load FASTA sequence
  seqs <- readDNAStringSet(fasta_path)
  gene_seq <- seqs[[1]]
  seq_length <- nchar(gene_seq)  # FIXED
  
  # Load exon data (GTF or BED)
  exons <- import(exon_file)
  
  # Create a range for the FASTA
  seq_range <- GRanges(seqnames = seqnames(exons)[1], ranges = IRanges(start = 1, end = seq_length))
  
  # Filter exons to those overlapping with the sequence
  trimmed_exons <- subsetByOverlaps(exons, seq_range)
  
  if (length(trimmed_exons) == 0) stop("No exons overlap with the FASTA sequence.")
  
  # Adjust exon coordinates relative to start of sequence
  trimmed_exons_rel <- shift(trimmed_exons, -start(seq_range) + 1)
  
  return(trimmed_exons_rel)
}


exons_trimmed <- trim_exons_to_fasta_range("ACADVL.fasta", "ACADVL_exons.gtf")

#######################


# Import full GTF
gtf <- import("acadvl_best_exons.gtf")

# Subset only exons that fall within known coordinates
acadvl_trimmed_exons <- gtf[start(gtf) >= 7219937 & end(gtf) <= 7225266]

# Save or use in plot
export(acadvl_trimmed_exons, "acadvl_trimmed_exons.gtf")




# Import full GTF
gtf <- import("acox2_best_exons.gtf")

# Subset only exons that fall within known coordinates
acox2_trimmed_exons <- gtf[start(gtf) >= 58505136 & end(gtf) <= 58537190]

# Save or use in plot
export(acox2_trimmed_exons, "acox2_trimmed_exons.gtf")


########################################################

plot_gene_with_multiple_motifs <- function(fasta_path, exon_bed_path, motif_patterns) {
  library(Biostrings)
  library(rtracklayer)
  
  # Load genomic sequence
  seqs <- readDNAStringSet(fasta_path)
  gene_seq <- seqs[[1]]
  seq_length <- nchar(as.character(gene_seq))
  
  # Load exon coordinates from BED or GTF
  exons <- import(exon_bed_path)
  
  # Define the gene region start from exon coordinates
  region_start <- start(range(exons))
  
  # Adjust exon positions to be relative to gene region
  exon_coords <- data.frame(
    start = start(exons) - region_start + 1,
    end = end(exons) - region_start + 1
  )
  
  
  # Convert coordinates to kilobases
  kb_seq_length <- seq_length / 1000
  kb_exon_coords <- exon_coords / 1000
  
  # Plot gene structure
  plot(c(0, kb_seq_length), c(0, 1), type = "n", axes = FALSE,
       xlab = "Genomic Position (kb)", ylab = "", main = "ACADVL Structure with Motifs",
       xlim = c(0, kb_seq_length))
  axis_ticks <- pretty(c(0, kb_seq_length))
  axis(1, at = axis_ticks, labels = paste0(axis_ticks))
  lines(c(0, kb_seq_length), c(0.5, 0.5), lwd = 2)
  
  for (i in 1:nrow(kb_exon_coords)) {
    rect(kb_exon_coords$start[i], 0.4, kb_exon_coords$end[i], 0.6, col = "black")
  }
  
  # Split motifs and assign colors
  motifs <- unlist(strsplit(motif_patterns, "\\|"))
  colors <- rainbow(length(motifs))
  
  # Plot motifs
  for (i in seq_along(motifs)) {
    motif <- motifs[i]
    matches <- matchPattern(motif, gene_seq)
    
    if (length(matches) == 0) {
      warning(paste("No matches for motif:", motif))
    } else {
      motif_positions <- start(matches)
      
      # Adjust motif positions relative to region start
      relative_motif_positions <- motif_positions - 1 + 1  # already relative to FASTA
      kb_motif_positions <- relative_motif_positions / 1000
      
      y_pos <- ifelse(i %% 2 == 1, 0.75, 0.25)
      arrows(
        x0 = kb_motif_positions,
        y0 = rep(y_pos, length(kb_motif_positions)),
        x1 = kb_motif_positions,
        y1 = rep(0.5, length(kb_motif_positions)),
        col = colors[i], lwd = 2, length = 0.2
      )
    }
  }
  
  # Add legend beside the plot
  par(xpd = TRUE)
  legend("topright", legend = motifs, col = colors, pch = 15, cex = 1.2, inset = c(0.67, -0.20), bty = "n",
         x.intersp = 0.5, y.intersp = 0.7)
  par(xpd = FALSE)
}



# === Example Usage ===
plot_gene_with_multiple_motifs(
  fasta_path = "ACADVL.fasta",
  exon_bed_path = "acadvl_trimmed_exons.gtf",
  motif_patterns = "TGGTGG|GGTTGG|TGTTG|TTGGA"  # Example motif  # Multiple motifs separated by '|'
)

plot_gene_with_multiple_motifs(
  fasta_path = "ACOX2.fa",
  exon_bed_path = "acox2_trimmed_exons.gtf",
  motif_patterns = "TGGTGG|GTTGGT|GTGTTG|GTTGGA"
)

