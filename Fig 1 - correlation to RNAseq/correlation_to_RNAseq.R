library(tidyverse)
library(openxlsx)


### inputs ###
geomx_data_file <- "~/NAS_data/szimmerman/WTA_manuscript/WTA_manuscript_code/data/mouse_CPA_CollapsedTargetCounts.xlsx"
rnaseq_data_file <- "~/NAS_data/szimmerman/WTA_manuscript/WTA_manuscript_code/data/Mouse_Int_RNASeq_NSTG.csv"
organism <- "mouse" # human or mouse
SampleID_column <- "SegmentDisplayName" # column with sample IDs
CellLine_column <- "Cell_line" # column with cell line IDs
ROIsize_column <- "ROI_size" # column with ROI sizes
RNAseq_expressed <- 1 # threshold in TPM for calling gene expressed in RNAseq

# plotting theme
my_theme <- theme_bw() + theme(panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank(),
                                axis.title.x = element_text(size = 8),
                                axis.text.x = element_text(color = "black", size = 8),
                                axis.text.y = element_text(color = "black", size = 8),
                                axis.title.y = element_text(size = 8),
                                legend.text = element_text(size = 8),
                                legend.title = element_text(size = 8))


## read in geomx data ##
counts <- read.xlsx(geomx_data_file, sheet = 3, startRow = 1, colNames = TRUE)
annot <- read.xlsx(geomx_data_file, sheet = 1, startRow = 1, colNames = TRUE)

# process geomx data
counts <- counts %>% pivot_longer(!TargetName, names_to = "Sample_ID", values_to = "count") %>% rename("Gene" = "TargetName")
annot <- annot %>% rename("cell_line" := !!CellLine_column, "ROI_size" := !!ROIsize_column, "Sample_ID" := !!SampleID_column) %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID, " ", ".")) %>%
  mutate(ROI_size = as.numeric(ROI_size))

counts <- inner_join(counts, annot, by = c("Sample_ID"))

# filter data - remove ROIs failing seq QC, glass and missing cell line ROIs
counts <- counts %>% filter(SequencingSaturation >= 50 & RawReads >= 5000 & 
                                            !(cell_line == "GLASS" | cell_line == "EML cl 1"))

## read in RNAseq data ##
if (organism == "human") {
  RNAseq <- read.delim(rnaseq_data_file, sep=',')
  RNAseq <- RNAseq %>% pivot_longer(!Gene, names_to = "cell_line", values_to = "TPM") %>%
    mutate(cell_line = str_replace_all(cell_line, "_[A-Z]*", ""))
  
  RNAseq[RNAseq$cell_line=='NCIH596', ]$cell_line <- 'H596' # rename one cell line to match geomx
  RNAseq <- RNAseq %>% filter(cell_line %in% unique(counts$cell_line)) # remove cell lines not in experiment
  counts_RNAseq <- inner_join(counts, RNAseq, by = c("Gene", "cell_line")) # to compare same cell lines
  counts_RNAseq_all <- inner_join(counts, RNAseq, by = c("Gene"), suffix = c("_WTA", "_RNAseq")) # to compare all cell lines by all
}

if (organism == "mouse") {
  RNAseq <- read.delim(rnaseq_data_file, sep = ",") %>% select(-GeneID) %>%
    pivot_longer(!Gene, names_to = "cell_line", values_to = "TPM")
  # fix cell line names
  counts[counts$cell_line=='Alpha TC1 cl 6',]$cell_line <- 'alphaTC'
  counts[counts$cell_line=='LL/2',]$cell_line <- 'LL2'
  counts[counts$cell_line=='209/MDCT',]$cell_line <- 'MDCT'
  counts[counts$cell_line=='B16-F10',]$cell_line <- 'B16F10'
  counts[counts$cell_line=='J774A.1',]$cell_line <- 'J774A'
  counts[counts$cell_line=='C8-D1A',]$cell_line <- 'C8D1A'
  counts[counts$cell_line=='NE-4C',]$cell_line <- 'NE_4C'
  counts[counts$cell_line=='N1E-115',]$cell_line <- 'N1E115'
  counts[counts$cell_line=='BW5147.3',]$cell_line <- 'BW51473'
  counts[counts$cell_line=='3T3',]$cell_line <- 'X3T3'
  RNAseq <- RNAseq %>% filter(cell_line %in% unique(counts$cell_line)) # remove cell lines not in experiment
  counts_RNAseq <- inner_join(counts, RNAseq, by = c("Gene", "cell_line")) # to compare same cell lines
  counts_RNAseq_all <- inner_join(counts, RNAseq, by = c("Gene"), suffix = c("_WTA", "_RNAseq")) # to compare all cell lines in DSP to all cell lines in RNAseq
}


## correlation to RNAseq (per ROI) ##
cor <- counts_RNAseq %>% group_by(Sample_ID, cell_line, ROI_size) %>% summarise(rho = cor(count, TPM, method = "spearman"))

# plot scatterplot of single ROI
cell_line_to_plot <- "MDCT"  ## MDCT, 200 um shown in figure
size_to_plot <- 200
ggplot(data = filter(counts_RNAseq, cell_line == cell_line_to_plot & ROI_size == size_to_plot), aes(x = log2(TPM), y = log2(count))) + 
  geom_point(alpha = 0.1) + 
  geom_label(data = filter(cor, cell_line == cell_line_to_plot & ROI_size == size_to_plot), aes(x = -Inf, y = Inf, label = paste("R = ", round(rho,2))), hjust = 0, vjust = 1, label.size = 0, inherit.aes = FALSE, size = 3.5) +
  ylab(bquote(log[2]*"(WTA count)")) + xlab(bquote(log[2]*"(RNA-seq TPM)")) +
  ggtitle(paste0(cell_line_to_plot, ", ", size_to_plot, " um")) +
  my_theme + theme(plot.title = element_text(size = 10))


## correlation to RNAseq - all x all ##
cor_all <- counts_RNAseq_all %>% group_by(Sample_ID, ROI_size, cell_line_WTA, cell_line_RNAseq) %>% 
  summarise(cor = cor(count, TPM, method = "spearman"))

cor_all <- cor_all %>% mutate(match = case_when(cell_line_WTA == cell_line_RNAseq ~ TRUE, TRUE ~ FALSE))

cor_all %>% 
  ggplot(aes(x = cell_line_WTA, y = cor, color = match)) + 
  geom_jitter(width = 0.15, size = 1) + 
  scale_color_manual(values = c("grey", "#3182bd"), labels = c("unmatched cell line", "matched cell line")) + 
  xlab("WTA cell line") + ylab("correlation to RNA-seq") + labs(color = "", ) +
  facet_grid(~ ROI_size, labeller = labeller(ROI_size = c("50" = "50 um", "200" = "200 um", "400" = "400 um")), scales = "free_x") + 
  my_theme + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) 

## correlation to RNAseq - expression bins ## 
p <- c(0, 0.25, 0.5, 0.75, 1) # quartiles
TPM_breaks <- RNAseq %>% filter(TPM >= RNAseq_expressed) %>% summarise(q = quantile(TPM, probs = p, names = F))
TPM_breaks <- as.numeric(TPM_breaks$q)

counts_RNAseq_all <- counts_RNAseq_all %>% mutate(TPM_bin = cut(TPM, TPM_breaks, labels = c("Q1", "Q2", "Q3", "Q4"))) 

# RNAseq all x all cor by bin
cor_bybin <- counts_RNAseq_all %>% group_by(Sample_ID, ROI_size, cell_line_WTA, cell_line_RNAseq, TPM_bin) %>% 
  summarise(cor = cor(count, TPM, method = "spearman"))
cor_bybin <- cor_bybin %>% mutate(match = case_when(cell_line_WTA == cell_line_RNAseq ~ TRUE, TRUE ~ FALSE))

quartile_to_plot <- "Q1"
cor_bybin %>% filter(TPM_bin == quartile_to_plot) %>% 
  ggplot(aes(x = cell_line_WTA, y = cor, color = match)) + 
  geom_jitter(width = 0.15, size = 1) + 
  scale_color_manual(values = c("grey", "#3182bd"), labels = c("unmatched cell line", "matched cell line")) + 
  xlab("WTA cell line") + ylab("correlation to RNA-seq") + 
  ggtitle(paste0(quartile_to_plot, " of expressed genes in RNAseq")) +
  labs(color = "") +
  facet_grid(~ ROI_size, labeller = labeller(ROI_size = c("50" = "50 um", "200" = "200 um", "400" = "400 um")), scales = "free_x") + 
  my_theme + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5))
 
