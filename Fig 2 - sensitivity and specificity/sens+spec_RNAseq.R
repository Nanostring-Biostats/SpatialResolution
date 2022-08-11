library(tidyverse)
library(openxlsx)
library(yardstick)

### inputs ###
geomx_data_file <- "~/NAS_data/szimmerman/WTA_manuscript/WTA_manuscript_code/data/Human_11CPA_CollapsedTargetCounts.xlsx"
rnaseq_data_file <- "~/NAS_data/szimmerman/WTA_manuscript/WTA_manuscript_code/data/Human_updated_RNAseq_ccle_TPM_GM.csv"
organism <- "human" # human or mouse
SampleID_column <- "SegmentDisplayName" # column with sample IDs in annot
CellLine_column <- "cell_line" # column with cell line IDs in annot
ROIsize_column <- "ROI_size" # column with ROI sizes in annot
RNAseq_expressed <- 1 # threshold in TPM for calling gene expressed in RNAseq
LOQ_column <- "LOQ.(Dev.Commercial.Human.WTA)" # column with LOQ in annot
LOQ_level <- 2 # LOQ level (# of standard deviations)
LOQ_floor <- 2 # minimum LOQ
RNAseq_expressed <- 1 # TPM level for calling gene expressed in RNAseq

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
annot <- annot %>% rename("cell_line" := !!CellLine_column, "ROI_size" := !!ROIsize_column, "Sample_ID" := !!SampleID_column, "LOQ" := !!LOQ_column) %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID, " ", ".")) %>%
  mutate(ROI_size = as.numeric(ROI_size))

counts <- inner_join(counts, annot, by = c("Sample_ID"))

# filter data - remove ROIs failing seq QC, glass and missing cell line ROIs
counts <- counts %>% filter(SequencingSaturation >= 50 & RawReads >= 5000 & 
                              !(cell_line == "GLASS" | cell_line == "EML cl 1"))

## set LOQ ##
counts <- counts %>% mutate(LOQ = case_when(LOQ < LOQ_floor ~ LOQ_floor, LOQ >= LOQ_floor ~ LOQ)) %>%
  mutate(DSP_det = case_when(count > LOQ ~ T, count <= LOQ ~ F))

# number of genes expressed per ROI
genes_expressed <- counts %>% 
  group_by(Sample_ID, cell_line, ROI_size) %>% 
  summarise(expressed = sum(DSP_det == T))

genes_expressed %>% 
  ggplot(aes(x = as.factor(ROI_size), y = expressed)) + 
  geom_boxplot() +
  ylim(0,11000) +
  xlab("ROI diameter (um)") + ylab("genes detected per ROI") +
  my_theme + theme(legend.title = element_blank())

## read in RNAseq data ##
if (organism == "human") {
  RNAseq <- read.delim(rnaseq_data_file, sep=',')
  RNAseq <- RNAseq %>% pivot_longer(!Gene, names_to = "cell_line", values_to = "TPM") %>%
    mutate(cell_line = str_replace_all(cell_line, "_[A-Z]*", ""))
  
  RNAseq[RNAseq$cell_line=='NCIH596', ]$cell_line <- 'H596' # rename one cell to match geomx
  RNAseq <- RNAseq %>% filter(cell_line %in% unique(counts$cell_line)) # remove cell lines not in experiment
  counts_RNAseq <- inner_join(counts, RNAseq, by = c("Gene", "cell_line")) # to compare same cell lines
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
}

# call detections
counts_RNAseq <- counts_RNAseq %>% mutate(RNAseq_det = case_when(TPM > RNAseq_expressed ~ T, TRUE ~ F),
                                    comp = case_when(DSP_det == T & RNAseq_det == T ~ "TP",
                                                     DSP_det == T & RNAseq_det == F ~ "FP",
                                                     DSP_det == F & RNAseq_det == T ~ "FN",
                                                     DSP_det == F & RNAseq_det == F ~ "TN"))

## ROC curves ##
options(yardstick.event_first = FALSE)
pred <- counts_RNAseq %>% group_by(ROI_size) %>% 
  roc_curve(truth = as.factor(RNAseq_det), count)

ggplot(pred, aes(x = 1 - specificity, y = sensitivity, color = as.factor(ROI_size))) +
  geom_line() + 
  labs(color = "ROI diameter (um)") +
  my_theme

## calculate sensivity and specificity ##
ss <- counts_RNAseq %>% group_by(ROI_size) %>% 
  summarise(TN = sum(comp == "TN"), TP = sum(comp == "TP"), FN = sum(comp == "FN"), FP = sum(comp == "FP"),
                                                                 sens = TP/(TP+FN), spec = 1-(FP/(FP+TN)))


