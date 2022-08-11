library(tidyverse)
library(ggplot2)
library(scales)
library(openxlsx)


### inputs ###
geomx_data_file <- "~/NAS_data/szimmerman/WTA_manuscript/WTA_manuscript_code/data/Human_23CPA_TargetCounts.xlsx"
rnascope_data_file <- "~/NAS_data/szimmerman/WTA_manuscript/WTA_manuscript_code/data/RNAscope_spot_counting.txt"
SampleID_column <- "Sample_ID" # column with sample IDs in annot
CellLine_column <- "cell_line" # column with cell line IDs in annot
ROIsize_column <- "ROI_size" # column with ROI sizes in annot
LOQ_level <- 2 # LOQ level (# of standard deviations)
LOQ_floor <- 2 # minimum LOQ
RNAscope_expressed <- 1 # transcripts/cell level in RNAscope for calling gene expressed

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
counts <- read.xlsx(geomx_data_file, sheet = 2, startRow = 1, colNames = TRUE)
annot <- read.xlsx(geomx_data_file, sheet = 1, startRow = 1, colNames = TRUE)

# process geomx data
counts <- counts %>% pivot_longer(!TargetName, names_to = "Sample_ID", values_to = "count") %>% rename("Gene" = "TargetName")
annot <- annot %>% rename("cell_line" := !!CellLine_column, "ROI_size" := !!ROIsize_column, "Sample_ID" := !!SampleID_column) %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID, "-", ".")) %>%
  mutate(ROI_size = as.numeric(ROI_size)) %>%
  mutate(cell_line = toupper(cell_line))

counts <- inner_join(counts, annot, by = c("Sample_ID"))

# filter data - remove ROIs failing seq QC and no template controls
counts <- counts %>% filter(SequencingSaturation >= 50 & RawReads >= 5000 & !(grepl("NTC", Sample_ID)))

## set LOQ ##
counts <- counts %>% mutate(LOQ = NegGeoMean_01 * (NegGeoSD_01^LOQ_level)) %>%
  mutate(LOQ = case_when(LOQ < LOQ_floor ~ LOQ_floor, LOQ >= LOQ_floor ~ LOQ)) %>%
  mutate(DSP_det = case_when(count > LOQ ~ T, count <= LOQ ~ F))


## read in RNAscope data ##
RNAscope <- read.delim("/home/rstudio/NAS_data/szimmerman/WTA_manuscript/WTA_manuscript_code/data/RNAscope_spot_counting.txt", sep = "\t") 

# process RNAscope data
RNAscope <- RNAscope %>%
  mutate(cell_line = toupper(cell_line)) %>%
  filter(Gene %in% unique(counts$Gene) & cell_line %in% unique(counts$cell_line)) 

# calculate mean transcripts per cell
RNAscope_sum <- RNAscope %>% group_by(cell_line, Gene) %>% summarise(mean_spots = mean(as.numeric(spots), na.rm = T), 
                                                                     sd_spots = sd(as.numeric(spots), na.rm = T), 
                                                                     CV = sd_spots/mean_spots, n_cells = n()) 


counts_RNAscope <- inner_join(counts, RNAscope_sum, by = c("Gene", "cell_line"))

# calculate transcripts per AOI estimate
counts_RNAscope <- counts_RNAscope %>% mutate(transcript_num = nuclei*mean_spots)

# call genes expressed
counts_RNAscope <- counts_RNAscope %>% mutate(RNAscope_detected = case_when(mean_spots >= 1 ~ T, TRUE ~ F))

counts_RNAscope <- counts_RNAscope %>% mutate(comp = case_when(DSP_det == T & RNAscope_detected == T ~ "TP",
                                                                DSP_det == T & RNAscope_detected == F ~ "FP",
                                                                DSP_det == F & RNAscope_detected == T ~ "FN",
                                                                DSP_det == F & RNAscope_detected == F ~ "TN"))


## sensitivity by transcripts per AOI ##
breaks_tps <- c(10,50,100,200,500,Inf)
counts_RNAscope <- counts_RNAscope %>% mutate(tps_bin = cut(transcript_num, breaks_tps))

tps_summary <- counts_RNAscope %>% group_by(tps_bin) %>% 
  summarise(TP = sum(comp == "TP"), TN = sum(comp == "TN"), FP = sum(comp == "FP"), FN = sum(comp == "FN"), n = n()) %>% 
  mutate(sens = TP/(FN+TP))

# plot
tps_summary %>% filter(tps_bin != "NA") %>% 
  ggplot(aes(x = tps_bin, y = sens)) + 
  geom_point(size = 2) +
  ylab("sensitivity") + xlab("transcripts per ROI") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_x_discrete(labels=c("(10,50]" = "10-50", "(50,100]" = "50-100", "(100,200]" = "100-200", "(200,500]" = "200-500", "(500,Inf]" = ">500")) + ylim(0,1)


## sensitivity by transcripts per cell ##
breaks_tpc <- c(1,2,5,10,20,Inf)
counts_RNAscope <- counts_RNAscope %>% mutate(tpc_bin = cut(mean_spots, breaks_tpc))
counts_RNAscope <- counts_RNAscope %>% group_by(ROI_size) %>% 
  mutate(mean_nuclei = mean(nuclei), size_nuclei = paste0(ROI_size, " ", "(", round(mean_nuclei,1), ")"))

tpc_summary <- counts_RNAscope %>% group_by(tpc_bin, ROI_size, size_nuclei) %>% 
  summarise(TP = sum(comp == "TP"), TN = sum(comp == "TN"), FP = sum(comp == "FP"), FN = sum(comp == "FN")) %>% 
  mutate(sens = TP/(FN+TP))

# plot
tpc_summary <- tpc_summary %>% ungroup() %>% mutate(size_nuclei = as.factor(size_nuclei)) %>% mutate(size_nuclei = fct_reorder(size_nuclei, ROI_size))

tpc_summary %>% 
  filter(tpc_bin != "NA") %>%
  ggplot(aes(x = tpc_bin, y = sens, color = size_nuclei)) + 
  geom_point(size = 2) +
  ylim(0,1) +
  xlab("number of transcripts per cell") + ylab("sensitivity") +
  labs(color = "ROI size (mean nuclei)") +
  scale_x_discrete(labels=c("(1,2]" = "1-2", "(2,5]" = "2-5", "(5,10]" = "5-10", "(10,20]" = "10-20", "(20,Inf]" = ">20")) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## specificity by transcripts per cell ##
counts_RNAscope_notexpressed <- counts_RNAscope %>% filter(mean_spots < 1)

specificity <- counts_RNAscope_notexpressed %>% group_by(ROI_size, size_nuclei) %>% 
  summarise(TP = sum(comp == "TP"), TN = sum(comp == "TN"), FP = sum(comp == "FP"), FN = sum(comp == "FN")) %>% 
  mutate(FPR = FP/(FP+TN))

# plot
specificity <- specificity %>% ungroup() %>% mutate(size_nuclei = as.factor(size_nuclei)) %>% mutate(size_nuclei = fct_reorder(size_nuclei, ROI_size))

specificity %>% 
  ggplot(aes(x = size_nuclei, y = 1-FPR)) + 
  geom_point(size = 2) +
  ylim(0,1) +
  xlab("AOI size (number of nuclei)") + ylab("specificity") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  
