setwd("/media/hkh/8TB/XUANTHANG/TCGA/METABRIC")

# Data generated from microarray

#1. Load METABRIC data and clinical
#Download 
metabric_clinical <- read.delim(file = "data_clinical_sample.txt", header = T)
metabric_data <- read.delim(file ="data_expression_median.txt", header = T )
#remove duplicate gene symbol
metabric_data <- metabric_data[!duplicated(metabric_data$Hugo_Symbol),]
#Add gene Symbol as rownames of metabric data, remove outliner data and change patient names
rownames(metabric_data) <- metabric_data$Hugo_Symbol
metabric_data <- metabric_data[,-c(1:2)]
colnames(metabric_data) <- gsub("\\.", "-", colnames(metabric_data))  #or [.]
# Convert patient as column and gene names as column
dataNorm <- tibble::rownames_to_column(as.data.frame(t(metabric_data)), "patient_barcode") 
head(dataNorm)[1:4,1:4]
# Remove other clinical column data that no need for analysis
colnames(metabric_clinical)
data_clinical <- metabric_clinical[,c("PATIENT_ID", "ER_STATUS", "HER2_STATUS", "PR_STATUS", "GRADE", "TUMOR_STAGE")]
names(data_clinical) <- c("patient_barcode", "ER_status_by_IHC", "HER2_status_by_IHC", "PR_status_by_IHC", "Grade", "Tumor_Stage")

# Remove outline clinical data if needed
# intersect(data_clinical$patient_barcode, dataNorm$patient_barcode)
# data_clinical <- subset(data_clinical, data_clinical$patient_barcode %in% dataNorm$patient_barcode)

# Merge normalize data and clinical data by dplyr tools
dataNorm <- left_join(dataNorm, data_clinical, by= "patient_barcode")

saveRDS(dataNorm, file = "METABRIC_data_expression_median.rds")

data <- readRDS(file = "METABRIC_data_expression_median.rds")
# Extract data from INO80 complex
# INO80.complex: "INO80", "INO80B", "INO80C", "INO80D", "INO80E", "MCRS1", "NFRKB", "RUVBL1", "RUVBL2", "ACTL6A", "ACTR5", "ACTP8", "TFPT", "UCHL5", "YY1")

INO80 <- dataNorm[, c("INO80", "ER_status_by_IHC", "HER2_status_by_IHC", "PR_status_by_IHC", "Grade", "Tumor_Stage")]

RUVBL1 <- dataNorm[, c("RUVBL1", "ER_status_by_IHC", "HER2_status_by_IHC", "PR_status_by_IHC", "Grade", "Tumor_Stage")]
RUVBL2 <- dataNorm[, c("RUVBL2", "ER_status_by_IHC", "HER2_status_by_IHC", "PR_status_by_IHC", "Grade", "Tumor_Stage")]
MCRS1 <- dataNorm[, c("MCRS1", "ER_status_by_IHC", "HER2_status_by_IHC", "PR_status_by_IHC", "Grade", "Tumor_Stage")]
YY1 <- dataNorm[, c("YY1", "ER_status_by_IHC", "HER2_status_by_IHC", "PR_status_by_IHC", "Grade", "Tumor_Stage")]

#Test expression of INO80 complex
compare_means(YY1 ~ ER_status_by_IHC, data = YY1,
              method = "t.test")
ggplot(YY1, aes(x=ER_status_by_IHC, y=YY1, fill=ER_status_by_IHC)) + 
  geom_boxplot() 



# Visualize data by ggplot package
library(ggpubr)
library(ggplot2)
library(dplyr)
library(forcats)
library(EnvStats)


# INO80 expression in ER IHC data
compare_means(INO80 ~ ER_status_by_IHC, data = INO80,
              method = "t.test")
ggplot(INO80, aes(x=ER_status_by_IHC, y=INO80, fill=ER_status_by_IHC)) + 
  geom_boxplot() 

my_comparisons <- list( c("Positive", "Negative")) # Visualize: Specify the comparisons you want
INO80 %>% dplyr::filter(!is.na(ER_status_by_IHC) & !is.na(INO80)) %>%
  group_by(ER_status_by_IHC) %>%
  ggplot(aes(x=ER_status_by_IHC, y=INO80, fill=ER_status_by_IHC)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01) +
  stat_n_text() + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + #defaults method = "wilcox.test"
  stat_compare_means(label.y = 9) 
# coord_cartesian(ylim = c(2.5,4)) 


# INO80 expression in IHC data
compare_means(INO80 ~ PR_status_by_IHC, data = INO80,
              method = "t.test")
ggplot(INO80, aes(x=PR_status_by_IHC, y=INO80, fill=PR_status_by_IHC)) + 
  geom_boxplot() 

my_comparisons <- list( c("Positive", "Negative")) # Visualize: Specify the comparisons you want
INO80 %>% dplyr::filter(!is.na(PR_status_by_IHC) & !is.na(INO80)) %>%
  group_by(PR_status_by_IHC) %>%
  ggplot(aes(x=PR_status_by_IHC, y=INO80, fill=PR_status_by_IHC)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01) +
  stat_n_text() + 
  # stat_compare_means(method = "t.test") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + #defaults method = "wilcox.test"
  stat_compare_means(label.y = 9) 
# coord_cartesian(ylim = c(2.5,4)) 



compare_means(INO80 ~ HER2_status_by_IHC, data = INO80,
              method = "t.test")
my_comparisons <- list( c("Positive", "Negative")) # Visualize: Specify the comparisons you want
INO80 %>% dplyr::filter(!is.na(HER2_status_by_IHC) & !is.na(INO80)) %>%
  group_by(HER2_status_by_IHC) %>%
  ggplot(aes(x=HER2_status_by_IHC, y=INO80, fill=HER2_status_by_IHC)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01) +
  stat_n_text() + 
  # stat_compare_means(method = "t.test") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + 
  stat_compare_means(label.y = 9) 
# coord_cartesian(ylim = c(2.5,4)) 



# INO80 expression in a correlation with Tumor grade/Stage
compare_means(INO80 ~ Tumor_Stage, data = INO80,
              method = "anova")
ggplot(INO80, aes(x=as.character(Grade), y=INO80, fill=Grade)) + 
  geom_boxplot() 
my_comparisons <- list(c("0", "1"), c("1", "2"), c("1", "3"), c("1", "4")) # Visualize: Specify the comparisons you want
INO80 %>% dplyr::filter(!is.na(Tumor_Stage) & !is.na(INO80)) %>%
  group_by(Tumor_Stage) %>%
  ggplot(aes(x=as.character(Tumor_Stage), y=INO80, fill=Tumor_Stage)) + 
  geom_boxplot() +
  stat_n_text() + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + #defaults method = "wilcox.test"
  stat_compare_means(label.y = 9) 


