setwd("/media/hkh/8TB/XUANTHANG/TCGA/HTseq_hg38_TCGA")
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(stringr)

query.exp.hg38 <- GDCquery(project = "TCGA-BRCA", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM",
                           sample.type = c("Primary Tumor","Solid Tissue Normal"))
GDCdownload(query.exp.hg38)
brca.data.exp <- GDCprepare(query = query.exp.hg38,
                     save = TRUE, 
                     save.filename = "brca.tcga_data.rds")
# brca.data <- readRDS("brca.tcga_data.rds")
# BRCA subtype
colnames(colData(brca.data))
plyr::count(colData(brca.data)$"paper_BRCA_Subtype_PAM50")
plyr::count(colData(brca.data)$"sample_type_id")
plyr::count(colData(brca.data)$"definition")

dataSubtype <- as.data.frame(colData(brca.data)[,c("patient", "paper_BRCA_Subtype_PAM50", "definition")])
dataSubtype <- tibble::rownames_to_column(dataSubtype, "bcr_patient_barcode")
plyr::count(dataSubtype$paper_BRCA_Subtype_PAM50)
plyr::count(dataSubtype$definition)

#Normalization data as gcContent method to adjust for GC-content effect on read counts
dataPrep <- TCGAanalyze_Preprocessing(object = brca.data.exp, 
                                      cor.cut = 0.6)

brca.dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                           geneInfo = geneInfoHT,
                                           method = "gcContent") 

##Transform data to log10 RPKM
dataNormlog <- as.data.frame(log10(brca.dataNorm + 1))
dataNormlog <- tibble::rownames_to_column(dataNormlog, "Ensembl_ID")


#####################################################################################################################################
#OR download from Xena
dataNormlog <- read.delim("TCGA-BRCA.htseq_fpkm.tsv", header = T) #No need to set rownames = 1st column because of GeneID
head(dataNormlog)[1:5, 1:5]

Gistic2_data <- read.delim("TCGA-BRCA.gistic.tsv", header =T)
head(Gistic2_data)[1:5,1:5]

#Rename string patient and gene names
dataNormlog$Ensembl_ID <- substr(dataNormlog$Ensembl_ID, 1, 15)
colnames(dataNormlog) <- substr(colnames(dataNormlog), 1, 15)
dataNormlog <- dataNormlog[!duplicated(colnames(dataNormlog))]
head(dataNormlog)[1:5, 1:5]

colnames(Gistic2_data) <- substr(colnames(Gistic2_data), 1, 15)
Gistic2_data$Gene.Symbol <- substr(Gistic2_data$Gene.Symbol, 1, 15)
head(Gistic2_data)[1:5, 1:5]
#Convert GeneID to symbol
library(EnsDb.Hsapiens.v86)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= dataNormlog$Ensembl_ID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

names(geneIDs1) <- c("SYMBOL","Ensembl_ID")
dataNormlog <- left_join(geneIDs1, dataNormlog, by="Ensembl_ID")

names(geneIDs1) <- c("SYMBOL","Gene.Symbol")
Gistic2_data <- left_join(geneIDs1, Gistic2_data, by="Gene.Symbol")

# Remove duplicate SYMBOL genes and set Ensemble_ID as rownames
dataNormlog <- dataNormlog[!duplicated(dataNormlog$SYMBOL),]
rownames(dataNormlog) <- dataNormlog$SYMBOL
head(dataNormlog)[1:5, 1:5]

Gistic2_data <- Gistic2_data[!duplicated(Gistic2_data$SYMBOL),]
rownames(Gistic2_data) <- Gistic2_data$SYMBOL
head(Gistic2_data)[1:5, 1:5]

# remove unneeded column
dataNormlog <- dataNormlog[,-c(1:2)]
Gistic2_data <- Gistic2_data[,-c(1:2)]

#Convert rownames to colnames and set patient as a row with bcr_patient_bacode
dataNormlog <- tibble::rownames_to_column(as.data.frame(t(dataNormlog)), "bcr_patient_barcode")
dataNormlog$bcr_patient_barcode <- gsub("\\.", "-", dataNormlog$bcr_patient_barcode) #or [.]
head(dataNormlog)[1:5, 1:5]

Gistic2_data <- tibble::rownames_to_column(as.data.frame(t(Gistic2_data)), "bcr_patient_barcode")
Gistic2_data$bcr_patient_barcode <- gsub("\\.", "-", Gistic2_data$bcr_patient_barcode) #or [.]
head(Gistic2_data)[1:5, 1:5]
INO80_gistic2 <- Gistic2_data[,c("bcr_patient_barcode", "INO80")]
names(INO80_gistic2) <- c("bcr_patient_barcode", "INO80_GISTIC2")
plyr::count(INO80_gistic2$INO80)

# BRCA clinical data
brca_subtype <- TCGAquery_subtype(tumor = "BRCA")
brca_subtype <- brca_subtype[,c(1,4:7,12)]
names(brca_subtype) <- c("bcr_patient_barcode","vital_status","days_to_birth",
                         "days_to_death","days_to_last_followup","BRCA_Subtype_PAM50")
brca_subtype$bcr_patient_barcode <- str_c(brca_subtype$bcr_patient_barcode, "-01")
# 1.2 get Clinical data (IHC data)
brca_clinical_query <- GDCquery(project = "TCGA-BRCA", 
                                data.category = "Clinical",
                                data.type = "Clinical Supplement", 
                                data.format = "BCR Biotab")
GDCdownload(brca_clinical_query)
brca_clinical_all <- GDCprepare(brca_clinical_query)
brca_IHC <- brca_clinical_all$clinical_patient_brca[,c("bcr_patient_barcode","er_status_by_ihc", 
                                                         "pr_status_by_ihc", "her2_status_by_ihc")][-c(1:2),]
names(brca_IHC) <- c("bcr_patient_barcode", "ER_status_by_IHC", "PR_status_by_IHC",  "HER2_status_by_IHC")
# Replace other data by NA value
# brca_ER_IHC <- brca_IHC[!brca_IHC$ER_status_by_IHC == "[Not Evaluated]" & !brca_IHC$ER_status_by_IHC == "Indeterminate",]
# brca_IHC$ER_status_by_IHC[brca_IHC$ER_status_by_IHC == "[Not Evaluated]"] <- NA
brca_IHC <- brca_IHC %>%
  mutate(ER_status_by_IHC=replace(ER_status_by_IHC, brca_IHC$ER_status_by_IHC == "[Not Evaluated]", NA)) %>%
  mutate(ER_status_by_IHC=replace(ER_status_by_IHC, brca_IHC$ER_status_by_IHC == "Indeterminate", NA)) %>%
  mutate(PR_status_by_IHC=replace(PR_status_by_IHC, brca_IHC$PR_status_by_IHC == "[Not Evaluated]", NA)) %>%
  mutate(PR_status_by_IHC=replace(PR_status_by_IHC, brca_IHC$PR_status_by_IHC == "Indeterminate", NA)) %>%
  mutate(PR_status_by_IHC=replace(PR_status_by_IHC, brca_IHC$PR_status_by_IHC == "[Not Available]", NA)) %>%
  mutate(PR_status_by_IHC=replace(PR_status_by_IHC, brca_IHC$PR_status_by_IHC == "Equivocal", NA)) %>%
  mutate(HER2_status_by_IHC=replace(HER2_status_by_IHC, brca_IHC$HER2_status_by_IHC == "[Not Evaluated]", NA)) %>%
  mutate(HER2_status_by_IHC=replace(HER2_status_by_IHC, brca_IHC$HER2_status_by_IHC == "Indeterminate", NA)) %>%
  mutate(HER2_status_by_IHC=replace(HER2_status_by_IHC, brca_IHC$HER2_status_by_IHC == "[Not Available]", NA)) %>%
  mutate(HER2_status_by_IHC=replace(HER2_status_by_IHC, brca_IHC$HER2_status_by_IHC == "Equivocal", NA)) %>%
  as.data.frame()

plyr::count(brca_IHC$ER_status_by_IHC)
plyr::count(brca_IHC$PR_status_by_IHC)
plyr::count(brca_IHC$HER2_status_by_IHC)

# Add string "-01" for cancer patient sample
brca_IHC$bcr_patient_barcode <- str_c(brca_IHC$bcr_patient_barcode, "-01")

#Merge data and clinical data based on 
dataNormlog <- left_join(dataNormlog, brca_subtype, by = "bcr_patient_barcode")
dataNormlog <- left_join(dataNormlog, brca_IHC, by = "bcr_patient_barcode")
dataNormlog <- left_join(dataNormlog, INO80_gistic2, by = "bcr_patient_barcode")

# Replace Normal_subtype tumor to Normal-like and NA to normal sample
# dataNormlog$BRCA_Subtype_PAM50[dataNormlog$BRCA_Subtype_PAM50=="Normal"] <- "Normal-like"
dataNormlog <- dataNormlog %>%
  mutate(BRCA_Subtype_PAM50 = replace(BRCA_Subtype_PAM50, BRCA_Subtype_PAM50 == "Normal", "Normal-like"))
plyr::count(dataNormlog$BRCA_Subtype_PAM50)

# dataNormlog$BRCA_Subtype_PAM50[dataNormlog$bcr_patient_barcode=="^*-11",] <- "Normal"
dataNormlog <- dataNormlog %>%
  mutate(BRCA_Subtype_PAM50 = ifelse(substr(dataNormlog$bcr_patient_barcode, nchar(dataNormlog$bcr_patient_barcode)-1, 
                                            nchar(dataNormlog$bcr_patient_barcode)) == "11" , "Normal", BRCA_Subtype_PAM50))
# check data
plyr::count(dataNormlog$BRCA_Subtype_PAM50)


saveRDS(dataNormlog, file = "dataNormlog_TCGA_BRCA.rds")
# "INO80", "INO80B", "INO80C", "INO80D", "INO80E", "MCRS1", "NFRKB", "RUVBL1", "RUVBL2", "ACTL6A", "ACTR5", "ACTP8", "TFPT", "UCHL5", "YY1"

# 1. INO80 expression in IHCdata
library(ggpubr)
library(ggplot2)
library(dplyr)
library(forcats)
library(EnvStats)
# Prepare INO80 Expression data
INO80 <- as.data.frame(dataNormlog[, c("INO80", "BRCA_Subtype_PAM50", "ER_status_by_IHC","PR_status_by_IHC", "HER2_status_by_IHC", "INO80_GISTIC2")])


# INO80 expression in IHC data
plyr::count(INO80$ER_status_by_IHC)
compare_means(INO80 ~ ER_status_by_IHC, data = INO80,
              method = "t.test")
my_comparisons <- list( c("Positive", "Negative")) # Visualize: Specify the comparisons you want
INO80 %>% dplyr::filter(!is.na(ER_status_by_IHC) & !is.na(INO80)) %>%
  group_by(ER_status_by_IHC) %>%
  ggplot(aes(x=ER_status_by_IHC, y=INO80, fill=ER_status_by_IHC)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.02) +
  stat_n_text() + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + #defaults method = "wilcox.test"
  stat_compare_means(label.y = 5) 
# coord_cartesian(ylim = c(2.5,4)) 


plyr::count(INO80$PR_status_by_IHC)
compare_means(INO80 ~ PR_status_by_IHC, data = INO80,
              method = "t.test")
my_comparisons <- list( c("Positive", "Negative")) # Visualize: Specify the comparisons you want
INO80 %>% dplyr::filter(!is.na(PR_status_by_IHC) & !is.na(INO80)) %>%
  group_by(PR_status_by_IHC) %>%
  ggplot(aes(x=PR_status_by_IHC, y=INO80, fill=PR_status_by_IHC)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.02) +
  stat_n_text() + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + #defaults method = "wilcox.test"
  stat_compare_means(label.y = 5) 
# coord_cartesian(ylim = c(2.5,4)) 

plyr::count(INO80$HER2_status_by_IHC)
compare_means(INO80 ~ HER2_status_by_IHC, data = INO80,
              method = "t.test")
my_comparisons <- list( c("Positive", "Negative")) # Visualize: Specify the comparisons you want
INO80 %>% dplyr::filter(!is.na(HER2_status_by_IHC) & !is.na(INO80)) %>%
  group_by(HER2_status_by_IHC) %>%
  ggplot(aes(x=HER2_status_by_IHC, y=INO80, fill=HER2_status_by_IHC)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.02) +
  stat_n_text() + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + 
  stat_compare_means(label.y = 5) 
# coord_cartesian(ylim = c(2.5,4)) 










# 2. INO80 expression in subtype
# Plotting
# Visualize data by ggplot package
library(ggpubr)
library(ggplot2)
library(dplyr)
library(forcats)
library(EnvStats)

head(INO80, 10)
# colored by the variable color
ggplot(INO80, aes(x=INO80, y=BRCA_Subtype_PAM50, fill=BRCA_Subtype_PAM50)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01) +
  stat_n_text()



# Compare Normal and breast cancer subtype expression
compare_means(INO80 ~ BRCA_Subtype_PAM50, data = INO80,
              method = "t.test")
ggplot(INO80, aes(x=BRCA_Subtype_PAM50, y=INO80, fill=BRCA_Subtype_PAM50)) + 
  geom_boxplot() 


my_comparisons <- list( c("Normal", "Basal"), c("Normal", "LumA"), c("Normal", "LumB"), 
                        c("Normal", "Her2"), c("Normal","Normal-like" )) # Visualize: Specify the comparisons you want
# Reorder following a precise order to arrange the column data 
INO80 %>% mutate(BRCA_Subtype_PAM50 = fct_relevel(BRCA_Subtype_PAM50, "Normal", "Basal", "LumA", "LumB", "Her2", "Normal-like")) %>%
  dplyr::filter(!is.na(BRCA_Subtype_PAM50) & !is.na(INO80)) %>%
  ggplot(aes(x=BRCA_Subtype_PAM50, y=INO80, fill=BRCA_Subtype_PAM50)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis = "y", stackdir = "center",
               fill = "white", binwidth = 0.01) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + 
  stat_compare_means(label.y = 6) +
  stat_n_text()
  
