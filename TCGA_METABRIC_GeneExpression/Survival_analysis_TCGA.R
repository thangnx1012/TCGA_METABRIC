setwd("/media/hkh/8TB/XUANTHANG/TCGA")
# BiocManager::install("RTCGA")
# installTCGA()
# BiocManager::install("RTCGA.clinical")
# BiocManager::install("RTCGA.rnaseq") #or download and install by hand

library(tidyverse)
library(RTCGA.clinical) # survival times
library(RTCGA.rnaseq) # genes' expression
# LIHC Liver Hepatocellular Carcinoma
# BRCA: Breast Cancer Invasive
# lung adenocarcinoma (LUAD) and squamous cell carcinoma (LUSC)
# SKCM: TCGA Skin Cutaneous Melanoma
TCGA.surv <- survivalTCGA(BRCA.clinical, LIHC.clinical, LUAD.clinical, LUSC.clinical, SKCM.clinical) 

TCGA.rnaseq <- expressionsTCGA(BRCA.rnaseq, LIHC.rnaseq, LUAD.rnaseq, LUSC.rnaseq, SKCM.rnaseq, 
                                    extract.cols = c("INO80|54617")) %>% 
                                    rename(cohort = dataset, INO80 = `INO80|54617`) %>%
                                    #filter(substr(bcr_patient_barcode, 14, 15) == "01") %>%  # only cancer samples
                                    mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12))
head(TCGA.surv)
# Joining survival times and INO80 gene’s expression.
TCGA.surv_rnaseq <- TCGA.surv %>%  left_join(TCGA.rnaseq, by = "bcr_patient_barcode")

head(TCGA.surv_rnaseq)
#Check patients have both clinical and expression information
table(TCGA.surv_rnaseq$cohort, useNA = "always")
# Remove NA patient
TCGA.surv_rnaseq <- TCGA.surv_rnaseq %>% filter(!is.na(cohort))
# The complete data used for further analysis is printed below
dim(TCGA.surv_rnaseq)
head(TCGA.surv_rnaseq)

library(survminer)
# Determining the optimal cutpoint for INO80 gene’s expression

TCGA.surv_rnaseq.cut <- surv_cutpoint(TCGA.surv_rnaseq,
                                           time = "times", event = "patient.vital_status",
                                           variables = c("INO80", "cohort"))
summary(TCGA.surv_rnaseq.cut)
#Plot the cutpoint for gene INO80
plot(TCGA.surv_rnaseq.cut, "INO80", palette = "npg")

# Categorize INO80 variable
TCGA.surv_rnaseq.cat <- surv_categorize(TCGA.surv_rnaseq.cut) 
# Fit and visualize Kaplan-Meier estimates of survival curves
# Plot by using RTCGA way
RTCGA::kmTCGA(TCGA.surv_rnaseq.cat, 
              explanatory.names = c("INO80", "cohort"),
              pval = TRUE, conf.int = TRUE, xlim = c(0,3000))
#Plot by using survminer way
library(survival)
fit <- survfit(Surv(times, patient.vital_status) ~ INO80 + cohort,
               data = TCGA.surv_rnaseq.cat)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for TRUE: show, FALSE Not show
  # point estimaes of survival curves.
  xlim = c(0,3000),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 1000,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)



