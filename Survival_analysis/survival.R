# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("RTCGA")
# BiocManager::install("RTCGA.clinical")
# BiocManager::install("RTCGA.mRNA")
# BiocManager::install("RTCGA.methylation")

# TCGA paclages
library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.mRNA)
library(RTCGA.methylation)

# Survival packages
library(survival)
library(survminer)
# Get your interesting dataset infor
infoTCGA()

# Section 1: Prepare data includes clinical information, methylation (Meth27), gene expression (Array)
# Clinical infor
?clinical
dim(BRCA.clinical)
head(BRCA.clinical[1:5])# too many infor just export to visual the data
# write.table(BRCA.clinical,"clinical.txt",sep="\t")
# Clinical gene expression by array
head(BRCA.methylation[1:5])

# Extract columns in clinical infors with the bar code of patients
clin <- survivalTCGA(BRCA.clinical,
                     extract.cols=c("admin.disease_code", "patient.gender","patient.stage_event.pathologic_stage","patient.stage_event.tnm_categories.pathologic_categories.pathologic_m",
                                    "patient.stage_event.tnm_categories.pathologic_categories.pathologic_n","patient.stage_event.tnm_categories.pathologic_categories.pathologic_t",
                                    "patient.age_at_initial_pathologic_diagnosis"))
head(clin)
dim(clin)

# Extract columns in gene expression of mRNA and the barcode of sample
BRCA.ex=BRCA.mRNA %>% 
  # then make it a tibble (nice printing while debugging)
  as_tibble() %>% 
  # then get just a few genes
  select(bcr_patient_barcode, PAX8, GATA3, ESR1) %>% 
  # then trim the barcode (see head(clin), and ?substr)
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% 
  # then join back to clinical data
  inner_join(clin, by="bcr_patient_barcode")
# Extract columns bar code and the cpg site
BRCA.meth=BRCA.methylation %>% 
  # then make it a tibble (nice printing while debugging)
  as_tibble() %>% 
  # then get just a few cpg
  select(bcr_patient_barcode, cg27653134, cg27654142) %>% 
  # then trim the barcode (see head(clin), and ?substr)
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12))  

head(clin)
head(BRCA.meth)
head(BRCA.ex)

BRCA.all=inner_join(x=BRCA.meth,y=BRCA.ex,by="bcr_patient_barcode")

head(BRCA.all)
# bcr_patient_barcode cg27653134 cg27654142   PAX8  GATA3   ESR1 times patient.vital_s… admin.disease_c… patient.gender
# <chr>                    <dbl>      <dbl>  <dbl>  <dbl>  <dbl> <dbl>            <dbl> <chr>            <chr>         
#   1 TCGA-A1-A0SD            0.691      0.0194 -0.542  2.87   3.08    437                0 brca             female        
# 2 TCGA-A2-A04N            0.420      0.0187  0.340  2.03   4.12   3153                0 brca             female        
# 3 TCGA-A2-A04P            0.196      0.0438 -1.35  -0.293  0.473   548                1 brca             female        
# 4 TCGA-A2-A04Q            0.575      0.0622  0.533 -0.803 -1.80   2179                0 brca             female        
# 5 TCGA-A2-A04T            0.152      0.0406  0.465 -2.32  -3.24   1950                0 brca             female        
# 6 TCGA-A2-A04U            0.0639     0.0125 -1.02  -3.76  -4.21    671                0 brca             female     
dim(BRCA.all)

# Section 2: Analysis
# Survival function
f <- Surv(BRCA.all$times, BRCA.all$patient.vital_status)
summary(f)
# Survival vs event
survfit(f~1) #Nothing
result1=summary(survfit(f~1))
head(result1)
result1
str(result1) # There is many factor inside can check later
# Check the hazard ratio h/h0=e^bx1
summary(survfit(Surv(times, patient.vital_status)~1, data=BRCA.all), time=seq(0, 1000, 100))
323/324
# Check for the standard error # standard deviation
sqrt(323/324*(1-323/324)/324)
# [1] 0.004376336
# CI=mean+/-z*SE
# 0.9956+0.004376336*1.96
# 0.9956-0.004376336*1.96


# Function vs one catagorical features (Sex)
# f1 <- survfit(Surv(times, patient.vital_status)~patient.gender, data=BRCA.all) #consider the bias
# f1

# Function vs one catagorical features (patient.age_at_initial_pathologic_diagnosis)
f2 <- survfit(Surv(times/365, patient.vital_status)~patient.age_at_initial_pathologic_diagnosis<65, data=BRCA.all) #consider the bias
f2
summary(f2) # it will return by 2 dataframe for each group
# KMC plot basic
plot(f2)
# install.packages("survminer")
library(survminer)

ggsurvplot(f2)

ggsurvplot(f2, conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("larger65", "less65"), legend.title="Sex",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Kaplan-Meier Curve for TCGA BRCA", 
           risk.table.height=.15)

ggsurvplot(f2, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("larger65", "less65"), legend.title="Sex",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Kaplan-Meier Curve for TCGA BRCA", 
           risk.table.height=.15)
# should not do f3
# rename to make it easier
# names(BRCA.all)=c("code","cg27653134","cg27654142","PAX8","GATA3","ESR1","times","status","code_disease","gender","patho_stage",
#                        "m_stage","n_stage","t_stage","age")

f3 <- survfit(Surv(times/365,status)~n_stage, data=BRCA.all) 
ggsurvplot(f3, conf.int=TRUE, pval=TRUE, risk.table=TRUE, legend.title="Sex",  
           title="Kaplan-Meier Curve for TCGA BRCA", 
           risk.table.height=.15)

# Univariate
summary(coxph(Surv(times/365,status)~n_stage, data=BRCA.all)) 
# Multivariate
summary(coxph(Surv(times/365,status)~cg27654142+cg27653134+GATA3+ESR1, data=BRCA.all)) 
# Combined with catagorical data
head(lung)
# What happened if we have catagorical data
fit <- coxph(Surv(time, status)~sex+agecat+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss, data=lung)
# Fit a Cox proportional hazards model
ggforest(fit, data = lung)
coxph(formula=Surv(time,Risk) ~ dm, data=df3)

