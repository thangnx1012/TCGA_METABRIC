setwd("/media/hkh/8TB/XUANTHANG/TCGA/CGDSR")

library(cgdsr)
mycgds <- CGDS("http://www.cbioportal.org/")
studies <- getCancerStudies(mycgds)

mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)

getCancerStudies(mycgds)[,c(1,2)]
getGeneticProfiles(mycgds,'brca_tcga')[,c(1:2)]
getCaseLists(mycgds,'brca_tcga')[,c(1:2)]
getProfileData(mycgds, "INO80", c("brca_tcga_gistic", "brca_tcga_mrna"), "brca_tcga_all")[c(1:4),]
getClinicalData(mycgds, "brca_tcga_all")[c(1:3),]

df <- getProfileData(mycgds, "INO80", c("brca_tcga_gistic", "brca_tcga_mrna"), "brca_tcga_all")
boxplot(df[,2] ~ df[,1], main="NF1 : CNA status vs mRNA expression", xlab="CNA status", ylab="mRNA expression", outpch = NA)
stripchart(df[,2] ~ df[,1], vertical=T, add=T, method="jitter",pch=1,col='red')

plot(mycgds, "brca_tcga", "INO80", c("brca_tcga_gistic", "brca_tcga_mrna"), "brca_tcga_all", skin = 'disc_cont')

df1 = getProfileData(mycgds, c("INO80","ESR1"), "brca_tcga_mrna" ,"brca_tcga_all")
plot(df1, main="INO80 and ESR1 mRNA expression", xlab="INO80 mRNA expression", ylab="ESR1 mRNA expression")
plot(mycgds, "brca_tcga", c("INO80","ESR1"), "brca_tcga_mrna" ,"brca_tcga_all")


df.pri = getProfileData(mycgds, "INO80", "prad_mskcc_mrna_median_Zscores", "prad_mskcc_primary")
df.met = getProfileData(mycgds, "INO80", "prad_mskcc_mrna_median_Zscores", "prad_mskcc_mets")


boxplot(list(t(df.pri),t(df.met)), main="INO80 expression in primary and metastatic tumors", xlab="Tumor type", ylab="INO80 mRNA expression",names=c('primary','metastatic'), outpch = NA)
stripchart(list(t(df.pri),t(df.met)), vertical=T, add=T, method="jitter",pch=1,col='red')
