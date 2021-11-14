#1. Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
setwd("/media/hkh/8TB/XUANTHANG/TCGA")


#2 Load TCGA data
GDCprojects = getGDCprojects()
head(GDCprojects[c("project_id", "name")])
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

query_TCGA <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

GDCdownload(query = query_TCGA)
tcga_data <- GDCprepare(query_TCGA)
colnames(colData(tcga_data))

gistic.query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Copy Number Variation",
                  data.type = "Gene Level Copy Number Scores",              
                  access="open")
GDCdownload(gistic.query)
gistic.data <- GDCprepare(gistic.query)




table(tcga_data@colData$vital_status)
table(tcga_data@colData$tumor_stage)
table(tcga_data@colData$definition)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$gender)
dim(assay(tcga_data))     # gene expression matrices.
head(assay(tcga_data)[,1:5]) # expression of first 6 genes and first 5 samples
head(rowData(tcga_data))    # ensembl id and gene id of the first 6 genes.

# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again
saveRDS(object = tcga_data,
        file = "brca.tcga_data.RDS",
        compress = FALSE)
#The data can be loaded with the following command
tcga_data <- readRDS(file = "brca.tcga_data.RDS")


#3 RNASeq Normalization
limma_pipeline = function(
  tcga_data,
  condition_variable,
  reference_group=NULL){
  
  design_factor = colData(tcga_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(tcga_data),
                samples=colData(tcga_data),
                genes=as.data.frame(rowData(tcga_data)))
  
  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number=30000, sort.by="p")
  
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 1000 most differentially expressed genes
    )
  )
}


limma_res <- limma_pipeline(tcga_data=tcga_data,
                           condition_variable="definition",
                           reference_group="Solid Tissue Normal")

#3.1.1 Visualization
plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  # include a legend for points
  legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
}

res_pca <- plot_PCA(limma_res$voomObj, "definition")

#4. Classification
# 4.1 Train and test paradigm
# Transpose and make it into a matrix object
d_mat <- as.matrix(t(limma_res$voomObj$E))

# As before, we want this to be a factor
d_resp <- as.factor(limma_res$voomObj$targets$definition)
# Divide data into training and testing set

# Set (random-number-generator) seed so that results are consistent between runs
set.seed(42)
train_ids <- createDataPartition(d_resp, p=0.75, list=FALSE)

x_train <- d_mat[train_ids, ]
x_test <- d_mat[-train_ids, ]

y_train <- d_resp[train_ids]
y_test <- d_resp[-train_ids]
# Train model on training dataset using cross-validation
res <- cv.glmnet(x = x_train, y = y_train,
                alpha = 0.5, family = "binomial")
# 4.2 Model Evaluation
# Test/Make prediction on test dataset
y_pred = predict(res, newx=x_test, type="class", s="lambda.min")
confusion_matrix = table(y_pred, y_test)

# Evaluation statistics
print(confusion_matrix)
# Getting genes that contribute for the prediction
res_coef = coef(res, s="lambda.min") # the "coef" function returns a sparse matrix
dim(res_coef)

head(res_coef) # in a sparse matrix the "." represents the value of zero
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef = res_coef[-1]
relevant_genes = names(res_coef) # get names of the (non-zero) variables.
length(relevant_genes) # number of selected genes

# get coefficients with non-zero values
res_coef = res_coef[res_coef[,1] != 0,]
# note how performing this operation changed the type of the variable
head(res_coef)
head(relevant_genes) # few select genes
head(limma_res$voomObj$genes)
relevant_gene_names = limma_res$voomObj$genes[relevant_genes,"external_gene_name"]
head(relevant_gene_names) # few select genes (with readable names now)
#4.4 Hierarchical clustering
# define the color palette for the plot
hmcol = colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)

# perform complete linkage clustering
clust = function(x) hclust(x, method="complete")
# use the inverse of correlation as distance.
dist = function(x) as.dist((1-cor(t(x)))/2)

# Show green color for genes that also show up in DE analysis
colorLimmaGenes = ifelse(
  # Given a vector of boolean values
  (relevant_genes %in% limma_res$topGenes$ensembl_gene_id),
  "green", # if true, return green for that value
  "white" # if false, return white for that value
)

# As you've seen a good looking heatmap involves a lot of parameters
gene_heatmap = heatmap.2(
  t(d_mat[,relevant_genes]),
  scale="row",          # scale the values for each gene (row)
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  col=hmcol,            # define the color map
  labRow=relevant_gene_names, # use gene names instead of ensembl annotation
  RowSideColors=colorLimmaGenes,
  labCol=FALSE,         # Not showing column labels
  ColSideColors=as.character(as.numeric(d_resp)), # Show colors for each response class
  dendrogram="both",    # Show dendrograms for both axis
  hclust = clust,       # Define hierarchical clustering method
  distfun = dist,       # Using correlation coefficient for distance function
  cexRow=.6,            # Resize row labels
  margins=c(1,5)        # Define margin spaces
)

# Extract the hierarchical cluster from heatmap to class "hclust"
hc <- as.hclust(gene_heatmap$rowDendrogram)

# Cut the tree into 2 groups, up-regulated in tumor and up-regulated in control
clusters <- cutree(hc, k=2)
table(clusters)
# selecting just a few columns so that its easier to visualize the table
gprofiler_cols = c("significant","p.value","overlap.size","term.id","term.name")

# make sure the URL uses https
set_base_url("https://biit.cs.ut.ee/gprofiler")

# Group 1, up in tumor
gprofiler(names(clusters[clusters %in% 1]))[, gprofiler_cols]
# Group 2, up in control
gprofiler(names(clusters[clusters %in% 2]))[, gprofiler_cols]







#5 Survival Analysis
# extract clinical data
clinical <- tcga_data@colData
dim(clinical)
# Data only interested in the "Primary solid Tumor" cases for survival
clin_df <- clinical[clinical$definition == "Primary solid Tumor",
                   c("patient",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "gender",
                     "tumor_stage")]
# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clin_df$deceased <- clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df$overall_survival <- ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# show first 10 samples
head(clin_df)

# single gene survival analysis
# let's extract the table of differential expression we got earlier
expr_df <- limma_res$topGenes
# print the first row, to see the gene name, the logFC value and the p-value
print(expr_df[expr_df$external_gene_name=="INO80", ])
# get the ensembl gene id of the first row
gene_id <- expr_df[expr_df$external_gene_name=="INO80", "ensembl_gene_id"]
# also get the common gene name of the first row
gene_name <- expr_df[expr_df$external_gene_name=="INO80", "external_gene_name"]

# visualize the gene expression distribution on the diseased samples (in black)
# versus the healthy samples (in red)
expr_diseased <- d_mat[rownames(clin_df), gene_id]
expr_healthy <- d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_id]

boxplot(expr_diseased, expr_healthy,
        names=c("Diseased", "Healthy"), main="Distribution of gene expression")
# get the expression values for the selected gene
clin_df$gene_value <- d_mat[rownames(clin_df), gene_id]

# find the median value of the gene and print it
median_value <- median(clin_df$gene_value)
print(median_value)

# divide patients in two groups, up and down regulated.
# if the patient expression is greater or equal to them median we put it
# among the "up-regulated", otherwise among the "down-regulated"
clin_df$gene <- ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

# we can fit a survival model, like we did in the previous section
fit <- survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)

# we can extract the survival p-value and print it
pval <- surv_pvalue(fit, data=clin_df)$pval
print(pval)
# and finally, we produce a Kaplan-Meier plot
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T, title=paste(gene_name))
