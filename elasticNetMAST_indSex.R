#!/bin/bash/env Rscript
rm(list=ls())
gc()

# Build TRS model for single iteration

#############################
#	Command line arguments 	#
#############################
args <- commandArgs(trailingOnly = TRUE)
cluster <- args[1]
sex <- args[2]
seed <- args[3]
fold <- args[4]
out_file <- args[5]
print(args)

#####################
#	Load libraries 	#
#####################
library(MAST)
library(Matrix)
library(plyr)
library(dplyr)
library(tibble)
library(lme4)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(glmnet)
library(gdata)
library(caret)
library(e1071)
library(muscat)
library(Hmisc)
# library(pROC)
# library(plotROC)

# Got an error for some cell clusters. See https://stackoverflow.com/questions/32826906/how-to-solve-protection-stack-overflow-issue-in-r-studio
options(expressions = 5e5)
options(mc.cores = detectCores()- 1)
print(packageVersion("MAST"))

#################
#	Load data 	#
#################
md <- read.csv("/.mounts/labs/awadallalab/scratch/ebader/CPTP_scRNAseq/ohs_phenotypes.csv", stringsAsFactors=FALSE, header = TRUE)
exp <- readRDS(paste0("/.mounts/labs/awadallalab/scratch/ebader/CPTP_scRNAseq/secondary_analysis/TRSnoGenetics/perInd/cluster",cluster,"_avgExpMat.Rds"))

# Data split - assigns each individual as training or test
sample_ids <- read.table(paste0("/.mounts/labs/awadallalab/scratch/ebader/CPTP_scRNAseq/secodary_analysis/TRSnoGenetics/perInd/dataSplits/sex/cluster",cluster,"/cluster",cluster,"dataSplit_",sex,"_seed",seed,"_fold",fold,".tsv"), sep = "\t", header = TRUE)


#````````````````````````````````````````````` DGE TESTING WITH MAST ```````````````````````````````````````````````````
	
#############################
#	Filter expression data 	# 
#############################
	# Include training samples only
	ids_train <- sample_ids %>% filter(dataset == "train" & scRNA.MISO.ID %in% dimnames(exp)[[2]])
	exp_train <- exp[ ,which(dimnames(exp)[[2]] %in% ids_train$scRNA.MISO.ID)]


#################
#	Covariates 	#
#################
	md$group <- as.factor(ifelse(md$RiskGroup %in% c(1,3), "healthy","unhealthy"))
	dge_cov <- md %>%
			filter(scRNA.MISO.ID %in% ids_train$scRNA.MISO.ID) %>%
	                select(scRNA.MISO.ID, batch,group,RiskGroup,age_visit,Sex)


##################################
# 	Filter lowly expressed genes #
##################################
	# Remove genes expressed in less than 10% of samples from either group
	healthy_ids <- dge_cov %>% filter(group == "healthy")
	s1 <- exp_train[ ,which(dimnames(exp_train)[[2]] %in% healthy_ids$scRNA.MISO.ID)]
	thresh_group1 <- ncol(s1)*0.1
	genes_group1 <- rowSums(s1 > 0) >= thresh_group1
	genes_group1 <- names(genes_group1[which(genes_group1 == TRUE)])

	unhealthy_ids <- dge_cov %>% filter(group == "unhealthy")
	s2 <- exp_train[ ,which(dimnames(exp_train)[[2]] %in% unhealthy_ids$scRNA.MISO.ID)]
	thresh_group2 <- ncol(s2)*0.1
	genes_group2 <- rowSums(s2 > 0) >= thresh_group2
	genes_group2 <- names(genes_group2[which(genes_group2 == TRUE)])

	genes_keep <- Reduce(intersect, list(genes_group1, genes_group2))
	print(length(genes_keep))
	exp_train <- exp_train[genes_keep, ]

################################
## MAST DGE modelling
################################
    # Create feature-level dataframe
    gene_ids <- read.table("/.mounts/labs/awadallalab/scratch/EpiCan/CPTP_scRNASeq/features.tsv", header = FALSE, sep = "\t", col.name = c("ensg_id", "hugo","feature_type"), stringsAsFactors = FALSE)
	gene_ids <- gene_ids %>% tibble::column_to_rownames('ensg_id')
    f = gene_ids[which(rownames(gene_ids) %in% dimnames(exp_train)[[1]]), ]
	
	# Create single-cell assay obj
	dge_cov <- dge_cov %>% arrange(scRNA.MISO.ID)
	exp_train <- log1p(exp_train)
	sca <- MAST::FromMatrix(exp_train, cData = dge_cov, fData = f)
	zlmGroup <- zlm(~ group + age_visit + batch, sca) # No need for random effects here bc there is no repeated measures
	zlmGroup.sum <- summary(zlmGroup, doLRT = "groupunhealthy")
    gene_ids <- gene_ids %>% tibble::rownames_to_column("primerid")
	
	# Combine with HUGO ids
	zlmGroupDt <- zlmGroup.sum$datatable
	zlmGroupDt <- merge(zlmGroupDt, gene_ids[ ,1:2], by = "primerid")

	# Extract hurdle model p-values and logFC coef
	fcHurdle <- merge(zlmGroupDt %>% filter(contrast == "groupunhealthy" & component == "H") %>% select(primerid, hugo, `Pr(>Chisq)`),zlmGroupDt %>% filter(contrast == "groupunhealthy" & component == "logFC") %>% select(primerid, hugo, coef, ci.hi, ci.lo),by = c("primerid","hugo"))

	# Genes that will be used in elastic net should have absolute logFC >= 3*sd of all estimates
	features_test <- fcHurdle %>% filter(abs(coef) >= 3*sd(fcHurdle$coef, na.rm = TRUE)) %>% pull(primerid)
 	print(length(features_test))


#````````````````````````````````````````````` ELASTIC NET REGRESSION ``````````````````````````````````````````````````````

	exp <- exp[features_test, which(dimnames(exp)[[2]] %in% sample_ids$scRNA.MISO.ID)]
	exp <- log1p(exp)

	# Get indices of training observations
		df <- as.data.frame(t(exp))
		df <- df %>% tibble::rownames_to_column("scRNA.MISO.ID")
		df_cov <- md %>%
	                	filter(scRNA.MISO.ID %in% sample_ids$scRNA.MISO.ID) %>%
	                	select(group, age_visit,scRNA.MISO.ID) %>%
						arrange(scRNA.MISO.ID)	# to ensure samples are in same order bc need to specify index for train and test data

		df <- merge(df_cov, df, by = "scRNA.MISO.ID")

	# Save row indices of train data to identify them for fitting elastic net
		train_idx <- which(df$scRNA.MISO.ID %in% ids_train$scRNA.MISO.ID)

#############
# 	Caret 	#
#############
	# Set training parameters
	tc <- trainControl(sampling = "up",
			   index = list(train_idx),
                           summaryFunction = twoClassSummary, 
                           classProbs = TRUE, 
                           savePredictions = TRUE)
	# Grid search to optimize lambda and alpha
	myGrid <- expand.grid(alpha = seq(0,1, length = 11), lambda = 10^seq(-3, 2, length = 50))
    mod <- caret::train(form = group ~ . -scRNA.MISO.ID ,
			    		data = df,
                        method = "glmnet", 
                        preProcess = c("center","scale"),
                        na.action = "na.omit", 
                        trControl = tc, 
                        tuneGrid = myGrid,
			    		metric = "ROC")

# If optimal lambda is 100 (max lambda) re-run with a new grid covering larger lambda values
    if(mod$bestTune$lambda == 100){
	    print("Optimal lambda is 100. Re-fitting the model with larger lambda values")
	    myGrid <- expand.grid(alpha = seq(0,1, length = 11), lambda = 10^seq(2,3, length = 20))
	    mod <- caret::train(form = group ~ . ,
                        	data = df[ ,!names(df) %in% c("Sex")],
                        	method = "glmnet",
                        	preProcess = c("center","scale"),
                        	na.action = "na.omit",
                        	trControl = tc,
                        	tuneGrid = myGrid)
    }


# Save the model
saveRDS(mod, file = paste0("../perInd/results/sex/cluster",cluster,"/elasticNet_",out_file,"_mod.Rds"))

#````````````````````````````````````````````` FIN ```````````````````````````````````````````````
	
