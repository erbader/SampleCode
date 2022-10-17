# Healthy aging Fig 1C
# Working locally

rm(list=ls())
gc()

#####################
# 	Load libraries 	#
#####################
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)

#################
# 	Load data 	#
#################
	metadata <- read.csv("/Users/ebader/Documents/sc-RNA-seq/ohs_phenotypes.csv", header = TRUE, stringsAsFactors=FALSE)

###################
# 	Format data   #
###################
	irsVars <- metadata %>% select("ParticipantId","BIOSPEC_HCT","BIOSPEC_WBC","BIOSPEC_PLT","BIOSPEC_MCV","BIOSPEC_MCHC","BIOSPEC_RDWCV", "RiskGroup","IRSnoAge")
	irsVars <- reshape2::melt(irsVars, id.vars=c("ParticipantId","RiskGroup","IRSnoAge"))

	# Add custom facet titles
	cbcLabels <- c("Hematocrit","WBC","Platelets","MCV","MCHC","RDWCV")
	names(cbcLabels) <- unique(irsVars$variable)

	plt <- ggplot(irsVars, aes(x=factor(RiskGroup), y=value, fill=factor(RiskGroup))) + 
		   geom_violin() + 
		   scale_fill_manual(name = "Risk\nGroup",
		                     values=c("#2EC4B6","#9F3AA3","#FF9F1C","#E71D36"),
	                       breaks = c(1,2,3,4),
	                       labels = c("YL","YH","AL","AH")) +
		   geom_boxplot(aes(x=factor(RiskGroup), y=value), width=0.1, fill="white") + 
		   facet_wrap(~variable, labeller=labeller(variable=cbcLabels), scales = "free", nrow = 2) +
		   labs(x="Risk Group", y=NULL) +
		   theme_bw() +
		   theme(strip.text.x = element_text(size=14),
		         axis.text.x = element_blank(),
		         axis.ticks.x = element_blank(),
		         axis.title.x = element_blank(),
		         axis.text.y = element_text(size = 12))

	
	ggsave(plt, file = "/Users/ebader/Documents/Healthy_aging/Figures/CBCvars_riskgroup_violin_Fig1C.png", width = 7, height = 5, device = "png")

#``````````````````````````````````````````````````````` FIN ``````````````````````````````````````````````````````````````````````````#