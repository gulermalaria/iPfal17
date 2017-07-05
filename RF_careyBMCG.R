
library(Biobase)
library(GEOquery)
library(stringr)
library(plyr)
library(randomForest)
library(e1071)
library(caret)
library(tidyverse)
# base_directory = getwd() 

set.seed(2017) # submission/publication year

setwd(paste(base_directory,'/data',sep = ''))

# STRINGS MUST BE FACTORS HERE
accession_data2 = read.table("mok_accession_data_removed_ifNoRNA_orIfGam_orIfNoPCH.csv", header = T, sep = ",")
accession_data = accession_data2[,4:48] # REMOVE DUPLICATE COLUMNS
rownames(accession_data) = accession_data[,1]
accession_data = accession_data[accession_data$Asexual_stage_hpi <= 10,] # EARLY RINGS ONLY

data = mutate(accession_data, # identify if resistant sample
              resistant1 = ifelse(is.na(parasite_clearance_halflife_hr) | is.na(K13_KPdBT_B),
                                  NA,
                                  'notNA'),
              resistantPCH = ifelse(parasite_clearance_halflife_hr >5,
                                    'yes',
                                    'no'),
              resistantMut = ifelse(K13_KPdBT_B >=2,
                                    'yes',
                                    'no'))
data = mutate(data,
              resistant = ifelse(is.na(resistant1),
                                 resistant1,
                                 ifelse(resistantPCH == resistantMut,
                                        resistantMut,
                                        NA)))
# remove extra variables that are essentially repeats or uninformative
data$resistant1 = NULL; data$resistantPCH = NULL; data$resistantMut = NULL        
data$parasite_clearance_halflife_hr = NULL; data$K13_KPdBT_B = NULL
data$GrpA_subtype = NULL;data$GrpB_subtype = NULL;data$GrpC_subtype = NULL
data$Patient_Code = NULL;data$Barcode = NULL; data$Field_site_Code = NULL;  
data$Hb_g_per_dL = NULL;data$log2_ratio_gam_per_total_parasite = NULL
data$SE_of_clearance = NULL; data$Clearance_rate_constant = NULL; 
data$PC99_h = NULL;data$PC95_h = NULL
data$PC90_h = NULL; data$PC50_h = NULL; data$Tlag_h= NULL
data$Intercept_tlag= NULL; data$R2 = NULL; data$Country = NULL;
data$ratio_gam_per_total_parasite = NULL;data$mean_gam_exp_ratio = NULL
data$K13_KPdBT_B = NULL; data$total_parasite_count = NULL #repetative with gam and asexual

data = subset(data, !is.na(data$resistant))

colnames(data)[which(colnames(data) == "crt_N326S")] = c("crt_mutation1")
colnames(data)[which(colnames(data) == "arps10_V127M")] = c("aprs_mutation")
colnames(data)[which(colnames(data) == "crt_I356T")] = c("crt_mutation2")
colnames(data)[which(colnames(data) == "fd_D193Y")] = c("fd_mutation")
colnames(data)[which(colnames(data) == "mdr2_T484I")] = c("mdr_mutation")

# begin by spliting datafile
indx <- sample(2, nrow(data), replace = T, prob = c(0.6,0.4))
training <- data[indx==1,] # training set
validation <- data[indx==2,] # validation dataset
# table(training$resistant)/nrow(training)
# table(validation$resistant)/nrow(validation) # similar distribution

training$resistant = as.factor(training$resistant)
# if class == factor, classifying tree, else regressor
# class(cross.sell.dev$resistant)

# MAKE FORMULA
# set up columns to evalutate and columns to avoid
varNames <- names(training)
varNames <- varNames[!varNames %in%  # Exclude IDs and Response variables
                       c("resistant","geo_accession")]

# add + sign between exploratory variables
varNames1 <- paste(varNames, collapse = "+")
# Add response variable and convert to a formula object
rf.form <- as.formula(paste("resistant", varNames1, sep = " ~ "))

# for (i in 1:ncol(training)) {
#   if (sum(is.na(training[,i]))>1) {
#     print(i) } }

training.rf <- randomForest(rf.form, data = training,
                            ntree=500, na.action = "na.omit", importance=T)
# plot(training.rf) # see how many trees are needed

varImpPlot(training.rf, # identify important variables
           sort = T,
           main="Variable Importance",
           n.var=10)

var.imp  = as.data.frame(training.rf$importance)
var.imp$Variables = row.names(var.imp) # make row names as columns

# set up for figure
var_importance <- data_frame(variable=rownames(var.imp),
                             importance=var.imp$MeanDecreaseAccuracy)

var_importance <- arrange(var_importance, desc(importance))
#write.table(var_importance,'variable_importance_metadata.csv',sep = ',', quote = F, col.names = T, row.names= T)

var_plot = var_importance
var_plot$variable = factor(var_plot$variable, levels = var_plot$variable[order(var_plot$importance)])

# Predicting response variable
training$predicted.response = predict(training.rf ,training)

# COnfusion matrix
confusionMatrix(data=training$predicted.response,
                reference=training$resistant,
                positive= 'yes')

# Predicting response variable
validation$predicted.response <- predict(training.rf,validation)

# Create Confusion Matrix
val = confusionMatrix(data = validation$predicted.response,
                reference = validation$resistant,
                positive='yes')

accur2 = paste('Accuracy', paste(round(val$overall[[1]],4)*100,'%',sep = ''), sep = ': ')
sens2 = paste('Sensitivity', paste(round(val$byClass[[1]],4)*100,'%',sep = ''), sep = ': ')
spec2 = paste('Specificity', paste(round(val$byClass[[2]],4)*100,'%',sep = ''), sep = ': ')

var_plot2 = var_plot

p2 = ggplot(var_plot2) + geom_point(aes(x=variable, y=importance*100)) + 
  xlab("Variable") + 
  ylab("Variable Importance") +
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=1),
        panel.spacing.x = unit(2, "lines"),
        plot.margin = unit(c(10,10,10,10),"mm")) +
  coord_flip() +
  annotate("rect", xmin = 3, xmax = 7, ymin = 3.22, ymax = 7.5,alpha = 1, fill = "white") +
  annotate("text", x = 6, y = 5.5, label = accur2, size = 4) +
  annotate("text", x = 5, y = 5.5, label = sens2, size = 4) +
  annotate("text", x = 4, y = 5.5, label = spec2, size = 4) 
# 
# tiff(filename = "variable_importance_paper_metadata1.tiff", width = 570, height = 525)
# plot(p)
# dev.off()

### WITH EXPRESSION DATA

# OPEN EXPRESSION DATA (PREPROCESSED AND NORMALIZED)
gset2 <- getGEO("GSE59097", GSEMatrix =TRUE,getGPL= F)
gset <- gset2[[1]] # CONVERT TO EXPRESSION SET
eset = exprs(gset)

# imput missing values
for (j in 1:(nrow(eset))) {
  imp = mean(eset[j,], na.rm = T)
  eset[j,is.na(eset[j,])] = imp
}

# correct probes
platform = read.table("GPL18893.txt", sep = "\t", header = T, stringsAsFactors = F)
#write.table(platform[,1:2], file = "GPL_mok.tsv", sep = "\t", col.names = F, row.names = F, quote=FALSE)

# CORRECTLY ANNOTE PROBES: BLAST ALL PROBE SEQUENCES TO PlasmoDB-28_Pfalciparum3D7_AnnotatedTranscripts ACQUIRED AUG 17 2016
# open blast results with evalue cutoff of 0.0001
setwd(paste(base_directory,'/data/blast_results_probes',sep = ''))
blast = read.table("blast_output_pf_e0001", sep = ",", header = T, stringsAsFactors = F)
blast1 = subset(blast, percent_identity >95) # remove if not perfect match
blast2 = subset(blast1, score > 100) # remove if low score
blast3 = subset(blast2, mismatch < 1) # remove if mismatches
blast4 = subset(blast3, gaps < 1) # remove if gaps
blast5 = blast4[,1:2] # remove extra columns (keep query and alignmnet)

# remove alternative splicing notation
select = blast5[which(grepl("\\.",blast5$alignment)),2]
blast5[which(grepl("\\.",blast5$alignment)),2] =substr(select,1,nchar(select)-2)

rownames(blast5) = NULL
blast6 = (unique(blast5))
blast_final = aggregate(alignment ~ query, blast6, paste, collapse = " ")
rownames(blast_final) = blast_final$query

rm(blast,blast1,blast2,blast3,blast4,blast5,blast6,platform,select)

setwd(paste(base_directory,'/data/',sep = ''))
platform = read.table("GPL18893.txt", sep = "\t", header = T, stringsAsFactors = F,
                      na.strings = c("NA",""))
platform = (subset(platform, !grepl("_Dd2",ID))) # remove genotyping probes
platform = (subset(platform, !grepl("_HB3",ID)))
platform = (subset(platform, !grepl("_3D7",ID)))
platform = (subset(platform, !grepl("_7G8",ID)))
platform = (subset(platform, !grepl("_D10",ID)))

rownames(platform) = platform$ID
platform_with_blast = merge(platform, blast_final, by = "row.names", all = T)

rm(platform)

# Correct probe labels
eset2 = merge(eset, platform_with_blast, by.x = "row.names", by.y = "Row.names", all = T)
eset3 = subset(eset2, !(is.na(eset2$alignment))) # remove if no alignment result (blast)
eset4 = (eset3[rowSums(is.na(eset3))<1000,]) # remove if no gene expression
eset4$Row.names = NULL; eset4$ORF = NULL; eset4$ORF_old = NULL; eset4$query = NULL;
eset4$SPOT_ID = NULL; eset4$SEQUENCE = NULL; eset4$CHROMOSOME = NULL; eset4$POISITION = NULL
eset4$ID = NULL
rownames(eset4) = make.names(eset4$alignment, unique=TRUE)
eset4$alignment = NULL

# flip
eset = t(eset4)

# STRINGS MUST BE FACTORS HERE
accession_data2 = read.table("mok_accession_data_removed_ifNoRNA_orIfGam_orIfNoPCH.csv", header = T, sep = ",")
accession_data = accession_data2[,4:48] # REMOVE DUPLICATE COLUMNS
rownames(accession_data) = accession_data[,1]
accession_data = accession_data[accession_data$Asexual_stage_hpi <= 10,] # EARLY RINGS ONLY

# MERGE ESET AND ACCESSION
data = merge(eset,accession_data[,c(22,36)], by = "row.names")
rownames(data) = data$Row.names; data$Row.names = NULL
colnames(data) <- lapply(colnames(data), function(x) { gsub("[.]", "_", x) })
colnames(data) <- lapply(colnames(data), function(x) { gsub("-", "_", x) })
data = mutate(data, 
              resistant1 = ifelse(is.na(parasite_clearance_halflife_hr) | is.na(K13_KPdBT_B),
                                  NA,
                                  'notNA'),
              resistantPCH = ifelse(parasite_clearance_halflife_hr >5,
                                    'yes',
                                    'no'),
              resistantMut = ifelse(K13_KPdBT_B >=2,
                                    'yes',
                                    'no'))
data = mutate(data,
              resistant = ifelse(is.na(resistant1),
                                 resistant1,
                                 ifelse(resistantPCH == resistantMut,
                                        resistantMut,
                                        NA)))

data$parasite_clearance_halflife_hr = NULL; data$K13_KPdBT_B = NULL; 
data$resistant1 = NULL; data$resistantPCH = NULL;data$resistantMut = NULL;

data = subset(data, !is.na(data$resistant))

# begin by spliting datafile
indx <- sample(2, nrow(data), replace = T, prob = c(0.6,0.4))
training <- data[indx==1,] # training set
validation <- data[indx==2,] # validation dataset
# table(training$resistant)/nrow(training)
# table(validation$resistant)/nrow(validation) # similar distribution

training$resistant = as.factor(training$resistant)
# if class == factor, classifying tree, else regressor
# class(cross.sell.dev$resistant)

# MAKE FORMULA
# set up columns to evalutate and columns to avoid
varNames <- names(training)
varNames <- varNames[!varNames %in%  # Exclude IDs and Response variables
                       c("resistant","geo_accession")]

# add + sign between exploratory variables
varNames1 <- paste(varNames, collapse = "+")
# Add response variable and convert to a formula object
rf.form <- as.formula(paste("resistant", varNames1, sep = " ~ "))

# for (i in 1:ncol(training)) {
#   if (sum(is.na(training[,i]))>1) {
#     print(i) } }

training.rf <- randomForest(rf.form, data = training,
                            ntree=500, na.action = "na.omit", importance=T)
# plot(training.rf) # see how many trees are needed

varImpPlot(training.rf, # identify important variables
           sort = T,
           main="Variable Importance",
           n.var=10)

library(ggplot2)
var.imp = as.data.frame(training.rf$importance)

# make row names as columns
var.imp$Variables = row.names(var.imp)
file_g = var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]
file_a = var.imp[order(var.imp$MeanDecreaseAccur,decreasing = T),]

var_importance <- data_frame(variable=rownames(var.imp),
                             importance=var.imp$MeanDecreaseAccuracy)
var_importance = mutate(var_importance,
                        variable = ifelse(nchar(variable)>13,
                                          paste(substr(variable,1,13),'and_others', sep = '_'),
                                          substr(variable,1,13)))
var_importance <- arrange(var_importance, desc(importance))
write.table(var_importance,'variable_importance_expression.csv',sep = ',', quote = F, col.names = T, row.names= T)


var_plot = var_importance[1:25,]

# Predicting response variable
training$predicted.response <- predict(training.rf ,training)

# COnfusion matrix
library(e1071)
library(caret)
confusionMatrix(data=training$predicted.response,
                reference=training$resistant,
                positive= 'yes')

# Predicting response variable
validation$predicted.response <- predict(training.rf ,validation)

# Create Confusion Matrix
val = confusionMatrix(data = validation$predicted.response,
                reference = validation$resistant,
                positive='yes')

accur = paste('Accuracy', paste(round(val$overall[[1]],4)*100,'%',sep = ''), sep = ': ')
sens = paste('Sensitivity', paste(round(val$byClass[[1]],4)*100,'%',sep = ''), sep = ': ')
spec = paste('Specificity', paste(round(val$byClass[[2]],4)*100,'%',sep = ''), sep = ': ')

var_plot$variable = as.character(var_plot$variable)
var_plot$variable[grepl('PF3D7_0222300',var_plot$variable)][1] = "PF3D7_0222300_and_others_1"
var_plot$variable[grepl('PF3D7_0222300',var_plot$variable)][2] =  "PF3D7_0222300_and_others_2"
var_plot$variable[grepl('PF3D7_0222300',var_plot$variable)][3] =  "PF3D7_0222300_and_others_3"
var_plot$variable[grepl('PF3D7_0222300',var_plot$variable)][4] =  "PF3D7_0222300_and_others_4"
var_plot$variable[grepl('PF3D7_0222300',var_plot$variable)][5] =  "PF3D7_0222300_and_others_5"

var_plot$variable = as.factor(var_plot$variable)
var_plot$variable = factor(var_plot$variable, levels = var_plot$variable[order(var_plot$importance)])

p1 = ggplot(var_plot) + geom_point(aes(x=variable, y=importance*100)) + 
  xlab("Variable") + 
  ylab("Variable Importance") +
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=1),
        panel.spacing.x = unit(2, "lines"),
        plot.margin = unit(c(10,10,10,10),"mm")) +
  coord_flip() +
  annotate("rect", xmin = 3, xmax = 7.5, ymin = 0.065, ymax = 0.1,alpha = 1, fill = "white") +
  annotate("text", x = 6.5, y = .0828, label = accur, size = 4) +
  annotate("text", x = 5.25, y = .0828, label = sens, size = 4) +
  annotate("text", x = 4, y = .0828, label = spec, size = 4)  


library(gridExtra)
g = arrangeGrob(p1,p2, nrow = 2)
ggsave(file = 'SF2_RF.pdf', g, width = 8, height = 11)
# add a and b to figure in illustrator