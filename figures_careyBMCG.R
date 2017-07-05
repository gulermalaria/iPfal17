# this file generates figures 1 and 2 and suppl. figures 1, 3 and 4 for Carey et al. BMC Sys Bio
library(ggplot2)


##### FIGURE 1: histogram of PCH to determine PCH cutoffs for locations
# requires functions (filter_accession_file) and files from DEGs_march172017
accession_cambodia = filter_accession_file("Cambodia")
accession_vietnam = filter_accession_file("Vietnam")
accession_cv = rbind(accession_cambodia,accession_vietnam)

accession_data_k = (accession_cv[,c(4,22,36)])
accession_data_k[(is.na(accession_data_k$K13_KPdBT)),2] = "missing"
accession_data_k[(accession_data_k$K13_KPdBT == 0),2] = "missing"
accession_data_k[(accession_data_k$K13_KPdBT == 1),2] = "reference allele"
accession_data_k[(accession_data_k$K13_KPdBT == 2),2] = "mutant allele"
accession_data_k[(accession_data_k$K13_KPdBT == 3),2] = "mixed"

hist_k_pch = ggplot(accession_data_k, aes(parasite_clearance_halflife_hr, fill = as.character(K13_KPdBT_B))) + 
  scale_fill_manual(values = c("grey","black", "red","blue")) +
  geom_histogram(binwidth = .25) + 
  geom_vline(aes(xintercept=5), color = "black", size = 1) + xlim(0,12) +
  facet_wrap(c("Country"), ncol=4) + 
  xlab("Parasite Clearance Halflife") + ylab("Number of samples") + 
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_blank(), strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=1),
        panel.grid.major.y = element_line(color = "grey"),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.margin.x = unit(2, "lines"),
        plot.margin = unit(c(10,10,10,10),"mm"),
        legend.title=element_blank(),
        legend.position = c(.88, 0.5),
        legend.text=element_text(size=12))

setwd(paste(base_directory,'/figures',sep=""))
tiff(filename = "Fig1ab_histogram_K13_pch.tiff", width = 650, height = 360)
plot(hist_k_pch)
dev.off()

### FIGURE 1: DENDROGRAM OF EXPRESSION PROFILES

change_res3 = function(sample_names_country,gset) {
  sml <- c()
  for (i in 1:nrow(sample_names_country)) {
    if (is.na(sample_names_country[i,3]) | is.na(sample_names_country[i,2])) {
      sml[i] = 'X' }
    else if (as.numeric(sample_names_country[i,3]) > 5) {
      if (as.numeric(sample_names_country[i,2]) >=2) {
        sml[i] = 'R' }
      else {sml[i] = 'X'} }
    else if (as.numeric(sample_names_country[i,3]) < 5){
      if (as.numeric(sample_names_country[i,2]) == 1) {
        sml[i] = 'S' }
      else {sml[i] = 'X'} }
    else {
      sml[i] = 'X' } }
  # eliminate samples marked as "X" (samples in gset that aren't R/S or aren't in cambodia)
  sel <- which(sml != "X")
  sml <- sml[sel]
  #sml <- paste("G", sml, sep="")    # set group names
  return(sml)
}


# load accession_file with data from all samples
accession_cambodia = filter_accession_file("Cambodia")
accession_vietnam = filter_accession_file("Vietnam")
accession_cv = rbind(accession_cambodia,accession_vietnam)

#both
sample_names_cv = merge(sample_names,accession_cv[,c(22,36)], by = "row.names", all = T)
rownames(sample_names_cv) = sample_names_cv$Row.names
sample_names_cv = sample_names_cv[,2:4]

# convert text resistance status to R/S CAMBODIA
gset_cv = change_res(sample_names_cv,gset)
sml_cv = change_res3(sample_names_cv,gset)

###### BOTH COUNTRIES

#all genes
eset = exprs(gset_cv)
eset[is.na(eset)] = 0
eset = t(eset)
distance = dist(eset)
clust = hclust(distance)
plot(clust)
library(ggdendro)
library(ggplot2)
# Build dendrogram object from hclust results
dend <- as.dendrogram(clust)
dend_data <- dendro_data(dend, type = "rectangle")
dend_data$labels[,3] = sml_cv # change sample names to R/S
seg = dend_data$segments
seg$y = seg$y/10
seg$yend = seg$yend/10
dend_data$segments = seg
p <- ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = dend_data$labels, aes(x=x, y=y, color = factor(label)), size=9, shape = 124) + #15
  scale_colour_manual(values = c("red","blue"), guide = F) +
  guides(fill = F) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(p)
setwd(paste(base_directory,'/figures',sep=""))
tiff(filename = "Fig1b_dendrogram.tiff", width = 640, height = 220)
plot(p)
dev.off()

### metabolic only
# screen for metabolic genes
eset = exprs(gset_cv)
eset[is.na(eset)] = 0
eset = merge(blast_final,eset, by = "row.names")
eset = merge(eset,dictionary, by.x = "alignment", by.y = "ORF")
eset$model_ORF = NULL
eset$alignment = NULL
eset$Row.names = NULL
eset$query = NULL
eset = t(eset)
distance = dist(eset)
clust = hclust(distance)
library(ggdendro)
# Build dendrogram object from hclust results
dend <- as.dendrogram(clust)
dend_data <- dendro_data(dend, type = "rectangle")
dend_data$labels[,3] = sml_cv # change sample names to R/S
seg = dend_data$segments
seg$y = seg$y/10
seg$yend = seg$yend/10
dend_data$segments = seg
p <- ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = dend_data$labels, aes(x=x, y=y, color = factor(label)), size=9, shape = 15) +
  scale_colour_manual(values = c("red","blue")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


### FIGURE 2:  BAR PLOT OF CURATION
# x = number
# y = group/subsystem
# color = subsystem (black, grey, white, patterned)
# for each subsystem, there will be # Rxns, # genes, # modifications, and # citations bar
setwd(paste(base_directory,'/data',sep="")) # make sure figure 2 prep has ';' not ',' in subsystems
sub = read.table("figure2_prep.txt", sep = "\t", header = T, stringsAsFactors = F)
sub = as.data.frame(sub[,1:5])

s <- strsplit(sub$Subsystems, split = ";")
sub2 = data.frame(Reaction = rep(sub$Reactions, sapply(s, length)), 
                  Subsystem = unlist(s),
                  References = rep(sub$References, sapply(s, length)),
                  Notes = rep(sub$Notes, sapply(s, length)),
                  Genes = rep(sub$Genes, sapply(s, length))
)
sub2$Subsystem = trimws(sub2$Subsystem)
sub1 = sub
sub = sub2
# make subsystem more abstract
for (i in 1:nrow(sub)) {
  sub$Subsystem[i] = strsplit(sub$Subsystem[i], split = " ")[[1]][1]
}
# add collumns
sub$modify = NA
sub$citation = NA
sub$gene_count = NA
sub$reaction_count = as.vector(rep(1, nrow(sub)))

sub$Subsystem <- lapply(sub$Subsystem, function(x) { gsub("exchange", "Exchange", x) })
sub$Subsystem <- lapply(sub$Subsystem, function(x) { gsub("Export", "Exchange", x) })
sub$Subsystem <- lapply(sub$Subsystem, function(x) { gsub("Biomass", "Exchange", x) })
sub$Subsystem <- lapply(sub$Subsystem, function(x) { gsub("Ion", "Others", x) })
sub$Subsystem <- lapply(sub$Subsystem, function(x) { gsub("Expression", "Exchange", x) })
sub$Subsystem <- lapply(sub$Subsystem, function(x) { gsub("Spontaneous", "Others", x) })
sub$Subsystem <- lapply(sub$Subsystem, function(x) { gsub("Gap-filled", "Others", x) })
sub$Subsystem <- lapply(sub$Subsystem, function(x) { gsub("\\bHemoglobin\\b", "Hemoglobin/Hemozoin", x) })
sub$Subsystem <- lapply(sub$Subsystem, function(x) { gsub("\\bHemozoin\\b", "Hemoglobin/Hemozoin", x) })
sub$Subsystem <- lapply(sub$Subsystem, function(x) { gsub("Hemoglobin/Hemoglobin/Hemozoin", "Hemoglobin/Hemozoin", x) })

mod = read.table("edits.csv", sep = ",", header = T, stringsAsFactors = F)
colnames(mod) = c("Modify",'Delete')

sub$References[sub$References == ''] = NA
sub$Notes[sub$Notes == ''] = NA

# biniarize modifications etc
for(i in 1:nrow(sub)) {
  if (as.logical(sum(grepl(sub$Reaction[i],mod$Modify)))) {
    sub$modify[i] = 1 }
  else {sub$modify[i] = 0 }
}
# sum(sub$modify) > length(mod$Modify) because of duplicate subsystem assigmnets
for(i in 1:nrow(sub)) {
  if (!is.na(sub$References[i])) {
    sub$citation[i] = 1 }
  else {sub$citation[i] = 0 }
}
for(i in 1:nrow(sub)) {
  if (is.na(sub$Genes[i])) {
    sub$gene_count[i] = 0 }
  else if ( nchar(as.character(sub$Gene[i])) < 2 ) {
    sub$gene_count[i] = 0 }
  else { sub$gene_count[i] = 1 }
}
sub = sub[,c(2,6:9)]

# make summary table
summary_table = data.frame(matrix(vector(), length(unique(sub$Subsystem)), 5,
                                  dimnames=list(c(), c("group", "modified", "citation","reactions","genes"))),
                           stringsAsFactors=F)
summary_table$group = as.vector(unique(sub$Subsystem))

# get count data from sub dataframe and add to summary_table
for (i in 1:nrow(summary_table)) {
  if (
    is.na(summary_table$group[[i]][1])
  ) { print(i) 
  }
  else { 
    summary_table$modified[i] = sum(sub[which(grepl(summary_table$group[[i]][1],sub$Subsystem)),2]);
    summary_table$citation[i] = sum(sub[which(grepl(summary_table$group[[i]][1],sub$Subsystem)),3]); 
    summary_table$reactions[i] = sum(sub[which(grepl(summary_table$group[[i]][1],sub$Subsystem)),5]);
    summary_table$genes[i] = sum(sub[which(grepl(summary_table$group[[i]][1],sub$Subsystem)),4]) 
  }
}
summary_table = na.omit(summary_table)
colnames(summary_table) = c("Subsystem",'modified for iPfal17','with citation',"total reactions","with gene annotation")
summary_table$Subsystem = unlist(summary_table$Subsystem)  #reshape for plotting

#calculate percent with gene anotations
sumt = summary_table[-c(7,8,11),]
sum(sumt$GeneAnnotations)/sum(sumt$TotalReactions)

#reshape data for plotting
library(reshape2)
df = melt(summary_table)

library(ggplot2)

df = melt(summary_table[,c(1:3,5)]) # remove total because that's in first half of plotting function
p = ggplot(data = summary_table, aes(x = Subsystem, y = `total reactions`)) +
  geom_bar(stat = "identity",width = .9, color = "black",fill = "black") + 
  scale_y_sqrt(minor_breaks = seq(0, 700, 50)) + coord_flip() +
  geom_bar(data = df, aes(x = Subsystem, y = value, fill = variable),
           stat = "identity", color = "black",position = "dodge") + 
  scale_fill_manual(values = c("grey80", "white","grey50")) + guides(fill = F) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=1),
        panel.grid.major.x = element_line(color = "grey"),
        panel.grid.minor.x = element_line(color = "grey"),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.margin.x = unit(2, "lines"),
        plot.margin = unit(c(10,10,10,10),"mm"),
        legend.title=element_blank(),
        legend.position = c(.82, 0.75),
        legend.text=element_text(size=8))  
p = p + guides(fill= guide_legend(reverse=T, element_text(size = 8))) 
#setwd(paste(base_directory,'/figures',sep=""))
ggsave("Fig2_curation_summary.pdf", p, width = 7, height = 5)

##### SUPPL. FIGURE 1 a&b: DOT PLOT OF FOLD CHANGE VALUES
setwd(paste(base_directory,'/data',sep=""))

cam = read.table("cambodia.csv", sep = ",", header = T)
#cam = subset(cam, cam$P.Value<0.05)
cam2 = cam %>% mutate(sig= P.Value<0.05)
cam2$logFC = 2^cam2$logFC
colnames(cam2)[2] = "FC"
cam3 = cam2[order(cam2$FC),]
library(ggplot2)
rownames(cam3) = NULL
cam3$model_ORF = rownames(cam3)
colnames(cam3)[3] = "gene"
cam3$gene = as.numeric(cam3$gene)
c_plot = ggplot(cam3, aes(x = gene, y = FC)) + 
  geom_point(aes(color = sig)) + 
  scale_color_manual(values = c("black","red"), guide=FALSE) +
  guides(fill = F) + ylim(.5,2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        panel.background = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=1),
        panel.grid.major.y = element_line(color = "grey"),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.margin.x = unit(2, "lines"),
        plot.margin = unit(c(10,10,10,10),"mm"))


setwd(paste(base_directory,'/figures',sep=""))
tiff(filename = "SF1a_cambodia_fc.tiff", width = 380, height = 340)
plot(c_plot)
dev.off()

setwd(paste(base_directory,'/data',sep=""))
vi = read.table("vietnam.csv", sep = ",", header = T)
#cam = subset(vi, vi$P.Value<0.05)
vi2 = vi %>% mutate(sig= P.Value<0.05)
vi2$logFC = 2^vi2$logFC
colnames(vi2)[2] = "FC"
vi3 = vi2[order(vi2$FC),]
library(ggplot2)
rownames(vi3) = NULL
vi3$model_ORF = rownames(vi3)
colnames(vi3)[3] = "gene"
vi3$gene = as.numeric(vi3$gene)
v_plot = ggplot(vi3, aes(x = gene, y = FC)) + 
  geom_point(aes(color = sig)) + 
  scale_color_manual(values = c("black","red"), guide=FALSE) +
  guides(fill = F) + ylim(.5,2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        panel.background = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=1),
        panel.grid.major.y = element_line(color = "grey"),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.margin.x = unit(2, "lines"),
        plot.margin = unit(c(10,10,10,10),"mm"))

setwd(paste(base_directory,'/figures',sep=""))
tiff(filename = "SF1b_vietnam_fc.tiff", width = 380, height = 340)
plot(v_plot)
dev.off()

##### SUPPL. FIGURE 3: gene state heatmap
# load file (gene, gene state, associated reactions, associated subsystems)
# gene, state, state, state, state, reactions (separated by spaces), subsystems (separated by spaces and no space in names)

setwd(paste(base_directory,'/data',sep=""))
file_c = read.csv("gene_states_c_feb2017.csv", sep = ',', header = T, stringsAsFactors = F)
file_v = read.csv("gene_states_v_feb2017.csv", sep = ',', header = T, stringsAsFactors = F)
file = merge(file_c,file_v, by.x = 'genes', by.y = 'genes')

file = as.data.frame(file)
colnames(file) = c("Gene", "Cambodia (resistant)","Cambodia (sensitive)",
                   "Vietnam (resistant)","Vietnam (sensitive)")
for (i in 1:nrow(file)) {
  if (as.numeric(file[i,2]) >0.1) {
    file[i,2] = 5 }
  if (as.numeric(file[i,4]) >0.1) {
    file[i,4] = 5 } }
for (i in 1:nrow(file)) {
  for (j in c(3,5)) {
    if (as.numeric(file[i,j]) >0.1) {
      file[i,j] = 1 } } }

library(ggplot2)
library(reshape2)

file3 = melt(file)
mybreaks = levels(as.factor(file3$value))
cols = c("white","blue","red")
ord2 = hclust( dist(file, method = "euclidean"), method = "ward.D" )$order
ord1 = hclust( dist(t(file[,2:5]), method = "euclidean"), method = "ward.D" )$order

file3$Gene = factor(file3$Gene, levels = file3$Gene[ord2])
levels(file3$variable) = gsub(" ", "\n", levels(file3$variable))
plot = ggplot(file3, aes(x = variable, y = Gene, fill = as.factor(value))) + 
  geom_tile(stat = "identity", color = "gray", size = .001) + 
  scale_fill_manual(values = cols,breaks = mybreaks)
p = plot + ggtitle(NULL) + labs(y = 'Gene States')+
  theme(axis.title.y = element_text(size = 36, 
                                    color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20, angle = 0, 
                                   hjust = .5, vjust = 0, 
                                   color = "black"),
        legend.position = "none",
        axis.ticks.y=element_blank(),
        plot.title = element_text( size = 14),
        panel.background = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.margin.x = unit(2, "lines"),
        plot.margin = unit(c(10,10,10,10),"mm")) 

setwd(paste(base_directory,'/figures',sep=""))
tiff(filename = "SF3_gene_states.tiff", width = 510, height = 900)
plot(p)
dev.off()


##### SUPPL. FIGURE 4: enrichment analysis on gene states
# load file (gene, gene state, associated reactions, associated subsystems)
# gene, state, state, state, state, reactions (separated by spaces), subsystems (separated by spaces and no space in names)
setwd(paste(base_directory,'/data',sep=""))
file_c = read.csv("gene_states_c.csv", sep = ',', header = T, stringsAsFactors = F)
file_v = read.csv("gene_states_v.csv", sep = ',', header = T, stringsAsFactors = F)
file = merge(file_c,file_v, by.x = 'genes', by.y = 'genes')
grs = read.table("GeneRxnSubsystems.txt", sep = "\t", header = T, stringsAsFactors = F)
  # MAKE SURE TO DO THE FOLLOWING: in xls file (output from matlab), replace '-' with ';' in Subsystems, replace '-' with '' in reactions, save as txt

file = as.data.frame(file)
colnames(file) = c("Genes","c_res","c_sens","v_res","v_sens")
rownames(file) = file$Genes
rownames(grs) = grs$Gene
file = merge(file, grs, by = "row.names")
file$Gene = NULL; file$Row.names = NULL; 
library(tidyverse) 
file = mutate(file,
              res_consensus = (c_res+v_res)/2,
              sens_consensus = (c_sens+v_sens)/2)
file$c_res = NULL; file$v_res = NULL; file$v_sens = NULL; file$c_sens = NULL
file$res_sens_dif = NULL; file$X = NULL
colnames(file)[1] = "gene"; colnames(file)[2] = "reactions"

# subset list for only on or only off genes
sens_on = as_data_frame(subset(file, as.numeric(file$sens_consensus) == 1))
res_on = as_data_frame(subset(file, as.numeric(file$res_consensus) == 1))
sens_off = as_data_frame(subset(file, as.numeric(file$sens_consensus) == 0))
res_off = as_data_frame(subset(file, as.numeric(file$res_consensus) == 0))
rm(file)
files_all = list(res_on,res_off,sens_on,sens_off)
#access each as files_all[[x]] x = 1:4

#load model stats
# reactions, subsystems
model = read.table("figure2_prep.txt", sep = "\t", header = T, stringsAsFactors = F)
model = as.data.frame(model)
model$References = NULL; model$Notes = NULL; model$Genes = NULL
colnames(model)[2] = "Old_system"
model$Subsystem <- lapply(model$Old_system, function(x) { 
  gsub(" ", "", x) })

# prep model doc # count elements in each subsystem
add1 = NA # model file
model$Subsystem[which(model$Reactions == '3MOBtmt')] = "Transport"
model$Old_system[which(model$Reactions == '3MOBtmt')] = "Transport"
model$Subsystem[which(model$Reactions == 'EX_lac_L(e)')] = "Exchange"
model$Old_system[which(model$Reactions == 'EX_lac_L(e)')] = "Exchange"
model_old = model
model$Subsystem <- lapply(model$Subsystem, function(x) { 
  gsub("Transport[a-zA-Z]*", "Transport ", x) })
model$Subsystem <- lapply(model$Subsystem, function(x) { 
  gsub("Export|Ion|Expression|ProteinProduction|LipidProduction", "Exchange", x) })
model$Subsystem <- lapply(model$Subsystem, function(x) { 
  gsub("Hemoglobindigestion|Hemozoinproduction", "Hemoglobin", x) })
model$Subsystem <- lapply(model$Subsystem, function(x) { 
  gsub("RedoxRedox", "Redox", x) })
model$Subsystem <- lapply(model$Subsystem, function(x) { 
  gsub("Others|Metabolism|AminoAcids|Carbohydrates|Cofactors|Lipids|Nucleotides", "", x) })
model$Subsystem <- lapply(model$Subsystem, function(x) { 
  gsub(";", " ", x) })
model$Subsystem <- lapply(model$Subsystem, function(x) { 
  gsub(",", " ", x) })
model$Subsystem <- lapply(model$Subsystem, function(x) { 
  gsub("  ", " ", x) })
model[which(model$Subsystem == ''),3] = 'Lipids'

for (i in 1:nrow(model)) {
  #print(i)
  split = strsplit(as.character(model$Subsystem[i]),' ')
  temp = model[i,]
  for (j in 1:length(split[[1]])) {
    temp$Subsystem = split[[1]][j]
    add1 = rbind(add1,temp)
  }
  model$Subsystem[i] = NA
}
check = add1[add1$Subsystem == '',]
# prev_check = ''
# for (i in 2:nrow(check)) {
#   reaction_check = check[i,1]
#   if (reaction_check == prev_check) {}
#   else {print(model[which(grepl(reaction_check,model$Reactions)),]) 
#   print(model_old[which(grepl(reaction_check,model_old$Reactions)),])}
#   prev_check = reaction_check
# }
model = na.omit(add1)

enrichment_analysis = function(input_file) {
  file = input_file
  
  file$Subsystems <- lapply(file$Subsystems, function(x) { 
    gsub(" ", "", x) })
  file$Subsystems <- lapply(file$Subsystems, function(x) { 
    gsub("Transport[a-zA-Z]*", "Transport ", x) })
  file$Subsystems <- lapply(file$Subsystems, function(x) { 
    gsub("Export|Ion|Expression|ProteinProduction|LipidProduction", "Exchange", x) })
  file$Subsystems <- lapply(file$Subsystems, function(x) { 
    gsub("RedoxRedox", "Redox", x) })
  file$Subsystems <- lapply(file$Subsystems, function(x) { 
    gsub("Others|Metabolism|AminoAcids|Carbohydrates|Cofactors|Lipids|Nucleotides", "", x) })
  file$Subsystems <- lapply(file$Subsystems, function(x) { 
    gsub("Hemoglobindigestion|Hemozoinproduction", "Hemoglobin", x) })
  file$Subsystems <- lapply(file$Subsystems, function(x) { 
    gsub(",", " ", x) })
  file$Subsystems <- lapply(file$Subsystems, function(x) { 
    gsub(";", " ", x) })
  file$Subsystems <- lapply(file$Subsystems, function(x) { 
    gsub("  ", " ", x) })
  file$Subsystems[which(file$Subsystems == '')] = 'Lipids'
  
  # separate if multiple subsystems in one row in both documents
  add1 = NA # subsytem file
  for (i in 1:nrow(file)) {
    split = strsplit(as.character(file$Subsystems[i]),' ')
    temp = file[i,]
    for (j in 1:length(split[[1]])) {
      temp$Subsystems = split[[1]][j]
      add1 = rbind(add1,temp)
    }
    file$Subsystems[i] = NA
  }
  
  add1$Subsystems[add1$Subsystems == ''] = NA
  file = na.omit(add1)
  file = unique( file ) #(remove duplicate rows)
  
  # fishers exact
  # set up matrix : subsystem and on list = 38, other and on list = 232,
  # not on list and  subsystem = 74, not on list and other = 3324
  # fisher.test(matrix(c(38,74,232,3324),nrow=2,ncol=2),alternative = "two.sided")
  
  # x = list length
  x = nrow(file)
  # w = model length
  w = nrow(model)
  fisher_result = with(file, table(Subsystems))
  for (i in 1:dim(with(file, table(Subsystems)))) {
    # y = subsystem and on list
    y = with(file, table(Subsystems))[i]
    print(y)
    # z = model and subsystem 
    z = with(model, table(Subsystem))[i]
    print(z)
    ft = fisher.test(matrix(c(y,x-y,z,w-z),nrow=2,ncol=2),
                     alternative = "two.sided")
    fisher_result[i] = ft$p.value }
  
  # multiple testing correction
  fisher_result_adj = p.adjust(fisher_result, method = "fdr", n = length(fisher_result))
  return(fisher_result_adj)
}

fisher_result_final = list(1,2,3,4)
for (i in 1:4) {
  fisher_result_final[[i]] = enrichment_analysis(files_all[[i]])
}

fisher_heatmap = merge(fisher_result_final[[1]], fisher_result_final[[2]], by = "row.names", all = T)
rownames(fisher_heatmap) = fisher_heatmap$Row.names; fisher_heatmap$Row.names = NULL
fisher_heatmap = merge(fisher_heatmap, fisher_result_final[[3]], by = "row.names", all = T)
rownames(fisher_heatmap) = fisher_heatmap$Row.names; fisher_heatmap$Row.names = NULL
fisher_heatmap = merge(fisher_heatmap, fisher_result_final[[4]], by = "row.names", all = T)
rownames(fisher_heatmap) = fisher_heatmap$Row.names; fisher_heatmap$Row.names = NULL
fisher_heatmap[is.na(fisher_heatmap)] = 1
colnames(fisher_heatmap) = c("res_on","res_off","sens_on","sens_off")

df = as.matrix(fisher_heatmap)
rownames(df)

library(ggplot2)
library(reshape2)

file2 = melt(df)
file2$value[file2$value <= 0.001] = 0.001
file2$value[(file2$value <= 0.01)&(file2$value > 0.001)] = 0.01
file2$value[(file2$value <= 0.05)&(file2$value > 0.01)] = 0.05
file2$value[(file2$value > 0.05)] = 1
colnames(file2) = c("Var1","Var2","value")
file2 = file2[which(grepl("res_on",file2$Var2)|grepl("sens_on",file2$Var2)),]
file2$Var2 = as.character(file2$Var2)
file2$Var2[which(grepl('res_on',file2$Var2))] = 'Resistant'
file2$Var2[which(grepl('sens_on',file2$Var2))] = 'Sensitive'
file2$Var2 = as.factor(file2$Var2)

file2$Var1 = as.character(file2$Var1)
file2$Var1[which(grepl('RedoxMitochondrialAntioxidantSystem',file2$Var1))] = 'MitochondiralRedox'
file2$Var1[which(grepl('PhosphatidyletanolaminePhosphatidylserine',file2$Var1))] = 'Phospholipids(PE/PS)'
file2$Var1[which(grepl('Phosphatidylcholine',file2$Var1))] = 'Phospholipids(PC)'
file2$Var1 = as.factor(file2$Var1)

plot = ggplot(file2, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = factor(value)),stat = "identity",  size = .001) + 
  scale_fill_manual(values = c("0.001" = "grey13","0.01" = "grey37","0.05"=  "grey55", "1"= "white")) + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, angle = 0, 
                                   hjust = 1, vjust = 0.5, 
                                   color = "black"),
        axis.text.x = element_text(size = 10, angle = 90, 
                                   hjust = 1, vjust = .5,
                                   color = "black"),
        plot.title = element_text( size = 14),
        panel.background = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.margin.x = unit(2, "lines"),
        plot.margin = unit(c(10,30,10,10),"mm"),
        legend.title=element_blank(), 
        legend.position = c(1.1, 0.4),
        legend.text=element_text(size=10),
        legend.key = element_rect(colour = 'black', size = 0.5, linetype='solid'))+ 
  coord_flip()  
setwd(paste(base_directory,'/figures',sep=""))
ggsave("SF4_enrichment.pdf", plot, width = 8.5, height = 2.8)
