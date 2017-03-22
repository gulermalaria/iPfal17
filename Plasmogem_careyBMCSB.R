library(tidyverse)
setwd(paste(base_directory,'/data',sep=""))
plasmo = read_csv('plasmogem.csv')
plasmo_e = subset(plasmo, grepl('essential',plasmo$phenotype))
plasmo_d = rbind(subset(plasmo, grepl('dispensable',plasmo$phenotype)),
                 subset(plasmo, grepl('fast',plasmo$phenotype)))

# write.table(plasmo_e,"plasmo_e.csv",sep = ',')
# write.table(plasmo_d,"plasmo_d.csv",sep = ',')

# take list of only genes
# plasmoDB
# search > new search > genes > annotations etc > geneIDs > 
# upload list
# add step > transform by orthology > Pfal 3D7

# above (essential) list converted to 3D7 gene IDs
essential3d7 = read_tsv('essential_plasmogem.txt')
colnames(essential3d7)[1] = "Genes"
# above (dispensible) list converted to 3D7 gene IDs
dispensible3d7 = read_tsv('dispensible_plasmogem.txt')
colnames(dispensible3d7)[1] = "Genes"

# run dictionary in DEGs_from_preprocessed script to get genes in model
dictionary = as.data.frame(dictionary)

# get model gene IDS added to essential list
essential_with_model = merge(essential3d7, dictionary, by.x = "Genes", by.y = "ORF")
plasmogem_essential_genes = essential_with_model$model_ORF #130
# get model gene IDS added to dispensible list
disp_with_model = merge(dispensible3d7, dictionary, by.x = "Genes", by.y = "ORF")
plasmogem_dispensible_genes = disp_with_model$model_ORF #130

# note dimensions decrease due to lack of gene ID mapping or nonmetabolic genes

#model predictions
predictions = read_csv('essential_Genes_feb2017.csv') #106

plasmogem_essential_genes = as_data_frame(plasmogem_essential_genes)
plasmogem_essential_genes = mutate(plasmogem_essential_genes,
                                   prediction = 'essential',
                                   our_prediction = NA,
                                   binary = NA)
colnames(plasmogem_essential_genes)[1] = 'gene'
for (i in 1:nrow(plasmogem_essential_genes)) {
  if (grepl(plasmogem_essential_genes$gene[i], predictions)) {
    plasmogem_essential_genes$our_prediction[i] = 'essential'
    plasmogem_essential_genes$binary[i] = 1
  }
  else {plasmogem_essential_genes$our_prediction[i] = 'dispensible'
        plasmogem_essential_genes$binary[i] = 0}
}

plasmogem_dispensible_genes = as_data_frame(plasmogem_dispensible_genes)
plasmogem_dispensible_genes = mutate(plasmogem_dispensible_genes,
                                     prediction = 'dispensible',
                                   our_prediction = NA,
                                   binary = NA)
colnames(plasmogem_dispensible_genes)[1] = c('gene')
for (i in 1:nrow(plasmogem_dispensible_genes)) {
  if (grepl(plasmogem_dispensible_genes$gene[i], predictions)) {
    plasmogem_dispensible_genes$our_prediction[i] = 'essential'
    plasmogem_dispensible_genes$binary[i] = 0
  }
  else {plasmogem_dispensible_genes$our_prediction[i] = 'dispensible'
        plasmogem_dispensible_genes$binary[i] = 1}
}
  
sum(plasmogem_dispensible_genes$binary)/nrow(plasmogem_dispensible_genes)
sum(plasmogem_essential_genes$binary)/nrow(plasmogem_essential_genes)
(sum(plasmogem_dispensible_genes$binary)+sum(plasmogem_essential_genes$binary))/
  (nrow(plasmogem_dispensible_genes)+nrow(plasmogem_essential_genes))

write.table(rbind(plasmogem_essential_genes[,1:3],plasmogem_dispensible_genes[,1:3]), 
            file='PlasmoGem_compared_to_predictions_Feb.csv',
            sep = ',',quote = F, row.names = F, col.names = T)
