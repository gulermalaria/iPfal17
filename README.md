# iPfal17
Code availability for the manuscript titled "Novel *Plasmodium falciparum* metabolic network reconstruction identifies shifts associated with clinical antimalarial resistance " by Carey *et al*.

The article is submitted BMC Genomics (DOI: 10.1186/s12864-017-3905-1) and available on BioRxiv (see here: http://www.biorxiv.org/content/early/2017/03/23/119941).


All data in subdirectory entitled data.

### model_curation_careyBMCG.m

This MATLAB script incorporates edits from 3 files (biomass_edits_final.xls, edits_additions_or_modifications_final.xls, and reaction_delete_edits_final.xls) into the Plata *et al.* reconstruction (plata_orig.xml) to generate our reconstruction *iPfal17* and genes_used_by_model.csv. 

### findUsedGenes.m 

This MATLAB function identifies genes in model associated with reactions. Necessary for test_exp_knockouts_careyBMCG.m. 

### test_exp_knockouts_careyBMCG.m 

This MATLAB script identifies essential genes and lethal enzyme inhibition in the model. It also compares to predictions in Plata et al. 2010. It requires the model from model_curation_careyBMCG.m and exp_predictions_KOs_wPlata.xls. It is used in generating Table 4 and Suppl. Table 7. It also generates essential_genes.csv for Plasmogem_careyBMCG.R.

### NEA_careyBMCG.m 

This MATLAB script implements Network Expression Analysis using MADE (Jensen and Papin, 2011) and iPfal17. Expression files (cambodia.csv and vietnam.csv) from DEGs_careyBMCG.Rmd and model from model_curation_careyBMCG.m are used. It writes the following files used in figures_careyBMCG.R:
	GeneRxnsSubsytems.xls
	figure2_prep.xls
	gene_states_c.csv
	gene_states_v.csv

### metabolic_tasks_careyBMCG.m 

This MATLAB script evaluates a model on 11 metabolic tasks (see Table  3). Some tasks are meant to fail. All are based on experimental literature or COBRA standards. It requires the model from model_curation_careyBMCG.m.

### DEGs_careyBMCG.Rmd

This R markdown script was used to calculate differential expression of ring-stage samples from GSE59097.   This script inputs preprocessed and normalized expression data to calculate differential expression. 
The following files are required:
	mok_accession_data_removed_ifNoRNA_orIfGam_orIfNoPCH.csv (from supplemental table in Mok *et al.* 2015)
	GPL18893.txt (acquired from GSE59097)
	blast_output_pf_e0001 (generated from blast_script_careyBMCG.rtf)
	genes_used_by_model.csv (acquired from model_curation_careyBMCG.m)
	Pfalciparum3D7_GeneAliases.csv (acquired from PlasmoDB.org)
The following files are generated:
	cambodia.csv
	vietnam.csv

### RF_careyBMCG.R 

This R script was used to generate random forest classifiers of metadata and expression data associated with GSE59097.   This script inputs preprocessed and normalized expression data, uses requires files mok_accession_data_removed_ifNoRNA_orIfGam_orIfNoPCH.csv, GPL18893.txt, blast_output_pf_e0001, and generates Suppl. Figure 2 and variable_importance_expression.csv. 

### Plasmogem_careyBMCG.R 

This R script is used to generate Suppl. Table 6. This script compares iPfal17 predictions (File essential_genes.csv) to *P. berghei* essentiality results generated by the PlasmoGEM team (in file plasmogem.csv which is subsetted into essential_plasmogem.txt and dispensible_plasmogem.txt).
Note: used dictionary from DEGs_careyBMCG.R

### figures_careyBMCG.R 

This R script generates Figures 1 and 2 and Suppl. Figures 1, 3 and 4 for Carey et al. BMC Sys Bio. To generate Figure 1, functions and files from DEGs_careyBMCG.R are needed. Figure 2 requires figure2_prep.txt and edits.csv. Suppl. Figure 1 requires cambodia.csv and vietnam.csv. Suppl. Figure 3 requires gene_states_c.csv and gene_states_v.csv. Suppl. Figure 4 requires GeneRxnSubsystems.txt and figure2_prep.txt.

### blast_script_careyBMCG.rtf

This text file lists commands to reannotate platform probe IDs. Probe sequences must first be reverse complemented. It will generate blast_output_pf_e0001.
