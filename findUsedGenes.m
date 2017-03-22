function [ gene_r genes_used] = findUsedGenes(model, genes)
% genes that are utilized
% remove genes not associated with any reactions
% [ gene_r genes_used] = findUsedGenes(model, genes)
%INPUTS
% model                 COBRA model structure
% genes                 string of single gene or cell array of multiple
%                       genes for which rxns are desired.
%

%OPUTPUTS
% new_model             new model without excess genes.  
% genes_r               genes removed
% 
% Adapted from FindRxnsFromGenes by Nathan Lewis 02/16/08/edited 04/05/09 
% (MUCH faster now -- NL)/edited 06/11/10 (yet even faster now -- NL)
% MAC 10/10/16


if ~iscell(genes)
    gene = genes;
    clear genes
    genes{1} = gene;
    clear gene
end

if iscell(genes{1})
    for i = 1:length(genes)
        gene(i) = genes{i};
    end
    clear genes
    genes = gene;
    clear gene
end

model2 = model;
model.genes = regexprep(model.genes,'-','_DASH_');
model.genes = regexprep(model.genes,'\.','_POINT_');
genes = regexprep(genes,'-','_DASH_');
genes = regexprep(genes,'\.','_POINT_');

%find where the genes are located in the rxnGeneMat
GeneID(1) = 0;
for j = 1:length(genes)
    Ind = find(~cellfun('isempty', regexp(model.genes,cat(2,'^',genes{j},'$'))));
    
            if ~isempty(Ind)
                GeneID(j) = Ind;
            end

end

if min(GeneID) == 0
    warning('A gene was not found in the model!')
    results = struct([]);
    if max(GeneID) ==0,results = struct([]);ListResults = {};
    return
    end
    Ind = find(GeneID==0);
    GeneID(Ind) = [];
    genes(Ind) = [];
end

gene_r = []; ID_to_keep = [];
matrikx = model2.rxnGeneMat;
for i = 1:length(GeneID)
    Ind_rxns = find(matrikx(:,GeneID(i))==1);
    % if Ind_rxns == empty matrix -> not associated with any reactions, so
    % remove columnn of rxnGeneMat
    if isempty(Ind_rxns)
        gene_r = [gene_r ; genes(GeneID(i))];
        model2.rxnGeneMat(:,GeneID(i)) = []
        ID_to_keep = [ID_to_keep ; (GeneID(i))];
    end
    disp(i)
end


model2.genes(ID_to_keep) = [];
genes_used = model2.genes


end