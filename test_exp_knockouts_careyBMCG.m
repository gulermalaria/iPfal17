model = pf_cobra
x = 0;
[ gene_r genes_used] = findUsedGenes(model, model.genes) 

%% single gene deletionsa
[grRatio_genes,grRateKO,~,delRxns,~] = singleGeneDeletion(model,'FBA', genes_used);
%[grRatio_genes_c,~,~,~,~] = singleGeneDeletion(pf_cobra,'FBA',genes_used);

newmodel_gene_deletions = table(genes_used, grRatio_genes,...
    'VariableNames',{'genes','KOgrowthRatio'});
%indx = find(single(newmodel_gene_deletions{:,2}) < x);
newmodel_gene_deletions2 = newmodel_gene_deletions(newmodel_gene_deletions{:,2} <=0.1,:);
writetable(newmodel_gene_deletions2,'essential_genes.csv')

%% single reaction deletions
[grRatio_rxns,~,~,~,~,~] = singleRxnDeletion(model,'FBA');

newmodel_rxn_deletions = table(model.rxns, grRatio_rxns,...
    'VariableNames',{'reactions','KOgrowthRatio'});

newmodel_rxn_deletions2 = newmodel_rxn_deletions(newmodel_rxn_deletions{:,2} <=0.1,:);

% interpreting rxn deletions
indx = find(single(newmodel_rxn_deletions{:,2}) <=0.1);
lethalKOrxn = newmodel_rxn_deletions(indx,:);
reaction_names = model.rxnNames(indx);
rxn_formulas = printRxnFormula(model,model.rxns(indx));
rxn_EC = model.rxnECNumbers(indx);
lethalKOrxn = [lethalKOrxn reaction_names rxn_formulas rxn_EC];
% check 

cd('C:\Users\mac9jc\Documents\MATLAB\curation\final');
file = (readtable('exp_predictions_KOs_wPlata.xls')); 
file2 = [file.Properties.VariableNames; table2cell(file)];
lethal = [file.Properties.VariableNames; file2(strcmpi('L',file2(:,4)),:)];
nonlethal = [file.Properties.VariableNames; file2(strcmpi('NL',file2(:,4)),:)];
[m,n] = size(lethal);
[q,r] = size(nonlethal);

for i = 2:m % skip header row
    if strcmp({''},lethal(i,2)); % if ec is blank, skip
        disp(i);
        disp('EC is blank');
    end
    s = sum(strcmpi(lethal(i,2),rxn_EC));
    if s == 0 % if lethal(i) EC has no match in model's list of lethal predictions
        if strcmp({''},lethal(i,11)); % and if there is not modification done in experiment
            disp(i);
            disp('Incorrect prediction');
            lethal{i,n} = 'Incorrect';
        else
            disp(i);
            disp('Conditional knockout, see below');
            lethal{i,n} = 'Condition KO';
        end
    elseif s == 1 % if there is one EC match 
         disp(i);
         disp('correct');
         lethal{i,n} = 'Correct';
    else disp(i) % otherwise display error
        disp('error');
    end
end

for i = 2:q
    if strcmp({''},nonlethal(i,2)); % if ec is blank, skip
        disp(i);
        disp('EC is blank');
    end
    s = sum(strcmpi(nonlethal(i,2),rxn_EC));
    if s == 1 % if nonlethal(i) EC has a match in model's list of lethal predictions
        if strcmp({''},nonlethal(i,11)) % and if there is not modification done in experiment
            disp(i);
            disp('Incorrect prediction');
            nonlethal{i,r} = 'Incorrect';
        else
            disp(i);
            disp('Conditional knockout, see below');
            nonlethal{i,r} = 'Conditional KO';
        end
    elseif s == 0 % if there are no EC matchs to KO list
         disp(i);
         disp('correct');
         nonlethal{i,r} = 'Correct';
    else disp(i) % otherwise display error
        disp('error');
    end
end

%% NOTES: {'Plasmepsin II'} and cytosolic lysyl-tRNA synthetase are lethal (try gene KO) but reaction isn't associated with any ECs, because there are many

%% conditonal deletions
opt = optimizeCbModel(model);
% carbamoyl-P synthase 6.3.5.5 (Gene PF13_0044)
% This gene governs two reactions that do the same thing with alt
% substrates
an = singleGeneDeletion(model,'FBA',{'PF13_0044'});
assert(an == 0,'KO is not lethal')
% LETHAL

%succinate dehydrogenase can act reversibly
opt = optimizeCbModel(model);
test1 = changeRxnBounds(model,'CYOR_u6m_mt',0,'b');
%test2 = changeRxnBounds(test1,'SUCD2_u6m_mt',0,'l');
test3 = addDemandReaction(test1,'q8[m]');
test4 = changeObjective(test3,'DM_q8[m]');
opt_no_hypo = optimizeCbModel(test4);


% ubiquinone synthesis
model_bounds = changeRxnBounds(model,'OMBZLM_mt',0,'b'); % 2nd to last step
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0 ,'KO is not lethal')
% NOT LETHAL

opt = optimizeCbModel(model);
% lactate dehydrogenase 1.1.1.27 (Genes PF13_0141, PF13_0144)
model_bounds = changeRxnBounds(model,'PYRt2r',0,'b'); % remove pyruvate import
model_bounds = changeRxnBounds(model_bounds,'LDH_L',0,'b'); % lactate dehydrogenase
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0 ,'KO is not lethal')
% LETHAL

% Ornithine decarboxylase/S-Adenosylmethionine decarboxylase
% 4.1.1.17/4.1.1.50 ( Gene 'PF10_0322')
model_bounds1 = changeRxnBounds(model,'PTRCt2',0,'u'); % remove putrescine import
model_bounds = changeRxnBounds(model_bounds1,{'ADMDC';'ORNDC'},0,'b'); % S-Adenosylmethionine decarboxylase reaction KO
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f ==0,'KO is not lethal')
% LETHAL

% Spermidine synthase  2.5.1.16	(Gene 'PF11_0301')
model_bounds1 = changeRxnBounds(model,'SPMDt2',0,'u'); % remove spermidine import 
model_bounds = changeRxnBounds(model_bounds1,{'SPMS2';'SPMS'},0,'b'); % Spermidine synthase KO
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f ==0,'KO is not lethal')
% LETHAL

% Adenosine deaminase	3.5.4.4	'PF10_0289' % experiment was done without
% any purines in media except MTA (methylthioadenosine)
model_bounds1 = changeRxnBounds(model,'HYXNt',0,'u'); % remove hypoxanthine import 
model_bounds1 = changeRxnBounds(model_bounds1,'INSt',0,'u'); % remove inosine import
model_bounds = changeRxnBounds(model_bounds1,{'MTAADA';'DADA';'ADA'},0,'b'); % Adenosine deaminase KO
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f ==0,'KO is not lethal')
% NON LETHAL ins[c] cna be generated via 'NTD11' using imp; not lethal in
% Plata either

% Purine nucleoside phosphorylase	2.4.2.1	'PFE0660c'
model_bounds1 = changeRxnBounds(model,{'HYXNt'},0,'u'); % remove hypoxanthine import 
model_bounds = changeRxnBounds(model_bounds1,{'PUNP1';'PUNP2';'PUNP3';'PUNP4';'PUNP5';'PUNP6';'PUNP7';'DURIPP';'PUNP8'},0,'b'); % Purine nucleoside phosphorylase KO
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0,'KO is not lethal')
% LETHAL


%ILETRS
model_bounds = changeRxnBounds(model,{'ILETRS'},0,'b'); % inhibit reaction in all compartments 
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f==0,'KO is not lethal')
% Lethal

%ribonucleotide reductase
model_bounds = changeRxnBounds(model,{'RNDR1','RNDR2','RNDR3','RNDR4','RNDR5'},0,'b'); % inhibit reaction in all compartments 
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0,'KO is not lethal')
%LETHAL
 
%%%%%%%%%%%%%%%5
% Glutathione Reductase	1.8.1.7	'PF14_0192'
model_bounds = changeRxnBounds(model,{'GTHOr_c','GTHOr_ap'},0,'b'); % inhibit reaction in all compartments 
model_bounds = changeRxnBounds(model,{'GTHOr_c'},0,'b');
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0,'KO is not lethal')
% Not Lethal

% Glutathione s-transferase
model_bounds = changeRxnBounds(model,{'GST12_c'},0,'b'); % inhibit reactionko_opt = optimizeCbModel(model_bounds);
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0,'KO is not lethal')
% NOT Lethal

% Glutathione synthase
model_bounds = changeRxnBounds(model,{'GTHS'},0,'b'); % inhibit reaction ko_opt = optimizeCbModel(model_bounds);
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0,'KO is not lethal')
% NOT Lethal


% alpha-ketoglutarate dehydrogenase	1.2.4.2	'PF08_0045'
model_bounds = changeRxnBounds(model,{'AKGDH_mt'},0,'b'); 
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0,'KO is not lethal')
% NL

% succinate-CoA ligase	6.2.1.4	'PF11_0097'
model_bounds = changeRxnBounds(model,{'SUCOAS1m_mt','SUCOAS_mt'},0,'b'); 
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0,'KO is not lethal')
% NL

% succinate dehydrogenase	1.3.5.1	
model_bounds = changeRxnBounds(model,{'SUCD2_u6m_mt'},0,'b'); 
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f ==0,'KO is not lethal')
% NL

% citrate synthase 2.3.3.1
model_bounds = changeRxnBounds(model,{'CS_mt'},0,'b'); 
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0,'KO is not lethal')
% NL

% isocitrate dehydrogenase 1.1.1.42
model_bounds = changeRxnBounds(model,{'ICDHyr_mt'},0,'b'); 
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0,'KO is not lethal')
% NL

% acetyl-coa carboxylase 6.4.1.2
model_bounds = changeRxnBounds(model,{'ACCOAC[ap]'},0,'b'); 
ko_opt = optimizeCbModel(model_bounds);
assert(ko_opt.f == 0,'KO is not lethal')
% NL

% acquaglyceroporin
[grRatio_genes,grRateKO,~,delRxns,~] = singleGeneDeletion(model,'FBA', {'PF11_0338'});
% non lethal

[grRatio_genes,grRateKO,~,delRxns,~] = singleGeneDeletion(model,'FBA', {'PFL1870c'});
% % NL 

