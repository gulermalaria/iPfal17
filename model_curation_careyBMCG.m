%% Setup environment
initCobraToolbox % this only needs to be run once
changeCobraSolver('gurobi5');
changeCobraSolver('gurobi5','MILP')
changeCobraSolver('gurobi5','LP')

%% load plata model
cd('C:\Users\mac9jc\Documents\MATLAB\curation')
plata = readCbModel('plata_orig.xml');
model = plata; % copy plata model for modifications

%% update apicoplast localization nomenclature (necessary for rewriting file at end)
% replace '_ap' with '[ap]'
model.mets = regexprep(model.mets,'_ap\>' ,'[ap]');
model.mets = regexprep(model.mets,'\[a\]\>' ,'[ap]');
model.rxns = regexprep(model.rxns,'_ap\>' ,'[ap]');
model.rxns = regexprep(model.rxns,'\[a\]\>' ,'[ap]');

%% deleting reactions
cd('C:\Users\mac9jc\Documents\MATLAB\curation\final')
file = table2cell(readtable('reaction_delete_edits_final.xls')); 
for j = 1: (length(file(:,1)))
    % check that row has reaction ID
    if length(file{j,1}) <= 1
        disp('no reaction name')
        break
    end
    
    % remove old version of reaction (if present)
    vec = file(j,:);
    new1 = removeRxns_old(model, vec{1}); 
    model = new1;
end
clearvars vec j new1 file

%% add/modify reactions 
file = table2cell(readtable('edits_additions_or_modifications_final.xls'));
for j = 1:(length(file(:,1)))
    % check that row has reaction ID
    if length(file{j,1}) <= 1
        disp('no reaction name')
        break
    end
    
    % make edit
    vec = file(j,:);
    new1 = removeRxns_old(model, vec{1}); % remove old version
    % add back correct reaction
    if vec{8} == 0 % reaction is not reversible
        revFlag = 'false'; 
    else % reaction is reversible
        revFlag = 'true'; 
    end
    lb = vec{9}; ub = vec{10};
    rxnName = vec{1}; metaboliteList = vec{3}; objCoeff = vec{11};
    subSystem = vec{7}; grRule = vec{4}; %geneNameList = vec(5);
    stoichCoeffList= {};
    if isnan(grRule) % if no genes, switch from nan to empty cell
        grRule = {};
    end
    new = addReaction(new1,rxnName,...
        metaboliteList,stoichCoeffList,revFlag,lb,ub,objCoeff,...
        subSystem,grRule); %,geneNameList,systNameList,checkDuplicate)
    model = new;
    indx = (strcmp(rxnName,model.rxns));
    reaction_ec = vec(13); reaction_notes = vec(14); reaction_ref = vec(15);
    model.rxnECNumbers(indx) = reaction_ec;
    model.rxnNotes(indx) = reaction_notes;
    model.rxnReferences(indx) = reaction_ref;
end

clearvars vec j new1 file indx reaction_ec reaction_notes new
clearvars ub lb revFlag rxnName stoichCoeffList objCoeff metaboliteList
clearvars  reaction_ref subSystem grRule

%% biomass modification
file = table2cell(readtable('biomass_edits_final.xls'));
for j = 1:length(file(:,1))
    % check that row has reaction ID
    if length(file{j,1}) <= 1
        disp('no reaction name')
        break
    end
    
    % make edit
    vec = file(j,:);
    new1 = removeRxns_old(model, vec{1}); % remove old version
    
    % add back correct reaction
    if vec{8} == 0 % reaction is not reversible
        revFlag = 'false'; 
    else % reaction is reversible
        revFlag = 'true'; 
    end
    lb = vec{9}; ub = vec{10};
    rxnName = vec{1}; metaboliteList = vec{3}; objCoeff = vec{11};
    subSystem = vec{7}; grRule = vec{4}; %geneNameList = vec(5);
    stoichCoeffList= {};
    if isnan(grRule) % if no genes, switch from nan to empty cell
        grRule = {};
    end
    new = addReaction(new1,rxnName,...
        metaboliteList,stoichCoeffList,revFlag,lb,ub,objCoeff,...
        subSystem,grRule); %,geneNameList,systNameList,checkDuplicate)
    model = new;
    indx = (strcmp(rxnName,model.rxns));
    reaction_ec = vec(13); reaction_notes = vec(14); reaction_ref = vec(15);
    model.rxnECNumbers(indx) = reaction_ec;
    model.rxnNotes(indx) = reaction_notes;
    model.rxnReferences(indx) = reaction_ref;
end

clearvars vec j new1 file indx reaction_ec reaction_notes new
clearvars ub lb revFlag rxnName stoichCoeffList objCoeff metaboliteList
clearvars  reaction_ref subSystem grRule

% check model numbers
[r4, ~] = size(model.rxns);
[g4, ~] = size(model.genes);
j = 0; %ECs
for i = 1:length(model.rxnECNumbers)
    if length(model.rxnECNumbers{i}) >1
        j = j+1;
    end
end
EC4 = j; % number of EC numbers in model

clearvars j vec i 

%% a couple additional edits easier here for formatting
model = addReaction(model,'EX_lac_L(e)',{'lac_L[e]'},(-1),false);
model.lb(strcmp('EX_pyr(e)',model.rxns)) = -1000;
model = changeRxnBounds(model,'SUCD2_u6m_mt',0,'l');
model.mets = regexprep(model.mets,'_ap\>' ,'[ap]');

%% check duplicate reactions and growth before saving
[model,removed1] = checkDuplicateRxn(model,1);
if ~isempty(removed1)
    warning('duplicate reactions, check removed1')
end
[model,removed2] = checkDuplicateRxn(model,2);
if ~isempty(removed2)
    warning('duplicate reactions, check removed2')
end

clearvars removed1 removed2

opt = optimizeCbModel(model,'max',0,0);
if opt.f < .1
    warning('fba of cobra <.1')
end
% [missing_mets,present_mets] = biomassPrecursorCheck(model);
% Doesnt work because can't add dm reactions for mets that have exchange
% reactions

%% convert to tiger model and check
addpath(genpath('C:/Users/mac9jc/Documents/MATLAB/tiger/'));
start_tiger
set_solver('gurobi');
set_solver_option('MaxTime',60*60); % max time a TIGER simulation is allowed to run
set_solver_option('IntFeasTol',1e-8); % cutoff number for interpreting a value as 0

pf_tiger = cobra_to_tiger(model, 'add_gpr','true','fast_gpr','false','reactions_only','false'); % make cobra model into tiger model
opt2 = fba(pf_tiger);
if opt2.val < .1
    warning('fba of tiger model <.1')
end

%% SAVE MODEL AS MATLAB
cd('C:\Users\mac9jc\Documents\MATLAB\curation\final')
save model_cobra_edited.mat model;
% Convert each COBRA model to a TIGER model % DO NOT DO THIS- PROMOTES ROUNDING ERRORS
% model_tiger_edited = cobra_to_tiger(model_cobra_edited);

%% SAVE AS XML 
writeCbModel(model,'sbml','2017_pf_model.xml',{'c','e','fv','m','ap'},...
   {'cytoplasm','extracellular','food_vacuole','mitochondria','apicoplast'} );

% %% save model without localization if needed for metdraw purposes
% model_wo_loc = model
% model_wo_loc.mets = regexprep(model_wo_loc.mets, '\[(m)\]', '\[c\]')
% model_wo_loc.mets = regexprep(model_wo_loc.mets, '\[(ap)\]', '\[c\]')
% model_wo_loc.metNames = regexprep(model_wo_loc.mets, '\[(m)\]', '\[c\]')
% model_wo_loc.mets = regexprep(model_wo_loc.mets, '\[(fv)\]', '\[c\]')
% model_wo_loc.metNames = regexprep(model_wo_loc.mets, '\[(fv)\]', '\[c\]')
% model_wo_loc.metNames = regexprep(model_wo_loc.mets, '\[(ap)\]', '\[c\]')
% model_wo_loc.subSystems = regexprep(model_wo_loc.subSystems,';','')
% model_wo_loc.subSystems = regexprep(model_wo_loc.subSystems,'\s(\w+)','')
% writeCbModel(model_wo_loc,'sbml','2017_pf_model_wo_loc.xml',{'c','e'},...
%    {'cytoplasm','extracellular'} );

 [~, genes_used_by_model] = findUsedGenes(model, model.genes);
 genes_used_by_model = table(genes_used_by_model);
 writetable(genes_used_by_model,'genes_used_by_model.csv') 
