
pf_cobra = model; % from model_curation_careyBMCSB

%% import expression data from individual countries for MADE 
gene_expression_cambodia = readtable('cambodia.csv','Delimiter',',');
gene_expression_vietnam = readtable('vietnam.csv','Delimiter',',');

%% Identify gene expression for integration
gene_in_pf_tiger_cambodia = ismember(gene_expression_cambodia.model_ORF, pf_tiger.genes);
gene_in_pf_tiger_vietnam = ismember(gene_expression_vietnam.model_ORF, pf_tiger.genes);
% 
% gene_left_out = ~ismember(pf_tiger.genes,gene_expression_cambodia.model_ORF);
% no_gene_expression_made_cambodia = pf_tiger.genes(gene_left_out,:);
% %     'COI' not in model
% %     'COXIII' not in model anymore
% %     'MAL13P1_56' % must swap . for _
% %     'PF14_0534' not in model anymore
% %     'PF07_0009' not in model anymore
% %     'PF10_0215' not in model anymore
% %     'PF13_0352' not in model anymore
% %     'PF14_0401' not in model anymore
% %     'PFL2510w' not in model anymore
% %     'PF3D7_0707100' not in microarray 
% %     'PF3D7_1438200' not in microarray 
% %     'PF3D7_0702800' not in microarray 
% %     'PF10_0409' not in microarray
% %     'mal_mito_2' not in microarray 
% %     'mal_mito_1' not in microarray 
% %     'MAL13P1.390' % must swap . for _

gene_expression_made_cambodia = gene_expression_cambodia(gene_in_pf_tiger_cambodia,:);
gene_expression_made_vietnam = gene_expression_vietnam(gene_in_pf_tiger_vietnam,:);
fdr_threshold = 0.5; % FDR threshold: above this, FC assumed to be none
clearvars gene_in_pf_tiger_cambodia gene_in_pf_tiger_vietnam 
clearvars gene_expression_cambodia gene_expression_vietnam

%% RUN MADE with growth threshold 80%

    tic
    pf_made_cambodia80 = made(...
        pf_tiger,... % model
        gene_expression_made_cambodia.logFC,... % fold change (log)
        gene_expression_made_cambodia.P_Value,... % p values
        'gene_names',gene_expression_made_cambodia.model_ORF,...
        'obj_frac', 0.8, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_vietnam80 = made(...
        pf_tiger,...
        gene_expression_made_vietnam.logFC,...
        gene_expression_made_vietnam.P_Value,...
        'gene_names',gene_expression_made_vietnam.model_ORF,...
        'obj_frac', 0.8, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    
%% RUN MADE with growth threshold 70%

    tic
    pf_made_cambodia70 = made(...
        pf_tiger,... % model
        gene_expression_made_cambodia.logFC,... % fold change (log)
        gene_expression_made_cambodia.P_Value,... % p values
        'gene_names',gene_expression_made_cambodia.model_ORF,...
        'obj_frac', 0.7, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_vietnam70 = made(...
        pf_tiger,...
        gene_expression_made_vietnam.logFC,...
        gene_expression_made_vietnam.P_Value,...
        'gene_names',gene_expression_made_vietnam.model_ORF,...
        'obj_frac', 0.7, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    
%% RUN MADE with growth threshold 60%

    tic
    pf_made_cambodia60 = made(...
        pf_tiger,... % model
        gene_expression_made_cambodia.logFC,... % fold change (log)
        gene_expression_made_cambodia.P_Value,... % p values
        'gene_names',gene_expression_made_cambodia.model_ORF,...
        'obj_frac', 0.6, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_vietnam60 = made(...
        pf_tiger,...
        gene_expression_made_vietnam.logFC,...
        gene_expression_made_vietnam.P_Value,...
        'gene_names',gene_expression_made_vietnam.model_ORF,...
        'obj_frac', 0.6, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
        
%% RUN MADE with growth threshold 50%

    tic
    pf_made_cambodia50 = made(...
        pf_tiger,... % model
        gene_expression_made_cambodia.logFC,... % fold change (log)
        gene_expression_made_cambodia.P_Value,... % p values
        'gene_names',gene_expression_made_cambodia.model_ORF,...
        'obj_frac', 0.5, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_vietnam50 = made(...
        pf_tiger,...
        gene_expression_made_vietnam.logFC,...
        gene_expression_made_vietnam.P_Value,...
        'gene_names',gene_expression_made_vietnam.model_ORF,...
        'obj_frac', 0.5, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    
%% RUN MADE with growth threshold 40%

    tic
    pf_made_cambodia40 = made(...
        pf_tiger,... % model
        gene_expression_made_cambodia.logFC,... % fold change (log)
        gene_expression_made_cambodia.P_Value,... % p values
        'gene_names',gene_expression_made_cambodia.model_ORF,...
        'obj_frac', 0.4, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_vietnam40 = made(...
        pf_tiger,...
        gene_expression_made_vietnam.logFC,...
        gene_expression_made_vietnam.P_Value,...
        'gene_names',gene_expression_made_vietnam.model_ORF,...
        'obj_frac', 0.4, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc    

    
%% RUN MADE with growth threshold 30%

    tic
    pf_made_cambodia30 = made(...
        pf_tiger,... % model
        gene_expression_made_cambodia.logFC,... % fold change (log)
        gene_expression_made_cambodia.P_Value,... % p values
        'gene_names',gene_expression_made_cambodia.model_ORF,...
        'obj_frac', 0.3, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_vietnam30 = made(...
        pf_tiger,...
        gene_expression_made_vietnam.logFC,...
        gene_expression_made_vietnam.P_Value,...
        'gene_names',gene_expression_made_vietnam.model_ORF,...
        'obj_frac', 0.3, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    
%% clear extra variables    
clearvars gene_expression_cambodia_n gene_expression_made_cambodia
clearvars gene_expression_made_vietnam gene_expression_vietnam_n
clearvars gene_left_out gene_in_pf_tiger_cambodia gene_in_pf_tiger_vietnam
clearvars fdr_threshold obj_value_desired no_gene_expression_made_cambodia
clearvars no_gene_expression_made_cambodia pf_obj_frac pf_obj_flux pf_obj_flux_result
clearvars required_growth 

%% confirm growth of each model
c_files = {pf_made_cambodia30;pf_made_cambodia40;pf_made_cambodia50;...
    pf_made_cambodia60;pf_made_cambodia70;pf_made_cambodia80};
v_files = {pf_made_vietnam30;pf_made_vietnam40;pf_made_vietnam50;...
    pf_made_vietnam60;pf_made_vietnam70;pf_made_vietnam80};

tiger_fba_results = cell(5,7,1); %preallocate results
tiger_fba_results{1,1} = 'Model';
tiger_fba_results{1,2} = '30'; tiger_fba_results{1,3} = '40';
tiger_fba_results{1,4} = '50'; tiger_fba_results{1,5} = '60';
tiger_fba_results{1,6} = '70'; tiger_fba_results{1,7} = '80';
tiger_fba_results{2,1} = 'c_res'; tiger_fba_results{3,1} = 'c_sens';
tiger_fba_results{4,1} = 'v_res'; tiger_fba_results{5,1} = 'v_sens';

for i = 1:6
    c_models = c_files{i};
    c_res = fba(c_models.models{1,1});
    c_sens = fba(c_models.models{1,2});
    v_models = v_files{i};
    v_res = fba(v_models.models{1,1});
    v_sens = fba(v_models.models{1,2});
    tiger_fba_results{2,(i+1)} = c_res.val;
    tiger_fba_results{3,(i+1)} = c_sens.val;
    tiger_fba_results{4,(i+1)} = v_res.val;
    tiger_fba_results{5,(i+1)} = v_sens.val;
end

clearvars v_models c_models c_res c_sens i v_res v_sens c_files v_files 

%% Save Cambodia MADE results for COBRA
c_files = {pf_made_cambodia30;pf_made_cambodia40;pf_made_cambodia50;...
pf_made_cambodia60;pf_made_cambodia70;pf_made_cambodia80};
v_files = {pf_made_vietnam30;pf_made_vietnam40;pf_made_vietnam50;...
pf_made_vietnam60;pf_made_vietnam70;pf_made_vietnam80};

gene_states = cell(2,6);
for i = 1:6
    c_models = c_files{i};
    v_models = v_files{i};
    % row 1 of cell arrays = all cambodia
    c_table = table(c_models.genes, c_models.gene_states(:,1),c_models.gene_states(:,2),'VariableNames',...
        {'genes','res_state','sens_state'});
    gene_states{1,i} = [c_table.Properties.VariableNames; table2cell(c_table)];
    % row 2 of cell arrays = all vietnam
    v_table = table(v_models.genes, v_models.gene_states(:,1),v_models.gene_states(:,2),'VariableNames',...
        {'genes','res_state','sens_state'});
    gene_states{2,i} = [v_table.Properties.VariableNames; table2cell(v_table)];
    % columns = growth thresholds
end
clearvars i

cobra_flux = cell(4,6); % each column = growth threshold
% rows = c_res, c_sens, v_res, v_sens
for i = 1:6 % for each threshold
    % identify genes to delete (if gene state == 0)
    res_index = (cell2mat(gene_states{1,i}(2:end,2)) == 0);
    res_index = vertcat(logical(0), res_index); % shift as first element is header
    c_deleted_res = gene_states{1,i}(res_index,1);
    res_index = (cell2mat(gene_states{1,i}(2:end,3)) == 0);
    res_index = vertcat(logical(0), res_index);
    c_deleted_sens = gene_states{1,i}(res_index,1);
    res_index = (cell2mat(gene_states{2,i}(2:end,2)) == 0);
    res_index = vertcat(logical(0), res_index);
    v_deleted_res = gene_states{2,i}(res_index,1);
    res_index = (cell2mat(gene_states{2,i}(2:end,3)) == 0);
    res_index = vertcat(logical(0), res_index);
    v_deleted_sens = gene_states{2,i}(res_index,1);
    
    % delete genes from model
    c_res = deleteModelGenes(pf_cobra,c_deleted_res);
    c_sens = deleteModelGenes(pf_cobra,c_deleted_sens);
    v_res = deleteModelGenes(pf_cobra,v_deleted_res);
    v_sens = deleteModelGenes(pf_cobra,v_deleted_sens);
    
    %predict flux once genes are deleted
    c_res_flux = optimizeCbModel(c_res);
    c_sens_flux = optimizeCbModel(c_sens);
    v_res_flux = optimizeCbModel(v_res);
    v_sens_flux = optimizeCbModel(v_sens);

    % store data
    cobra_flux{1,i} = {c_res c_res_flux.f};
    cobra_flux{2,i} = {c_sens c_sens_flux.f};
    cobra_flux{3,i} = {v_res v_res_flux.f};
    cobra_flux{4,i} = {v_sens v_sens_flux.f};
end

for i = 1:6
    for j = 1:4
    flux = cobra_flux{j,i}{2};
        if (flux < .1)
            disp(i)
            disp(j)
            warning('Cobra flux == 0')
        end
        if (isempty(flux))
            disp(i)
            disp(j)
            warning('Infeasible solution')
        end
    end
end
clearvars  c_deleted_res c_deleted_sens v_deleted_res v_deleted_sens c_res
clearvars v_res v_sens c_sens c_files c_table i c_res_flux c_models c_sens_flux
clearvars c_res_flux_result res_index v_files v_models v_table
clearvars v_sens_flux v_res_flux
clearvars pf_made_*

%% Rxn KO studies
rxn_KO = cell(4,6); % each column = growth threshold
rxn_infeas = cell(4,6);
% rows = c_res, c_sens, v_res, v_sens
for i = 6%1:6
    disp(i)
    for j = 1:4
        disp(j)
        model = cobra_flux{j,i}{1};
        [grRatio,~,~,~,~,~] = singleRxnDeletion(model,'FBA',model.rxns);
        lethalKO = model.rxns(grRatio < 0.1);
        infeas = model.rxns(isnan(grRatio));
        rxn_KO{j,i} = lethalKO;
        rxn_infeas{j,i} = infeas;
    end
end
rxn_KO_80 = rxn_KO(:,6);
rxn_Infeas_80 = rxn_infeas(:,6);
[m1, ~] = size(rxn_Infeas_80{1});[m2, ~] = size(rxn_Infeas_80{2});
[m3, ~] = size(rxn_Infeas_80{3});[m4, ~] = size(rxn_Infeas_80{4});
if (m1+m2+m3+m4)>0
    warning('infeasible KOs')
end
clearvars m1 m2 m3 m4 rxn_infeas

% 80 consensus res essential
consensus_res80 = intersect(rxn_KO_80{1},rxn_KO_80{3});
% 80 consensus sens essential
consensus_sens80 = intersect(rxn_KO_80{2},rxn_KO_80{4});
% 80 unique res essential
unique_res80 = setdiff(consensus_res80,consensus_sens80);
% 80 unique sens essential
unique_sens80 = setdiff(consensus_sens80,consensus_res80);
%essential to all
all_80 = intersect(consensus_sens80,consensus_res80);

%% get reactions EC and formula for lethal rxn KOs
res = cell(length(unique_res80),5); res(:,1)= unique_res80;
sens = cell(length(unique_sens80),5); sens(:,1) = unique_sens80;
all = cell(length(all_80),5); all(:,1) = all_80;
model1 = cobra_flux{1,6}{1}; % 80 threshold cambodia resistant model
for i = 1:length(unique_res80)
    res{i,2} = printRxnFormula(model1,unique_res80(i));
    indx = strcmp(unique_res80(i), model1.rxns);
    res{i,3} = model1.rxnECNumbers(indx);
    res{i,4} = model1.subSystems(indx);
    genes = findGenesFromRxns(model1,unique_res80(i));
    gene1 = '';
    for j = 1:length(genes{1,1})
        gene1 = strcat(gene1,{' '},genes{1,1}(j));
    end
    res{i,5} = gene1;
end
model2 = cobra_flux{2,6}{1}; % 80 threshold cambodia sens model
for i = 1:length(unique_sens80)
    sens{i,2} = printRxnFormula(model2,unique_sens80(i));
    indx = strcmp(unique_sens80(i), model2.rxns);
    sens{i,3} = model2.rxnECNumbers(indx);
    sens{i,4} = model2.subSystems(indx);
    genes = findGenesFromRxns(model2,unique_sens80(i));
    gene1 = '';
    for j = 1:length(genes{1,1})
        gene1 = strcat(gene1,{' '},genes{1,1}(j));
    end
    sens{i,5} = gene1;
end
model = pf_cobra;
% use original model
for i = 1:length(all_80)
    all{i,2} = printRxnFormula(model,all_80(i));
    indx = strcmp(all_80(i), model.rxns);
    all{i,3} = model.rxnECNumbers(indx);
    all{i,4} = model.subSystems(indx);    
    genes = findGenesFromRxns(model,all_80(i));
    gene1 = '';
    for j = 1:length(genes{1,1})
        gene1 = strcat(gene1,{' '},genes{1,1}(j));
    end
    all{i,5} = gene1;
end
clearvars i flux gene1 indx infeas model1 model2 opt2 j EC4 g4 r4
% all, res, sens SUPPLE T 5, T5 & T6
[m, ~] = size(res)
for i = 1:m
    if iscell(res{i,5})
        res(i,5) = res{i,5};
    end
end
[m, ~] = size(sens)
for i = 1:m
    if iscell(sens{i,5})
        sens(i,5) = sens{i,5};
    end
end
[m, ~] = size(all)
for i = 1:m
    if iscell(all{i,5})
        all(i,5) = all{i,5};
    end
end
res = cell2table(res,'VariableNames',{'Reactions','Formula','EC','Subsystems','Genes'});
writetable(res,'uniqueResEssentialRxns.xls') 
sens = cell2table(sens,'VariableNames',{'Reactions','Formula','EC','Subsystems','Genes'});
writetable(sens,'uniqueSensEssentialRxns.xls') 
all = cell2table(all,'VariableNames',{'Reactions','Formula','EC','Subsystems','Genes'});
writetable(all,'consensus_allEssentialRxns.xls') 
clearvars all res sens

%% fva
c_res = cobra_flux{1,6}{1};
v_res = cobra_flux{3,6}{1};
c_sens = cobra_flux{2,6}{1};
v_sens = cobra_flux{4,6}{1};
c_r_opt = optimizeCbModel(c_res);
v_r_opt = optimizeCbModel(v_res);
c_s_opt = optimizeCbModel(c_sens);
v_s_opt = optimizeCbModel(v_sens);


fva_results = struct('c_res_min',[],'c_res_max',[],...
    'c_sens_min',[],'c_sens_max',[],...
    'v_res_min',[],'v_res_max',[],...
    'v_sens_min',[],'v_sens_max',[]);
[fva_results.c_res_min,fva_results.c_res_max] = fluxVariability(c_res);
[fva_results.c_sens_min,fva_results.c_sens_max] = fluxVariability(c_sens);
[fva_results.v_res_min,fva_results.v_res_max] = fluxVariability(v_res);
[fva_results.v_sens_min,fva_results.v_sens_max] = fluxVariability(v_sens);
fva_results.reactions = pf_cobra.rxns;
fva_results2 = struct2table(fva_results);
% use in figures
writetable(fva_results2,'fva_results.xls') 
clearvars all res sens

clearvars c_res v_res v_res_flux v_s_opt c_sens v_sens v_sens_flux v_s_opt
clearvars v_r_opt fva_results i j c_r_opt c_res_flux c_s_opt c_sens_flux
clearvars minFlux maxFlux

%% enrichment prep
model = pf_cobra;
gene_states_80 = [gene_states{1,6},gene_states{2,6}(:,2:3)];
gene_states_80 = gene_states_80(2:end,:);
[m, ~] = size(gene_states_80);
[ gene_r, genes_used] = findUsedGenes(model, model.genes);
react_sub = cell(m,3,1);
for j = 1:m % for each gene
    react_sub{j,1} = gene_states_80{j,1};
    gene = react_sub{j,1};
    if ~ismember(gene,genes_used)
        continue
    end
    reactions = struct2cell(findRxnsFromGenes(model,gene));
    [q,~] = size(reactions{1});
    if q==0
        continue
    end
    reactions_list = [];
    subsystem_list = [];
    for p = 1:q
        reactions_list= strcat(reactions_list,' -',reactions{1}{p,1});
        subsystem_list = strcat(subsystem_list,' -',reactions{1}{p,3});
    end
    react_sub{j,2} = reactions_list;
    react_sub{j,3} = subsystem_list;
end
react_sub = cell2table(react_sub,'VariableNames',{'Gene','Reactions','Subsystems'});
writetable(react_sub,'GeneRxnSubsystems.xls') %  SEE R FILE FOR COPY PASTE MODIFICATIONS
clearvars p q reactions_list subsystem_list react_sub reactions gene genes_used
clearvars m n gene_r gene_states_80

%% figure 2 prep
model = pf_cobra; rxns = model.rxns; sub = model.subSystems; 
ref = model.rxnReferences; not = model.rxnNotes;
[m,n] = size(rxns); g_rxn = cell(m,1,1);

for j = 1:m % for each gene
    disp(j)
    an = findGenesFromRxns(model,rxns{j});
    [q,~] = size(an{1});
    if q==0
        g_rxn{j,1} = [];
        continue
    else
        gene_list = [];
        for p = 1:q
            gene_list= strcat(gene_list,' -',an{1}{p,1});
        end
        g_rxn{j,1} = gene_list;
    end 
end
react_sub = table(rxns,sub,ref,not,g_rxn,'VariableNames',{'Reactions',...
    'Subsystems','References','Notes','Genes'});
writetable(react_sub,'figure2_prep.xls') %  replace '-' with spaces
clearvars not ref q gene_list p g_rxn rxns sub m n

%% get reactions EC and formula for flux enrichment
model = pf_cobra; %c_res, c_sens, v_res, v_sens
c_res_flux = optimizeCbModel(cobra_flux{1,6}{1,1}); c_res_flux = ...
    c_res_flux.x;
c_sens_flux = optimizeCbModel(cobra_flux{2,6}{1,1}); c_sens_flux = ...
    c_sens_flux.x;
v_res_flux = optimizeCbModel(cobra_flux{3,6}{1,1}); v_res_flux = ...
    v_res_flux.x;
v_sens_flux = optimizeCbModel(cobra_flux{4,6}{1,1}); v_sens_flux = ...
    v_sens_flux.x;
fluxes = table(model.rxns, c_res_flux,c_sens_flux,v_res_flux,v_sens_flux,...
    model.subSystems,'VariableNames',...
    {'Rxn','c_res','c_sens','v_res','v_sens','subsystems'});

% add subsystem to flux
writetable(fluxes,'fluxes_paper.csv') % DONT replace ';'

%% gene states
gene_states 
    % row 1 of cell arrays = all cambodia
    % row 2 of cell arrays = all vietnam
    % columns = growth thresholds
gene_states_print_c = gene_states{1,6}
gene_states_print_v = gene_states{2,6}
writetable(cell2table(gene_states_print_c),'gene_states_c.csv') 
    % delete first row in these files (nonsense header)
writetable(cell2table(gene_states_print_v),'gene_states_v.csv') 

