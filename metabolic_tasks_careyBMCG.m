
% this script runs a series of metabolic tasks on a cobra model structure
% if a task fails, the model is unable to perform behaviors observed in
% vitro or in vivo. the tasks include: KO reactions (lethal KOs observed
% experimentally should be recapitulated in silico), ATP production (model
% should not produce ATP if no import is permitted), purine necessity
% (except hypoxathine, purines are not supplemneted in vitro, but are present in vivo, the
% parasite should grow with and without purines), sugar requirement (while
% it is unclear in the literature, Pf likely cannot grow unless glucose and
% or fructose are present), production or import of riboflavin, nicotinamide, pyridoxine, and thiamine

opt = optimizeCbModel(model);
lethal = 0;
frac = opt.f*.9; 

%% Do antimetabolites of riboflavin, nicotinamide, pyridoxine, and thiamine 
% inhibit growth? (Nutritional requirements of Plasmodium falciparum 
% in culture. II. Effects of antimetabolites in a semi-defined medium.
% Geary TG, Divo AA, Jensen JB

% ribo
% findRxnsFromMets(model,'ribflv[c]')
ribo = removeRxns(model,{'RBFK';'RIBFLVt2';'ACP1'});
opt = optimizeCbModel(ribo);
if ~(opt.f < frac)
    warning('Model grows with riboflavin antimetabolite')
end
% PASSES

% nico
% findRxnsFromMets(model,'ncam[c]')
nico = removeRxns(model,{'NNAM';'EX_nicotinamide2'});
opt = optimizeCbModel(nico);
if ~(opt.f < frac)
    warning('Model grows with nicotinamide antimetabolite')
end
% Fails, but there are many ways to make NAD+ and NADP+
% off target affects of ncam antimetabolite, because NAD is similiar in
% structure?

% pyro
% findRxnsFromMets(model,'pydxn[c]')
pyro = removeRxns(model,{'PYDXNK';'PDXPP';'PYDXNtr'});
opt = optimizeCbModel(pyro);
if ~(opt.f < frac)
    warning('Model grows with pyridoxine antimetabolite')
end
% Fails, maybe antimetabolite has off target effects?

% thi
% findRxnsFromMets(model,'thm[c]')
thi = removeRxns(model,{'TMDPK';'THMP';'THMDP';'THMt3'});
opt = optimizeCbModel(thi);
if ~(opt.f < frac)
    warning('Model grows with thiamine antimetabolite')
end
% FAILS

%% Can model produce ATP if no exchange is allowed
exchange_rxns = model.rxns(findExcRxns(model));
atp_pre = changeRxnBounds(model,exchange_rxns,0,'u');
atp = addDemandReaction(atp_pre,'atp[c]');
atp = changeObjective(atp,'DM_atp[c]');
opt = optimizeCbModel(atp);
if ~(opt.f < frac)
    warning('Model produces ATP when all import reactions are blocked')
end
% PASSES

%% Media purine composition
% Can parasite grow when hypoxanthine is the only purine imported? 
pur = changeRxnBounds(model,{'MTAADA';'ADEt';'dIMP_t';'INSt';'ADNt';...
    'GUAt';'DADNt4';'DGSNt';'DINt';'GSNt';'PAPt';'XANt'},0,'u');
opt = optimizeCbModel(pur);
if ~(opt.f > (lethal))
    warning('Model does not grow in in vitro conditions with hypoxanthine as the sole purine')
end
% PASSES

pur1 = changeRxnBounds(pur,'HYXNt',0,'b');
opt_no_hypo = optimizeCbModel(pur1);
if ~(opt_no_hypo.f < (frac))
    warning('One purine (hypoxanthine) is not necessary for growth')
end
% Passes

% Fails if MTAADA is not included in the RxnBounds list becuase of MTAADA 
% generation of hypoxanthine from spermidine (which is not in RPMI 1640 
% media, and is in low concentrations in human serum; thus, we've removed 
% it from our test 

% What if only Adenine is available?
pur = changeRxnBounds(model,{'MTAADA';'HYXNt';'dIMP_t';'INSt';'ADNt';...
    'GUAt';'DADNt4';'DGSNt';'DINt';'GSNt';'PAPt';'XANt'},0,'u');
opt = optimizeCbModel(pur);
if ~(opt.f > (lethal))
    warning('Model does not grow in in vitro conditions with adenine as the sole purine')
end
pur1 = changeRxnBounds(pur,'ADEt',0,'b');
opt_no_hypo = optimizeCbModel(pur1);
if ~(opt_no_hypo.f == (lethal))
    warning('Model does not grow in in vitro conditions with adenine as the sole purine')
end
% PASSES

% What if only Guanine is available?
pur = changeRxnBounds(model,{'MTAADA';'HYXNt';'ADEt';'dIMP_t';'INSt';'ADNt';...
    'DADNt4';'DGSNt';'DINt';'GSNt';'PAPt';'XANt'},0,'u');
opt = optimizeCbModel(pur);
if ~(opt.f > (lethal))
    warning('Model does not grow in in vitro conditions with guanine as the sole purine')
end
% FAILS
pur1 = changeRxnBounds(pur,'GUAt',0,'b');
opt_no_hypo = optimizeCbModel(pur1);
if ~(opt_no_hypo.f == (lethal))
    warning('Model does not grow in in vitro conditions with guanine as the sole purine')
end
% PASSES

% What if only Inosine is available?
pur = changeRxnBounds(model,{'MTAADA';'HYXNt';'ADEt';'dIMP_t';'ADNt';...
    'GUAt';'DADNt4';'DGSNt';'DINt';'GSNt';'PAPt';'XANt'},0,'u');
opt = optimizeCbModel(pur);
if ~(opt.f > (lethal))
    warning('Model does not grow in in vitro conditions with Inosine as the sole purine')
end
pur1 = changeRxnBounds(pur,'INSt',0,'b');
opt_no_hypo = optimizeCbModel(pur1);
if ~(opt_no_hypo.f == (lethal))
    warning('One purine (Inosine) is not necessary for growth')
end
% PASSES

% What if only Adenosine is available?
pur = changeRxnBounds(model,{'MTAADA';'HYXNt';'ADEt';'dIMP_t';'INSt';...
    'GUAt';'DADNt4';'DGSNt';'DINt';'GSNt';'PAPt';'XANt'},0,'u');
opt = optimizeCbModel(pur);
if ~(opt.f > (lethal))
    warning('Model does not grow in in vitro conditions with Adenosine as the sole purine')
end
pur1 = changeRxnBounds(pur,'ADNt',0,'b');
opt_no_hypo = optimizeCbModel(pur1);
if ~(opt_no_hypo.f == (lethal))
    warning('One purine (Adenosine) is not necessary for growth')
end
% PASSES

% What if only Guanosine is available?
pur = changeRxnBounds(model,{'MTAADA';'HYXNt';'ADEt';'dIMP_t';'INSt';...
    'GUAt';'DADNt4';'DGSNt';'DINt';'ADNt';'PAPt';'XANt'},0,'u');
opt = optimizeCbModel(pur);
if ~(opt.f > (lethal))
    warning('Model does not grow in in vitro conditions with Guanosine as the sole purine')
end
% FAILS
pur1 = changeRxnBounds(pur,'GSNt',0,'b'); 
opt_no_hypo = optimizeCbModel(pur1);
if ~(opt_no_hypo.f == (lethal))
    warning('One purine (Guanosine) is not necessary for growth')
end
% PASSES

% isopentenyl pyrophosphate (the product of the isoprenoid pathway) is
% essential ( http://dx.doi.org/10.1371/journal.pbio.1001138 )
% eliminate all other ap products?
ipp = removeRxns(model,'DMPPS_ap'); % remove original reaction that prpoduces ipp (dmpp_ap)
indx = strfind(model.mets,'[ap]');
indx = find(~cellfun('isempty', indx));
ap_mets = model.mets(indx);
apicoplast_rxns = findRxnsFromMets(ipp,ap_mets);
ipp = removeRxns(ipp,apicoplast_rxns);
opt = optimizeCbModel(ipp);
ipp_dmpp = addReaction(ipp,'dmpp_supplementation','dmpp[c] ->');
opt_dmpp = optimizeCbModel(ipp_dmpp);
if ~(opt_dmpp.f > (lethal))
    warning('Model fails to grow with no apicoplast and ipp supplementation')
end
% FAILS

% Can model grow on glucose as sole sugar source?
%sugars = maltose (malt) fructose (f6p, fdp, fru) mannose (gdpmann, man1p
%man6p, man) glucose (g1p, g6p, glc_D, udpg)  
test = changeRxnBounds(model,'GLCt1r',0,'u');
opt = optimizeCbModel(test);
if ~(opt.f == (lethal))
    warning('Alternative sugars can replace glucose as sole sugar source (incorrectly)')
end
% FAILS, many reversible reactions that can produce downstream mets

% Can model grow without isoleucine in media?
ile = changeRxnBounds(model,'ILEt2r',0,'u');
opt = optimizeCbModel(ile);
if ~(opt.f == (lethal))
    warning('Model grows without isoleucine')
end
% Passes

% Can model growth without p-aminobenzoic acid?
opt = optimizeCbModel(removeRxns(model,'DHPS2'));
if ~(opt.f == (lethal))
    warning('Model grows without p-aminobenzoic acid')
end
% PASSES

% Is growth reduced without tyrosine supplementation?
tyr = changeRxnBounds(model,'TYRt2r',0,'u');
opt = optimizeCbModel(tyr);
if ~(opt.f > (frac))
    warning('tyrosine required for growth')
end
if ~(opt.f < frac)
    warning('Growth is not reduced without tyrosine')
end
% fails

% Is growth reduced without methionine supplementation?
met = changeRxnBounds(model,{'METt2r';'METLEUex'},0,'u');
opt = optimizeCbModel(met);
if ~(opt.f > frac)
    warning('methionine required for growth')
end
if ~(opt.f < frac)
    warning('Growth is not reduced without methionine')
end
% fails

% Is growth reduced without proline supplementation?
pro = changeRxnBounds(model,'PROt2r',0,'u');
opt = optimizeCbModel(pro);
if ~(opt.f > frac)
    warning('Proline required for growth')
end
if ~(opt.f < frac)
    warning('Growth is not reduced without proline')
end
% fails

% Is growth reduced without glutamate supplementation?
glutamate = changeRxnBounds(model,'GLUt2r',0,'u');
opt = optimizeCbModel(glutamate);
if ~(opt.f > frac)
    warning('glutamate required for growth')
end
if ~(opt.f < frac)
    warning('Growth is not reduced without glutamate')
end
% fails

% Is growth reduced without glutamine supplementation? 
glutamine = changeRxnBounds(model,'GLNt2r',0,'u');
glutamine = changeRxnBounds(model,{'THRGLNexR';'SERGLNexR';'ALAGLNexR';'CYSGLUexR'},0,'l');
opt = optimizeCbModel(glutamine);
if ~(opt.f > frac)
    warning('glutamine required for growth')
end
if ~(opt.f < frac)
    warning('Growth is not reduced without glutamine')
end
% fails

%% Leak test
% ensure no metabolite can be produced from nothing
% set lower bounds of all exchange and sink reactions to 0, optimize for a demand
% reaction for each metabolite

ex_rx = model.rxns(findExcRxns(model));
model_mod = changeRxnBounds(model,ex_rx,0,'l');
vec = model_mod.mets;

model_mod = changeRxnBounds(model_mod,{'trdrd_exp','EX_hb','fldox_exp',...
    'acp_exp','EX_pyr(e)'},0,'u');
for i = 1:length(model_mod.mets)
    %disp(i)
    met_test = model_mod.mets(i);
    [model2, newRxn] = addDemandReaction( model_mod,met_test);
    if sum(strcmp(newRxn,model2.rxns)) == 0
        vec{i} = 0; % BUT CHECK IN NEXT STEP
    else model2 = changeObjective(model2,newRxn);
        opt = optimizeCbModel(model2);
        if opt.f <= 0.000001 
            vec{i} = 0;
        else vec{i} = 1;
        end
    end
end
leaky_mets = model_mod.mets(logical(cell2mat(vec)))
% 
% t = table(model2.rxns(abs(opt.x)>0),...
% printRxnFormula(model2,model2.rxns(abs(opt.x)>0)),opt.x(abs(opt.x)>0))

vec2 = ex_rx;
for i = 1:length(ex_rx)
    model2 = changeObjective(model_mod,ex_rx(i));
    opt = optimizeCbModel(model2);
    if opt.f <= 0 
        vec2{i} = 0;
    else vec2{i} = 1;
    end
end
leakymets2 = ex_rx(logical(cell2mat(vec2)))
% NO LEAKY METS

