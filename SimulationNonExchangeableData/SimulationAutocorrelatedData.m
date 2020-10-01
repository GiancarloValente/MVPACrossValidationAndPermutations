addmexlibraries;

clear variables



Params = struct('nSamplesperRun',10,'nVoxels',100,'nRuns',10,'correlationWidth',5,...
    'correlationS',[],'designType','random','crossValidation','kfold','nFolds',2,...
    'repeatCrossValidation', true,'nRepetitionsCrossValidation',20,'nPerm',1e3,...
    'permutationType','beforesplit','nDatasets',2000);

Params.nSamplesPerClass1    = Params.nSamplesperRun*Params.nRuns/2;
Params.nSamplesPerClass2    = Params.nSamplesperRun*Params.nRuns/2;

Srange = [.2 .5 .8 1.2];

for iS  = numel(Srange):-1:1
    Params.correlationS = Srange(iS);
    
    pvals(iS)=  getpvaluesScenarios(Params);
    
end




save('ResultsCorrelatedSamplesScenariosPermutationsAndRandomization.mat')

function pvals = getpvaluesScenarios(Params)

%  permutations
Params.randomizationTest = false;

% Scenario A: completely random design
fixedDesign                 = false;
Params.designType           = 'random';
pvals.pA                    = getfalsepositiverate(Params,fixedDesign);

% Scenario B: design is fixed across subjects, it is blocked, using
% permutations
fixedDesign                 = true;
Params.designType           = 'blocked';
pvals.pB                    = getfalsepositiverate(Params,fixedDesign);

% Scenario C: design is fixed across subjects, it is alternate
fixedDesign                 = true;
Params.designType           = 'alternate';
pvals.pC                    = getfalsepositiverate(Params,fixedDesign);


%  randomization
Params.randomizationTest = true;
% Scenario B: design is fixed across subjects, it is blocked, using
% randomizations
fixedDesign                 = true;
Params.designType           = 'blocked';
pvals.pBrandomization       = getfalsepositiverate(Params,fixedDesign);


% Scenario C: design is fixed across subjects, it is alternate, using
% randomizations
fixedDesign                 = true;
Params.designType           = 'alternate';
pvals.pCrandomization       = getfalsepositiverate(Params,fixedDesign);

end






function [p]                    = getfalsepositiverate(Params,fixedDesign)

p                               = zeros(1,Params.nDatasets);
l                               = generateLabels(Params);

for idrep                       = 1:Params.nDatasets
    x                           = generateautocorrelateddata(Params);
    
    if ~fixedDesign
        l                       = generateLabels(Params);
    end
    
    [Splits]                    = generatesplits(Params,l);
    
    
    [err,errPerm]               = classifywithpermutationsmultirunLIBSVM(x,l,Splits,Params);
    p(idrep)                    = (sum(sum(errPerm,2) <= sum(err))+1)./(Params.nPerm + 1);
end
end