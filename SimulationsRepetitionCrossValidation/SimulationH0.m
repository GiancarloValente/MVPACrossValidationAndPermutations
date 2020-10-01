function SimulationH0_v1(algname)

addmexlibraries

% initialize the seed of the random number generator
rng(10717);


% select the right number of cores you want to use, if you want the defaults, comment from lines 10 to 15 and let Matlab take care of this once you call a parfor
p1 = gcp('nocreate');
if isempty(p1)
    parpool('AllCores',36);
    p1 = gcp;
    p1.IdleTimeout          = 1e8;
end


effectrange                 = 0; % if effectrange is zero, then we are under H0
Folds                       = [2 5 10];
nDatasets                   =  1000;

PermType                    = {'beforesplit','withinsplit','withinsplittrain','withinsplittest'};

p                           = zeros(nDatasets,numel(Folds),numel(PermType));
prep                        = zeros(nDatasets,numel(Folds),numel(PermType));
tvec                        = zeros(nDatasets,1);
counter                     = 0;
errvect                     = cell(nDatasets,numel(Folds),numel(PermType));

t0                          = clock;
Params                      = initializeparams('differenceMagnitude',effectrange,...
    'repeatCrossValidation',1,'nRepetitionsCrossValidation',20,...
    'nSamplesPerClass1',40,'nSamplesPerClass2',40,'nDatasets',nDatasets,'nPerm',1e3,...
    'permutationType',PermType,'nRuns',10,'nFolds',5);
for iDataset                = 1:Params.nDatasets
    t1=clock;
    [x,l]                   = generatedatamultirun(Params);
    
    for iFold           = 1:numel(Folds)
        if Folds(iFold) == Params.nRuns
            repXVal     = false;
            numrepXVal  = 1;
        else
            repXVal     = true;
            numrepXVal  = Params.nRepetitionsCrossValidation;
        end
        Params1            = Params;
        Params1.repeatCrossValidation = repXVal;
        Params1.nRepetitionsCrossValidation = numrepXVal;
        Params1.nFolds  = Folds(iFold);
        for iPerm = 1:numel(PermType)
            Params1.permutationType  = PermType{iPerm};
        
        Splits                  = generatesplits(Params1,l);
        switch lower(algname)
            case 'svm'
                [err,errPerm]       = classifywithpermutationsmultirunLIBSVM(x,l,Splits,Params1);
            case 'gnb'
                [err,errPerm]       = classifywithpermutationsmultirunGNB_parallel(x,l,Splits,Params1);
            case 'lrl2'
                [err,errPerm]       = classifywithpermutationsmultirunLiblinear(x,l,Splits,Params1);
        end
        
        
        
        prep(iDataset,iFold,iPerm) = (sum(sum(errPerm,2) <= sum(err))+1)./(Params.nPerm + 1);
        temp                = [Splits.repetitionidx]==1;
        p(iDataset,iFold,iPerm)   = ...
            (sum(sum(errPerm(:,temp),2) <= sum(err(temp)))+1)./(Params.nPerm + 1);
        
        errvect{iDataset,iFold,iPerm} = err;
        
        
        
        
        
        end
    end
    
    t2=clock;
    counter                 = counter + 1;
    
    tvec(counter)          = etime(t2,t1);
    fprintf('Done, permutation time nr %d, iteration %d in %2.4f seconds, total time = %5.2f seconds\n',iPerm,iDataset,etime(t2,t1), etime(t2,t0));
    fprintf('Estimated time to completion: %4.2f seconds\n',etime(t2,t0)/counter*(numel(tvec)-counter));
    
end


nameSave    = sprintf('Results%sH0_%dsamples_%druns_%drepCV_%diterations.mat', upper(algname),Params.nSamplesPerClass1+Params.nSamplesPerClass2,...
    Params.nRuns,Params.nRepetitionsCrossValidation,Params.nDatasets);
save(nameSave,'p','prep','effectrange','Folds','errvect');



