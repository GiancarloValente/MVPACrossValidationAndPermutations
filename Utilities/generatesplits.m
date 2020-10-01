function [Splits]               = generatesplits(Params,l)
% dividing the available runs in training-testing. If multiple runs are
% available, then the cross-validation is done on a run-level

nSamplesClass1                  = Params.nSamplesPerClass1;
nSamplesClass2                  = Params.nSamplesPerClass2;
nSamplesTotal                   = nSamplesClass1 + nSamplesClass2;
nRuns                           = Params.nRuns;

% if (strcmpi(Params.crossValidation,'kfold') && nRuns == Params.nFolds)
%     Params.crossValidation      = 'lro';
% end




switch lower(Params.crossValidation)
    case 'lro'
        Splits(nRuns)           = struct('train',[],'test',[],'repetitionidx',[]);
        for irun                = 1:nRuns
            Splits(irun).test   = find(l(:,2) == irun);
            Splits(irun).train  = setdiff(1:nSamplesTotal,Splits(irun).test);
            Splits(irun).runtest= irun;
            Splits(irun).runtrain ...
                                = setdiff(1:Params.nRuns,irun);
        end
        
    case 'kfold'
        nFolds                  = Params.nFolds;
        nRunsperFold            = floor(nRuns./nFolds);
        if nRunsperFold*nFolds  ~= nRuns
            disp('warning: the number of runs is not a multiple of the number of cv folds');
            disp('the last cv fold will contain more runs');
        end
        if Params.repeatCrossValidation 
            nCVRep              = Params.nRepetitionsCrossValidation;
        else
            nCVRep              = 1;
        end
        Splits(nFolds*nCVRep)   = struct('train',[],'test',[],'repetitionidx',[]);
        for iCVRep              = 1:nCVRep
            runidx              = randperm(nRuns);
            for iFold           = 1:nFolds
                if iFold        < nFolds
                   runsTest     = runidx(nRunsperFold*(iFold-1)+1  : nRunsperFold*iFold);
                else
                   runsTest     = runidx(nRunsperFold*(nFolds-1)+1 : end);
                end
                runsTrain        = setdiff(runidx,runsTest);
                Splits(nFolds*(iCVRep-1)+iFold).test ...
                                = find(ismember(l(:,2),runsTest));
                Splits(nFolds*(iCVRep-1)+iFold).train ...
                                = find(ismember(l(:,2),runsTrain));  
                Splits(nFolds*(iCVRep-1)+iFold).repetitionidx ...
                                = iCVRep;
                Splits(nFolds*(iCVRep-1)+iFold).runtest ...
                                = runsTest;
                Splits(nFolds*(iCVRep-1)+iFold).runtrain ...
                                = runsTrain;
                            
                 
            
            end
        end
        
    case 'holdout'
        nRunsTest                = floor(nRuns*Params.percentageHoldOut);
        
        if Params.repeatHoldout 
            nHORep              = Params.nRepetitionsHoldout;
        else
            nHORep              = 1;
        end
        Splits(nHORep)          = struct('train',[],'test',[],'repetitionidx',[]);
        for iHoldouot           = 1:nHORep
            runsTest            = randperm(nRuns,nRunsTest);
            runsTrain           = setdiff(1:nRuns,runsTest);
            Splits(iHoldouot).test ...
                                = find(ismember(l(:,2),runsTest));
            Splits(iHoldouot).train ...
                                = find(ismember(l(:,2),runsTrain));
            Splits(iHoldouot).repetitionidx ...
                                = iHoldouot;
            Splits(iHoldouot).runtest ...
                                = runsTest;
            Splits(iHoldouot).runtrain ...
                                = runsTrain;
        end
        
        
end






