function SH1_RunVariability_LargeNumberOfRepetitions_ErrorStats

addmexlibraries

rng(10717);

p1 = gcp('nocreate');
if isempty(p1)
    parpool('AllCores',36);
    p1 = gcp;
    p1.IdleTimeout          = 1e8;
end


effectrange                 = .07:.07:.5;
% choosing only the fourth, and zooming in with more repetitions of CV
effectrange                 = effectrange(4);
Folds                       = [2 5 10 20];
nDatasets                   =  1;

PermType                    = 'beforesplit';

p                           = zeros(nDatasets,numel(Folds),numel(effectrange));
prep                        = zeros(nDatasets,numel(Folds),numel(effectrange));
tvec                        = zeros(nDatasets*numel(effectrange),1);
counter                     = 0;
errvect                     = cell(nDatasets,numel(Folds),numel(effectrange));
errPermAllvect              = cell(nDatasets,numel(Folds),numel(effectrange));
splitAllvect                = cell(nDatasets,numel(Folds),numel(effectrange));

t0                          = clock;
for iEffectsize         = 1:numel(effectrange)
    Params                      = initializeparams('differenceMagnitude',effectrange(iEffectsize),...
        'runVariance',.3,...
        'repeatCrossValidation',1,'nRepetitionsCrossValidation',20,...
        'nSamplesPerClass1',40,'nSamplesPerClass2',40,'nDatasets',nDatasets,'nPerm',1e6,...
        'permutationType',PermType,'nRuns',20,'nFolds',5);
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
            
            
            Splits                  = generatesplits(Params1,l);

            [err,errPerm,~,ePall]   = classifywithpermutationsmultirunLIBSVM_detailederrors(x,l,Splits,Params1);
   
                            

            
            prep(iDataset,iFold,iEffectsize) = (sum(sum(errPerm,2) <= sum(err))+1)./(Params.nPerm + 1);
            temp                = [Splits.repetitionidx]==1;
            p(iDataset,iFold,iEffectsize)   = ...
                (sum(sum(errPerm(:,temp),2) <= sum(err(temp)))+1)./(Params.nPerm + 1);
  

            
            
            errvect{iDataset,iFold,iEffectsize} = err;
            

            
                                            
                                            
            conditionsvector                    = [Splits.test];
            repetitionsvector                   = [Splits.repetitionidx];
            conditionsvectorunwrapped           = conditionsvector + (repetitionsvector-1)*(Params1.nSamplesPerClass1 + Params1.nSamplesPerClass2);
            errorvectorP_detailed               = zeros([numel(conditionsvectorunwrapped) Params1.nPerm]);
            errorvectorP_detailed(conditionsvectorunwrapped,:) ...
                                                = reshape(ePall,[numel(conditionsvectorunwrapped),Params1.nPerm]);
                                            
            errPermAllvect{iDataset,iFold,iEffectsize}  ...
                                                = errorvectorP_detailed;
            splitAllvect{iDataset,iFold,iEffectsize} ...
                                                = Splits;
            
        end
        
        t2=clock;
        counter                 = counter + 1;
        
        tvec(counter)          = etime(t2,t1);
        fprintf('Done, effectsize nr %d, iteration %d in %2.4f seconds, total time = %5.2f seconds\n',iEffectsize,iDataset,etime(t2,t1), etime(t2,t0));
        fprintf('Estimated time to completion: %4.2f seconds\n',etime(t2,t0)/counter*(numel(tvec)-counter));
        
    end
    
    
end
nameSave    = sprintf('Results%sH1_Detailed_%dsamples_%druns_RunVariance%.1f_%drepCV_%diterations.mat','SVM',Params.nSamplesPerClass1+Params.nSamplesPerClass2,...
    Params.nRuns,Params.runVariance,Params.nRepetitionsCrossValidation,Params.nDatasets);
save(nameSave,'p','prep','errPermAllvect','splitAllvect','l','effectrange','Folds','errvect','-v7.3');

