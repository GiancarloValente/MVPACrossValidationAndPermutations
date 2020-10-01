function [E,Er,p,prep]  = comparexvalconnectome_iteration(features,labels,EVNames,class1,class2,varargin)

numberOfSubjectsPerBatch        = 20;
Folds                           = [2 5 10 20];
numberOfPermutations            = 50;
numvoxused                      =  50;
vox_reduction                   = 'random' ;

if numel(varargin)>=1
    numberOfPermutations        = varargin{1};
end
if numel(varargin)>=2
    numvoxused                  = varargin{2};
end
if numel(varargin)>=3
    vox_reduction              = varargin{3};
end
if numel(varargin)>=4      
    numberOfSubjectsPerBatch    = varargin{4};
end
if numel(varargin)>=5      
    Folds                       = varargin{5};
end



indexClass1                     = find(ismember(EVNames,class1));
indexClass2                     = find(ismember(EVNames,class2));
goodvertices                    = ~isnan(features(1,:)); %remove nans


% randomizing the order of subjects
totsubjects         = numel(unique(labels(:,2)));
tmp                 = randperm(totsubjects);
idx                 = (tmp-1)*numel(EVNames)+ [1:numel(EVNames)]';
idx                 = idx(:);



subjforvoxelselection            = mod(totsubjects,numberOfSubjectsPerBatch);
if subjforvoxelselection < 5 
    subjforvoxelselection = subjforvoxelselection + numberOfSubjectsPerBatch;
end



tmp1                            = features(idx(indexClass1:numel(EVNames):numel(EVNames)*subjforvoxelselection),goodvertices);
tmp2                            = features(idx(indexClass2:numel(EVNames):numel(EVNames)*subjforvoxelselection),goodvertices);
[~,ptemp]                       = ttest2(tmp1,tmp2);

indicesgoodvertices             = find(goodvertices);
tmp                             = false(size(goodvertices));


switch lower(vox_reduction)
    case 'random'
        
        nonsigvox                       = find(ptemp > .8);
        indvoxchosen                    = nonsigvox(randperm(numel(nonsigvox),numvoxused));
                
%       sigvox                          = find(ptemp<0.05);
%        indvoxchosen                    = [sigvox(randperm(numel(sigvox),0)) nonsigvox(randperm(numel(nonsigvox),numvoxused))];
    

    case 'univariate'
        [~,ranking]                     = sort(ptemp,'descend');
        indvoxchosen                    = ranking(1:numvoxused);
    case 'rfe'
        [mdl]                           = fitcsvm(zscore([tmp1;tmp2]),[ones(size(tmp1,1),1); -ones(size(tmp2,1),1)]);
        A                               = zeros(size(tmp1,1)+size(tmp2,1),1);
        A(mdl.IsSupportVector)          = mdl.Alpha;
        M                               = (A.*[ones(size(tmp1,1),1); -ones(size(tmp2,1),1)])'*zscore([tmp1;tmp2]);
        [~,ranking]                     = sort(abs(M),'ascend');
        indvoxchosen                    = ranking(1:numvoxused);

end
tmp(indicesgoodvertices(indvoxchosen))  = goodvertices(indicesgoodvertices(indvoxchosen));
goodvertices                            = tmp;
featuresred                             = features(idx((subjforvoxelselection*numel(EVNames) + 1):end),goodvertices);
labelsred                               = labels(idx((subjforvoxelselection*numel(EVNames) + 1):end),:);


numberOfSubjects                = size(featuresred,1)/numel(EVNames);
numberOfBatches                 = floor(numberOfSubjects/numberOfSubjectsPerBatch);

indexSubjects                   = reshape(randperm(numberOfSubjects,numberOfBatches*numberOfSubjectsPerBatch),numberOfBatches,numberOfSubjectsPerBatch);
indexSubjects                   = sort(indexSubjects,2);




Params                      = initializeparams('differenceMagnitude',nan,...
    'repeatCrossValidation',1,'nRepetitionsCrossValidation',20,...
    'nSamplesPerClass1',numberOfSubjectsPerBatch,'nSamplesPerClass2',numberOfSubjectsPerBatch,'nDatasets',nan,...
    'nPerm',numberOfPermutations,'permutationType','beforesplit','nRuns',numberOfSubjectsPerBatch,'nFolds',2);
Er                  = zeros(numberOfBatches,numel(Folds));
E                   = zeros(numberOfBatches,numel(Folds));

prep                = zeros(numberOfBatches,numel(Folds));
p                   = zeros(numberOfBatches,numel(Folds));

fprintf('\nAnalyzing batch  %.0f of %2.0f',0,numberOfBatches);
for indBatch = 1:numberOfBatches
    fprintf(repmat('\b',1,8));
    fprintf('%2.0f of %2.0f',indBatch,numberOfBatches);
    thisIndexSubjects   = indexSubjects(indBatch,:);
    thisfeaturesIndex   = (thisIndexSubjects(:)-1)*numel(EVNames) + [1:numel(EVNames)];
    thisfeaturesIndex   = thisfeaturesIndex';
    thisfeaturesIndex   = thisfeaturesIndex(:);
    x                   = featuresred(thisfeaturesIndex,:);
    l                   = labelsred(thisfeaturesIndex,:);
    
    x                   = [x(indexClass1:numel(EVNames):end,:); x(indexClass2:numel(EVNames):end,:)];
    l                   = [l(indexClass1:numel(EVNames):end,:); l(indexClass2:numel(EVNames):end,:)];
    [~,~,l(:,2)]        = unique(l(:,2));
    ltemp               = l;
    l(ltemp(:,1) == indexClass1,1) = 1;
    l(ltemp(:,1) == indexClass2,1) = 2;
    
    
    
    
    for iFold           = 1:numel(Folds)
        if Folds(iFold) == numberOfSubjectsPerBatch
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
        
        
        Splits              = generatesplits(Params1,l);
        [err,errPerm]       = classifywithpermutationsmultirunLIBSVM(x,l,Splits,Params1);
        
        Er(indBatch,iFold)   = sum(err);
        E(indBatch,iFold)   = sum(err(1:Params1.nFolds));
        
        prep(indBatch,iFold)          = (sum(sum(errPerm,2) <= sum(err))+1)./(Params.nPerm + 1);
        p(indBatch,iFold)             = (sum(sum(errPerm(:,1:Params1.nFolds),2) <= sum(err (1:Params1.nFolds)))+1)./(Params.nPerm + 1);
    end
end
fprintf('\n');
