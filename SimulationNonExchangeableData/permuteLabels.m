function  [lp]                  = permuteLabels(l,Params,Splits)
% shadowing the original version to implement randomization rather than permutation  

nSplits                         = numel(Splits);
runSizes                        = histcounts(l(:,2),[1:Params.nRuns+1]);
if all(runSizes == runSizes(1))
    equalrunsize                = true;
else
    equalrunsize                = false;
end


if equalrunsize
    selectrun                   = zeros(runSizes(1),Params.nRuns);
    for irun                    = 1:Params.nRuns
        selectrun(:,irun)       = find(l(:,2)==irun);
    end
    
        
    switch lower(Params.permutationType)
        case 'beforesplit'
           
            if Params.randomizationTest
                tmp = binornd(1,.5,[1 Params.nRuns]);
                lp = l;
                for idrun = 1:Params.nRuns
                    if tmp(idrun)
                        lp(lp(:,2)==idrun,1) = 3-lp(lp(:,2)==idrun,1) ;
                    end
                end
                
            else
                
                [~,randsort]        = sort(rand(runSizes(1),Params.nRuns));
                randsort            = randsort + [0:runSizes(1):(Params.nRuns-1)*runSizes(1)];
                ltemp(selectrun,1)  = l(selectrun(randsort),1);
                lp                  = repmat(ltemp,[1 1 nSplits]);
            end
           
        case 'withinsplit'
            labelsonly          = l(:,1);
            ltemp               = zeros(size(l,1),numel(Splits));
            [~,randsort]        = sort(rand(runSizes(1),Params.nRuns.*numel(Splits)));
            randsort            = randsort + repmat([0:runSizes(1):(Params.nRuns-1)*runSizes(1)],1,numel(Splits));
            temp                = reshape(selectrun(randsort),[Params.nRuns*runSizes(1) numel(Splits)]);
            ltemp(repmat(selectrun(:),[1 numel(Splits)])+[0:numel(labelsonly):(numel(Splits)-1)*numel(labelsonly)] ) ...
                                = labelsonly(temp);
            lp                  = zeros([size(l) numel(Splits)]);
            lp(:,1,:)           = ltemp;
            lp(:,2,:)           = repmat(l(:,2),[1 numel(Splits)]);
           
            

        case 'withinsplittrain'
            nrunsperSplit       = numel([Splits(1).runtest]); 
            nrunsperrep         = numel(unique([Splits([Splits.repetitionidx]==1).runtest]));
            nreps               = max([Splits.repetitionidx]);

            testruns            = [Splits.runtest];
            testruns            = reshape(testruns,[nrunsperSplit,numel(testruns)./nrunsperSplit]);
            
            temp                = 0:Params.nRuns:(Params.nRuns*numel(Splits))-1;
            testruns            = testruns + temp;
            testruns            = testruns(:)';
            
            
            labelsonly          = l(:,1);
            ltemp               = zeros(size(l,1),numel(Splits));
            [~,randsort]        = sort(rand(runSizes(1),Params.nRuns.*numel(Splits)));
            randsort(:,testruns)= repmat([1:runSizes(1)]',1,numel(testruns));
            
            
            randsort            = randsort + repmat([0:runSizes(1):(Params.nRuns-1)*runSizes(1)],1,numel(Splits));
            temp                = reshape(selectrun(randsort),[Params.nRuns*runSizes(1) numel(Splits)]);
            ltemp(repmat(selectrun(:),[1 numel(Splits)])+[0:numel(labelsonly):(numel(Splits)-1)*numel(labelsonly)] ) ...
                                = labelsonly(temp);
            lp                  = zeros([size(l) numel(Splits)]);
            lp(:,1,:)           = ltemp;
            lp(:,2,:)           = repmat(l(:,2),[1 numel(Splits)]);
            
            
            
        case 'withinsplittest'
            nrunsperSplit       = numel([Splits(1).runtrain]); 
%             nrunsperrep         = numel(unique([Splits([Splits.repetitionidx]==1).runtrain]));
%             nreps               = max([Splits.repetitionidx]);

            trainruns           = [Splits.runtrain];
            trainruns            = reshape(trainruns,[nrunsperSplit,numel(trainruns)./nrunsperSplit]);
            
            temp                = 0:Params.nRuns:(Params.nRuns*numel(Splits))-1;
            trainruns            = trainruns + temp;
            trainruns            = trainruns(:)';
            
            
            labelsonly          = l(:,1);
            ltemp               = zeros(size(l,1),numel(Splits));
            [~,randsort]        = sort(rand(runSizes(1),Params.nRuns.*numel(Splits)));
            randsort(:,trainruns)= repmat([1:runSizes(1)]',1,numel(trainruns));
            
            
            randsort            = randsort + repmat([0:runSizes(1):(Params.nRuns-1)*runSizes(1)],1,numel(Splits));
            temp                = reshape(selectrun(randsort),[Params.nRuns*runSizes(1) numel(Splits)]);
            ltemp(repmat(selectrun(:),[1 numel(Splits)])+[0:numel(labelsonly):(numel(Splits)-1)*numel(labelsonly)] ) ...
                                = labelsonly(temp);
            lp                  = zeros([size(l) numel(Splits)]);
            lp(:,1,:)           = ltemp;
            lp(:,2,:)           = repmat(l(:,2),[1 numel(Splits)]);
            
            
    end
    
    
    
else
    
    
    switch lower(Params.permutationType)
        case 'beforesplit'
            ltemp               = l;
            
            for indrun          = 1:Params.nRuns
                trrun           = l(:,2)==indrun;
                lrun            = l(trrun,1);
                ltemp(trrun,1)  = lrun(randperm(size(lrun,1)));
            end
            lp                  = repmat(ltemp,[1 1 nSplits]);
        case 'withinsplit'
            lp                  = zeros([size(l) nSplits]);
            for iSplit          = 1:numel(Splits)
                ltemp           = l;
                for indrun      = 1:Params.nRuns
                    trrun       = l(:,2)==indrun;
                    lrun        = l(trrun,1);
                    ltemp(trrun,1)  ...
                        = lrun(randperm(size(lrun,1)));
                end
                lp(:,:,iSplit)  = ltemp;
            end
        case 'withinsplittrain'
            lp                  = zeros([size(l) nSplits]);
            for iSplit          = 1:numel(Splits)
                ltemp           = l;
                trainRuns       = unique(l(Splits(iSplit).train,2));
                for indrun      = trainRuns(:)'
                    trrun       = l(:,2)==indrun;
                    lrun        = l(trrun,1);
                    ltemp(trrun,1)  ...
                        = lrun(randperm(size(lrun,1)));
                end
                lp(:,:,iSplit)  = ltemp;
            end
        case 'withinsplittest'
            lp                  = zeros([size(l) nSplits]);
            for iSplit          = 1:numel(Splits)
                ltemp           = l;
                testRuns        = unique(l(Splits(iSplit).test,2));
                for indrun      = testRuns(:)'
                    trrun       = l(:,2)==indrun;
                    lrun        = l(trrun,1);
                    ltemp(trrun,1)  ...
                        = lrun(randperm(size(lrun,1)));
                end
                lp(:,:,iSplit)  = ltemp;
            end
            
            
    end
end
