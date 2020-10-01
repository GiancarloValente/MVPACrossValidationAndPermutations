clear variables
load(fullfile('C:\Users\giancarlo.valente\Documents\Research\PermutationsClassifiers\PermutationsClassifiersCode\SimulationsRepetitionCrossValidation',...
    'ResultsSVMH1_80samples_20runs_RunVariance0.3_20repCV_1000iterations.mat'));

rangeRepeatedCV = 0:5:80*20;
rangeCV         = 0:80;

errHistRepCVH1    = zeros(numel(rangeRepeatedCV),numel(Folds),numel(effectrange));
errHistRepCVKernelH1 = zeros(size(errHistRepCVH1));
errHistCVH1    = zeros(numel(rangeCV),numel(Folds),numel(effectrange));
errHistCVKernelH1 = zeros(size(errHistCVH1));
meanerrRepCVH1     = zeros(numel(Folds),numel(effectrange));
meanerrCVH1     = zeros(numel(Folds),numel(effectrange));

for iFold = 1:numel(Folds)
    sum1fold = @(x) sum(x(1:Folds(iFold)));
    for ieff = 1:numel(effectrange)
        thiserr = errvect(:,iFold,ieff);
        errHistRepCVH1(:,iFold,ieff) = hist(cellfun(@sum,thiserr),rangeRepeatedCV);
        errHistCVH1(:,iFold,ieff)    = hist(cellfun(sum1fold,thiserr),rangeCV);
        errHistRepCVKernelH1(:,iFold,ieff) = ksdensity(cellfun(@sum,thiserr),rangeRepeatedCV);
        errHistCVKernelH1(:,iFold,ieff) = ksdensity(cellfun(sum1fold,thiserr),rangeCV);
        meanerrRepCVH1(iFold,ieff)   = mean(cellfun(@sum,thiserr))/max(rangeRepeatedCV);
        meanerrCVH1(iFold,ieff)   = mean(cellfun(sum1fold,thiserr))/max(rangeCV);

    end
end

load(fullfile('C:\Users\giancarlo.valente\Documents\Research\PermutationsClassifiers\PermutationsClassifiersCode\CovarianceDecomposition',...
    'ResultsSVMH1_Detailed_80samples_20runs_RunVariance0.3_20repCV_1iterations_1e6Perm'),'errPermAllvect');

errHistRepCVH0    = zeros(numel(rangeRepeatedCV),numel(Folds));
errHistRepCVKernelH0 = zeros(size(errHistRepCVH0));
errHistCVH0    = zeros(numel(rangeCV),numel(Folds));
errHistCVKernelH0 = zeros(size(errHistCVH0));
meanerrRepCVH0 = zeros(1,numel(Folds));
meanerrCVH0  = zeros(1,numel(Folds));
for iFold = 1:numel(Folds)
    thiserr = sum(errPermAllvect{iFold});
    sum1fold = @(x) sum(x(1:Folds(iFold)));
    
    errHistRepCVH0(:,iFold) = hist(thiserr,rangeRepeatedCV);
    errHistCVH0(:,iFold)    = hist(sum(errPermAllvect{iFold}(1:80,:)),rangeCV);
    errHistRepCVKernelH0(:,iFold) = ksdensity(thiserr,rangeRepeatedCV);
    errHistCVKernelH0(:,iFold) = ksdensity(sum(errPermAllvect{iFold}(1:80,:)),rangeCV);
    meanerrRepCVH0(iFold)   = mean(thiserr)/max(rangeRepeatedCV);
    meanerrCVH0(iFold)   = mean(sum(errPermAllvect{iFold}(1:80,:)))/max(rangeCV);
end


%%
ieff = 5;
colors = [1 0 0; 0 0 1; 0 1 0; .3 .3 .3];
figure(1);
clf;
set(gcf,'units','centimeters','position',[2 2 25 20],'color',[1 1 1]);

axes('units','normalized','position',[.05 .05 .9 .4]);
hold all

for idx = 1:3
    tmp = errHistRepCVKernelH1(:,idx,ieff);
    p1(idx) = plot(rangeRepeatedCV/max(rangeRepeatedCV),tmp./sum(tmp),...
        'color',colors(idx,:),'linewidth',2);
end

for idx = 1:3
    tmp = errHistRepCVKernelH0(:,idx);
    plot(rangeRepeatedCV/max(rangeRepeatedCV),tmp./sum(tmp),...
        'color',colors(idx,:),'linestyle','--','linewidth',2);
end
xlim([.15 .65])
for idx = 1:3
    tmp = errHistRepCVKernelH0(:,idx)./sum(errHistRepCVKernelH0(:,idx));
    thr = rangeRepeatedCV(find(cumsum(tmp)<.05,1,'last'))/max(rangeRepeatedCV);
    plot([thr thr],ylim,'linewidth',2,'linestyle','-.','color',colors(idx,:));
end
plot([0.5 0.5],ylim,'linewidth',2,'color','k')

legend(p1,{'Split-Half','5-Fold','10-Fold'},'Fontsize',12);
title('Repeated Cross-Validation','Fontsize',20);


axes('units','normalized','position',[.05 .55 .9 .4]);
hold all
for idx = 1:4
    tmp = errHistCVKernelH1(:,idx,ieff);
    p2(idx) = plot(rangeCV/max(rangeCV),tmp./sum(tmp),...
        'color',colors(idx,:),'linewidth',2);
end
for idx = 1:4
    tmp = errHistCVKernelH0(:,idx);
    plot(rangeCV/max(rangeCV),tmp./sum(tmp),...
        'color',colors(idx,:),'linestyle','--','linewidth',2);
end

for idx = 1:4
    tmp = errHistCVKernelH0(:,idx)./sum(errHistCVKernelH0(:,idx));
    thr = rangeCV(find(cumsum(tmp)<.05,1,'last'))/max(rangeCV);
    plot([thr thr],ylim,'linewidth',2,'linestyle',':','color',colors(idx,:));
end
plot([0.5 0.5],ylim,'linewidth',2,'color','k')
xlim([.15 .65])
legend(p2,{'Split-Half','5-Fold','10-Fold','Leave-one-out'},'Fontsize',12);
title('Cross-Validation','Fontsize',20);

 set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
 print(gcf,'-dpdf',fullfile(pwd,'ComparisonNullAndAlternativeDistributions'));


