load('ResultsSVMH1_Detailed_80samples_20runs_RunVariance0.3_20repCV_1iterations_1e6Perm.mat');


clear R
numsamples = 80; numrep = 20; 
numruns = max(l(:,2));
v   = @(x)x(~isnan(x(:)));
%% 

for idx = 1:numel(splitAllvect)
    Splits                              = splitAllvect{idx};
    conditionsvector                    = [Splits.test];
    repetitionsvector                   = [Splits.repetitionidx];
    repetitionsvector                   = ones(size(conditionsvector,1),1)*repetitionsvector;
    conditionsvectorunwrapped           = conditionsvector + (repetitionsvector-1)*numsamples;
    
    e                                   = errPermAllvect{idx};
    C                                   = cov(e');
    Cnan                                = C + diag(nan*ones(size(C,1),1));
    runindices                          = zeros(size(conditionsvectorunwrapped));
    lrep                                = repmat(l,[max(repetitionsvector(:)),1]);
    runindices(:)                       = lrep(conditionsvectorunwrapped(:),2);
    runindices                          = runindices+(repetitionsvector-1)*numruns;
    
    for idrun                           = max(runindices(:)):-1:1
        [~,temp]                        = find(runindices==idrun);
        rep                             = floor((temp(1)-1)/Folds(idx))+1;
        runsinrep                       = ((rep-1)*numruns +  1):rep*numruns;
        runsinrepsamefold               = setdiff(unique(runindices(:,temp)),idrun);
        runsinrepdifffold               = setdiff(runsinrep,unique(runindices(:,temp)));
        
        thismod                         = mod(idrun,numruns);
        if thismod == 0
            thismod = numruns; % 1-based quirk...
        end
        runsameacrossreps               = thismod+[0:max(repetitionsvector(:))-1]*numruns;
        runsameacrossreps               = setdiff(runsameacrossreps,idrun);
        runsdiffacrossreps              = setdiff(1:max(runindices(:)),runsinrep);
        runsdiffacrossreps              = setdiff(runsdiffacrossreps,runsameacrossreps);
        
        
        [a0]                            = ismember(runindices,idrun);
        [arunssamefold]                 = ismember(runindices,runsinrepsamefold);
        [arunsdifffold]                 = ismember(runindices,runsinrepdifffold);
        [arundiffrep]                   = ismember(runindices,runsameacrossreps);
        [arunsdiffrunsdiffrep]          = ismember(runindices,runsdiffacrossreps);
        [arunsdiffrundiffreprunintraining] = and(arunsdiffrunsdiffrep,~any(arundiffrep));
        [arunsdiffrundiffreprunnotintraining] = and(arunsdiffrunsdiffrep,any(arundiffrep));

        
        
        
        R(idx).C                                = C;
        R(idx).Cwithinrun(idrun,:)               = (v(Cnan(a0,a0)));
        R(idx).Cwithinfold(idrun,:)              = (v(Cnan(a0,arunssamefold)));
        R(idx).Cbetweenfolds(idrun,:)            = (v(Cnan(a0,arunsdifffold)));
        R(idx).Crunacrossreps(idrun,:)           =  (v(Cnan(a0,arundiffrep)));
        R(idx).Crundiffacrossreps(idrun,:)       = (v(Cnan(a0,arunsdiffrunsdiffrep)));
        R(idx).Crundiffacrossrepsrunintraining(idrun,:)       = (v(Cnan(a0,arunsdiffrundiffreprunintraining)));
        R(idx).Crundiffacrossrepsrunnotintraining(idrun,:)       = (v(Cnan(a0,arunsdiffrundiffreprunnotintraining)));
        
        
    end
    R(idx).SumC                            = sum(C(:));
    R(idx).SumCoffDiag                     = nansum(Cnan(:));
    if idx < 4
        tmp     = zeros(1,numrep);
        for idrep = 1:numrep
            tmp(idrep) = nansum(v(Cnan((1:numsamples)+(idrep-1)*numsamples,(1:numsamples)+(idrep-1)*numsamples)));
        end
        R(idx).SumCoffdiagwithinrep      = tmp;
    end
    
    
end
save('CovarianceDecomposition.mat','R') % the calculation in this cell is very heavy, you can run it only once, and then if you need it again start from the next command

%%

 load('CovarianceDecomposition.mat','R')

fnames = fieldnames(R); fnames([1 9:end]) = [];

cvnames = {'Split Half','Five Fold','ten fold','LRO'};
v   = @(x)x(~isnan(x(:)));




Cpartition  = zeros(3,7);
Cpartitionmeans = zeros(3,7);
Cpartitionmeansperrun = zeros(3,7,400);

Cpartitionnumitems = zeros(3,7);


for idx = 1:3
    for idnames= 1:7
        Cpartition(idx,idnames) = sum(v(R(idx).(fnames{idnames})));
        Cpartitionmeans(idx,idnames) = mean(v(R(idx).(fnames{idnames})));
        Cpartitionmeansperrun(idx,idnames,:) = mean(R(idx).(fnames{idnames}),2);
        Cpartitionnumitems(idx,idnames) = numel(v(R(idx).(fnames{idnames})));
    end
end
%%
dirdata  = pwd;

figure(2); clf
set(gcf,'color',[1 1 1],'units','centimeters','position',[5 2 27 15]);

CovValues         =[ [sum(diag(R(1).C)) sum(diag(R(1).C)) sum(diag(R(1).C))]' Cpartition(:,1:5) ...
    [R(1:3).SumC]'];
CovValues         = [CovValues(:,1) sum(CovValues(:,[2 3 4]),2) sum(CovValues(:,[5 6]),2) CovValues(:,end)];



axes('position',[.1 .15 .86 .75]);
hold all
maxall = @(x)max(x(:));
for ix = numel(R)-1:-1:1
    tmp(ix) = max(structfun(maxall,R(ix)));
end
ylimval = [0 maxall(tmp)*1.05];
    

h  = bar([1:size(CovValues,2)],CovValues');
h(1).FaceColor = [1 0 0]; 
h(2).FaceColor = [0 0 1]; 
h(3).FaceColor = [0 1 0]; 


set(gca,'Xtick',1:size(CovValues,2),'XtickLabel',{'Main Diagonal','Within Repetition','Between Repetitions','Total'},...
    'XTickLabelRotation',0,'fontsize',14);
hold all
ylim(ylimval);
ys = ylim;

legend(h,{'Split Half','Five-Fold','Ten-Fold'},'location','northwest','Fontsize',14);

xlim([0.5 4.5])

set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
print(gcf,'-dpdf',fullfile(dirdata,'CovarianceDecomposition'));


%%

figure(22); clf
set(gcf,'color',[1 1 1],'units','centimeters','position',[5 2 27 15]);

CovValues         =[ [sum(diag(R(1).C)) sum(diag(R(1).C)) sum(diag(R(1).C))]' Cpartition(:,1:5) ...
    [R(1:3).SumC]'];


axes('position',[.05 .2 .8 .7]);
hold all
maxall = @(x)max(x(:));
for ix = numel(R)-1:-1:1
    tmp(ix) = max(structfun(maxall,R(ix)));
end
ylimval = [0 maxall(tmp)*1.05];
    
patch([ 0.5 1.5 1.5 0.5],[max(ylimval)*ones(1,2) 0*ones(1,2)],[1 0.7 0.5],'FaceAlpha',.4);
patch([ 1.5 4.5 4.5 1.5],[max(ylimval)*ones(1,2) 0*ones(1,2)],[.5 1 1],'FaceAlpha',.4);
patch([ 4.5 6.5 6.5 4.5],[max(ylimval)*ones(1,2) 0*ones(1,2)],[1 1 .5],'FaceAlpha',.4);
patch([ 6.5 7.5 7.5 6.5],[max(ylimval)*ones(1,2) 0*ones(1,2)],[.7 .7 .7],'FaceAlpha',.4);



h  = bar([1:size(CovValues,2)],CovValues');
h(1).FaceColor = [1 0 0]; 
h(2).FaceColor = [0 0 1]; 
h(3).FaceColor = [0 1 0]; 


set(gca,'XtickLabel',{'Diagonal','Within Run','Within Fold','Between Folds','Same Run (reps)','Different Run (reps)','Total'},...
    'XTickLabelRotation',45','fontsize',11);
hold all
% ylim(ylim*1.2);
ylim(ylimval);
ys = ylim;
plot([1.5 1.5],ys,'linewidth',.5,'linestyle','--','Color',[.3 .3 .3]);
plot([4.5 4.5],ys,'linewidth',.5,'linestyle','--','Color',[.3 .3 .3]);
plot([6.5 6.5],ys,'linewidth',.5,'linestyle','--','Color',[.3 .3 .3]);
annotation('textbox','position',[.18 .79 .28 .1],'String','Within CV repetition',...
    'horizontalalignment','center','VerticalAlignment','middle','Fontsize',13,'linestyle','none');
annotation('textbox','position',[.515 .79 .2 .1],'String','Between CV repetitions',...
    'horizontalalignment','center','VerticalAlignment','middle','Fontsize',13,'linestyle','none');
legend(h,{'Split Half','Five-Fold','Ten-Fold'},'position',[.87 .75 .1 .1],'Fontsize',14);

xlim([0.5 7.5])

set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
print(gcf,'-dpdf',fullfile(dirdata,'CovarianceDecompositionDetailed'));

%%

figure(3); clf
set(gcf,'color',[1 1 1],'units','centimeters','position',[5 2 14 18]);

h1 = axes('position',[.1 .7 .85 .25]);

tmp         =Cpartition(:,[6 7 5]) ;
    h  = bar([1:size(tmp,2)], tmp');
h(1).FaceColor = [1 0 0]; 
h(2).FaceColor = [0 0 1]; 
h(3).FaceColor = [0 1 0]; 

ylabel('Correlation Values');




h2 = axes('position',[.1 .4 .85 .25]);

tmp         =Cpartitionnumitems(:,[6 7 5]) ;
    h  = bar([1:size(tmp,2)], tmp');
h(1).FaceColor = [1 0 0]; 
h(2).FaceColor = [0 0 1]; 
h(3).FaceColor = [0 1 0]; 

ylabel('Number of items');




h3 = axes('position',[.1 .1 .85 .25]);

tmp         =Cpartitionmeans(:,[6 7 5]) ;
    h  = bar([1:size(tmp,2)], tmp');
h(1).FaceColor = [1 0 0]; 
h(2).FaceColor = [0 0 1]; 
h(3).FaceColor = [0 1 0]; 
set(gca,'XtickLabel',{'Run in Training','Run not in training','Total'},...
    'XTickLabelRotation',0,'fontsize',14);
ylabel('Mean Correlation Values','fontsize',11);

 set(gcf,'PaperOrientation','portrait','PaperPositionMode','auto');
 print(gcf,'-dpdf',fullfile(dirdata,'CovarianceAcrossrepetitionsDecomposition'));
 
 
 
figure(1); clf
set(gcf,'color',[1 1 1],'units','centimeters','position',[2 2 25 15]);

tmp = [ .25*80*ones(4,1)  [mean(R(1).SumCoffdiagwithinrep(:)) mean(R(2).SumCoffdiagwithinrep(:)) mean(R(3).SumCoffdiagwithinrep(:)) R(4).SumCoffDiag]'];
tmp = [tmp sum(tmp,2)];
h  = bar([1:size(tmp,2)], tmp');
h(1).FaceColor = [1 0 0]; 
h(2).FaceColor = [0 0 1]; 
h(3).FaceColor = [0 1 0]; 
h(4).FaceColor = [.5 .5 .5];
    
ylabel('Covariance');
set(gca,'xtick',[1 2 3],'Xticklabel',{'Main Diagonal','Off Diagonal','Total'},'Fontsize',16);

%%
figure(12); clf
set(gcf,'color',[1 1 1],'units','centimeters','position',[2 2 25 15]);

f1 = @(x) mean(sum(reshape(sum(x,2),[20,size(x,1)/20]),1));

tmp = [ .25*80*ones(4,1)   [f1(R(1).Cwithinrun) f1(R(2).Cwithinrun)  f1(R(2).Cwithinrun) sum(R(4).Cwithinrun(:))]'    ...
     [f1(R(1).Cwithinfold) f1(R(2).Cwithinfold)  f1(R(2).Cwithinfold) sum(R(4).Cwithinfold(:))]'...
     [f1(R(1).Cbetweenfolds) f1(R(2).Cbetweenfolds)  f1(R(2).Cbetweenfolds) sum(R(4).Cbetweenfolds(:))]'];
tmp = [tmp nansum(tmp,2)];
h  = bar([1:size(tmp,2)], tmp');
h(1).FaceColor = [1 0 0]; 
h(2).FaceColor = [0 0 1]; 
h(3).FaceColor = [0 1 0]; 
h(4).FaceColor = [.5 .5 .5];
    
ylabel('Covariance');
set(gca,'xtick',[1 2 3 4 5],'Xticklabel',{'Main Diagonal','Within run','Within Fold','Between Folds','Total'},'Fontsize',16);

legend('Split Half','Five-Fold','Ten-Fold','Leave-one-out','location','northwest','Fontsize',16);


set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
print(gcf,'-dpdf',fullfile(dirdata,'CovarianceDecompositionNoRepetitions'));
