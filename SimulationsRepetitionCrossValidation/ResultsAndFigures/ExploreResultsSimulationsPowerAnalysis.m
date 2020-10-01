clear variables
close all
algused        = {'SVM','GNB','LRL2'};
% filename       = @(alg) ['Results' alg 'H1_80samples_20runs_20repCV_1000iterations.mat'];
varianceRuns   = {'0.0','0.3','0.5'};

idVariance     = 3; % change this to 1 or 2 to plot the results of the other scenarios!
filename       = @(alg) ['Results' alg 'H1_80samples_20runs_RunVariance' varianceRuns{idVariance} '_20repCV_1000iterations.mat'];

load( filename(algused{1}),'Folds','effectrange');


thused                      = [0.05];

ptheta                      = zeros(numel(Folds),numel(effectrange),numel(algused),numel(thused));
pthetarep                   = zeros(size(ptheta));
ci                          = zeros(numel(Folds),numel(effectrange),numel(algused),numel(thused),2);
cirep                       = zeros(size(ci));
ecv                         = zeros(numel(Folds),numel(effectrange),numel(algused),1000);
ecvrep                      = zeros(numel(Folds),numel(effectrange),numel(algused),1000);

for iAlg = numel(algused):-1:1
    load( filename(algused{iAlg}));
    
    for iFold = 1:numel(Folds)
        
        for ith = numel(thused):-1:1
            [ptheta(iFold,:,iAlg,ith),ci(iFold,:,iAlg,ith,:)] = binofit(sum(squeeze(p(:,iFold,:))<thused(ith)),size(p,1),.05);
            [pthetarep(iFold,:,iAlg,ith),cirep(iFold,:,iAlg,ith,:)] = binofit(sum(squeeze(prep(:,iFold,:))<thused(ith)),size(prep,1),.05);
        end
        
        for ieff = 1:numel(effectrange)
            thiserr = errvect(:,iFold,ieff);
            
            for idrep  = 1:numel(thiserr)
                ecvrep(iFold,ieff,iAlg,idrep) = sum(thiserr{idrep})/(80*(20*(iFold<=3) + 1*(iFold==4)));
                ecv(iFold,ieff,iAlg,idrep) = sum(thiserr{idrep}(1:Folds(iFold)))/80;
                
            end


        end
    end
    
end
%%

FoldName   = {'Split half','5-Fold','10-Fold','LRO','Split half REPEATED','5-Fold REPEATED',' 10-Fold REPEATED'};
axXsep  = .05;
axXoffsetleft  = .08;
axXoffsetright = 0.03;
axwidth    = (1 - axXoffsetleft -axXoffsetright-axXsep*(numel(thused)-1))/numel(thused);

axYoffsetbottom = 0.12;
axYoffsettop    = .12;
axYsep      = 0.02;
axheight   = (1 - axYoffsettop - axYoffsetbottom);
axXpos     = axXoffsetleft + (axwidth+axXsep)*[0:numel(thused)-1];
axYpos     =flip( axYoffsetbottom + (axheight+axYsep)*[0:0]);
v = @(x)x(:);





colors = [1 0 0; 0 0 1; 0 1 0; .3 .3 .3];
for iAlg = 1 : numel(algused)
    figure(iAlg), clf;
    set(gcf,'color',[1 1 1],'units','centimeters','position',[5 5 25 15]);
    
    for ith = 1:numel(thused)
        
        h0  = axes('position',[axXpos(ith) axYpos(1) axwidth axheight],'box','off');
        hold all
        
        for iFold = 1:numel(Folds)
            
            
            thisp = (ptheta(iFold,:,iAlg,ith));
            thisciL = ci(iFold,:,iAlg,ith,1);
            thisciU = ci(iFold,:,iAlg,ith,2);
            
            
            
            plot(effectrange,thisp,'linewidth',2,'color',colors(iFold,:).*.9,'linestyle','-.');
            hf = fill([effectrange flip(effectrange)],[thisciU flip(thisciL)],colors(iFold,:).*.9,'FaceAlpha',.05);
            hf.EdgeColor = [1 1 1];
            
            
            
            if iFold < numel(Folds)
                thisp = (pthetarep(iFold,:,iAlg,ith));
                thisciL = cirep(iFold,:,iAlg,ith,1);
                thisciU = cirep(iFold,:,iAlg,ith,2);
                
                
                
                plot(effectrange,thisp,'linewidth',2,'color',colors(iFold,:),'linestyle','-');
                hf = fill([effectrange flip(effectrange)],[thisciU flip(thisciL)],colors(iFold,:),'FaceAlpha',.05);
                hf.EdgeColor = [1 1 1];
            end
            
            
        end
        
        
        plot(xlim,(thused(ith))*[1 1],'color',[0 0 0],'linewidth',2);
        hl = findobj(gca,'type','line');
        hl(1) = [];
        legend(flip(hl),flip(FoldName([4,7,3,6,2,5,1])),'location','northwest','Fontsize',12)
        xlabel('\textbf{Injected Class Difference}','interpreter','latex','Fontsize',15);
        ylabel('\textbf{Power}','interpreter','latex','Fontsize',15);
        grid on
        grid minor         
    end
    title(['\textbf{ Power of cross-validation schemes, ' algused{iAlg} '}'],'Fontsize',18,'Fontweight','bold','interpreter','latex');
    namesave = ['H1_RunVariance' varianceRuns{idVariance} '_' algused{iAlg} '.pdf'];
    set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');

     print(gcf,'-dpdf',fullfile(pwd,namesave));

end



%%

ecvmean     = mean(ecv,4);
ecvsterr    = std(ecv,0,4)./sqrt(1000);
ecvrepmean = mean(ecvrep,4);
ecvrepstderr= std(ecvrep,0,4)./sqrt(1000);



FoldName   = {'Split half','5-Fold','10-Fold','LRO','Split half REPEATED','5-Fold REPEATED',' 10-Fold REPEATED'};
axXsep  = .05;
axXoffsetleft  = .08;
axXoffsetright = 0.03;
axwidth    = (1 - axXoffsetleft -axXoffsetright-axXsep*(numel(thused)-1))/numel(thused);

axYoffsetbottom = 0.12;
axYoffsettop    = .12;
axYsep      = 0.02;
axheight   = (1 - axYoffsettop - axYoffsetbottom);
axXpos     = axXoffsetleft + (axwidth+axXsep)*[0:numel(thused)-1];
axYpos     =flip( axYoffsetbottom + (axheight+axYsep)*[0:0]);
v = @(x)x(:);





colors = [1 0 0; 0 0 1; 0 1 0; .3 .3 .3];
for iAlg = 1: numel(algused)
    figure(iAlg+4), clf;
    set(gcf,'color',[1 1 1],'units','centimeters','position',[5 5 25 15]);
    
    for ith = 1
        
        h0  = axes('position',[axXpos(ith) axYpos(1) axwidth axheight],'box','off');
        hold all
        
        for iFold = 1:numel(Folds)
            
            
%             thisp = (ptheta(iFold,:,iAlg,ith));
%             thisciL = ci(iFold,:,iAlg,ith,1);
%             thisciU = ci(iFold,:,iAlg,ith,2);
            
            
            thiserr = ecvmean(iFold,:,iAlg);
            thiserrL = ecvmean(iFold,:,iAlg)-ecvsterr(iFold,:,iAlg);
            thiserrU = ecvmean(iFold,:,iAlg)+ecvsterr(iFold,:,iAlg);
                        
            
            plot(effectrange,thiserr,'linewidth',2,'color',colors(iFold,:).*.9,'linestyle','-.'); 
            %'marker','d','markerfacecolor',colors(iFold,:),'MarkerEdgeColor',[0 0 0],'markersize',5);
            hf = fill([effectrange flip(effectrange)],[thiserrU flip(thiserrL)],colors(iFold,:).*.9,'FaceAlpha',.05);
            hf.EdgeColor = [1 1 1];
            
            
            
            if iFold < numel(Folds)
                thisp = (pthetarep(iFold,:,iAlg,ith));
                thisciL = cirep(iFold,:,iAlg,ith,1);
                thisciU = cirep(iFold,:,iAlg,ith,2);
                
                thiserr = ecvrepmean(iFold,:,iAlg);
                thiserrL = ecvrepmean(iFold,:,iAlg)-ecvrepstderr(iFold,:,iAlg);
                thiserrU = ecvrepmean(iFold,:,iAlg)+ecvrepstderr(iFold,:,iAlg);
                
                
                plot(effectrange,thiserr,'linewidth',2,'color',colors(iFold,:),'linestyle','-');
                hf = fill([effectrange flip(effectrange)],[thiserrU flip(thiserrL)],colors(iFold,:),'FaceAlpha',.05);
                hf.EdgeColor = [1 1 1];
            end
            
            
        end
        
        
%         plot(xlim,(.5*[1 1]),'color',[0 0 0],'linewidth',2);
        ylim([.25 .5])
        hl = findobj(gca,'type','line');
        
        legend(flip(hl),flip(FoldName([4,7,3,6,2,5,1])),'location','southwest','Fontsize',12)
        xlabel('\textbf{Injected Class Difference}','interpreter','latex','Fontsize',15);
        ylabel('\textbf{Average cross-validation error}','interpreter','latex','Fontsize',15);
        grid on
        grid minor         
    end
    title(['\textbf{ Errors of cross-validation schemes, ' algused{iAlg} '}'],'Fontsize',18,'Fontweight','bold','interpreter','latex');
    namesave = ['H1Errors_' varianceRuns{idVariance} '_' algused{iAlg},'.pdf'];
    set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');

      print(gcf,'-dpdf',fullfile(pwd,namesave));

end



