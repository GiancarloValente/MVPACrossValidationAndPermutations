clear variables
algused        = {'SVM','GNB','LRL2'};
% % filename       = @(alg) ['Results' alg 'H1_80samples_20runs_100repCV_1000iterations.mat'];

filename       = @(alg) ['Results' alg 'H1_80samples_20runs_RunVariance0.3_100repCV_1000iterations.mat'];

load( filename(algused{1}),'Folds','effectrange','ppartial');
thused         = [0.05];
colors = [1 0 0; 0 0 1; 0 1 0; .3 .3 .3];

numrep  = size(ppartial,3);

for iAlg =1: numel(algused)
    load( filename(algused{iAlg}));
    
    powerincrementalmean        = zeros([size(ppartial,2) size(ppartial,3) numel(thused)]);
    powerincrementalci          = zeros([size(ppartial,2) size(ppartial,3) 2 numel(thused)]);
    for ith = 1: numel(thused)
        for iFold = 1:size(ppartial,2)-1
        [powerincrementalmean(iFold,:,ith),powerincrementalci(iFold,:,:,ith)] = binofit(squeeze(sum(ppartial(:,iFold,:)<thused(ith),1)),size(p,1),.05);
        end
    end
    
    
    figure(iAlg), clf;
    set(gcf,'color',[1 1 1],'units','centimeters','position',[5 5 25 10]);
    h0  = axes('position',[.08 .15 .85 .7],'box','off');
     hold all
    ith = 1;
     for iFold = 1:3
         
         
         
          plot(1:numrep,powerincrementalmean(iFold,:,ith),'linewidth',1,'color',colors(iFold,:).*.9,'linestyle','-');
          hf = fill([1:numrep flip(1:numrep)],[powerincrementalci(iFold,:,2,ith) flip(powerincrementalci(iFold,:,1,ith))],colors(iFold,:).*.9,'FaceAlpha',.05);
          hf.EdgeColor = [1 1 1];
    
          
     end
    
     xlabel('\textbf{Number of Cross-validation repetitions}','interpreter','latex','Fontsize',15);
     ylabel('\textbf{Power}','interpreter','latex','Fontsize',15);
     hl = findobj(gca,'type','line');
     legend(flip(hl),{'Split Half','5-fold','10-Fold'},'location','southeast','Fontsize',12)
     
     title(['\textbf{Effect of Repetition on Power, ' algused{iAlg} '}'],'Fontsize',18,'Fontweight','bold','interpreter','latex');
    namesave = ['H1_' algused{iAlg} '_100repetitions'];
    set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');

    print(gcf,'-dpdf',fullfile(pwd,namesave));

end
    
   