clear variables
load ResultsCorrelatedSamplesScenariosPermutationsAndRandomization.mat      


corrplot_x = -Params.correlationWidth:Params.correlationWidth;

for idS = numel(Srange):-1:1
    corr_shape = normpdf(-Params.correlationWidth:Params.correlationWidth,0,Srange(idS));
    corr_shape = corr_shape./max(corr_shape);
    corr_shape_all(idS,:) = corr_shape;
    
%     
%     C(:,:,idS) = full(spdiags([ones(Params.nSamplesperRun,1)*corr_shape],-Params.correlationWidth:Params.correlationWidth,...
%     Params.nSamplesperRun,Params.nSamplesperRun));
    
end

%%
figure(1), clf
set(gcf,'units','centimeters','position',[2 2 27 18],'color',[1 1 1]);
set(groot, 'defaultAxesTickLabelInterpreter','latex');
y_start = .2;
y_end   = .1;
y_s     = 0.05;
y  = (1-(y_start+y_end))/4;
y_height = y - y_s;

h0 = annotation('textbox','position',[.0 .9 .25 .06],'linestyle','none','horizontalAlignment','center',...
    'verticalAlignment','middle','Fontsize',15,'string','Correlation between consecutive samples',...
    'interpreter','Latex')
h0 = annotation('textbox','position',[.25 .9 .7 .06],'linestyle','none','horizontalAlignment','center',...
    'verticalAlignment','middle','Fontsize',20,'string','Type-I Error rate by design and permutation',...
    'interpreter','Latex')
for idS = 1:4
    axes('position',[.05,y_start + y*(4-idS),.15,y_height]);
    plot(-Params.correlationWidth:Params.correlationWidth,corr_shape_all(idS,:),...
        'linewidth',1,'marker','o','markersize',5,'markerfacecolor',[0 0 1]);
    grid on
    set(gca,'Fontsize',12);
    if idS == 4
        xlabel('Distance','Interpreter','Latex');
    else
        set(gca,'XTick',[]) ;
    end
%     ylabel('Correlation')
    
    
    axes('position',[.25,y_start + y*(4-idS),.7,y_height]);
    [ptheta, ci] = binofit(sum([pvals(idS).pA(:) pvals(idS).pB(:)...
        pvals(idS).pBrandomization(:) pvals(idS).pC(:)...
        pvals(idS).pCrandomization(:)]<0.05),numel(pvals(idS).pA));
    
    thisp = ptheta;
    thisciL = ci(:,1)';
    thisciU = ci(:,2)';
                
    Pbar    = [thisciL; thisp-thisciL; thisciU-thisp];
    ylimval = [0 max(ci(:))*1.05];
    
    
    bhandle= bar(1:5, [Pbar'],'stacked');
%     bhandle=  get(gca,'children');
    
    bhandle(1).FaceColor = [1 1 1];
    bhandle(1).EdgeColor = [1 1 1];
%     bhandle(2).FaceColor = [1 1 1]*.95;
    bhandle(2).FaceColor = 'flat';
    bhandle(2).CData     = [1 .6 0.4; 1 .3 .8 ; 1 .6 .4; 1 .3 .8;  1 .6 .4];
    bhandle(2).FaceAlpha = .4;
    bhandle(2).EdgeColor = [.3 .3 .3];
    
%     bhandle(3).FaceColor = [1 1 1]*.95;
    bhandle(3).FaceColor = 'flat';
    bhandle(3).CData     = [1 .6 0.4; 1 .3 .8 ; 1 .6 .4; 1 .3 .8;  1 .6 .4];
    bhandle(3).FaceAlpha = .4;

    bhandle(3).EdgeColor = [.3 .3 .3];
    
    bhandle(2).LineWidth = 1.5;
    bhandle(3).LineWidth = 1.5;
    ylim([0 max(ci(:))*1.05]);
    xlim([.3 5.7])
    set(gca,'Fontsize',14)
    
    if idS == 4
        xTick = 1:5;
        set(gca,'xtick',xTick)
        yTick = get(gca,'ytick');
        set(gca,'xticklabel',[])
        
        xTickLabel = {{'\emph{Random} Design';'\textbf{Unrestricted} Perm'},...
            {'\emph{Blocked} Design';'\textbf{Unrestricted} Perm'},...
            {'\emph{Blocked} Design';'\textbf{Restricted} Perm'},...
            {'\emph{Alternate} Design';'\textbf{Unrestricted} Perm'},...
            {'\emph{Alternate} Design';'\textbf{Restricted} Perm'}};
        for k = 1:length(xTick)
            text(xTick(k),yTick(1)-0.55*(yTick(end)-yTick(1)),xTickLabel{k},...
                'HorizontalAlignment','center','Fontsize',12,'Rotation',25,'Interpreter','latex')
        end
        
%         set(gca,'Xticklabels',{'Random Design \newline Unrestr Perm ','Blocked Perm','Blocked Rand',...
%             'Alternate Perm','Alternate Rand'});
%         set(gca,'XTickLabelRotation',0,'Fontsize',12);
        
    else
        set(gca,'XTick',[]) ;
    end

    hold all

    plot(xlim,[0.05 .05],'linewidth',1.5,'linestyle','--','color',[0 0 0 ]);
 
end
 set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
print(gcf,'-dpdf','SimulationsDesignsCorrelatedSamples.pdf');



%% this is play only, remove it!
return
load ResultsCorrelatedSamplesScenarios_RandomDesignFixedAcrossSubjects.mat

pall = reshape([pvals.pB],[numel(pvals(1).pB) numel(pvals)]);
[p,ci] = binofit(sum(pall<0.05),numel(pvals(1).pB));
figure(2), clf
set(gcf,'units','centimeters','position',[2 2 27 18],'color',[1 1 1]);
[~,tmp] = sort(p);
plot(1:numel(pvals),p(tmp),'linewidth',2);
hold all
fill([1:numel(pvals) flip(1:numel(pvals))],[(ci(tmp,2)')  flip(ci(tmp,1)')],[1 0 0 ],'FaceAlpha', .3)

