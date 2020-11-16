

clear variables
algused        = {'SVM','GNB','LRL2'};
filename       = @(alg) ['Results' alg 'H0_80samples_20runs_20repCV_1000iterations.mat'];
load( filename(algused{1}),'Folds');
saveplots                   = true;

thused                      = [0.01 0.05 .1];
PermType                    = {'beforesplit','withinsplit','withinsplittrain','withinsplittest'};

ptheta                      = zeros(numel(Folds),numel(PermType),numel(algused),numel(thused));
pthetarep                   = zeros(size(ptheta));
ci                          = zeros([size(ptheta) 2 ]);
cirep                       = zeros(size(ci));

for iAlg = numel(algused):-1:1
    load( filename(algused{iAlg}));
    
    for iFold = 1:numel(Folds)
        
        for ith = numel(thused):-1:1
            [ptheta(iFold,:,iAlg,ith),ci(iFold,:,ith,iAlg,:)] = binofit(sum(squeeze(p(:,iFold,:))<thused(ith)),size(p,1),.05);
            [pthetarep(iFold,:,iAlg,ith),cirep(iFold,:,ith,iAlg,:)] = binofit(sum(squeeze(prep(:,iFold,:))<thused(ith)),size(prep,1),.05);
        end
    end
    
end
%%

for idfigure = 1:2
    
    Permscheme = {'Before','Within','Within (Train)','Within (test)'};
    if idfigure ==1 
        FoldName   = {'Split half','5-Fold','10 Fold','LRO'};
    else
        FoldName   = {'Split half','5-Fold','10 Fold'};
        Folds      = Folds(1:3);
    end
        
    axXsep  = .04;
    axXoffset  = .08;
    axwidth    = (1 - axXoffset -axXsep*3)/3;
    
    axYoffsetbottom = 0.12;
    axYoffsettop    = .12;
    axYsep      = 0.02;
    if idfigure==1
        axheight   = (1 - axYoffsettop - axYoffsetbottom -axYsep*3)/4;
    else
        axheight   = (1 - axYoffsettop - axYoffsetbottom -axYsep*2)/3;
    end
    
    if idfigure == 1
        estvalue = ptheta;
        confint = ci;
    else
        estvalue = pthetarep;
        confint = cirep;
    end
    
    axXpos     = axXoffset + (axwidth+axXsep)*[0:2];
    if idfigure == 1
        axYpos     =flip( axYoffsetbottom + (axheight+axYsep)*[0:3]);
    else
        axYpos     =flip( axYoffsetbottom + (axheight+axYsep)*[0:2]);
    end

        
    v = @(x)x(:);
    for iAlg = 1:numel(algused)
        figure(iAlg+3*(idfigure-1)), clf;
        set(gcf,'color',[1 1 1],'units','centimeters','position',[5 5 25 15]);
        for iFold = 1:numel(Folds)
            
            for ith = 1:numel(thused)
                thisp = (estvalue(iFold,:,iAlg,ith));
                thisciL = confint(iFold,:,ith,iAlg,1);
                thisciU = confint(iFold,:,ith,iAlg,2);
                
                Pbar    = [thisciL; thisp-thisciL; thisciU-thisp];
                
                
                %                 subplot(numel(Folds),numel(thused),(iFold-1)*numel(thused)+ith);
                
                h0  = axes('position',[axXpos(ith) axYpos(iFold) axwidth axheight],'box','off');
                bar(1:numel(Permscheme), [Pbar'],'stacked');
                bhandle=  get(gca,'children');
                bhandle(end).FaceColor = [1 1 1];
                bhandle(end).EdgeColor = [1 1 1];
                bhandle(1).FaceColor = [1 1 1]*.95;
                bhandle(1).EdgeColor = [.3 .3 .3];
                bhandle(2).FaceColor = [1 1 1]*.95;
                bhandle(2).EdgeColor = [.3 .3 .3];
                
                bhandle(2).LineWidth = 1;
                bhandle(1).LineWidth = 1;
                ylim([0 max(v(confint(:,:,ith,:,2)))*1.05]);
                xlim([0.5 4.5]);
                set(gca,'Fontsize',10)
                
                if iFold == numel(Folds)
                    set(gca,'Xticklabels',[Permscheme]);
                    set(gca,'XTickLabelRotation',35);
                else
                    set(gca,'XTick',[]) ;
                end
                hold all
                plot(xlim,(thused(ith))*[1 1],'color',[0 0 0],'linewidth',2);
                
                if iFold == 1
                    ht = title(sprintf('${\\alpha} = \\textbf{%1.2f}$',thused(ith)),'Interpreter','latex','fontsize',14);
                end
                
                if ith == 1
                    ylabel(['\textbf{' FoldName{iFold} '}'],'fontsize',14,'interpreter','latex','Fontweight','bold')
                end
                
                if ith == 2 && iFold == 1
                    temp = ylim;
                    if idfigure==1
                        textann     = ['Cross-validation Type-I error rate, ' algused{iAlg}];
                    else
                        textann     = ['Repeated Cross-validation Type-I error rate, ' algused{iAlg}];
                    end
                    ha = annotation('textbox',[.15 .95 .65 .05],'String',textann);
                    ha.FontSize = 20;
                    ha.FontWeight = 'bold';
                    ha.HorizontalAlignment = 'center';
                    ha.Interpreter ='latex';
                    ha.LineStyle = 'none';
                end
            end
            
        end
        
        if saveplots
            set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
            if idfigure==1
                namesave = ['H0_' algused{iAlg} '_CV'];

            else
                namesave = ['H0_' algused{iAlg} '_RepeatedCV'];
            end
            print(gcf,'-dpdf',fullfile(pwd,namesave));
            
        end
    end
end
close all

%%






