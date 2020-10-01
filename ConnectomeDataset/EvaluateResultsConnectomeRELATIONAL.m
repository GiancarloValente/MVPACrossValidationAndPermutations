clear variables
saveplots = false;

dirdata = pwd;
load(fullfile(dirdata,'ResultsRELATIONAL_1000perm_100iter.mat'));
e       = zeros(2,4,size(res,2),size(res,1));
pow     = zeros(2,4,size(res,2),size(res,1));

for ind = 1:size(res,1)
    for ind1 = 1:size(res,2)
        
        e(1,:,ind1,ind)      = mean(res(ind,ind1).E)/40;
        e(2,:,ind1,ind)      = mean(res(ind,ind1).Er)/40./[20 20 20 1];
        pow(1,:,ind1,ind)    = mean(res(ind,ind1).p < 0.05);
        pow(2,:,ind1,ind)    = mean(res(ind,ind1).prep<0.05);
    end
end




%%
figure(1);
clf
set(gcf,'color',[1 1 1],'units','centimeters','position',[5 2 20 10]);

axXsep  = .05;
axXoffset  = .08;
axXoffsetright = .03;
axwidth    = (1 - axXoffset- axXoffsetright -axXsep*2)/3;


axYoffsetbottom = 0.15;
axYoffsettop    = .12;
axYsep      = 0.02;
axheight   = (1 - axYoffsettop - axYoffsetbottom -axYsep*0)/1;

axXpos     = axXoffset + (axwidth+axXsep)*[0:2];
axYpos     =flip( axYoffsetbottom + (axheight+axYsep)*0);
colorofbars = [  0.8500    0.3250    0.0980;  0    0.4470    0.7410];
titles    = {'Scenario I','Scenario II','Scenario III'};
for ind = 1:size(res,1)
    
    h0      = axes('position',[axXpos(ind) axYpos axwidth axheight],'box','off');
    thise   = e(:,:,:,ind);
    meane   = mean(thise,3);
    thisq   = quantile(bootstrp(1e4,@mean,shiftdim(thise,2)),[.025 .975]);
    
    x       = sort([1:3:11 2:3:8]);
    hold all
    for idx = 1:numel(x)
        bar(x(idx),meane(idx),'barwidth',.9,'facecolor',colorofbars(double(mod(idx,2)==1)+1,:));
        errorbar(x(idx),meane(idx),meane(idx)-thisq(1,idx),thisq(2,idx)-meane(idx),'LineWidth',2,...
            'color',[0 0 0]);
    end
    set(gca,'fontsize',8,'Xtick',sort(x(:)),...
        'xticklabels',{'SplitHalf','SplitHalf Rep','5-Fold','5-Fold Rep','10-Fold','10-Fold rep','LOO'},...
        'XTickLabelRotation',35);
    ylim([0 0.55]);
    plot(xlim,[0.5 .5],'k--');
    if ind==1
        ylabel('Average Error rate');
    end
    title(titles{ind})
    
    
end
if saveplots

 set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
 print(gcf,'-dpdf',fullfile(dirdata,'ErrorsConnectomeRelational'));
end


figure(2);
clf
set(gcf,'color',[1 1 1],'units','centimeters','position',[5 2 20 10]);


for ind = 1:size(res,1)
    
    h0      = axes('position',[axXpos(ind) axYpos axwidth axheight],'box','off');
    thisp   = pow(:,:,:,ind);
    meanp   = mean(thisp,3);
    thisq   = quantile(bootstrp(1e4,@mean,shiftdim(thisp,2)),[.025 .975]);
    
    x       = sort([1:3:11 2:3:8]);
    hold all
    for idx = 1:numel(x)
        bar(x(idx),meanp(idx),'barwidth',.9,'facecolor',colorofbars(double(mod(idx,2)==1)+1,:));
        errorbar(x(idx),meanp(idx),meanp(idx)-thisq(1,idx),thisq(2,idx)-meanp(idx),'LineWidth',2,...
            'color',[0 0 0]);
    end
    set(gca,'fontsize',8,'Xtick',sort(x(:)),...
        'xticklabels',{'SplitHalf','SplitHalf Rep','5-Fold','5-Fold Rep','10-Fold','10-Fold rep','LOO'},...
        'XTickLabelRotation',35);
    ylim([0 1]);
    plot(xlim,[0.05 .05],'k--');
    
    
     if ind==1
        ylabel('Percentage significant datasets');
    end
    title(titles{ind});
end
if saveplots

 set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
 print(gcf,'-dpdf',fullfile(dirdata,'PowerConnectomeRelational'));
end


%%
figure(3);
clf
set(gcf,'color',[1 1 1],'units','centimeters','position',[1 1 20 20]);

logit= @(x)log(x./(1-x));
invlogit = @(x)1./(1+exp(-x));
nrows = 7;
axXsep  = .01;
axXoffset  = .02;
axXoffsetright = .05;
axwidth    = (1 - axXoffset- axXoffsetright -axXsep*nrows)/nrows;


axYoffsetbottom = 0.03;
axYoffsettop    = .05;
axYsep      = 0.02;
axheight   = (1 - axYoffsettop - axYoffsetbottom -axYsep*nrows)/nrows;

axXpos     = axXoffset + (axwidth+axXsep)*[0:nrows-1];
axYpos     =flip( axYoffsetbottom + (axheight+axYsep)*[0:nrows-1]);

tmp         = [[ones(4,1); 2*(ones(3,1))] [1:4 1:3]',];


schemenames ={'Split Half','5-Fold','10-Fold','LOO','Split Half Rep','5-Fold Rep','10-Fold rep'};




for idx1 = 1:nrows
    for idx2 = [1:(idx1-1)]
        h0      = axes('position',[axXpos(idx1) axYpos(idx2) axwidth axheight],'box','off');
        % adding here very tiny displacements to aid visualization. Power is first transformed with logit transformation, noise with  very small standard deviation is added, and afterwards brought back to 0-1 with an inverse logit transform 
        plot(invlogit(logit(squeeze((pow(tmp(idx1,1),tmp(idx1,2),:,:))))+.05*randn(100,3)),...
            invlogit(logit(squeeze((pow(tmp(idx2,1),tmp(idx2,2),:,:))))+.05*randn(100,3)),'.')
        hold all
        plot([0 1],[0 1]);
        xlim([0 1]); ylim([0 1]);
        axis square
        set(h0,'XTick',[0 1],'Ytick',[0 1],'FontSize',6,'YAxisLocation','right');
        if idx1 == 5 || idx2 == 5
            h0.Color = [.95 1 .95];
        end


        if idx2==1
            title(schemenames{idx1},'Fontsize',10);
        end
        if idx1 == nrows
            ylabel(schemenames{idx2},'Fontweight','bold','Fontsize',10)

        end
    end
    h0      = axes('position',[axXpos(idx1) axYpos(idx1) axwidth axheight],'box','on');
    hold all
    for idx3 = 1:size(pow,4)
        histogram(pow(tmp(idx1,1),tmp(idx1,2),:,idx3),linspace(0,1,31));
    end
    set(h0,'Ytick',[0],'FontSize',6);
    axis square
    if idx1 == nrows
        set(gca,'YAxisLocation','right')
        ylabel(schemenames{idx2+1},'Fontweight','bold','Fontsize',10)
     
    end

end



if saveplots
    set(gcf,'InvertHardCopy','off')
 set(gcf,'PaperOrientation','portrait','PaperPositionMode','auto');
 print(gcf,'-dpdf',fullfile(dirdata,'CorrPlotPowerConnectomeRelational'));
end

