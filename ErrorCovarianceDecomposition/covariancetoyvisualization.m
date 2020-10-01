clear

xvalcode  = zeros(12,12);

split1 = [1 2 3; 4 5 6];
split2 = [7 8 10; 9 11 12];

level_main_diag = 1;
level_within_testing_fold = 2;
level_between_testing_fold = 3;


same_run_across_splits = 4;
level_across_splits = 5;

for idsp = 1:2
    for idx = split1(idsp,:)
        for idy = split1(idsp,:)
            xvalcode(idx,idy) = level_within_testing_fold;
        end
    end
end
for idx = split1(1,:)
    for idy = split1(2,:)
        xvalcode(idx,idy)=  level_between_testing_fold;
        xvalcode(idy,idx) = xvalcode(idx,idy);
    end
end

for idsp = 1:2
    for idx = split2(idsp,:)
        for idy = split2(idsp,:)
            xvalcode(idx,idy) = level_within_testing_fold;
        end
    end
end
for idx = split2(1,:)
    for idy = split2(2,:)
        xvalcode(idx,idy)=  level_between_testing_fold;
        xvalcode(idy,idx) = xvalcode(idx,idy);
    end
end

for idx = 1:size(xvalcode,2)
    xvalcode(idx,idx) = level_main_diag;
end

xvalcode(1:6,7:end) = level_across_splits;
xvalcode(7:end,1:6) = level_across_splits;

for idx  = 1:6
    xvalcode(idx,idx+6) = same_run_across_splits;
    xvalcode(idx+6,idx) = same_run_across_splits;
end




%%
figure(1), clf
set(gcf,'units','centimeters','position',[2 2 20 20],'Color',[1 1 1]);
hold all
h = axes('units','normalized','position',[.05 .05 .9 .9]);
h.YAxis.Visible = 'off';
h.XAxis.Visible = 'off';


colorcode = [1 1 0; 1 .5 0; 1 .2 .2; 0.5 .5 1; .5 .7 .2];
diffstep = 0.008;
xgrid = linspace(0,1,13); xstep  = diff(xgrid(1:2))-diffstep;
ygrid = linspace(1,0,13); ystep  = diff(ygrid([2 1]))-diffstep;
for idx = 1:12
    for idy = 1:12
        thispos = [xgrid(idx) ygrid(idy) xstep ystep];
        drawCurvedRectangle(thispos,colorcode(xvalcode(idx,idy),:));
    end
end
axis square




xlimval = xlim;
ylimval = ylim;

positions = linspace(-.1,1.25,6); 
strings = {['\textsl{Within} iteration' newline '\textbf{Same Unit}'],...
    ['\textsl{Within} iteration' newline '\textbf{Same Fold}'],...
    ['\textsl{Within} iteration' newline '\textbf{Different Fold}'],...
    ['\textsl{Between} iterations' newline '\textbf{Same Unit}'],...
    ['\textsl{Between} iterations' newline '\textbf{Different Units}']};
for idx = 1:5
    rectangle('position',[positions(idx) -.05 xstep xstep],...% 1.23
        'FaceColor',colorcode(idx,:),'EdgeColor',colorcode(idx,:)*.5);
    t = text(positions(idx)+.03,-.1, strings{idx} ,'Clipping','off','Fontsize',12,...% 1.18
        'Interpreter','latex','fontweight','bold',...
        'HorizontalAlignment','center'); 
end

% rectangle('position',[0 1.12 xstep xstep],'FaceColor',colorcode(1,:),'EdgeColor',colorcode(1,:)*.5);
% text(xstep+.01,1.16,'Main \\ Diagonal','fontsize',12,'interpreter','latex')
axis square
xlim(xlimval+[-0.15 +0.15]);
ylim(ylimval+[-0.15 +0.15]);


annotation('doublearrow',[.1 .1],[.5 .78]) % 0.12
h2 = text(-0.12,.7,'\textbf{Iteration 1}','interpreter','latex','fontsize',18,'Fontweight','bold'); %-.1
set(h2,'rotation',90)

annotation('doublearrow',[.1 .1],[.5 .78]-.3)
h2 = text(-0.12,.2,'\textbf{Iteration 2}','interpreter','latex','fontsize',18,'Fontweight','bold'); %-.1
set(h2,'rotation',90)


annotation('doublearrow',[.15 .49],[.83 .83]) % .813
text(.17,1.18,'\textbf{Iteration 1}','interpreter','latex','fontsize',18,'Fontweight','bold'); % 1.17

annotation('doublearrow',[.51 .84],[.83 .83])
text(.65,1.18,'\textbf{Iteration 2}','interpreter','latex','fontsize',18,'Fontweight','bold');

tmp = linspace(0.03,.945,12);
runcode = repmat([1:6],1,2);
for idx = 1:12
    text(tmp(idx),1.11,num2str(runcode(idx)),'Fontsize',15,'interpreter','latex')
end

tmp = linspace(1.04,0.12,12);
for idx = 1:12
    text(-.04,tmp(idx),num2str(runcode(idx)),'Fontsize',15,'interpreter','latex')
end

annotation('textbox',[.08 .895 .2 .03],'String','\textsl{Iteration 1}:','interpreter','latex', 'VerticalAlignment','Middle',...
    'horizontalalignment','center','fontsize',18,'linestyle','None');

annotation('textbox',[.25 .875 .23 .05],'String',{ ...
    '\textsl{Fold 1}: 1 2 3'; '\textsl{Fold 2}: 4 5 6'},'interpreter','latex', 'VerticalAlignment','Middle',...
    'horizontalalignment','center','fontsize',16,'linestyle','None');

annotation('textbox',[.53 .895 .2 .03],'String','\textsl{Iteration 2}:','interpreter','latex', 'VerticalAlignment','Middle',...
    'horizontalalignment','center','fontsize',18,'linestyle','None');

annotation('textbox',[.70 .875 .23 .05],'String',{ ...
    '\textsl{Fold 1}: 1 2 4'; '\textsl{Fold 2}: 3 5 6'},'interpreter','latex', 'VerticalAlignment','Middle',...
    'horizontalalignment','center','fontsize',16,'linestyle','None');


annotation('textbox',[.0 .945 1 .05],'String', ...
    '\textbf{Covariance Decomposition of repeated Cross-Validation}','interpreter','latex',...
    'VerticalAlignment','Middle',...
    'horizontalalignment','center','fontsize',19,'linestyle','None');

dirdata = pwd;
set(gcf,'PaperOrientation','portrait','PaperPositionMode','auto');
print(gcf,'-dpdf',fullfile(dirdata,'CovarianceVisualization.pdf'));

function drawCurvedRectangle(position,color)

rectangle('position',position,'Curvature',0,'FaceColor',color,'EdgeColor',[.3 .3  .3],'linewidth',3.5);

end