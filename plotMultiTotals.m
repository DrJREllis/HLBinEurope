function [] = plotMultiTotals(citMat,ocean,rsdArray,comArray,...
    NrsdFin,NcomFin,IrsdFin,IcomFin,ErsdFin,EcomFin,CrsdFin,CcomFin,numInf,citInf,cellInf,...
    rsdInfTime,comInfTime,rsdVecTime,comVecTime,tMax,figname)


% Easting = (2635:3770);
% Northing = (1468:2468);

Nreps = size(numInf,2);

citNums2 = find(citMat>0);

n = 128; % dark grey => decrease n; light grey => increase n; 
tmp = [n n n]/255; % grey color
mymap = [tmp;jet(255)];     % my colormap

figure(1)
tiledlayout(1,2);
nexttile

f=zeros(size(citMat));
f(citMat>0)=100*sum(NrsdFin+NcomFin>0,3)/Nreps; f(ocean~=1) = nan;
f(citMat>0) = max(f(citMat>0),100*0.006);
f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
tmp = [0.5 0.5 0.5]; % map color
tmp2 = [0.7 0.7 0.7]; % citrus color
mymap = [tmp; tmp2; jet(255)];     % my colormap
colormap(mymap)
c=colorbar;
clim([0 100])
hold on
xticks([]); yticks([]);
plot(800:1000,ones(201,1)*900,'linewidth',3,'color','black')
plot(ones(31,1)*800,885:915,'linewidth',2,'color','black')
plot(ones(31,1)*1000,885:915,'linewidth',2,'color','black')
text(900,850,'200 km','HorizontalAlignment','center','interpreter','latex','FontSize',15)
hold off
% xlim([min(Ix1)-20, max(Ix1)+20])
% ylim([min(gridY-Iy1)-20, max(gridY-Iy1)+20])
title('Probability of infestation','interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;
ylabel(c, '\%','interpreter','latex','FontSize',20);


nexttile

f=zeros(size(citMat));
f(citMat>0)=mean(min(rsdVecTime,comVecTime),2,'omitnan')/365; f(ocean~=1) = nan;
f(citMat>0) = max(f(citMat>0),6*0.005);
f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
tmp = [0.5 0.5 0.5]; % map color
tmp2 = [0.7 0.7 0.7]; % citrus color
mymap = [tmp; tmp2; jet(255)];     % my colormap
colormap(mymap)
c=colorbar;
clim([0 6])
hold on
xticks([]); yticks([]);
plot(800:1000,ones(201,1)*900,'linewidth',3,'color','black')
plot(ones(31,1)*800,885:915,'linewidth',2,'color','black')
plot(ones(31,1)*1000,885:915,'linewidth',2,'color','black')
text(900,850,'200 km','HorizontalAlignment','center','interpreter','latex','FontSize',15)
hold off
% xlim([min(Ix1)-20, max(Ix1)+20])
% ylim([min(gridY-Iy1)-20, max(gridY-Iy1)+20])
title('Mean time of infestation','interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;
 ylabel(c, 'Years','interpreter','latex','FontSize',20);



set(1,'paperunits','centimeters');
set(1,'papersize',[28 12]);
set(1,'paperposition',[-2 0 30 12]);
% set(1,'papersize',[35 34]);
% set(1,'paperposition',[0 -5 35 38]);


print(1,'-dpdf',append('Figures/',figname,'_VectorTotals.pdf'));
print(1,'-dtiff',append('Figures/',figname,'_VectorTotals.tiff'));
print(1,'-deps',append('Figures/',figname,'_VectorTotals.eps'));
savefig(append('Figures/',figname,'_VectorTotals.fig'));

close 

    %%

figure(1)
tiledlayout(1,2);
nexttile

f=zeros(size(citMat));
f(citMat>0)=100*sum((IcomFin+IrsdFin+CcomFin+CrsdFin+EcomFin+ErsdFin)>0,3)/Nreps; f(isnan(f))=0;
f(citMat>0) = max(f(citMat>0),100*0.0006);
f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
tmp = [0.5 0.5 0.5]; % map color
tmp2 = [0.7 0.7 0.7]; % citrus color
mymap = [tmp; tmp2; jet(2550)];     % my colormap
colormap(mymap)
c=colorbar;
clim([0 100])
hold on
xticks([]); yticks([]);
plot(800:1000,ones(201,1)*900,'linewidth',3,'color','black')
plot(ones(31,1)*800,885:915,'linewidth',2,'color','black')
plot(ones(31,1)*1000,885:915,'linewidth',2,'color','black')
text(900,850,'200 km','HorizontalAlignment','center','interpreter','latex','FontSize',15)
hold off
% xlim([min(Ix1)-10, max(Ix1)+10])
% ylim([min(gridY-Iy1)-10, max(gridY-Iy1)+10])
title('Probability of HLB infection','interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;
ylabel(c, '\%','interpreter','latex','FontSize',20);


nexttile

f=zeros(size(citMat));
f(citMat>0)=mean(min(rsdInfTime,comInfTime),2,'omitnan')/365; f(ocean~=1) = nan;
f(citMat>0) = max(f(citMat>0),(tMax/365)*0.0005);
f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
tmp = [0.5 0.5 0.5]; % map color
tmp2 = [0.7 0.7 0.7]; % citrus color
colormap(mymap)
c=colorbar;
clim([0 tMax/365])
hold on
xticks([]); yticks([]);
% plot(800:1000,ones(201,1)*900,'linewidth',3,'color','black')
% plot(ones(31,1)*800,885:915,'linewidth',2,'color','black')
% plot(ones(31,1)*1000,885:915,'linewidth',2,'color','black')
% text(900,850,'200 km','HorizontalAlignment','center','interpreter','latex','FontSize',15)
hold off
% xlim([min(Ix1)-20, max(Ix1)+20])
% ylim([min(gridY-Iy1)-20, max(gridY-Iy1)+20])
title('Mean time of infection','interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;
 ylabel(c, 'Years','interpreter','latex','FontSize',20);


set(1,'paperunits','centimeters');
set(1,'papersize',[28 12]);
set(1,'paperposition',[-2 0 30 12]);
% set(1,'papersize',[35 34]);
% set(1,'paperposition',[0 -5 35 38]);


print(1,'-dpdf',append('Figures/',figname,'_InfTotals.pdf'));
print(1,'-dtiff',append('Figures/',figname,'_InfTotals.tiff'));
print(1,'-deps',append('Figures/',figname,'_InfTotals.eps'));
savefig(append('Figures/',figname,'_InfTotals.fig'));

close

%%

figure(1)

tInt = (1:tMax)/365;
% numInf(1,:)=zeros(1,size(numInf,2));
% citInf(1,:)=zeros(1,size(citInf,2));
% cellInf(1,:)=zeros(1,size(cellInf,2));
for k=1:Nreps
while sum(isnan(numInf(:,k)),'all')>0
    numInf(isnan(numInf(:,k)),k)=numInf(find(isnan(numInf(:,k)))-1,k);
end
while sum(isnan(cellInf(:,k)),'all')>0
    cellInf(isnan(cellInf(:,k)),k)=cellInf(find(isnan(cellInf(:,k)))-1,k);
end
while sum(isnan(citInf(:,k)),'all')>0
    citInf(isnan(citInf(:,k)),k)=citInf(find(isnan(citInf(:,k)))-1,k);
end
end
% meanInfested = mean(numInfested,2);
% ind = find(tInt'<=20 & ~isnan(meanInfested));
% numInfested=sort(numInfested,'descend');

hold on
patch([tInt fliplr(tInt)], 100*[prctile(numInf(1:tMax,:),2.5,2)' flipud(prctile(numInf(1:tMax,:),97.5,2))']/sum((rsdArray>0)+(comArray>0)),...
    'r','FaceAlpha',0.3,'EdgeColor','none')
patch([tInt fliplr(tInt)], 100*[prctile(cellInf(1:tMax,:),2.5,2)' flipud(prctile(cellInf(1:tMax,:),97.5,2))']/sum(rsdArray+comArray>0),...
    [0.9290 0.6940 0.1250],'FaceAlpha',0.3,'EdgeColor','none')
patch([tInt fliplr(tInt)], 100*[prctile(citInf(1:tMax,:),2.5,2)' flipud(prctile(citInf(1:tMax,:),97.5,2))']/sum((rsdArray)+(comArray)),...
    [0 0.4470 0.7410],'FaceAlpha',0.3,'EdgeColor','none')


plot(tInt,100*median(numInf(1:tMax,:),2)/sum((rsdArray>0)+(comArray>0)),'linewidth',3,'color','r')
plot(tInt,100*median(cellInf(1:tMax,:),2)/sum(rsdArray+comArray>0),'linewidth',3,'color',[0.9290 0.6940 0.1250])
plot(tInt,100*median(citInf(1:tMax,:),2)/sum(rsdArray+comArray),'linewidth',3,'color',[0 0.4470 0.7410])

% plot(tInt,numInf(1:tMax,:)/sum((rsdArray>0)+(comArray>0)),'linewidth',1,'color',[0 0.4470 0.7410 0.2])
% plot(tInt,cellInf(1:tMax,:)/sum(rsdArray+comArray>0),'linewidth',1,'color',[0.9290 0.6940 0.1250 0.2])

hold off
% title('Total spread','interpreter','latex','FontSize',20)
% xlim([0 t(tStep)/365])
% ylim([0 0.3])
hl=legend('Vector','Infected cells','Infected citrus','Location','northwest');
% legend('Location','eastoutside')
set(hl,'interpreter','latex','FontSize',18);
ax = gca;
ax.FontSize = 18;
ylabel('\%','interpreter','latex','FontSize',20)
xlabel('Years','interpreter','latex','FontSize',20);

% annotation('textarrow',[0.275 0.275],[0.668 0.718],'String','N')
% annotation('textarrow',[0.275 0.275],[0.195 0.245],'String','N')



set(1,'paperunits','centimeters');
set(1,'papersize',[16 12]);
set(1,'paperposition',[0 0 16 12]);
% set(1,'papersize',[35 34]);
% set(1,'paperposition',[0 -5 35 38]);


print(1,'-dpdf',append('Figures/',figname,'_DPGs.pdf'));
print(1,'-dtiff',append('Figures/',figname,'_DPGs.tiff'));
print(1,'-deps',append('Figures/',figname,'_DPGs.eps'));
savefig(append('Figures/',figname,'_DPGs.fig'));