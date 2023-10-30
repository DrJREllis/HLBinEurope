function [] = plotSingleVec(k,citMat,spyr,nYrs,realRsdCit,realComCit,adjNrsd,adjNcom,ocean,Arsd,Nrsd,Acom,Ncom,rsdVecTime,comVecTime,figname)

figure(1)
tiledlayout(1,3);
nexttile
f=zeros(size(citMat));
f(citMat>0)=Nrsd(:,end).*adjNrsd+Ncom(:,end).*adjNcom; 
f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
n = 128; % dark grey => decrease n; light grey => increase n; 
tmp = [n n n]/255; % grey color
mymap = [tmp;jet(255)];     % my colormap
colormap(mymap);
c = colorbar;

ylabel(c,'Density','interpreter','latex','FontSize',20);

xticks([]); yticks([]);

hold off
title('Final vector density','interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;

nexttile
f=zeros(size(citMat));
f(citMat>0)=min(rsdVecTime,comVecTime)/365; f(isnan(f))=0; f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
n = 128; % dark grey => decrease n; light grey => increase n; 
tmp = [n n n]/255; % grey color
mymap = [tmp;jet(255)];     % my colormap
colormap(mymap);
c=colorbar;

ylabel(c,'Year','interpreter','latex','FontSize',20);

xticks([]); yticks([]);


hold off
title('Arrival time','interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;

nexttile

Arsd(Arsd>1)=1; Acom(Acom>1)=1;
sumArsd = sum(Arsd.*realRsdCit);
while sum(isnan(sumArsd),'all')>0
    sumArsd(isnan(sumArsd))=sumArsd(find(isnan(sumArsd))-1);
end
sumNrsd = sum(Nrsd.*realRsdCit);
while sum(isnan(sumNrsd),'all')>0
    sumNrsd(isnan(sumNrsd))=sumNrsd(find(isnan(sumNrsd))-1);
end
sumAcom = sum(Acom.*realComCit);
while sum(isnan(sumAcom),'all')>0
    sumAcom(isnan(sumAcom))=sumAcom(find(isnan(sumAcom))-1);
end
sumNcom = sum(Ncom.*realComCit);
while sum(isnan(sumNcom),'all')>0
    sumNcom(isnan(sumNcom))=sumNcom(find(isnan(sumNcom))-1);
end
totCit=sum(realRsdCit)+sum(realComCit);

plot(0:1/spyr:nYrs,sumArsd/totCit,'--','linewidth',3,'color',[0.4660 0.6740 0.1880])
hold on
plot(0:1/spyr:nYrs,sumNrsd/totCit,'--','linewidth',3,'color',[0.3010 0.7450 0.9330])

plot(0:1/spyr:nYrs,sumAcom/totCit,'linewidth',3,'color',[0.4660 0.6740 0.1880])
plot(0:1/spyr:nYrs,sumNcom/totCit,'linewidth',3,'color',[0.3010 0.7450 0.9330])

ax = gca;
ax.FontSize = 18;

% hl=legend('Arsd','Nrsd','Acom','Ncom','Location','northwest');
hl=legend('','','$E^V$','$N^V$','Location','northeast');
set(hl,'interpreter','latex','FontSize',16);

ylabel('Proportion of cells infested','interpreter','latex','FontSize',20);
xlabel('Time (years)','interpreter','latex','FontSize',20);

hold off
% title('Vector population over time','interpreter','latex','FontSize',20)


set(1,'paperunits','centimeters');
set(1,'papersize',[40 12]);
set(1,'paperposition',[-2 0 44 12]);
% set(1,'papersize',[35 34]);
% set(1,'paperposition',[0 -5 35 38]);

figname=append(figname,'_',num2str(k));

    print(1,'-dpdf',append('Figures/SingleSims/',figname,'_VectorMaps.pdf'));
    print(1,'-dtiff',append('Figures/SingleSims/',figname,'_VectorMaps.tiff'));
    print(1,'-deps',append('Figures/SingleSims/',figname,'_VectorMaps.eps'));
    savefig(append('Figures/SingleSims/',figname,'_VectorMaps.fig'));
