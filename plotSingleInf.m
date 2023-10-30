function [] = plotSingleInf(k,citMat,spyr,nYrs,rsdArray,comArray,ocean,Ersd,Ecom,Crsd,Ccom,Irsd,Icom,Rrsd,Rcom,rsdInfTime,comInfTime,figname)

figure(1)
tiledlayout(1,3);
nexttile
f=zeros(size(citMat));
f(citMat>0)=(Ersd(:,end)+Ecom(:,end)+Crsd(:,end)+Ccom(:,end)+Irsd(:,end)+Icom(:,end))/100; 
f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
n = 128; % dark grey => decrease n; light grey => increase n; 
tmp = [n n n]/255; % grey color
mymap = [tmp;jet(255)];     % my colormap
colormap(mymap);
c=colorbar;
hold on 

ylabel(c,'Density','interpreter','latex','FontSize',20);
xticks([]); yticks([]);
% plot(800:1000,ones(201,1)*900,'linewidth',3,'color','black')
% plot(ones(31,1)*800,885:915,'linewidth',2,'color','black')
% plot(ones(31,1)*1000,885:915,'linewidth',2,'color','black')
% text(900,850,'200 km','HorizontalAlignment','center','interpreter','latex','FontSize',15)

ax = gca;
ax.FontSize = 18;

hold off
title('Final infected citrus','interpreter','latex','FontSize',20)

nexttile
f=zeros(size(citMat));
f(citMat>0)=min(rsdInfTime,comInfTime)/365; f(isnan(f))=0; f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
n = 128; % dark grey => decrease n; light grey => increase n; 
tmp = [n n n]/255; % grey color
mymap = [tmp;jet(255)];     % my colormap
colormap(mymap);
c=colorbar;
ylabel(c,'Year','interpreter','latex','FontSize',20);
hold on 

xticks([]); yticks([]);
% plot(800:1000,ones(201,1)*900,'linewidth',3,'color','black')
% plot(ones(31,1)*800,885:915,'linewidth',2,'color','black')
% plot(ones(31,1)*1000,885:915,'linewidth',2,'color','black')
% text(900,850,'200 km','HorizontalAlignment','center','interpreter','latex','FontSize',15)

ax = gca;
ax.FontSize = 18;

hold off
title('Infection time','interpreter','latex','FontSize',20)

nexttile

sumErsd = sum(Ersd);
while sum(isnan(sumErsd),'all')>0
    sumErsd(isnan(sumErsd))=sumErsd(find(isnan(sumErsd))-1);
end
sumCrsd = sum(Crsd);
while sum(isnan(sumCrsd),'all')>0
    sumCrsd(isnan(sumCrsd))=sumCrsd(find(isnan(sumCrsd))-1);
end
sumIrsd = sum(Irsd);
while sum(isnan(sumIrsd),'all')>0
    sumIrsd(isnan(sumIrsd))=sumIrsd(find(isnan(sumIrsd))-1);
end
sumRrsd = sum(Rrsd);
while sum(isnan(sumRrsd),'all')>0
    sumRrsd(isnan(sumRrsd))=sumRrsd(find(isnan(sumRrsd))-1);
end
sumEcom = sum(Ecom);
while sum(isnan(sumEcom),'all')>0
    sumEcom(isnan(sumEcom))=sumEcom(find(isnan(sumEcom))-1);
end
sumCcom = sum(Ccom);
while sum(isnan(sumCcom),'all')>0
    sumCcom(isnan(sumCcom))=sumCcom(find(isnan(sumCcom))-1);
end
sumIcom = sum(Icom);
while sum(isnan(sumIcom),'all')>0
    sumIcom(isnan(sumIcom))=sumIcom(find(isnan(sumIcom))-1);
end
sumRcom = sum(Rcom);
while sum(isnan(sumRcom),'all')>0
    sumRcom(isnan(sumRcom))=sumRcom(find(isnan(sumRcom))-1);
end
totCit=0.1*sum(rsdArray)+sum(comArray);
sumScom = totCit - (sumEcom+sumCcom+sumIcom+sumRcom);

% plot(0:1/spyr:nYrs,0.1*sumErsd/totCit,'--','linewidth',3,'color',[0.4660 0.6740 0.1880])
% hold on
% plot(0:1/spyr:nYrs,0.1*sumCrsd/totCit,'--','linewidth',3,'color',[0.3010 0.7450 0.9330])
% plot(0:1/spyr:nYrs,0.1*sumIrsd/totCit,'--','linewidth',3,'color',[0.6350 0.0780 0.1840])
% plot(0:1/spyr:nYrs,0.1*sumRrsd/totCit,'--','linewidth',3,'color',[0.10 0.10 0.10])

% plot(0:1/spyr:nYrs,sumScom/totCit,'linewidth',3,'color',[0.4660 0.6740 0.1880])
% hold on
% plot(0:1/spyr:nYrs,sumEcom/totCit,'linewidth',3,'color',[0.3010 0.7450 0.9330])
% plot(0:1/spyr:nYrs,sumCcom/totCit,'linewidth',3,'color',[0.9290 0.6940 0.1250])
% plot(0:1/spyr:nYrs,sumIcom/totCit,'linewidth',3,'color',[0.6350 0.0780 0.1840])
% plot(0:1/spyr:nYrs,sumRcom/totCit,'linewidth',3,'color',[0.10 0.10 0.10])

a=area([sumScom' sumEcom' sumCcom' sumIcom' sumRcom']/totCit,'linewidth',1);
a(1).FaceColor=[0.4660 0.6740 0.1880];
a(2).FaceColor=[0.3010 0.7450 0.9330];
a(3).FaceColor=[0.9290 0.6940 0.1250];
a(4).FaceColor=[0.6350 0.0780 0.1840];
a(5).FaceColor=[0.10 0.10 0.10];
ylim([0 1])
xlim([1 length(sumScom)])
xticks([1:5*spyr:nYrs*spyr+1])
xticklabels([0:5:nYrs])

% plot(1/size(detection,1):nYrs/size(detection,1):nYrs,(sum(detection(:,:),2)>0)...
%     *ceil(100*max([max(sumEcom),max(sumCcom),max(sumIcom)])/totCit)/100 ...
%      -ceil(100*max([max(sumEcom),max(sumCcom),max(sumIcom)])/totCit)/200,'*','markersize',3,'color','r')
% ylim([0 ceil(100*max([max(sumEcom),max(sumCcom),max(sumIcom)])/totCit)/100])
ylim([0 1])

ax = gca;
ax.FontSize = 18;

% hl=legend(''E rsd','C rsd','I rsd','R rsd','E com','C com','I com','R com','Detection','Location','northwest');
hl=legend('$S^H$','$E^H$','$C^H$','$I^H$','$R^H$','Location','southwest');
set(hl,'interpreter','latex','FontSize',14);


ylabel('Proportion of total citrus infected','interpreter','latex','FontSize',20);
xlabel('Time (years)','interpreter','latex','FontSize',20);

hold off
% title('Infection over time','interpreter','latex','FontSize',20)

set(1,'paperunits','centimeters');
set(1,'papersize',[40 12]);
set(1,'paperposition',[-2 0 44 12]);
% set(1,'papersize',[35 34]);
% set(1,'paperposition',[0 -5 35 38]);

figname=append(figname,'_',num2str(k))

    print(1,'-dpdf',append('Figures/SingleSims/',figname,'_InfectionMaps.pdf'));
    print(1,'-dtiff',append('Figures/SingleSims/',figname,'_InfectionMaps.tiff'));
    print(1,'-deps',append('Figures/SingleSims/',figname,'_InfectionMaps.eps'));
    savefig(append('Figures/SingleSims/',figname,'_InfectionMaps.fig'));
