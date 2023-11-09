function [] = plotStoryboard(k,citMat,spyr,nYrs,rsdArray,comArray,ocean,Ersd,Ecom,Crsd,Ccom,Irsd,Icom,Rrsd,Rcom,figname)

figure(1)
tiledlayout(3,2);

T = size(Icom,2);

for j=1:5
    nexttile
    t = ceil(j*(T-1)/5)+1;
    f=zeros(size(citMat));
    f(citMat>0)=(Ersd(:,t)+Ecom(:,t)+Crsd(:,t)+Ccom(:,t)+Irsd(:,t)+Icom(:,t))/100; 
    if sum(isnan(f))>0
        disp("ERROR")
    end
    f(citMat>0) = max(f(citMat>0),(ceil(10*max((Ersd(:,end)+Ecom(:,end)+Crsd(:,end)+Ccom(:,end)+Irsd(:,end)+Icom(:,end))/100))/10)*0.005);
    f(ocean~=1) = nan;
    h=imagesc(f(:,end:-1:1)');
    set(h, 'AlphaData', h.CData>=0)
    tmp = [0.5 0.5 0.5]; % map color
    tmp2 = [0.7 0.7 0.7]; % citrus color
    mymap = [tmp; tmp2; jet(255)];     % my colormap
    colormap(mymap);
    c=colorbar;
    clim([0 ceil(10*max((Ersd(:,end)+Ecom(:,end)+Crsd(:,end)+Ccom(:,end)+Irsd(:,end)+Icom(:,end))/100))/10])
    ylabel(c,'Infected citrus density','interpreter','latex','FontSize',20);
    title(['Year ',num2str(nYrs*j/5)],'interpreter','latex','FontSize',20)
    
    xticks([]); yticks([]);
    
    ax = gca;
    ax.FontSize = 18;

end


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

yyaxis right 
a=area([sumRcom' sumIcom' sumCcom' sumEcom' sumScom']/totCit,'linewidth',1);
a(5).FaceColor=[0.4660 0.6740 0.1880];
a(4).FaceColor=[0.3010 0.7450 0.9330];
a(3).FaceColor=[0.9290 0.6940 0.1250];
a(2).FaceColor=[0.6350 0.0780 0.1840];
a(1).FaceColor=[0.10 0.10 0.10];
ylim([0 1.1-ceil((sumScom(end)/totCit)*10)/10])
xlim([1 length(sumScom)])
xticks([1:2*spyr:nYrs*spyr+1])
xticklabels([0:2:nYrs])

% plot(1/size(detection,1):nYrs/size(detection,1):nYrs,(sum(detection(:,:),2)>0)...
%     *ceil(100*max([max(sumEcom),max(sumCcom),max(sumIcom)])/totCit)/100 ...
%      -ceil(100*max([max(sumEcom),max(sumCcom),max(sumIcom)])/totCit)/200,'*','markersize',3,'color','r')
% ylim([0 ceil(100*max([max(sumEcom),max(sumCcom),max(sumIcom)])/totCit)/100])

ax = gca;
ax.FontSize = 18;

% hl=legend(''E rsd','C rsd','I rsd','R rsd','E com','C com','I com','R com','Detection','Location','northwest');
hl=legend(a(5:-1:1),'$S^H$','$E^H$','$C^H$','$I^H$','$R^H$','Location','Northwest');
set(hl,'interpreter','latex','FontSize',14);


ylabel('Proportion of total citrus infected','interpreter','latex','FontSize',20);
xlabel('Time (years)','interpreter','latex','FontSize',20);

ax.YAxis(2).Color = [0 0 0];
ax.YAxis(1).Color = [0 0 0];
yyaxis left
yticks([])
hold off
% title('Infection over time','interpreter','latex','FontSize',20)

set(1,'paperunits','centimeters');
set(1,'papersize',[28 32]);
set(1,'paperposition',[-2.5 -1.5 30 34]);
% set(1,'papersize',[35 34]);
% set(1,'paperposition',[0 -5 35 38]);

figname=append(figname,'_',num2str(k))

    print(1,'-dpdf',append('Figures/SingleSims/',figname,'_Storyboard.pdf'));
    print(1,'-dtiff',append('Figures/SingleSims/',figname,'_Storyboard.tiff'));
    print(1,'-deps',append('Figures/SingleSims/',figname,'_Storyboard.eps'));
    savefig(append('Figures/SingleSims/',figname,'_Storyboard.fig'));
