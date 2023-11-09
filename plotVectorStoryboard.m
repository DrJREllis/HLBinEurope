function [] = plotVectorStoryboard(k,citMat,spyr,nYrs,rsdArray,comArray,ocean,Nrsd,Ncom,adjNrsd,adjNcom,figname)

figure(1)
tiledlayout(3,2);

T = size(Ncom,2);
tArray = [1,2,3,5,10]*spyr+1;

for j=1:5
    nexttile
%     t = ceil(j*(T-1)/5)+1;
    t=tArray(j);
    f=zeros(size(citMat));
    f(citMat>0)=(Nrsd(:,t).*adjNrsd+Ncom(:,t).*adjNcom); 
    if sum(isnan(f))>0
        disp("ERROR")
    end
   f(citMat>0) = max(f(citMat>0),(ceil(100*max((Nrsd(:,end).*adjNrsd+Ncom(:,end).*adjNcom)))/100)*0.000005);
    f(ocean~=1) = nan;
    h=imagesc(f(:,end:-1:1)');
    set(h, 'AlphaData', h.CData>=0)
    tmp = [0.5 0.5 0.5]; % map color
    tmp2 = [0.7 0.7 0.7]; % citrus color
    mymap = [tmp; tmp2; jet(255000)];     % my colormap
    colormap(mymap);
    c=colorbar;
    clim([0 ceil(100*max((Nrsd(:,end).*adjNrsd+Ncom(:,end).*adjNcom)))/100])
    ylabel(c,'Vector density','interpreter','latex','FontSize',20);
    title(['Year ',num2str((tArray(j)-1)/spyr)],'interpreter','latex','FontSize',20)
    
    xticks([]); yticks([]);
    
    ax = gca;
    ax.FontSize = 18;

end


nexttile

sumNrsd = sum(Nrsd);
while sum(isnan(sumNrsd),'all')>0
    sumNrsd(isnan(sumNrsd))=sumNrsd(find(isnan(sumNrsd))-1);
end
sumNcom = sum(Ncom);
while sum(isnan(sumNcom),'all')>0
    sumNcom(isnan(sumNcom))=sumNcom(find(isnan(sumNcom))-1);
end

plot((sumNcom+sumNrsd)./(sum(rsdArray>0)+sum(comArray>0)),'LineWidth',3,'color',"#A2142F")
ylim([0 1])
yyaxis right 
ylim([0 1])
xlim([1 length(sumNcom)])
xticks([1:2*spyr:nYrs*spyr+1])
xticklabels([0:2:nYrs])

% plot(1/size(detection,1):nYrs/size(detection,1):nYrs,(sum(detection(:,:),2)>0)...
%     *ceil(100*max([max(sumEcom),max(sumCcom),max(sumIcom)])/totCit)/100 ...
%      -ceil(100*max([max(sumEcom),max(sumCcom),max(sumIcom)])/totCit)/200,'*','markersize',3,'color','r')
% ylim([0 ceil(100*max([max(sumEcom),max(sumCcom),max(sumIcom)])/totCit)/100])

ax = gca;
ax.FontSize = 18;


ylabel('Proportion of total citrus infested','interpreter','latex','FontSize',20);
xlabel('Time (years)','interpreter','latex','FontSize',20);

ax.YAxis(2).Color = [0 0 0];
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

    print(1,'-dpdf',append('Figures/SingleSims/',figname,'_VectorStoryboard.pdf'));
    print(1,'-dtiff',append('Figures/SingleSims/',figname,'_VectorStoryboard.tiff'));
    print(1,'-deps',append('Figures/SingleSims/',figname,'_VectorStoryboard.eps'));
    savefig(append('Figures/SingleSims/',figname,'_VectorStoryboard.fig'));
