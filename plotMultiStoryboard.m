function [] = plotMultiStoryboard(citMat,ocean,spyr,rsdArray,comArray,...
    NrsdFin,NcomFin,IrsdFin,IcomFin,ErsdFin,EcomFin,CrsdFin,CcomFin,numInf,citInf,cellInf,...
    rsdInfTime,comInfTime,rsdVecTime,comVecTime,tMax,figname)


% Easting = (2635:3770);
% Northing = (1468:2468);

Nreps = size(numInf,2);

citNums2 = find(citMat>0);



figure(1)
tiledlayout(3,2);

tArray = [1,2,3,5,10]*spyr+1;

for j=1:5
    nexttile
%     t = ceil(j*(T-1)/5)+1;
    t=tArray(j);
f=zeros(size(citMat));
f(citMat>0)=100*sum(NrsdFin(:,t,:)+NcomFin(:,t,:)>0,3)/Nreps; f(ocean~=1) = nan;
f(citMat>0) = max(f(citMat>0),max(max(f))*0.006);
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
    title(['Year ',num2str((tArray(j)-1)/spyr)],'interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;
ylabel(c, 'Probability of infestation (\%)','interpreter','latex','FontSize',20);

end

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
set(1,'papersize',[28 32]);
set(1,'paperposition',[-2.5 -1.5 30 34]);
% set(1,'papersize',[35 34]);
% set(1,'paperposition',[0 -5 35 38]);


print(1,'-dpdf',append('Figures/',figname,'_VectorTotalsStoryboard.pdf'));
print(1,'-dtiff',append('Figures/',figname,'_VectorTotalsStoryboard.tiff'));
print(1,'-deps',append('Figures/',figname,'_VectorTotalsStoryboard.eps'));
savefig(append('Figures/',figname,'_VectorTotalsStoryboard.fig'));

close 

    %%

figure(1)
tiledlayout(3,2);

tArray = [1,2,3,5,10]*spyr+1;

for j=1:5
    nexttile
    t=tArray(j);
f=zeros(size(citMat));
f(citMat>0)=100*sum((IcomFin(:,t,:)+IrsdFin(:,t,:)+CcomFin(:,t,:)+CrsdFin(:,t,:)+EcomFin(:,t,:)+ErsdFin(:,t,:))>0,3)/Nreps; f(isnan(f))=0;
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
xticks([]); yticks([]);
title(['Year ',num2str((tArray(j)-1)/spyr)],'interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;
ylabel(c, 'probability of infection (\%)' ,'interpreter','latex','FontSize',20);
end

nexttile

f=zeros(size(citMat));
f(citMat>0)=mean(min(rsdInfTime,comInfTime),2,'omitnan')/365; f(ocean~=1) = nan;
f(citMat>0) = max(f(citMat>0),(tMax/365)*0.005);
f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
tmp = [0.5 0.5 0.5]; % map color
tmp2 = [0.7 0.7 0.7]; % citrus color
mymap = [tmp; tmp2; jet(255)];     % my colormap
colormap(mymap)
c=colorbar;
clim([0 tMax/365])
hold on
xticks([]); yticks([]);
plot(800:1000,ones(201,1)*900,'linewidth',3,'color','black')
plot(ones(31,1)*800,885:915,'linewidth',2,'color','black')
plot(ones(31,1)*1000,885:915,'linewidth',2,'color','black')
text(900,850,'200 km','HorizontalAlignment','center','interpreter','latex','FontSize',15)
hold off
% xlim([min(Ix1)-20, max(Ix1)+20])
% ylim([min(gridY-Iy1)-20, max(gridY-Iy1)+20])
title('Mean time of infection','interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;
 ylabel(c, 'Year','interpreter','latex','FontSize',20);


set(1,'paperunits','centimeters');
set(1,'papersize',[28 32]);
set(1,'paperposition',[-2.5 -1.5 30 34]);
% set(1,'papersize',[35 34]);
% set(1,'paperposition',[0 -5 35 38]);


print(1,'-dpdf',append('Figures/',figname,'_InfTotalsStoryboard.pdf'));
print(1,'-dtiff',append('Figures/',figname,'_InfTotalsStoryboard.tiff'));
print(1,'-deps',append('Figures/',figname,'_InfTotalsStoryboard.eps'));
savefig(append('Figures/',figname,'_InfTotalsStoryboard.fig'));

close
