function [] = plotMaps(citMat,ocean,rsdArray,comArray,figname)

figure(1)
tiledlayout(1,3);
nexttile
f=zeros(size(citMat));
f(citMat>0)=comArray+rsdArray; f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
n = 128; % dark grey => decrease n; light grey => increase n; 
tmp = [n n n]/255; % grey color
mymap = [tmp;jet(255)];     % my colormap
colormap(mymap);
colorbar
hold on 

xticks([]); yticks([]);
plot(800:1000,ones(201,1)*900,'linewidth',3,'color','black')
plot(ones(31,1)*800,885:915,'linewidth',2,'color','black')
plot(ones(31,1)*1000,885:915,'linewidth',2,'color','black')
text(900,850,'200 km','HorizontalAlignment','center','interpreter','latex','FontSize',15)


hold off
title('Total citrus','interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;

nexttile
f=zeros(size(citMat));
f(citMat>0)=rsdArray; f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
n = 128; % dark grey => decrease n; light grey => increase n; 
tmp = [n n n]/255; % grey color
mymap = [tmp;jet(255)];     % my colormap
colormap(mymap);
colorbar
hold on 

xticks([]); yticks([]);
plot(800:1000,ones(201,1)*900,'linewidth',3,'color','black')
plot(ones(31,1)*800,885:915,'linewidth',2,'color','black')
plot(ones(31,1)*1000,885:915,'linewidth',2,'color','black')
text(900,850,'200 km','HorizontalAlignment','center','interpreter','latex','FontSize',15)


hold off
title('Residential citrus','interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;

nexttile
f=zeros(size(citMat));
f(citMat>0)=comArray; f(ocean~=1) = nan;
h=imagesc(f(:,end:-1:1)');
set(h, 'AlphaData', h.CData>=0)
n = 128; % dark grey => decrease n; light grey => increase n; 
tmp = [n n n]/255; % grey color
mymap = [tmp;jet(255)];     % my colormap
colormap(mymap);
colorbar
hold on 

xticks([]); yticks([]);
plot(800:1000,ones(201,1)*900,'linewidth',3,'color','black')
plot(ones(31,1)*800,885:915,'linewidth',2,'color','black')
plot(ones(31,1)*1000,885:915,'linewidth',2,'color','black')
text(900,850,'200 km','HorizontalAlignment','center','interpreter','latex','FontSize',15)


hold off
title('Commercial citrus','interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;


set(1,'paperunits','centimeters');
set(1,'papersize',[36 12]);
set(1,'paperposition',[-2 0 40 12]);
% set(1,'papersize',[35 34]);
% set(1,'paperposition',[0 -5 35 38]);


    print(1,'-dpdf',append('Figures/',figname,'_CitrusMaps.pdf'));
    print(1,'-dtiff',append('Figures/',figname,'_CitrusMaps.tiff'));
    print(1,'-deps',append('Figures/',figname,'_CitrusMaps.eps'));
    savefig(append('Figures/',figname,'_CitrusMaps.fig'));
