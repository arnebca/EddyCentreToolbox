%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ShiptrackDemo.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute grid

clearvars

x = linspace(-50,50,40); y = linspace(-50,50,40);
[X,Y] = meshgrid(x,y);
delta = 40; alpha = 7;

vel = 0.5*10^-1.*sqrt(X.^2+Y.^2).*exp(-(sqrt(X.^2+Y.^2)./delta).^alpha);

xcom = (-(Y./X)./sqrt(1+(Y./X).^2)).*sign(X).*vel;
ycom = (1./sqrt(1+(Y./X).^2)).*sign(X).*vel;

R = 6371*10^3; % radius in km
reflon = -50.5; reflat = 14.5;
newlon = reflon + (X*10^3./R) * (180/pi)/(cos(reflat*pi/180));
newlat = reflat + (Y*10^3./R) * (180/pi);

%% Crossed

figure(1),clf
set(gcf,'color','w','units','normalized','position',[0.05 0.1 0.8 0.65])
subplot(121), hold on
contourf(newlon,newlat,vel,32,'linestyle','none'),colormap(flipud(cmocean('tempo'))),cl = colorbar; title(cl,'[m/s]')
quiver(newlon,newlat,xcom,ycom,'color','k')
plot(reflon,reflat,'p','markerfacecolor','r','markeredgecolor','k','markersize',16)
set(gca,'linewidth',1.2,'fontsize',16),xlabel('longitude'),ylabel('latitude')
text(0.05,0.95,'crossed','units','normalized','fontsize',24,'fontweight','bold','color','w')
axis equal

diag_nr = 5;

sectionlon = diag(newlon,diag_nr); sectionlat = diag(newlat,diag_nr);

plot(sectionlon,sectionlat,'color','r','linewidth',2)

sectionu = diag(xcom,diag_nr); sectionv = diag(ycom,diag_nr);

[xd_full,yd_full] = lonlat2xy(sectionlon,...
        sectionlat,-50.5,14.5);

E1radmean = 29.72;

subplot(122),hold on
% mean eddy radius
l3 = viscircles([0 0],E1radmean,'color',[0.7 0.7 0.7],'linewidth',2);

% velocity vectors
l1 = quiver(xd_full'*10^-3,yd_full'*10^-3,sectionu,sectionv,'color','black','linewidth',1.5);

% eddy centre
l2 = plot(0,0,'p','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',12);
grid on,xlabel('eastward distance [km]'),ylabel('northward distance [km]'), axis equal
axrange = 80;
xlim([-axrange axrange]),ylim([-axrange axrange])
set(gca,'box','on','linewidth',1.3,'tickdir','out','fontsize',14,'xtick',-axrange:20:axrange,'ytick',-axrange:20:axrange)

% print converted eddy centre to longitude, latitude; mean radius
text(0.05,0.15,strcat('lon=',num2str(-50.5),',',{' '},'lat=',num2str(14.5)),'units','normalized','fontsize',16,'fontweight','bold')
text(0.05,0.05,strcat('mean radius=',num2str(round(E1radmean,1)),'km'),'units','normalized','fontsize',16,'fontweight','bold')

legend([l1,l2,l3],'velocity vector','eddy centre','eddy radius','location','northeast')

%% Half-crossed

figure(2),clf
set(gcf,'color','w','units','normalized','position',[0.05 0.1 0.8 0.65])
subplot(121), hold on
contourf(newlon,newlat,vel,32,'linestyle','none'),colormap(flipud(cmocean('tempo'))),cl = colorbar; title(cl,'[m/s]')
quiver(newlon,newlat,xcom,ycom,'color','k')
plot(reflon,reflat,'p','markerfacecolor','r','markeredgecolor','k','markersize',16)
set(gca,'linewidth',1.2,'fontsize',16),xlabel('longitude'),ylabel('latitude')
text(0.05,0.95,'half-crossed','units','normalized','fontsize',24,'fontweight','bold','color','w')
axis equal

diag_nr = -5;

sectionlon = diag(newlon,diag_nr); sectionlat = diag(newlat,diag_nr);

plot(sectionlon(17:end),sectionlat(17:end),'color','r','linewidth',2)

sectionu = diag(xcom,diag_nr); sectionv = diag(ycom,diag_nr);

[xd_full,yd_full] = lonlat2xy(sectionlon,...
        sectionlat,-50.5,14.5);

E1radmean = 30.6177;

subplot(122),hold on
% mean eddy radius
l3 = viscircles([0 0],E1radmean,'color',[0.7 0.7 0.7],'linewidth',2);

% velocity vectors
l1 = quiver(xd_full(17:end)'*10^-3,yd_full(17:end)'*10^-3,sectionu(17:end),sectionv(17:end),'color','black','linewidth',1.5);

% eddy centre
l2 = plot(0,0,'p','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',12);
grid on,xlabel('eastward distance [km]'),ylabel('northward distance [km]'), axis equal
axrange = 80;
xlim([-axrange axrange]),ylim([-axrange axrange])
set(gca,'box','on','linewidth',1.3,'tickdir','out','fontsize',14,'xtick',-axrange:20:axrange,'ytick',-axrange:20:axrange)

% print converted eddy centre to longitude, latitude; mean radius
text(0.05,0.15,strcat('lon=',num2str(-50.5),',',{' '},'lat=',num2str(14.5)),'units','normalized','fontsize',16,'fontweight','bold')
text(0.05,0.05,strcat('mean radius=',num2str(round(E1radmean,1)),'km'),'units','normalized','fontsize',16,'fontweight','bold')

legend([l1,l2,l3],'velocity vector','eddy centre','eddy radius','location','northeast')

%% Within

figure(3),clf
set(gcf,'color','w','units','normalized','position',[0.05 0.1 0.8 0.65])
subplot(121), hold on
contourf(newlon,newlat,vel,32,'linestyle','none'),colormap(flipud(cmocean('tempo'))),cl = colorbar; title(cl,'[m/s]')
quiver(newlon,newlat,xcom,ycom,'color','k')
plot(reflon,reflat,'p','markerfacecolor','r','markeredgecolor','k','markersize',16)
set(gca,'linewidth',1.2,'fontsize',16),xlabel('longitude'),ylabel('latitude')
text(0.05,0.95,'within','units','normalized','fontsize',24,'fontweight','bold','color','w')
axis equal

diag_nr = 2;

sectionlon = diag(newlon,diag_nr); sectionlat = diag(newlat,diag_nr);

plot(sectionlon(13:end-13),sectionlat(13:end-13),'color','r','linewidth',2)

sectionu = diag(xcom,diag_nr); sectionv = diag(ycom,diag_nr);

[xd_full,yd_full] = lonlat2xy(sectionlon,...
        sectionlat,-50.5,14.5);

subplot(122),hold on

% velocity vectors
l1 = quiver(xd_full(13:end-13)'*10^-3,yd_full(13:end-13)'*10^-3,sectionu(13:end-13),sectionv(13:end-13),'color','black','linewidth',1.5);

% eddy centre
l2 = plot(0,0,'p','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',12);
grid on,xlabel('eastward distance [km]'),ylabel('northward distance [km]'), axis equal
axrange = 80;
xlim([-axrange axrange]),ylim([-axrange axrange])
set(gca,'box','on','linewidth',1.3,'tickdir','out','fontsize',14,'xtick',-axrange:20:axrange,'ytick',-axrange:20:axrange)

% print converted eddy centre to longitude, latitude; mean radius
text(0.05,0.15,strcat('lon=',num2str(-50.5),',',{' '},'lat=',num2str(14.5)),'units','normalized','fontsize',16,'fontweight','bold')

legend([l1,l2],'velocity vector','eddy centre','location','northeast')

%% Manual

figure(4),clf
set(gcf,'color','w','units','normalized','position',[0.05 0.1 0.8 0.65])
subplot(121), hold on
contourf(newlon,newlat,vel,32,'linestyle','none'),colormap(flipud(cmocean('tempo'))),cl = colorbar; title(cl,'[m/s]')
quiver(newlon,newlat,xcom,ycom,'color','k')
plot(reflon,reflat,'p','markerfacecolor','r','markeredgecolor','k','markersize',16)
set(gca,'linewidth',1.2,'fontsize',16),xlabel('longitude'),ylabel('latitude')
text(0.05,0.95,'manual','units','normalized','fontsize',24,'fontweight','bold','color','w')
axis equal

sectionlon1 = diag(newlon,5); sectionlat1 = diag(newlat,5);
sectionu1 = diag(xcom,5); sectionv1 = diag(ycom,5);

sectionlon2 = newlon(20,:); sectionlat2 = newlat(20,:);
sectionu2 = xcom(20,:); sectionv2 = ycom(20,:);

% lon
sectionlon(1:length(sectionlon1(1:20))) = sectionlon1(1:20);
sectionlon(length(sectionlon(1:20))+1:length(sectionlon(1:20))+length(sectionlon2(1:24))) = fliplr(sectionlon2(1:24));
% lat
sectionlat(1:length(sectionlat1(1:20))) = sectionlat1(1:20);
sectionlat(length(sectionlat(1:20))+1:length(sectionlat(1:20))+length(sectionlat2(1:24))) = fliplr(sectionlat2(1:24));
% u
sectionu(1:length(sectionu1(1:20))) = sectionu1(1:20);
sectionu(length(sectionu(1:20))+1:length(sectionu(1:20))+length(sectionu2(1:24))) = fliplr(sectionu2(1:24));
% v
sectionv(1:length(sectionv1(1:20))) = sectionv1(1:20);
sectionv(length(sectionv(1:20))+1:length(sectionv(1:20))+length(sectionv2(1:24))) = fliplr(sectionv2(1:24));

plot(sectionlon,sectionlat,'color','r','linewidth',2)

xd_full =   [  -37.2287  -34.6577  -32.0872  -29.5173  -26.9478  -24.3789  -21.8105  -19.2426  -16.6752  -14.1083  -11.5420   -8.9762   -6.4109   -3.8461   -1.2819    1.2818    3.8450    6.4076...
    8.9697   11.5313    8.9688    6.4063    3.8438    1.2813   -1.2813   -3.8438   -6.4063   -8.9688  -11.5313  -14.0938  -16.6563  -19.2188  -21.7813  -24.3438  -26.9063  -29.4688...
  -32.0314  -34.5939  -37.1564  -39.7189  -42.2814  -44.8439  -47.4064  -49.9689];

yd_full =   [  -49.9663  -47.4039  -44.8416  -42.2792  -39.7168  -37.1544  -34.5921  -32.0297  -29.4673  -26.9049  -24.3426  -21.7802  -19.2178  -16.6554  -14.0931  -11.5307   -8.9683   -6.4059...
   -3.8436   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812...
   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812   -1.2812];

subplot(122),hold on

% velocity vectors
l1 = quiver(xd_full',yd_full',sectionu,sectionv,'color','black','linewidth',1.5);

% eddy centre
l2 = plot(0,0,'p','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',12);
grid on,xlabel('eastward distance [km]'),ylabel('northward distance [km]'), axis equal
axrange = 80;
xlim([-axrange axrange]),ylim([-axrange axrange])
set(gca,'box','on','linewidth',1.3,'tickdir','out','fontsize',14,'xtick',-axrange:20:axrange,'ytick',-axrange:20:axrange)

% print converted eddy centre to longitude, latitude; mean radius
text(0.05,0.15,strcat('lon=',num2str(-50.5),',',{' '},'lat=',num2str(14.5)),'units','normalized','fontsize',16,'fontweight','bold')

legend([l1,l2],'velocity vector','eddy centre','location','northeast')