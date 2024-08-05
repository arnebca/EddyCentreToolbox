function [output,xyuv] = eddyreconstruction(lonvec_all,latvec_all,uvec_all,vvec_all,timevec_all,...
    gridding,avg_int,depthvec,depth_int,ind,inddel,smooth,nr,subnr,lonlatbound,AVISO,...
    name,xaxis,save,projectdir,GN,r,shiptrack,x0,y0,ref,rx,NumOfIt,lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   eddyreconstruction - Eddy centre visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Arne Bendinger
%
%   Questions and bug reports to
%   'abendinger@geomar.de' or 'arne.bendinger@web.de'
%
%   Version 1.2 (29/06/2020), Author: Arne Bendinger
%   --> visualise eddy radius not using viscircles
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DESCRIPTION
%
%   This function illustrates and computes the eddy centre for a section of
%   interest or randomly distributed samples.
%   The function accepts both raw and already gridded/smoothed data.
%   If using non-gridded data, the function will project the values on a
%   regular grid making use of Gaussian weights.
%   If desired load gridded satellite altimetry data products in such way
%   as described in EddyCentreToolbox\Example\RunSADCP_example.m
%
%   The main computation is done by 'eddycentre_GN.m'. See the documention
%   for more information, i.e. the mathematical theory and the input
%   parameters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUT
%
%       lonvec_all  - Longitude from the entire cruise [1,N]
%
%       latvec_all  - Latitude from the entire cruise [1,N]
%
%       uvec_all    - Vector or array of zonal velocity in m/s [M,N].
%
%       vvec_all    - Vector or array of meridional velocity in m/s [M ,N]
%
%       timevec_all - Matlab time (mtime) from the entire cruise [1,N].
%                     Leave empty if not available.
%
%       gridding    - Choose "gridding" for along-track gridding or
%                     "nogridding" for data that is already gridded and/or
%                     smoothed.
%
%       avg_int     - Along-track average interval in km
%
%       depthvec    - Depth vector in m
%
%       depth_int   - Depth index for depth interval/level of interst.
%                     Velocities are averaged for given depth interval.
%
%       ind         - Specify indices for section of interest [1,N]. The
%                     section of interest is then extracted from
%                     lonvec_all, latvec_all, uvec_all, vvec_all, etc.
%
%       inddel      - Indices along section of interest to be deleted (e.g.
%                     when ship was on station or off-track). In a future
%                     version, CTD stations might be automatically removed.
%
%       smooth      - Influnce and cut-off radius [xrad xcut yrad ycut]
%
%       nr          - Number of section (for plotting reasons)
%
%       subnr       - Subnumber of section (for plotting reasons and in
%                     case more than one eddy is along section of interest)
%
%       lonlatbound - Western, eastern longitude and southern, northern
%                     latitude boundary for m_map [lon1 lon2 lat1 lat2]
%
%       AVISO       - AVISO SLA covering the entire cruise. Leave empty if
%                     not available.
%
%       name        - Name of research cruise or campaign
%
%       xaxis       - Choose 'lonaxis' or 'lataxis' for the lateral
%                     section's x-axis label
%
%       save        - Choose 'yes' or 'no' for saving option
%
%       projectdir  - Directory path for saved figures
%
%       GN          - Choose 'yes' for optimisation or 'no' for no
%                     optimisation
%
%       r           - Addtionally crop section of interest by specifying
%                     the amount of indices to be neglected from start and
%                     end of section [r1 r2]. Leave empty for no extra
%                     cropping (default). This option might be removed in a
%                     future version.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SEE 'eddycentre_GN.m' DOCUMENTATION FOR ADDITIONAL INPUT
%
%       shiptrack,x0,y0,ref,rx,NumOfIt,lambda
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   OUTPUT
%
%       output - Optimisation output from 'eddycentre_GN.m'
%       xyuv   - Gridded longitude,latitude and zonal,meridional velocity
%                for the section of interest
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare input

if isempty(r)
    r = [1 0];
end

lon_old = lonvec_all; lat_old = latvec_all;

lonvec_all(inddel) = [];
latvec_all(inddel) = [];
uvec_all(:,inddel) = [];
vvec_all(:,inddel) = [];
if ~isempty(AVISO) && ~isempty(timevec_all)
    timevec_all(inddel) = [];
end

if ~isempty(AVISO) && ~isempty(timevec_all)
    t = cellstr(datestr(timevec_all));
end

if isempty(ind)
    ind = 1:length(lonvec_all);
end

lon = lonvec_all(ind);
lat = latvec_all(ind);
u = uvec_all(:,ind);
v = vvec_all(:,ind);
if ~isempty(AVISO) && ~isempty(timevec_all)
    t = cellstr(datestr(timevec_all(ind)));
end

%% Plot section of interest

figure(1), clf, hold on
set(gcf,'color','w','units','normalized','position',[0.05 0.1 0.9 0.75])

subplot('position',[0.05 0.1 0.3 0.8]),hold on
m_proj('mercator','longitudes',[lonlatbound(1) lonlatbound(2)],'latitudes',[lonlatbound(3) lonlatbound(4)]);
m_gshhs_l('patch',[.7 .7 .7]);
m_grid('tickdir','in','linewidth',1.2,'fontsize',20);

m_plot(lon_old,lat_old,'o','color','k','markersize',1)
m_plot(lonvec_all,latvec_all,'o','color','k','markersize',1)
m_plot(lon,lat,'o','color','r','markersize',3)

text(0.05,0.95,strcat(name,{' '},'section',num2str(nr)),'units','normalized','fontweight','bold','fontsize',20)

if ~(size(v,1) == 1)
    subplot('position',[0.45 0.1 0.5 0.35]),hold on
    contourf(v,16,'linestyle','none'),caxis([-0.5 0.5]),colormap(cmocean('balance')),axis ij
    set(gca,'box','on','layer','top','linewidth',1.2,'fontsize',16),xlabel('index'),ylabel('depth index'),cl = colorbar; title(cl,'v [m/s]')
    text(0.05,0.15,'meridional','fontweight','bold','fontsize',20,'units','normalized')
    subplot('position',[0.45 0.55 0.5 0.35]),hold on
    contourf(u,16,'linestyle','none'),caxis([-0.5 0.5]),colormap(cmocean('balance')),axis ij
    set(gca,'box','on','layer','top','linewidth',1.2,'fontsize',16),xlabel('index'),ylabel('depth index'),cl = colorbar; title(cl,'u [m/s]')
    text(0.05,0.15,'zonal','fontweight','bold','fontsize',20,'units','normalized')
    if ~isempty(AVISO) && ~isempty(timevec_all)
        title(strcat(t(1),{' '},'-',{' '},t(end)))
    end
end
%% Grid and average velocity data

if strcmp(gridding,'gridding')
    
    for i = 1:length(lon)-1
        dist(1) = 0;
        dist(i+1) = dist(i) + sw_dist([lat(i) lat(i+1)],[lon(i) lon(i+1)],'km');
    end
    
    distnew = 0:avg_int:round(dist(end));
    
    for i = 1:length(distnew)-1
        ind_toaverage{i} = find(dist >= distnew(i) & dist < distnew(i+1));
        unew(:,i) = nanmean(u(:,ind_toaverage{i}),2);
        vnew(:,i) = nanmean(v(:,ind_toaverage{i}),2);
        lonquiver(i) = nanmean(lon(ind_toaverage{i}),2);
        latquiver(i) = nanmean(lat(ind_toaverage{i}),2);
    end
    
    [X,Y] = meshgrid(distnew(1:end-1),depthvec);
    
    Uinterp = obana3(unew,X,Y,X,Y,smooth(1),smooth(2),smooth(3),smooth(4));
    Vinterp = obana3(vnew,X,Y,X,Y,smooth(1),smooth(2),smooth(3),smooth(4));
    
elseif strcmp(gridding,'nogridding')
    
    lonquiver = lon;
    latquiver = lat;
    Uinterp = u;
    Vinterp = v;
    
end

%% Plot altimeter and SADCP data

figure(2), clf, hold on
set(gcf,'color','w','units','normalized','position',[0.05 0.1 0.9 0.75])

subplot('position',[0.05 0.1 0.3 0.8]),hold on
m_proj('mercator','longitudes',[lonlatbound(1) lonlatbound(2)],'latitudes',[lonlatbound(3) lonlatbound(4)]);
m_grid('tickdir','in','linewidth',1.2,'fontsize',20,'layer','top');
m_plot(lonvec_all,latvec_all,'.','color',[0.8 0.8 0.8],'markersize',3)


% traw = cellstr(datestr(timevec_all));

% AVISO
if ~isempty(AVISO) && ~isempty(timevec_all)
    
    % corresponding ADT data to start of the section
    for i = 1:length(AVISO.date)
        if strcmp(t{1}(1:11),AVISO.date{i}(1:11))
            %         if strcmp(traw{1}(1:11),AVISO.date{i}(1:11))
            AVISOeq = i;
            break
        end
    end
    
    % Plot sea level anomaly (SLA). SLA range is set to [-0.2 0.2]. Change
    % if desired. Instead of SLA, absolute dynamic topography can also be
    % plotted. To do this, comment rows 237:239 and uncomment 242:244.
    
    % AVISO SLA
    m_contourf(AVISO.LONAVISO',AVISO.LATAVISO',AVISO.SLA(:,:,AVISOeq)*100,16,'linestyle','none','linewidth',1.5);
    cl = colorbar('fontsize',14); ylabel(cl,'SLA [cm]')
    caxis([-20 20]) % SLA
    
    % AVISO ADT
    %     m_contourf(AVISO.LONAVISO',AVISO.LATAVISO',AVISO.ADT(:,:,AVISOeq)*100,16,'linestyle','none','linewidth',1.5);
    %     cl = colorbar('fontsize',14); ylabel(cl,'ADT [cm]')
    %     caxis([]) % manually specify ADT range
    
    % AVISO velocities
    HP = m_vec(1.5,AVISO.LONAVISO',AVISO.LATAVISO',AVISO.UGOSA(:,:,AVISOeq),AVISO.VGOSA(:,:,AVISOeq),'shaftwidth',0.5,'headlength',4,'edgeclip','on');
    set(HP,'facecolor',[0.3 0.3 0.3],'edgecolor',[0.3 0.3 0.3])
    %     set(HP,'facecolor','k','edgecolor','k')
    
end

% coastline
if strcmp(save,'yes')
    m_gshhs_l('patch',[.5 .5 .5]);
else
    m_gshhs_l('patch',[.5 .5 .5]);
end

if length(depth_int) == 1
    
    if strcmp(gridding,'gridding')
        
        if avg_int <= 2
            m_vec(1.5,lonquiver(1:3:end),latquiver(1:3:end),Uinterp(depth_int,1:3:end),Vinterp(depth_int,1:3:end));
        else
            m_vec(1.5,lonquiver,latquiver,Uinterp(depth_int,:),Vinterp(depth_int,:));
        end
    elseif strcmp(gridding,'nogridding')
        m_vec(1.5,lonquiver,latquiver,Uinterp(depth_int,:),Vinterp(depth_int,:));
    else
        error('Unknown option for gridding. Choose "gridding" or "nogridding".')
    end
    
    if isempty(subnr)
        title(strcat(name,{' '},'section',num2str(nr),{' '},num2str(depthvec(depth_int)),'m'),'fontsize',20)
    else
        title(strcat(name,{' '},'section',num2str(nr),'.',num2str(subnr),{' '},num2str(depthvec(depth_int)),'m'),'fontsize',20)
    end
    
else
    
    if strcmp(gridding,'gridding')
        
        if avg_int <= 2
            m_vec(1.5,lonquiver(1:3:end),latquiver(1:3:end),nanmean(Uinterp(depth_int,1:3:end)),nanmean(Vinterp(depth_int,1:3:end)));
        else
            m_vec(1.5,lonquiver,latquiver,nanmean(Uinterp(depth_int,:)),nanmean(Vinterp(depth_int,:)));
        end
        
    elseif strcmp(gridding,'nogridding')
        m_vec(1.5,lonquiver,latquiver,nanmean(Uinterp(depth_int,:)),nanmean(Vinterp(depth_int,:)));
    else
        error('Unknown option for gridding. Choose "gridding" or "nogridding".')
    end
    
    if isempty(subnr)
        title(strcat(name,{' '},'section',num2str(nr),{' '},'mean',{' '},num2str(depthvec(depth_int(1))),'-',num2str(depthvec(depth_int(end))),'m'),'fontsize',20)
    else
        title(strcat(name,{' '},'section',num2str(nr),'.',num2str(subnr),{' '},'mean',{' '},num2str(depthvec(depth_int(1))),'-',num2str(depthvec(depth_int(end))),'m'),'fontsize',20)
    end
end

% legend
m_vec(1.5,lonlatbound(2)-1.5,lonlatbound(4)-0.35,0.5,0)
m_text(lonlatbound(2)-1.5,lonlatbound(4)-0.25,'0.5 m/s','fontsize',14,'fontweight','bold')

if strcmp(xaxis,'lonaxis')
    [XAXIS,DEPTHVEC] = meshgrid(lonquiver,depthvec);
elseif strcmp(xaxis,'lataxis')
    [XAXIS,DEPTHVEC] = meshgrid(latquiver,depthvec);
end

if ~(size(v,1) == 1)
    subplot('position',[0.45 0.1 0.5 0.35]),hold on
    contourf(XAXIS,DEPTHVEC,Vinterp,-1:0.02:1,'linestyle','none'),caxis([-0.5 0.5]),colormap(cmocean('balance')),axis ij
    [C,h] = contour(XAXIS,DEPTHVEC,Vinterp,-1:0.1:-0.1,'color','k','linestyle','--');
    clabel(C,h,'color','k','fontsize',8,'fontweight','bold')
    [C,h] = contour(XAXIS,DEPTHVEC,Vinterp,0.1:0.1:1,'color','k','linestyle','-');
    clabel(C,h,'color','k','fontsize',8,'fontweight','bold')
    [C,h] = contour(XAXIS,DEPTHVEC,Vinterp,[0 0],'color','k','linewidth',1.5);
    clabel(C,h,'color','k','fontsize',8,'fontweight','bold')
    set(gca,'box','on','layer','top','linewidth',1.2,'fontsize',16),ylim([0 max(depthvec)])
    
    plot(XAXIS(1,:),ones(1,length(XAXIS(1,:)))*depthvec(depth_int(1)),'--','color','k','linewidth',1.5)
    if length(depth_int) > 1
        plot(XAXIS(1,:),ones(1,length(XAXIS(1,:)))*depthvec(depth_int(end)),'--','color','k','linewidth',1.5)
    end
    
    if strcmp(xaxis,'lonaxis')
        xlabel('longitude')
    elseif strcmp(xaxis,'lataxis')
        xlabel('latitude')
    end
    
    ylabel('depth [m]'),cl = colorbar; ylabel(cl,'v [m/s]')
    text(0.05,0.15,'meridional','fontweight','bold','fontsize',20,'units','normalized')
    
    subplot('position',[0.45 0.55 0.5 0.35]),hold on
    contourf(XAXIS,DEPTHVEC,Uinterp,-1:0.02:1,'linestyle','none'),caxis([-0.5 0.5]),colormap(cmocean('balance')),axis ij
    [C,h] = contour(XAXIS,DEPTHVEC,Uinterp,-1:0.1:-0.1,'color','k','linestyle','--');
    clabel(C,h,'color','k','fontsize',8,'fontweight','bold')
    [C,h] = contour(XAXIS,DEPTHVEC,Uinterp,0.1:0.1:1,'color','k','linestyle','-');
    clabel(C,h,'color','k','fontsize',8,'fontweight','bold')
    [C,h] = contour(XAXIS,DEPTHVEC,Uinterp,[0 0],'color','k','linewidth',1.5);
    clabel(C,h,'color','k','fontsize',8,'fontweight','bold')
    set(gca,'box','on','layer','top','linewidth',1.2,'fontsize',16),ylim([0 max(depthvec)])
    
    plot(XAXIS(1,:),ones(1,length(XAXIS(1,:)))*depthvec(depth_int(1)),'--','color','k','linewidth',1.5)
    if length(depth_int) > 1
        plot(XAXIS(1,:),ones(1,length(XAXIS(1,:)))*depthvec(depth_int(end)),'--','color','k','linewidth',1.5)
    end
    
    if strcmp(xaxis,'lonaxis')
        xlabel('longitude')
    elseif strcmp(xaxis,'lataxis')
        xlabel('latitude')
    end
    
    ylabel('depth [m]'),cl = colorbar; ylabel(cl,'u [m/s]')
    text(0.05,0.15,'zonal','fontweight','bold','fontsize',20,'units','normalized')
    
    if ~isempty(AVISO) && ~isempty(timevec_all)
        title(strcat(t(1),{' '},'-',{' '},t(end)))
    end
end
%% save

if strcmp(GN,'no') && strcmp(save,'yes')
    
    if strcmp(save,'yes')
        
        if length(depth_int) == 1
            if isempty(subnr)
                export_fig(strcat(projectdir,strcat(name,'_','section',num2str(nr),'_smoothed_',num2str(depthvec(depth_int)),'m')),'-png','-r50');
            else
                export_fig(strcat(projectdir,strcat(name,'_','section',num2str(nr),'.',num2str(subnr),'_smoothed_',num2str(depthvec(depth_int)),'m')),'-png','-r50');
            end
        else
            if isempty(subnr)
                export_fig(strcat(projectdir,strcat(name,'_','section',num2str(nr),'_smoothed_','mean',num2str(depthvec(depth_int(1))),'_',num2str(depthvec(depth_int(end))),'m')),'-png','-r50');
            else
                export_fig(strcat(projectdir,strcat(name,'_','section',num2str(nr),'.',num2str(subnr),'_smoothed_','mean',num2str(depthvec(depth_int(1))),'_',num2str(depthvec(depth_int(end))),'m')),'-png','-r50');
            end
        end
        
    end
    
end

%% Gauss-Newton optimisation

if strcmp(GN,'yes')
    
    range_int = r(1):length(lonquiver)-r(2);
    
%     lon = lonquiver(range_int);
%     lat = latquiver(range_int);
    
    
    if length(depth_int) == 1
%         [ output ] = eddycentre_GN(lon,lat,Uinterp(depth_int,range_int),Vinterp(depth_int,range_int),shiptrack,x0,y0,ref,rx,NumOfIt,lambda);
        [ output ] = eddycentre_GN(lonquiver,latquiver,Uinterp(depth_int,:),Vinterp(depth_int,:),shiptrack,x0,y0,ref,rx,NumOfIt,lambda,range_int);
    else
%         [ output ] = eddycentre_GN(lon,lat,nanmean(Uinterp(depth_int,range_int)),nanmean(Vinterp(depth_int,range_int)),shiptrack,x0,y0,ref,rx,NumOfIt,lambda);
        [ output ] = eddycentre_GN(lonquiver,latquiver,nanmean(Uinterp(depth_int,:)),nanmean(Vinterp(depth_int,:)),shiptrack,x0,y0,ref,rx,NumOfIt,lambda,range_int);
    end
    
    %% Save Gauss-Newton plot
    
    if strcmp(save,'yes')
        
        if length(depth_int) == 1
            if isempty(subnr)
                export_fig(strcat(projectdir,strcat(name,'_','GNsection',num2str(nr),'_smoothed_',num2str(depthvec(depth_int)),'m','_','GaussNewton')),'-png','-r50');
            else
                export_fig(strcat(projectdir,strcat(name,'_','GNsection',num2str(nr),'.',num2str(subnr),'_smoothed_',num2str(depthvec(depth_int)),'m','_','GaussNewton')),'-png','-r50');
            end
        else
            if isempty(subnr)
                export_fig(strcat(projectdir,strcat(name,'_','GNsection',num2str(nr),'_smoothed_','mean',num2str(depthvec(depth_int(1))),'_',num2str(depthvec(depth_int(end))),'m','_','GaussNewton')),'-png','-r50');
            else
                export_fig(strcat(projectdir,strcat(name,'_','GNsection',num2str(nr),'.',num2str(subnr),'_smoothed_','mean',num2str(depthvec(depth_int(1))),'_',num2str(depthvec(depth_int(end))),'m','_','GaussNewton')),'-png','-r50');
            end
        end
        
    end
    
    %% Include Gauss-Newton results in map projection
    
    close figure 1
    
    figure(1), clf, hold on
    set(gcf,'color','w','units','normalized','position',[0.05 0.1 0.9 0.75])
    
    subplot('position',[0.05 0.1 0.3 0.8]),hold on
    m_proj('mercator','longitudes',[lonlatbound(1) lonlatbound(2)],'latitudes',[lonlatbound(3) lonlatbound(4)]);
    m_grid('tickdir','in','linewidth',1.2,'fontsize',20,'layer','top');
    m_plot(lonvec_all,latvec_all,'.','color',[0.8 0.8 0.8],'markersize',3)
    
    % AVISO
    if ~isempty(AVISO) && ~isempty(timevec_all)
        
        % corresponding ADT data to start of the section
        for i = 1:length(AVISO.date)
            if strcmp(t{1}(1:11),AVISO.date{i}(1:11))
                AVISOeq = i;
                break
            end
        end
        
        % Plot sea level anomaly (SLA). SLA range is set to [-0.2 0.2]. Change
        % if desired. Instead of SLA, absolute dynamic topography can also be
        % plotted. To do this, comment rows 237:239 and uncomment 242:244.
        
        % AVISO SLA
        m_contourf(AVISO.LONAVISO',AVISO.LATAVISO',AVISO.SLA(:,:,AVISOeq)*100,16,'linestyle','none','linewidth',1.5);
        cl = colorbar('fontsize',14); ylabel(cl,'SLA [cm]')
        caxis([-20 20]) % SLA
        
        % AVISO ADT
        %     m_contourf(AVISO.LONAVISO',AVISO.LATAVISO',AVISO.ADT(:,:,AVISOeq)*100,16,'linestyle','none','linewidth',1.5);
        %     cl = colorbar('fontsize',14); ylabel(cl,'ADT [cm]')
        %     caxis([]) % manually specify ADT range
        
        % AVISO velocities
        HP = m_vec(1.5,AVISO.LONAVISO',AVISO.LATAVISO',AVISO.UGOSA(:,:,AVISOeq),AVISO.VGOSA(:,:,AVISOeq),'shaftwidth',0.5,'headlength',4,'edgeclip','on');
        set(HP,'facecolor',[0.3 0.3 0.3],'edgecolor',[0.3 0.3 0.3])
        %     set(HP,'facecolor','k','edgecolor','k')
        
        
    end
    
    % coastline
    if strcmp(save,'yes')
        m_gshhs_l('patch',[.5 .5 .5]);
    else
        m_gshhs_l('patch',[.5 .5 .5]);
    end
    
    if length(depth_int) == 1
        
        if strcmp(gridding,'gridding')
            
            if avg_int <= 2
                m_vec(1.5,lonquiver(1:3:end),latquiver(1:3:end),Uinterp(depth_int,1:3:end),Vinterp(depth_int,1:3:end));
            else
                m_vec(1.5,lonquiver,latquiver,Uinterp(depth_int,:),Vinterp(depth_int,:));
            end
        elseif strcmp(gridding,'nogridding')
            m_vec(1.5,lonquiver,latquiver,Uinterp(depth_int,:),Vinterp(depth_int,:));
        else
            error('Unknown option for gridding. Choose "gridding" or "nogridding".')
        end
        
        if isempty(subnr)
            title(strcat(name,{' '},'section',num2str(nr),{' '},num2str(depthvec(depth_int)),'m'),'fontsize',20)
        else
            title(strcat(name,{' '},'section',num2str(nr),'.',num2str(subnr),{' '},num2str(depthvec(depth_int)),'m'),'fontsize',20)
        end
        
    else
        
        if strcmp(gridding,'gridding')
            
            if avg_int <= 2
                m_vec(1.5,lonquiver(1:3:end),latquiver(1:3:end),nanmean(Uinterp(depth_int,1:3:end)),nanmean(Vinterp(depth_int,1:3:end)));
            else
                m_vec(1.5,lonquiver,latquiver,nanmean(Uinterp(depth_int,:)),nanmean(Vinterp(depth_int,:)));
            end
            
        elseif strcmp(gridding,'nogridding')
            m_vec(1.5,lonquiver,latquiver,nanmean(Uinterp(depth_int,:)),nanmean(Vinterp(depth_int,:)));
        else
            error('Unknown option for gridding. Choose "gridding" or "nogridding".')
        end
        
        if isempty(subnr)
            title(strcat(name,{' '},'section',num2str(nr),{' '},'mean',{' '},num2str(depthvec(depth_int(1))),'-',num2str(depthvec(depth_int(end))),'m'),'fontsize',20)
        else
            title(strcat(name,{' '},'section',num2str(nr),'.',num2str(subnr),{' '},'mean',{' '},num2str(depthvec(depth_int(1))),'-',num2str(depthvec(depth_int(end))),'m'),'fontsize',20)
        end
    end
    
    % plot eddy boundary
%     if strcmp(shiptrack,'crossed') || strcmp(shiptrack,'half-crossed')
    if ~strcmp(shiptrack,'within')
        [newx, newy] = m_ll2xy([output.eddycentre_lonlat(1) output.eddycentre_lonlat(1) + output.radvel.radx],...
            [output.eddycentre_lonlat(2) output.eddycentre_lonlat(2) + output.radvel.rady]);
        rcircle = sqrt((newx(1) - newx(2))^2 + (newy(1) - newy(2))^2);
        rectangle('position',[newx(1)-rcircle newy(1)-rcircle rcircle*2 rcircle*2],'curvature',[1,1],'linewidth',3.5,'edgecolor','k')
        rectangle('position',[newx(1)-rcircle newy(1)-rcircle rcircle*2 rcircle*2],'curvature',[1,1],'linewidth',2.5,'edgecolor','r')
    end

    
    text(0.05,0.125,strcat('lon=',num2str(round(output.eddycentre_lonlat(1),2)),',',{' '},'lat=', ...
        num2str(round(output.eddycentre_lonlat(2),2))),'units','normalized','fontsize',18,'fontweight','bold')
    
    if strcmp(shiptrack,'crossed') || strcmp(shiptrack,'random')
        text(0.05,0.05,strcat('mean radius =',{' '},num2str(round(mean(output.radius),1)),'km'),...
            'units','normalized','fontweight','bold','fontsize',18)
    elseif strcmp(shiptrack,'half-crossed')
        text(0.05,0.05,strcat('radius =',{' '},num2str(round(mean(output.radius),1)),'km'),...
            'units','normalized','fontweight','bold','fontsize',18)
    end
    
%     m_plot(output.eddycentre_lonlat(1),output.eddycentre_lonlat(2),'p','markeredgecolor','k','markerfacecolor','r','markersize',8)
    
    % legend
    m_vec(1.5,lonlatbound(2)-1.5,lonlatbound(4)-0.35,0.5,0)
    m_text(lonlatbound(2)-1.5,lonlatbound(4)-0.25,'0.5 m/s','fontsize',14,'fontweight','bold')
    
    if strcmp(xaxis,'lonaxis')
        [XAXIS,DEPTHVEC] = meshgrid(lonquiver,depthvec);
    elseif strcmp(xaxis,'lataxis')
        [XAXIS,DEPTHVEC] = meshgrid(latquiver,depthvec);
    end
    
    if ~(size(v,1) == 1) && ~strcmp(shiptrack,'random')
        
        subplot('position',[0.45 0.1 0.5 0.35]),hold on
        contourf(XAXIS,DEPTHVEC,Vinterp,-1:0.02:1,'linestyle','none'),caxis([-0.5 0.5]),colormap(cmocean('balance')),axis ij
        [C,h] = contour(XAXIS,DEPTHVEC,Vinterp,-1:0.1:-0.1,'color','k','linestyle','--');
        clabel(C,h,'color','k','fontsize',8,'fontweight','bold')
        [C,h] = contour(XAXIS,DEPTHVEC,Vinterp,0.1:0.1:1,'color','k','linestyle','-');
        clabel(C,h,'color','k','fontsize',8,'fontweight','bold')
        [C,h] = contour(XAXIS,DEPTHVEC,Vinterp,[0 0],'color','k','linewidth',1.5);
        clabel(C,h,'color','k','fontsize',8,'fontweight','bold')        
        set(gca,'box','on','layer','top','linewidth',1.2,'fontsize',16),ylim([0 max(depthvec)])
        
        plot(XAXIS(1,:),ones(1,length(XAXIS(1,:)))*depthvec(depth_int(1)),'--','color','k','linewidth',1.5)
        if length(depth_int) > 1
            plot(XAXIS(1,:),ones(1,length(XAXIS(1,:)))*depthvec(depth_int(end)),'--','color','k','linewidth',1.5)
        end
        
        if strcmp(xaxis,'lonaxis')
            xlabel('longitude')
        elseif strcmp(xaxis,'lataxis')
            xlabel('latitude')
        end
        
        ylabel('depth [m]'),cl = colorbar; ylabel(cl,'v [m/s]')
        text(0.05,0.15,'meridional','fontweight','bold','fontsize',20,'units','normalized')
        
        subplot('position',[0.45 0.55 0.5 0.35]),hold on
        contourf(XAXIS,DEPTHVEC,Uinterp,-1:0.02:1,'linestyle','none'),caxis([-0.5 0.5]),colormap(cmocean('balance')),axis ij
        [C,h] = contour(XAXIS,DEPTHVEC,Uinterp,-1:0.1:-0.1,'color','k','linestyle','--');
        clabel(C,h,'color','k','fontsize',8,'fontweight','bold')
        [C,h] = contour(XAXIS,DEPTHVEC,Uinterp,0.1:0.1:1,'color','k','linestyle','-');
        clabel(C,h,'color','k','fontsize',8,'fontweight','bold')
        [C,h] = contour(XAXIS,DEPTHVEC,Uinterp,[0 0],'color','k','linewidth',1.5);
        clabel(C,h,'color','k','fontsize',8,'fontweight','bold')        
        set(gca,'box','on','layer','top','linewidth',1.2,'fontsize',16),ylim([0 max(depthvec)])
        
        plot(XAXIS(1,:),ones(1,length(XAXIS(1,:)))*depthvec(depth_int(1)),'--','color','k','linewidth',1.5)
        if length(depth_int) > 1
            plot(XAXIS(1,:),ones(1,length(XAXIS(1,:)))*depthvec(depth_int(end)),'--','color','k','linewidth',1.5)
        end
        
        if strcmp(xaxis,'lonaxis')
            xlabel('longitude')
        elseif strcmp(xaxis,'lataxis')
            xlabel('latitude')
        end
        
        ylabel('depth [m]'),cl = colorbar; ylabel(cl,'u [m/s]')
        text(0.05,0.15,'zonal','fontweight','bold','fontsize',20,'units','normalized')
        if ~isempty(AVISO) && ~isempty(timevec_all)
            title(strcat(t(1),{' '},'-',{' '},t(end)))
        end
    end
    
    %% save
    
    if strcmp(save,'yes')
        
        if length(depth_int) == 1
            if isempty(subnr)
                export_fig(strcat(projectdir,strcat(name,'_','section',num2str(nr),'_smoothed_',num2str(depthvec(depth_int)),'m','_','GaussNewton')),'-png','-r50');
            else
                export_fig(strcat(projectdir,strcat(name,'_','section',num2str(nr),'.',num2str(subnr),'_smoothed_',num2str(depthvec(depth_int)),'m','_','GaussNewton')),'-png','-r50');
            end
        else
            if isempty(subnr)
                export_fig(strcat(projectdir,strcat(name,'_','section',num2str(nr),'_smoothed_','mean',num2str(depthvec(depth_int(1))),'_',num2str(depthvec(depth_int(end))),'m','_','GaussNewton')),'-png','-r50');
            else
                export_fig(strcat(projectdir,strcat(name,'_','section',num2str(nr),'.',num2str(subnr),'_smoothed_','mean',num2str(depthvec(depth_int(1))),'_',num2str(depthvec(depth_int(end))),'m','_','GaussNewton')),'-png','-r50');
            end
        end
        
    end
    
elseif strcmp(GN,'no')
    output = 'no optimisation performed';
    disp('No optimisation performed')
else
    error('Unknown option for "GN". Choose "yes" for optimisation or "no" for no optimisiation.')
end

xyuv.lonquiver = lonquiver; xyuv.latquiver = latquiver;
xyuv.Uinterp = Uinterp; xyuv.Vinterp = Vinterp;

end
