function [ output ] = eddycentre_GN(lon,lat,u,v,shiptrack,x0,y0,ref,rx,NumOfIt,lambda,range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   eddycentre_GN - Gauss-Newton optimisation eddy centre determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Arne Bendinger
%
%   Questions and bug reports to
%   'abendinger@geomar.de' or 'arne.bendinger@web.de'
%   
%
%   Version 1.8 (03/05/2021)
%   --> Inclusion of range variable
%   Version 1.7 (27/03/2021)
%   --> Include inverse power law for outer ring velocity structure
%   Version 1.6 (01/12/2020)
%   --> 
%   Version 1.5 (11/09/2020)
%   --> Graphical adjustments
%   Version 1.4 (20/07/2020)
%   --> remove input option 'rotation_angle' and 'origin'. No longer needed
%   Version 1.3 (29/06/2020)
%   --> corrected normalised position and velocity vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DESCRIPTION
%
%   This function determines the centre of an eddy based on a single
%   track of horizontal velocity or randomly distributed velocity sampled.
%   Assuming an axisymmetric and nontranslating vortex, the algorithm
%   relies on an optimisation of the azimuthal velocity component in the
%   framework of a cyclindrical coordinate system.
%   The optimisation, in turn, is performed by a nonlinear least-squares
%   Gauss-Newton method which minimises an objective function by computing
%   a descent or search direction in an
%   iterative procedure. A damping factor is implemented, known as the
%   Armijo rule or sufficient decrease condition, that prevents the line
%   search from taking overly long steps.
%   In addition, an optional penalty term is introduced (L2 regularisation
%   or ridge regression) to the residual sum of squares. It may avoid
%   overfitting issues. In practice, it turned out to be useful when the
%   ship transect is either really close to the center or exactly passing
%   through it.
%
%   The algorithm may strongly be dependent on the initialisation values.
%   Therefore, it is suggested to create a cluster of different start
%   values to ensure the convergence to a local minimum or global minimum
%   associated with the lowest residual sum of squares.
%
%   Since the velocity data may be noisy near the eddy rim due to
%   interactions with the surrounding flow, velocity samples from the
%   beginning and end of the track can be exluded in an iterative procedure
%   until the residual sum of squares converges in the ideal case or just
%   to take the mean of several estimates.
%
%   Input variables and independent parameters are all scaled in the range
%   [-1 1] with the original sign of the vector values being preserved.
%
%   Call ShiptrackDemo() for an illustration and demo version.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUT
%
%       lon,lat,u,v - longitude, latitude, zonal velocity [m/s], meridional
%                     velocity [m/s](the vectors must be all the same size
%                     and either all column or row vectors).
%
%       shiptrack   - Choose the type of transect through the eddy. For
%                     options 'crossed' and 'half-crossed' the tracks do
%                     not necessarily need to be perfectly fitted to the
%                     eddy velocity structure, i.e. velocities from the
%                     outer ring may be included.
%                     Choose type 'within' when the track lies within the
%                     inner core featuring no velocity maxima or when the
%                     data is noisy such that the two biggest velocity
%                     maxima are not associated with corresponding eddy
%                     boundary.
%
%                     1) 'crossed': ship track entirely crosses eddy
%                     featuring a velocity maximum at each eddy flank
%                     2) 'half-crossed': ship track does not entirely cross
%                     eddy featuring just one velocity maximum
%                     3) 'within': ship track is within inner core, not
%                     featuring any velocity maximum that is associated
%                     with an eddy boundary
%                     4) 'random': random ship track or velocity samples.
%                     Only determines an eddy centre. Further
%                     characteristics such as radius might be implemented
%                     in a future version.
%
%       ref         - Reference longitude,latitude for normalisation. Leave
%                     empty for taking the mean lon,lat as default or
%                     specify otherwise [lon lat].
%
%       x0,y0       - Scalar or vector of initialisation values in eastward
%                     and northward direction in [m] with respect to the
%                     reference longitude,latitude
%
%       rx          - Number of velocity samples that are iteratively
%                     excluded from the start and end of transect. All
%                     individual eddy estimates are averaged for the final
%                     results. Leave empty for default value rx = 0 for
%                     which no velocity samples are excluded.
%
%       NumOfIt     - Number of iterations. Leave empty for default value
%                     NumOfIt = 10. Note that the Gauss-Newton rate of
%                     convergence is approximately quadratic.
%
%       lambda      - Factor of penalty term. Leave empty for default value
%                     lambda = 0 for which the residual sum of squares is
%                     not penalised. Consider increasing the value for
%                     lambda when transect crosses (or is close to) the
%                     eddy centre to avoid overfitting. Large values tend
%                     to introduce underfitting.
%
%       (origin      - No longer included
%
%                     Define the origin of the coordinate system.
%
%                     1) 'eddycentre': the final eddy centre is taken as
%                     the coordinate system origin. The eddycentre with the
%                     lowest sum of squares is displayed
%                     2) 'reflonlat': the defined reference longitude,
%                     latitude (ref) is taken as the coordinate system
%                     origin. It displays all eddy centre estimate for the
%                     given cluster of initial values and highlights the
%                     one with lowest sum of squares.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   OUTPUT
%
%       output. ...
%
%       xvec,yvec          - Position vectors [m] relative to reference
%                            point
%
%       eddycentre_lonlat  - Eddy centre longitude,latitude
%
%       radius             - Radius [km] estimate for each given radial
%                            section
%
%       azimuthal_vel      - Maximum azimuthal velocity [m/s] for each
%                            given radial section
%
%       vorticity          - Inner core vorticity [1/s]
%
%       Ro                 - Rossby number of gradient flow
%
%       lambda             - Outer ring exponential decay scale [km] for 
%                            each given radial section:
%                            vtheta = Vmax*exp(-(r-Rmax)/lambda)
%
%       power              - Outer ring inverse power law exponent for each
%                            given radial section:
%                            vtheta = Vmax*(Rmax/r)^n
%
%       n0                 - Maximum sea surface heigth signal [cm] for
%                            each given radial section
%
%       radvel.radx,rady   - x,y components of the distance between eddy
%                            centre and maximum azimuthal velocity
%                            [degree]. Needed for plotting radius in m_map
%                            projection (r = sqrt(radx^2+rady^2))
%
%       centripetal        - Centripetal acceleration [m/s^2]
%
%       radial1,2          - Radius along ship track (separated in two
%                            radial sections [km]
%
%       azimuthal1,2       - Azimuthal velocity structure along ship track
%                            (separated in two radial sections) [m/s]
%
%       sumofsquares       - Residual sum of squares as a function of
%                            iteration for cluster that is associated with
%                            best eddy centre guess.
%
%       minsumofsquares    - Minimum sum of squares for each cluster
%
%       eddycentre_cluster - All eddy centre estimates for given cluster
%
%       mincluster         - Minimum sum of squares, cluster index, and
%                            corresponding eddy centre components [m] for
%                            each iterative section as given by rx
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   REQUIRED FUNCTIONS
%
%       sw_dist.m           - Calculate distance and angle between velocity
%                             samples (from m_map toolbox)
%       findpeaks.m         - Find velocity maxima associated with eddy
%                             boundary
%       lonlat2xy.m         - Convert longitude,latitude to Cartesian
%                             coordinates
%       normalise.m         - Normalse data in the range [-1 1] with
%                             preserving sign
%       rescale_diff.m      - Rescale normalised data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   EXAMPLE
%
%       Create cluster of starting values in [km] with respect to reference
%       longitude, latitude
%       x0 = linspace(-5,5,21)*10^3;
%       y0 = linspace(-5,5,21)*10^3;
%
%       [output] = eddycentre_GN(lon,lat,u,v,'crossed',x0,y0,[],[],[],[])
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% illustration and demo of shiptrack types
if nargin==0
    ShiptrackDemo
    return
end

%% vector size and default values

% The input lon,lat,u,v is the data for the whole section. The automatic
% data extraction is performed on the data that is set by range

% transpose if not a column vector
if iscolumn(lon)
    lon = lon';
    lat = lat';
    u = u';
    v = v';
end

lon_vis = lon;
lat_vis = lat;
u_vis = u;
v_vis = v;

lon = lon(range);
lat = lat(range);
u = u(range);
v = v(range);

% check wheter all input vector are of same size
if ~isequal(length(lon),length(lat),length(u),length(v))
    error('Input vectors lon,lat,u,v must all be of same size')
end

% reference point for normalisation
if isempty(ref)
    reflon = nanmean(lon); reflat = nanmean(lat);
%     reflon = nanmean(lon_GN); reflat = nanmean(lat_GN);
else
    reflon = ref(1);
    reflat = ref(2);
end

% number of samples exluded from ship track
if isempty(rx)
    rx = 0;
end

% number of iterations
if isempty(NumOfIt)
    NumOfIt = 10;
end

% penalty term factor
if isempty(lambda)
    lambda = 0;
end

% create longitude, latitude cluster
[LONCLUSTER,LATCLUSTER] = meshgrid(x0,y0);

%
if strcmp(shiptrack,'crossed')
    
    [pks,locs] = findpeaks(sqrt(u.^2+v.^2));
%     [pks,locs] = findpeaks(sqrt(u_GN.^2+v_GN.^2));
    
    if length(locs) >= 2
        
        [velpeak,velind] = sort(pks,'descend'); vel2peak = velpeak(1:2);
        rangeind = locs(min(velind(1:2))):locs(max(velind(1:2)));
        
        figure(1),clf
        set(gcf,'color','w','units','normalized','position',[0.1 0.2 0.8 0.6])
        
        subplot(121),hold on
        quiver(lon,lat,u,v,1); hold on, p1 = plot(reflon,reflat,'x','markersize',12,'linewidth',3);
%         quiver(lon_GN,lat_GN,u_GN,v_GN,1); hold on, p1 = plot(reflon,reflat,'x','markersize',12,'linewidth',3);
        legend(p1,'cluster reference point'),axis equal
        set(gca,'box','on','linewidth',1.2,'fontsize',14),grid on,xlabel('longitude'),ylabel('latitude')
        
        subplot(122),hold on
        co = get(gca,'colororder');
        area(rangeind,sqrt(u(rangeind).^2+v(rangeind).^2),'facecolor',[.7 .7 .7],'facealpha',0.5)
%         area(rangeind,sqrt(u_GN(rangeind).^2+v_GN(rangeind).^2),'facecolor',[.7 .7 .7],'facealpha',0.5)
        p1 = plot(sqrt(u.^2+v.^2),'o-','markerfacecolor',co(1,:),'markeredgecolor','k','markersize',6);
%         p1 = plot(sqrt(u_GN.^2+v_GN.^2),'o-','markerfacecolor',co(1,:),'markeredgecolor','k','markersize',6);
        p2 = plot(u,'o-','markerfacecolor',co(2,:),'markeredgecolor','k','markersize',5); p3 = plot(v,'o-','markerfacecolor',co(3,:),'markeredgecolor','k','markersize',5);
%         p2 = plot(u_GN,'o-','markerfacecolor',co(2,:),'markeredgecolor','k','markersize',5); p3 = plot(v_GN,'o-','markerfacecolor',co(3,:),'markeredgecolor','k','markersize',5);
        p4 = plot(rangeind, sqrt(u(rangeind).^2+v(rangeind).^2),'o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',8);
%         p4 = plot(rangeind, sqrt(u_GN(rangeind).^2+v_GN(rangeind).^2),'o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',8);
        plot(rangeind, u(rangeind),'o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',8);
        plot(rangeind, v(rangeind),'o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',8);
        line([rangeind(1) rangeind(1)],[sqrt(u(rangeind(1)).^2+v(rangeind(1)).^2) 0],'color','k','linestyle','--','linewidth',1.5)
        line([rangeind(end) rangeind(end)],[sqrt(u(rangeind(end)).^2+v(rangeind(end)).^2) 0],'color','k','linestyle','--','linewidth',1.5)
        p5 = plot([locs(min(velind(1:2))) locs(max(velind(1:2)))],[pks(min(velind(1:2))) pks(max(velind(1:2)))],'v','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',14);
        legend([p1,p2,p3,p4,p5],'speed','u','v','to be fitted','speed max','location','best')
        title('eddy entirely crossed')
        set(gca,'box','on','linewidth',1.2,'fontsize',14),grid on, ylabel('[m/s]'),xlabel('index')
        line([1 length(lon)],[0 0],'color','k','linewidth',2),xlim([1 length(lon)])
        
        if length(locs) > 2
            disp('More than two velocity maxima found. The two highest velocity maximum will be associated with eddy boundary. Correct? Press any key to continue...')
            pause
        elseif length(locs) == 2
            disp('Two velocity maxima found. Press any key to continue...')
            pause
        end
        
    elseif length(locs) == 1
        error('Only one velocity maximum found. Consider the shiptrack options "half-crossed", "within", or "random".')
    end
    
elseif strcmp(shiptrack,'half-crossed')
    
    [pks,locs] = findpeaks(sqrt(u.^2+v.^2));
    
    [velpeak,velind] = sort(pks,'descend'); vel2peak = velpeak(1);
    rangeind = locs(velind(1)):length(lon);
    
    figure(1),clf
    set(gcf,'color','w','units','normalized','position',[0.1 0.2 0.8 0.6])
    
    subplot(121),hold on
    quiver(lon,lat,u,v,1); hold on, p1 = plot(reflon,reflat,'x','markersize',12,'linewidth',3);
    legend(p1,'cluster reference point'),axis equal
    set(gca,'box','on','linewidth',1.2,'fontsize',14),grid on,xlabel('longitude'),ylabel('latitude')
    
    subplot(122),hold on
    co = get(gca,'colororder');
    area(rangeind,sqrt(u(rangeind).^2+v(rangeind).^2),'facecolor',[.7 .7 .7],'facealpha',0.5)
    p1 = plot(sqrt(u.^2+v.^2),'o-','markerfacecolor',co(1,:),'markeredgecolor','k','markersize',6);
    p2 = plot(u,'o-','markerfacecolor',co(2,:),'markeredgecolor','k','markersize',5); p3 = plot(v,'o-','markerfacecolor',co(3,:),'markeredgecolor','k','markersize',5);
    p4 = plot(rangeind(1):rangeind(end), sqrt(u(rangeind(1):rangeind(end)).^2+v(rangeind(1):rangeind(end)).^2),'o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',8);
    plot(rangeind, u(rangeind),'o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',8);
    plot(rangeind, v(rangeind),'o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',8);
    line([rangeind(1) rangeind(1)],[vel2peak 0],'color','k','linestyle','--','linewidth',1.5)
    line([rangeind(end) rangeind(end)],[vel2peak 0],'color','k','linestyle','--','linewidth',1.5)
    p5 = plot(rangeind(1),vel2peak ,'v','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',14);
    legend([p1,p2,p3,p4,p5],'speed','u','v','to be fitted','speed max','location','best')
    title('eddy not entirely crossed')
    set(gca,'box','on','linewidth',1.2,'fontsize',14),grid on, ylabel('[m/s]'),xlabel('index')
    line([1 length(lon)],[0 0],'color','k','linewidth',2),xlim([1 length(lon)])
    
    if length(locs) > 1
        disp('More than one velocity maximum found. Consider choosing shiptrack option "crossed" or press any key to continue.')
                pause
    end
    
    q = input('Does the section head from the outer ring to the inner core? (y/n) >','s');
    if strcmp(q,'y')
        
    elseif strcmp(q,'n')
        
        lon = fliplr(lon); lat = fliplr(lat); u = fliplr(u); v = fliplr(v);
        
        [pks,locs] = findpeaks(sqrt(u.^2+v.^2));
        
        [velpeak,velind] = sort(pks,'descend'); vel2peak = velpeak(1);
        rangeind = locs(velind(1)):length(lon);
        
        figure(1),clf
        set(gcf,'color','w','units','normalized','position',[0.2 0.2 0.6 0.4])
        
        subplot(121),hold on
        quiver(lon,lat,u,v,1); hold on, p1 = plot(reflon,reflat,'x','markersize',12,'linewidth',3);
        legend(p1,'cluster reference point'),axis equal
        set(gca,'box','on','linewidth',1.2,'fontsize',14),grid on,xlabel('longitude'),ylabel('latitude')
        
        subplot(122),hold on
        co = get(gca,'colororder');
        area(rangeind,sqrt(u(rangeind).^2+v(rangeind).^2),'facecolor',[.7 .7 .7],'facealpha',0.5)
        p1 = plot(sqrt(u.^2+v.^2),'o-','markerfacecolor',co(1,:),'markeredgecolor','k','markersize',6);
        p2 = plot(u,'o-','markerfacecolor',co(2,:),'markeredgecolor','k','markersize',5); p3 = plot(v,'o-','markerfacecolor',co(3,:),'markeredgecolor','k','markersize',5);
        p4 = plot(locs(velind(1)):length(lon), sqrt(u(locs(velind(1)):length(lon)).^2+v(locs(velind(1)):length(lon)).^2),'o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',8);
        plot(rangeind, u(rangeind),'o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',8);
        plot(rangeind, v(rangeind),'o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',8);
        line([rangeind(1) rangeind(1)],[vel2peak 0],'color','k','linestyle','--','linewidth',1.5)
        line([rangeind(end) rangeind(end)],[vel2peak 0],'color','k','linestyle','--','linewidth',1.5)
        p5 = plot(rangeind(1),vel2peak ,'v','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','k','markersize',14);
        legend([p1,p2,p3,p4,p5],'speed','u','v','to be fitted','speed max','location','southeast')
        title('eddy not entirely crossed')
        set(gca,'box','on','linewidth',1.2,'fontsize',14),grid on, ylabel('[m/s]'),xlabel('index')
        
    else
        error('Please type "y" for yes or "n" for no')
    end
    
elseif strcmp(shiptrack,'within')
    
    figure(1),clf
    set(gcf,'color','w','units','normalized','position',[0.2 0.2 0.6 0.4])
    
    subplot(121)
    q1 = quiver(lon,lat,u,v); hold on, p1 = plot(reflon,reflat,'x','markersize',12);
    grid on, legend(p1,'cluster reference point'), xlabel('lon'),ylabel('lat')
    
    subplot(122)
    plot(sqrt(u.^2+v.^2),'o-'), grid on, ylabel('m/s'),hold on,xlabel('index')
    
    rangeind = 1:length(lon);
    
    disp('The section lies within the eddy core. Optimise for the eddy centre. Press any key to continue...')
    pause
    
elseif strcmp(shiptrack,'random')
    
    rangeind = 1:length(lon);
    
else
    error('Unknown option for "shiptrack". Choose "crossed", "half-crossed", or "within".')
end

%% rotate ship track coordinates and velocities

% original input data
% lon_original = lon; lat_original = lat;
% u_original = u; v_original = v;

% if ~isempty(rotation_angle)
%     [lonr,latr,ur,vr] = rotate_uvxy(lon,lat,u,v,rotation_angle);
%     lon = lonr; lat = latr;
%     u = ur; v = vr;
% end

%% position vector for full ship track
[xd_full,yd_full] = lonlat2xy(lon,lat,reflon,reflat);
[xd_vis,yd_vis] = lonlat2xy(lon_vis,lat_vis,reflon,reflat);

% [xd_full_original,yd_full_original] = lonlat2xy(lon_original,lat_original,reflon,reflat);

%% optimisation

for cluster = 1:length(x0)*length(y0)
    
    for section_range = 1:rx+1
        
        range = cell(1,rx+1);
        
        if strcmp(shiptrack,'half-crossed')
            range{section_range} = rangeind(section_range):length(lon);
        else
            range{section_range} = rangeind(section_range):rangeind(end-(section_range-1));
        end
        
        %% cluster/section counter
        clc
        disp(strcat('cluster #',num2str(cluster),',','section #',num2str(section_range)))
        disp(strcat('progress:',num2str(round(100*cluster/(length(x0)*length(y0)),0)),'%'))
        
        %% longitude,latitude of velocity samples
        lontrackn{section_range} = lon(range{section_range});
        lattrackn{section_range} = lat(range{section_range});
        
        %% velocity samples
        un{section_range} = u(range{section_range});
        vn{section_range} = v(range{section_range});
        
        %% normalise velocity vectors in the range [-1,1]
        %         [un_norm{section_range},min_un{section_range},max_un{section_range}] = normalised_diff(un{section_range});
        %         [vn_norm{section_range},min_vn{section_range},max_vn{section_range}] = normalised_diff(vn{section_range});
        [un_norm{section_range},vn_norm{section_range},~,~] = normalise(un{section_range},vn{section_range});
        
        %% convert pairs of longitude,latitude into cartesian coordinates (x,y)
        [xd{section_range},yd{section_range}] = lonlat2xy(lontrackn{section_range},...
            lattrackn{section_range},reflon,reflat);
        
        %% normalise position vectors in the range [-1,1]
        %         [xd_norm{section_range},min_xd{section_range},max_xd{section_range}] = normalised_diff(xd{section_range});
        %         [yd_norm{section_range},min_yd{section_range},max_yd{section_range}] = normalised_diff(yd{section_range});
        
        %% normalise initial eddy center estimate in the range [-1,1]
        %         [xc_norm{section_range},min_xc{section_range},max_xc{section_range}] = normalised_diff([LONCLUSTER(cluster) xd{section_range}]);
        %         [yc_norm{section_range},min_yc{section_range},max_yc{section_range}] = normalised_diff([LATCLUSTER(cluster) yd{section_range}]);
        
        %% normalise position vectors and eddy centre estimate in the range [-1 1]
        [xcyc{section_range},xdyd{section_range},min_range{section_range},max_range{section_range}] = normalise([LONCLUSTER(cluster) LATCLUSTER(cluster)],[xd{section_range} yd{section_range}]);
        
        %% storage for eddy center estimates
        %         xyc_sum{cluster}{section_range}(1) = xc_norm{section_range}(1); % eddy centre x-component
        %         xyc_sum{cluster}{section_range}(2) = yc_norm{section_range}(1); % eddy centre y-component
        
        xyc_sum{cluster}{section_range}(1) = xcyc{section_range}(1); % eddy centre x-component
        xyc_sum{cluster}{section_range}(2) = xcyc{section_range}(2); % eddy centre y-component
        
        %% gradient descent/search direction vector
        p(1) = 0; p(2) = 0;
        
        %% residual vector
        
        uu = un_norm{section_range}; vv = vn_norm{section_range};
        %         xx = xc_norm{section_range}(2:end); yy = yc_norm{section_range}(2:end);
        xx = xdyd{section_range}(1:length(xd{section_range})); yy = xdyd{section_range}(length(xd{section_range})+1:end);
        
        for iteration = 1:NumOfIt
            
            xxc = xyc_sum{cluster}{section_range}(iteration,1);
            yyc = xyc_sum{cluster}{section_range}(iteration,2);
            
            F{cluster}{section_range}(iteration,:) = -uu.*cos(atan((yy-yyc)./(xx-xxc))) - vv.*sin(atan((yy-yyc)./(xx-xxc)));
            
            % sum of squared residuals
            sumofsquares{cluster}(section_range,iteration) = sum(F{cluster}{section_range}(iteration,:).^2) + lambda.*(xxc^2 + yyc^2);
            
            %% Jacobi matrix
            
            % dF/dxxc
            J{iteration}(:,1) = (uu.*(yy - yyc).^2)./((xx - xxc).^3.*((yy - yyc).^2./(xx - xxc).^2 + 1).^(3/2)) + ...
                (vv.*(yy - yyc).^3)./((xx - xxc).^4.*((yy - yyc).^2./(xx - xxc).^2 + 1).^(3/2)) - ...
                (vv.*(yy - yyc))./((xx - xxc).^2.*((yy - yyc).^2./(xx - xxc).^2 + 1).^(1/2));
            
            % dF/dyyc
            J{iteration}(:,2) = vv./((xx - xxc).*((yy - yyc).^2./(xx - xxc).^2 + 1).^(1/2)) - (uu.*(2.*yy - 2.*yyc))./(2.*(xx - xxc).^2.*((yy - yyc).^2./(xx - xxc).^2 + 1).^(3/2)) - (vv.*(2.*yy - 2.*yyc).*(yy - yyc))./(2.*(xx - xxc).^3.*((yy - yyc).^2./(xx - xxc).^2 + 1).^(3/2));
            
            
            %% determine search direction (solving set of normal equations)
            p(iteration+1,:) = (transpose(J{iteration})*J{iteration})\(-transpose(J{iteration})*F{cluster}{section_range}(iteration,:)');
            
            %% Armijo rule
            
            alpha = 1; % maximum step length
            tau = 0.5; % iterative step length factor
            c1 = 0.5; % Armijo constant
            
            for j = 1:10
                
                xyc_sum_armijo{cluster}{section_range}(iteration,:) = alpha*p(iteration+1,:) + xyc_sum{cluster}{section_range}(iteration,:);
                
                xxc_a = xyc_sum_armijo{cluster}{section_range}(iteration,1);
                yyc_a = xyc_sum_armijo{cluster}{section_range}(iteration,2);
                
                Farmijo{cluster}{section_range}(iteration,:) = -uu.*cos(atan((yy-yyc_a)./(xx-xxc_a))) - vv.*sin(atan((yy-yyc_a)./(xx-xxc_a)));
                
                sumofsquares_armijo{cluster}(section_range,iteration) = sum(Farmijo{cluster}{section_range}(iteration,:).^2) + lambda.*(xxc_a^2 + yyc_a^2);
                
                if sumofsquares_armijo{cluster}(section_range,iteration) <= sumofsquares{cluster}(section_range,iteration) + c1*alpha*(J{iteration}'*F{cluster}{section_range}(iteration,:)')'*p(iteration+1,:)'
                    break
                else
                    alpha = alpha*tau;
                end
                
            end
            
            %% update eddy center estimate using appropriate step length
            xyc_sum{cluster}{section_range}(iteration+1,:) = alpha*p(iteration+1,:) + xyc_sum{cluster}{section_range}(iteration,:);
            
        end
        
        %% rescale normalised data
        
        % eddy centre in [km]
        [eddycentre{cluster}{section_range}(:,1)] = rescale_diff(xyc_sum{cluster}{section_range}(:,1),min_range{section_range},max_range{section_range});
        [eddycentre{cluster}{section_range}(:,2)] = rescale_diff(xyc_sum{cluster}{section_range}(:,2),min_range{section_range},max_range{section_range});
        
        % theta
        theta_rescale{cluster}{section_range} = atand((yd{section_range}-eddycentre{cluster}{section_range}(end,2))./(xd{section_range}-eddycentre{cluster}{section_range}(end,1)));
        
        % azimuthal, radial velocity
        Vtheta{cluster}{section_range} = -un{section_range}.*sind(theta_rescale{cluster}{section_range}) + vn{section_range}.*cosd(theta_rescale{cluster}{section_range});
        Vr{cluster}{section_range} = un{section_range}.*cosd(theta_rescale{cluster}{section_range}) + vn{section_range}.*sind(theta_rescale{cluster}{section_range});
        
        %% minimum sum of squares
        
        minsumofsquares(section_range,cluster) = sumofsquares{cluster}(section_range,end);
        
        %% clear search direction and Jacobian
        
        clear p J
        
    end % section_range
    
end % cluster

if strcmp(shiptrack,'crossed')
    disp('Eddy entirely crossed. Optimisation satisfying? Derived characteristics are based on two radial sections.')
elseif strcmp(shiptrack,'half-crossed')
    disp('Eddy not entirely clossed. Optimisation satisfying? Derived characteristics are based on one radial section.')
    if strcmp(q,'n')
        disp('Input vectors (lon,lat,u,v) were transposed prior to optimisation')
    end
elseif strcmp(shiptrack,'within')
    disp('Section is within eddy core. Optimisation satisfying? Eddy radius and other characteristics cannot be computed.')
elseif strcmp(shiptrack,'random')
    disp('Optimisation satisfying? Derived characteristics may be included in a future version.')
end

%% Find best eddy centre estimates based on the minimum sum of squares

[mincluster(:,1),mincluster(:,2)] = min(minsumofsquares,[],2);

for section_range = 1:rx+1
    mincluster(section_range,3:4) = eddycentre{mincluster(section_range,2)}{section_range}(end,:);
end

eddystat(1,1:2) = [mean(mincluster(1:rx+1,3)) mean(mincluster(1:rx+1,4))];
if rx > 0
    eddystat(2,1:2) = [std(mincluster(1:rx+1,3)) std(mincluster(1:rx+1,4))];
end

%% Split each section in two separate radial sections

if ~strcmp(shiptrack,'random')
    radtrack = sqrt((xd_full-eddystat(1,1)).^2 + (yd_full-eddystat(1,2)).^2);
    [indvalue,indmin] = min(radtrack);
    
    radtrack_vis = sqrt((xd_vis-eddystat(1,1)).^2 + (yd_vis-eddystat(1,2)).^2);
    
    % theta
    theta_full = atand((yd_full-eddystat(1,2))./(xd_full-eddystat(1,1)));
    theta_vis = atand((yd_vis-eddystat(1,2))./(xd_vis-eddystat(1,1)));
    
    % azimuthal, radial velocity
    Vtheta_full = -u.*sind(theta_full) + v.*cosd(theta_full);
    Vtheta_vis = -u_vis.*sind(theta_vis) + v_vis.*cosd(theta_vis);
    
    %% Derive eddy characteristics
    
    % radius and azimuthal velocity for each radial section
    
    radial1 = fliplr(radtrack(1:indmin-1))*10^-3;
    radial2 = radtrack(indmin+1:end)*10^-3;
    azimuthal1 = fliplr(abs(Vtheta_full(1:indmin-1)));
    azimuthal2 = abs(Vtheta_full(indmin+1:end));
    
    % extended radius and azimuthal velocity for each radial section
    radial1_vis = fliplr(radtrack_vis(1:find(radtrack_vis == indvalue)-1))*10^-3;
    radial2_vis = radtrack_vis(find(radtrack_vis == indvalue)+1:end)*10^-3;
    azimuthal1_vis = fliplr(abs(Vtheta_vis(1:find(radtrack_vis == indvalue)-1)));
    azimuthal2_vis = abs(Vtheta_vis(find(radtrack_vis == indvalue)+1:end));
    
    % longitude, latitude for each radial section
    
    lon1 = fliplr(lon(1:indmin-1));
    lon2 = lon(indmin+1:end);
    lat1 = fliplr(lat(1:indmin-1));
    lat2 = lat(indmin+1:end);

    % longitude, latitude for each radial section
    lon1_vis = fliplr(lon_vis(1:find(radtrack_vis == indvalue)-1));
    lon2_vis = lon_vis(find(radtrack_vis == indvalue)+1:end);
    lat1_vis = fliplr(lat_vis(1:find(radtrack_vis == indvalue)-1));
    lat2_vis = lat_vis(find(radtrack_vis == indvalue)+1:end);
    
    % gravitational acceleration
    g = 9.81;
    
    % maximum azimuthal velocity and speed-based radius for first velocity maximum
    if strcmp(shiptrack,'crossed') || strcmp(shiptrack,'half-crossed')
        
        [pks1,locs1] = findpeaks(azimuthal1);
        [a1,b1] = sort(pks1,'descend');
        
        if ~isempty(pks1)
            
            % maximum azimuhtal velocity corresponding with maximum radius
            %         [a1,b1] = max(azimuthal1); % could raise problems with radial
            %         sections that are barely distinguishable from the surrounding
            % Rmax in km
            E1rad(1) = radial1(locs1(b1(1)));
            % Vmax in m/s
            E1vel(1) = a1(1);
            % Sense of rotation
            Velsign = fliplr(Vtheta_full(1:indmin-1));
            E1rot(1) = sign(Velsign(locs1(b1(1))));
            % Relative vorticity (zeta) in 1/s
            E1vort(1) = 2*(E1vel(1)/(E1rad(1)*10^3));
            % Coriolis parameter at rim in 1/s
            coriolis_rim(1) = 2*(2*pi/86400)*sind(lat1(locs1(b1(1))));
            % Rossby number of gradient flow
            E1Ro(1) = E1vel(1)/(coriolis_rim(1)*E1rad(1)*10^3);
            % Outer ring exponential decay scale in km: lambda
            % --> vtheta = Vmax*exp(-(r-Rmax)/lambda)
%             fexp1 = fit(radial1(locs1(b1(1)):end)'-radial1(locs1(b1(1))),azimuthal1(locs1(b1(1)):end)','exp1');
            fexp1 = fit(radial1_vis(locs1(b1(1)):end)'-radial1_vis(locs1(b1(1))),azimuthal1_vis(locs1(b1(1)):end)','exp1');
            E1lambda(1) = abs(1/fexp1.b);
            % Outer ring inverse power law: n
            % --> vtheta = Vmax*(Rmax/r)^n
%             fpow1 = fit(radial1_vis(locs1(b1(1)))./radial1_vis(locs1(b1(1)):end)',azimuthal1_vis(locs1(b1(1)):end)','power1');
            fpow1 = fit(radial1_vis(locs1(b1(1)))./radial1_vis(locs1(b1(1)):end)',azimuthal1_vis(locs1(b1(1)):end)','power1');
            E1power(1) = fpow1.b;
            % Sea surface height signal (n0) in cm
            coriolis_in(1) = 2*(2*pi/86400)*sind(lat1(1)); % Inner ring Coriolis parameter
%             coriolis_out(1) = 2*(2*pi/86400)*sind(lat1(end)); % Outer ring Coriolis parameter
            coriolis_out(1) = 2*(2*pi/86400)*sind(lat1_vis(end)); % Outer ring Coriolis parameter
            E1n0(1) = (E1vel(1)/(2*g))*(E1vel(1) + coriolis_in(1)*E1rad(1)*10^3 + 2*E1lambda(1)*10^3*coriolis_out(1))*10^2;
            % Centrifugal accerelation in m/s^2
            E1centripetal(1) = E1vel(1)^2/(E1rad(1)*10^3);
            
        else
            
            % maximum azimuhtal velocity corresponding with maximum radius
            [a1,b1] = max(azimuthal1); 
            % Rmax in km
            E1rad(1) = radial1(locs1(b1(1)));
            % Vmax in m/s
            E1vel(1) = a1(1);
            % Relative vorticity (zeta) in 1/s
            E1vort(1) = 2*(E1vel(1)/(E1rad(1)*10^3));
            % Coriolis parameter at rim in 1/s
            coriolis_rim(1) = 2*(2*pi/86400)*sind(lat1(locs1(b1(1))));
            % Rossby number of gradient flow
            E1Ro(1) = E1vel(1)/(coriolis_rim(1)*E1rad(1)*10^3);
%             % Outer ring decay scale (lambda) in km
%             f1 = fit(radial1(locs1(b1(1)):end)'-radial1(locs1(b1(1))),azimuthal1(locs1(b1(1)):end)','exp1');
%             E1lambda(1) = abs(1/f1.b);
%             % Sea surface height signal (n0) in cm
%             coriolis_in(1) = 2*(2*pi/86400)*sind(lat1(1)); % Inner ring Coriolis parameter
%             coriolis_out(1) = 2*(2*pi/86400)*sind(lat1(end)); % Outer ring Coriolis parameter
%             E1n0(1) = (E1vel(1)/(2*g))*(E1vel(1) + coriolis_in(1)*E1rad(1)*10^3 + 2*E1lambda(1)*10^3*coriolis_out(1))*10^2;
            % Centrifugal accerelation in m/s^2
            E1centripetal(1) = E1vel(1)^2/(E1rad(1)*10^3);
            
            disp('No velocity maximum found along first radial section. Outer ring decay is not computed. Use radius estimate with caution.')
            
        end       
        
    end
    
    % consider second velocity maximum if eddy was entirely crossed
    if strcmp(shiptrack,'crossed')
        
        [pks2,locs2] = findpeaks(azimuthal2);
        [a2,b2] = sort(pks2,'descend');
        
        if ~isempty(pks2)
            
            % maximum azimuhtal velocity corresponding with maximum radius
            %             [a2,b2] = max(azimuthal2); % could raise problems with radial
            %             sections that are barely distinguishable from the surrounding
            % Rmax in km
            E1rad(2) = radial2(locs2(b2(1)));
            % Vmax in m/s
            E1vel(2) = a2(1);   
            % Relative vorticity (zeta) in 1/s
            E1vort(2) = 2*(E1vel(2)/(E1rad(2)*10^3));
            % Coriolis parameter at rim in 1/s
            coriolis_rim(2) = 2*(2*pi/86400)*sind(lat2(locs2(b2(1))));
            % Rossby number of gradient flow
            E1Ro(2) = E1vel(2)/(coriolis_rim(2)*E1rad(2)*10^3);
            % Outer ring decay scale (lambda) in km
%             fexp2 = fit(radial2(locs2(b2(1)):end)'-radial2(locs2(b2(1))),azimuthal2(locs2(b2(1)):end)','exp1');
            fexp2 = fit(radial2_vis(locs2(b2(1)):end)'-radial2_vis(locs2(b2(1))),azimuthal2_vis(locs2(b2(1)):end)','exp1');
            E1lambda(2) = abs(1/fexp2.b);
            % Outer ring inverse power law: n
            % --> vtheta = Vmax*(Rmax/r)^n
%             fpow2 = fit(radial2(locs2(b2(1)))./radial2(locs2(b2(1)):end)',azimuthal2(locs2(b2(1)):end)','power1');
            fpow2 = fit(radial2_vis(locs2(b2(1)))./radial2_vis(locs2(b2(1)):end)',azimuthal2_vis(locs2(b2(1)):end)','power1');
            E1power(2) = fpow2.b;
            % Sea surface height signal (n0) in cm
            coriolis_in(2) = 2*(2*pi/86400)*sind(lat2(1)); % Inner ring Coriolis parameter
%             coriolis_out(2) = 2*(2*pi/86400)*sind(lat2(end)); % Outer ring Coriolis parameter
            coriolis_out(2) = 2*(2*pi/86400)*sind(lat2_vis(end)); % Outer ring Coriolis parameter            
            E1n0(2) = (E1vel(2)/(2*g))*(E1vel(2) + coriolis_in(2)*E1rad(2)*10^3 + 2*E1lambda(2)*10^3*coriolis_out(2))*10^2;
            % Centrifugal accerelation in m/s^2
            E1centripetal(2) = E1vel(2)^2/(E1rad(2)*10^3);
            
        else
            
            % maximum azimuhtal velocity corresponding with maximum radius
            [a2,b2] = max(azimuthal2);
            % Rmax in km
            E1rad(2) = radial2(b2);
            % Vmax in m/s
            E1vel(2) = a2;
            % Relative vorticity (zeta) in 1/s
            E1vort(2) = 2*(E1vel(2)/(E1rad(2)*10^3));
            % Coriolis parameter at rim in 1/s
            coriolis_rim(2) = 2*(2*pi/86400)*sind(lat2(b2));
            % Rossby number of gradient flow
            E1Ro(2) = E1vel(2)/(coriolis_rim(2)*E1rad(2)*10^3);
%             % Outer ring decay scale (lambda) in km
%             f2 = fit(radial2(locs2(b2(1)):end)'-radial2(locs2(b2(1))),azimuthal2(locs2(b2(1)):end)','exp1');
%             E1lambda(2) = abs(1/f2.b);
%             % Sea surface height signal (n0) in cm
%             coriolis_in(2) = 2*(2*pi/86400)*sind(lat2(1)); % Inner ring Coriolis parameter
%             coriolis_out(2) = 2*(2*pi/86400)*sind(lat2(end)); % Outer ring Coriolis parameter
%             E1n0(2) = (E1vel(2)/(2*g))*(E1vel(2) + coriolis_in(2)*E1rad(2)*10^3 + 2*E1lambda(2)*10^3*coriolis_out(2))*10^2;
            % Centrifugal accerelation in m/s^2
            E1centripetal(2) = E1vel(2)^2/(E1rad(2)*10^3);
            
            disp('No velocity maximum found along second radial section. Sea surface height signal and outer ring decay is not computed. Use radius estimate with caution.')
            
        end
    end
    
    % mean radius, azimuthal velocity, vorticity, Rossby number
    if strcmp(shiptrack,'crossed') || strcmp(shiptrack,'half-crossed')
        E1radmean = nanmean(E1rad); %E1radstd = nanstd(E1rad); % Rmax
    end

elseif strcmp(shiptrack,'random')
    
    radtrack = sqrt((xd_full-eddystat(1,1)).^2 + (yd_full-eddystat(1,2)).^2);
    
    % theta
    theta_full = atand((yd_full-eddystat(1,2))./(xd_full-eddystat(1,1)));
    
    % azimuthal, radial velocity
    Vtheta_full = -u.*sind(theta_full) + v.*cosd(theta_full);
    
%     radiusgrid = 0:10:round(max(radtrack),-1);
    
%     for i = 1:length(radiusgrid)-1
%         indgrid{i} = find(radtrack > radiusgrid(i) & radtrack < radiusgrid(i+1));
%         Vtheta_grid(i) = nanmedian(abs(Vtheta_full(indgrid{i})));
%     end
    
%     [pks1,locs1] = findpeaks(Vtheta_grid);
%     [a1,b1] = sort(pks1,'descend');
    
%     E1rad = radiusgrid(locs1(b1(1)));
%     E1vel = Vtheta_grid(locs1(b1(1)));
   
    [radtrack_sorted,sorted_ind] = sort(radtrack*10^-3);
    
    p1 = polyfit(radtrack_sorted,abs(Vtheta_full(sorted_ind)),4);
    pfit1 = polyval(p1,radtrack_sorted);
%     plot(radtrack_sorted,pfit1,'color','k')

    [pks1,locs1] = findpeaks(pfit1);
    [a1,b1] = sort(pks1,'descend');
    
    E1rad = radtrack_sorted(locs1(b1(1)));
    E1vel = pfit1(locs1(b1(1)));
    
    E1velmean = E1vel; E1radmean = E1rad;

end

%% Express eddy centre back to longitude, latitude

R = 6371*10^3; % radius in km

newlon = reflon + (eddystat(1,1)/R) * (180/pi)/(cos(reflat*pi/180));
newlat = reflat + (eddystat(1,2)/R) * (180/pi);

eddycentre_lonlatrad = [newlon newlat];

%% Eddy radius expressed in degrees

% determine the mean radius considering both velocity maxima
if strcmp(shiptrack,'crossed')
    radx = mean(abs(newlon - [(newlon + ((xd{1}(1)-eddystat(1,1))/R) * (180/pi)/(cos(newlat*pi/180))) (newlon + ((xd{1}(end)-eddystat(1,1))/R) * (180/pi)/(cos(newlat*pi/180)))]));
    rady = mean(abs(newlat - [(newlat + ((yd{1}(1)-eddystat(1,2))/R) * (180/pi)) (newlat + ((yd{1}(end)-eddystat(1,2))/R) * (180/pi))]));
elseif strcmp(shiptrack,'half-crossed')
    radx = mean(abs(newlon - [(newlon + ((xd{1}(1)-eddystat(1,1))/R) * (180/pi)/(cos(newlat*pi/180)))]));
    rady = mean(abs(newlat - [(newlat + ((yd{1}(1)-eddystat(1,2))/R) * (180/pi))]));
elseif strcmp(shiptrack,'random')
    radx = mean(abs(newlon - [(newlon + ((E1radmean*10^3)/R) * (180/pi)/(cos(newlat*pi/180)))]));
    rady = mean(abs(newlat - [(newlat + (0/R) * (180/pi)) (newlat + ((yd{1}(end)-eddystat(1,2))/R) * (180/pi))]));
end

%% Plot result

% if strcmp(origin,'eddycentre')

figure(2),clf,hold on
set(gcf,'color','w','units','normalized','position',[0.05 0.1 0.9 0.75])

% eddy center map
subplot('position',[0.01 0.2 0.5 0.7]),hold on

% longitude, latitude cluster
plot((LONCLUSTER - eddystat(1,1))*10^-3,(LATCLUSTER - eddystat(1,2))*10^-3,'.','color',[0.7 0.7 0.7]);
l4 = plot(x0(1),y0(1),'.','color',[0.7 0.7 0.7]);

% mean eddy radius
if ~strcmp(shiptrack,'within')
    th = linspace(0,2*pi,100); xcirc = E1radmean*cos(th) + 0; ycirc = E1radmean*sin(th) + 0;
    plot(xcirc,ycirc,'color',[.7 .7 .7],'linewidth',2)
end

% velocity vectors
% l1 = quiver((xd_full-eddystat(1,1))*10^-3,(yd_full-eddystat(1,2))*10^-3,u,v,1.5,'color','black','linewidth',1.5);
l1 = quiver((xd_vis-eddystat(1,1))*10^-3,(yd_vis-eddystat(1,2))*10^-3,u_vis,v_vis,1.5,'color','black','linewidth',1.5);

% eddy centre
l3 = plot(0,0,'p','markerfacecolor','k','markeredgecolor','k','markersize',12);
grid on,xlabel('eastward distance [km]'),ylabel('northward distance [km]'), axis equal
% axrange = round(max(abs((sqrt(xd_full.^2 + yd_full.^2)-sqrt(eddystat(1,1).^2 + eddystat(1,2).^2))*10^-3))*1.5,-1);
axrange = round(max(abs((sqrt(xd_vis.^2 + yd_vis.^2)-sqrt(eddystat(1,1).^2 + eddystat(1,2).^2))*10^-3))*1.1,-1);
xlim([-axrange axrange]),ylim([-axrange axrange])
if axrange <= 100
    set(gca,'box','on','linewidth',1.3,'tickdir','out','fontsize',14,'xtick',-axrange:10:axrange,'ytick',-axrange:10:axrange)   
elseif axrange <= 200 && axrange > 100
    set(gca,'box','on','linewidth',1.3,'tickdir','out','fontsize',14,'xtick',-axrange:20:axrange,'ytick',-axrange:20:axrange)
elseif axrange > 200
    set(gca,'box','on','linewidth',1.3,'tickdir','out','fontsize',14,'xtick',-axrange:50:axrange,'ytick',-axrange:50:axrange)
end

% print converted eddy centre to longitude, latitude; mean radius
text(0.05,0.15,strcat('lon=',num2str(round(newlon,2)),',',{' '},'lat=',num2str(round(newlat,2))),'units','normalized','fontsize',16,'fontweight','bold')
if strcmp(shiptrack,'crossed') || strcmp(shiptrack,'random') 
    text(0.05,0.05,strcat('mean radius=',num2str(round(E1radmean,1)),'km'),'units','normalized','fontsize',16,'fontweight','bold')
elseif strcmp(shiptrack,'half-crossed')
    text(0.05,0.05,strcat('radius=',num2str(round(E1radmean,1)),'km'),'units','normalized','fontsize',16,'fontweight','bold')
end

if rx == 0
    l2 = plot(eddycentre{mincluster(1,2)}{1}(:,1)*10^-3 - eddycentre{mincluster(1,2)}{1}(end,1)*10^-3,...
        eddycentre{mincluster(1,2)}{1}(:,2)*10^-3 - eddycentre{mincluster(1,2)}{1}(end,2)*10^-3,'o-','color','r','linewidth',1.2);
    legend([l1,l2,l3,l4],'velocity vector',sprintf('eddy centre\niteration steps'),'eddy centre','initial value cluster','location','northeast')
else
    legend([l1,l3,l4],'velocity vector','mean eddy centre','initial value cluster','location','northeast')
end

plot(0,0,'p','markerfacecolor','k','markeredgecolor','k','markersize',12);

if ~strcmp(shiptrack,'random')
    subplot('position',[0.55 0.62 0.3 0.29]),hold on
    
    if length(radial1) > 1
%         plot(radial1,azimuthal1,'o-','markeredgecolor','k','markerfacecolor',[0.7 0.7 0.7],'color','k')
        plot(radial1_vis,azimuthal1_vis,'o-','markeredgecolor','k','markerfacecolor',[0.7 0.7 0.7],'color','k')
    end
    
    if length(radial2) > 1
%         plot(radial2,azimuthal2,'o-','markeredgecolor','k','markerfacecolor',[0.7 0.7 0.7],'color','k')
        plot(radial2_vis,azimuthal2_vis,'o-','markeredgecolor','k','markerfacecolor',[0.7 0.7 0.7],'color','k')
    end
    
    set(gca,'box','on','linewidth',1.3,'tickdir','out','fontsize',14), grid on
    xlim([0 max([radial1_vis radial2_vis])])
    xlabel('radius [km]'),ylabel('azimuthal velocity [m/s]')
    
elseif strcmp(shiptrack,'random')
    
    subplot('position',[0.55 0.62 0.3 0.29]),hold on
    
    plot(radtrack*10^-3,abs(Vtheta_full),'o','markeredgecolor','k','markerfacecolor',[0.7 0.7 0.7],'color','k')

    set(gca,'box','on','linewidth',1.3,'tickdir','out','fontsize',14), grid on
%     xlim([0 max([radial1 radial2])])
    xlabel('radius [km]'),ylabel('azimuthal velocity [m/s]')

end

subplot('position',[0.55 0.2 0.3 0.29]), hold on
for section_range = 1:rx+1
    plot(1:NumOfIt,sumofsquares{mincluster(section_range,2)}(section_range,:),'o-','markeredgecolor','k','markerfacecolor',[0.7 0.7 0.7],'color','k')
end
xlabel('iteration'),ylabel('sum of squares')
text(0.35,0.8,strcat('$$\min f(x) = \sum_{j=1}^{m} r_j^2(x)$$ =',{' '},num2str(round(mincluster(1,1),3))),'Interpreter','latex',...
    'units','normalized','FontSize',16);
set(gca,'box','on','linewidth',1.3,'tickdir','out','fontsize',14), grid on
xlim([1 NumOfIt])

%% Output variable

% position vectors
output.xvec = xd_full; output.yvec = yd_full;

% eddy centre estimate in degrees
output.eddycentre_lonlat = eddycentre_lonlatrad;

% derived eddy characteristics
if strcmp(shiptrack,'crossed') || strcmp(shiptrack,'half-crossed')
    output.radius = E1rad;
    output.azimuthal_vel = E1vel;
    output.rotation = E1rot;
    output.vorticity = E1vort;
    output.Ro = E1Ro;
    output.lambda = E1lambda;
    output.power = E1power;
    output.n0 = E1n0;
    output.centripetal = E1centripetal;
    output.radvel.radx = radx; output.radvel.rady = rady;
    
    if length(radial1) > 1
        output.radial1 = radial1;
        output.azimuthal1 = azimuthal1;
    end
    
    if length(radial2) > 1
        output.radial2 = radial2;
        output.azimuthal2 = azimuthal2;
    end
    
elseif strcmp(shiptrack,'random')
    
    output.radius = E1rad;
    output.azimuthal_vel = E1vel;
    output.radvel.radx = radx; output.radvel.rady = rady;
    
end

% Reduction of residual sum of squares for
output.sumofsquares = sumofsquares{mincluster(1,2)};

% minimum sum of squares for each cluster index
output.minsumofsquares = minsumofsquares;

% all eddy centre estimate for given cluster
output.eddycentre_cluster = eddycentre;

% cluster index and eddy center estimate with lowest sum of squares
output.mincluster = mincluster;

end
