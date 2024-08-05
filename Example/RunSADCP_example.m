%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EddyReconstructionExample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load SADCP and AVISO data
load sadcp_example

AVISO = nc_getall_nostruct('dataset-duacs-nrt-global-merged-allsat-phy-l4_1560507452500.nc');

AVISO.sla(AVISO.sla < -2*10^5) = NaN;
AVISO.adt(AVISO.adt < -2*10^5) = NaN;
AVISO.ugos(AVISO.ugos < -2*10^5) = NaN;
AVISO.vgos(AVISO.vgos < -2*10^5) = NaN;
AVISO.ugosa(AVISO.ugosa < -2*10^5) = NaN;
AVISO.vgosa(AVISO.vgosa < -2*10^5) = NaN;

[AVISO.LON, AVISO.SLA, AVISO.UGOS, AVISO.VGOS, AVISO.UGOSA, AVISO.VGOSA, AVISO.ADT] = ...
    lon360to180(AVISO.longitude,'sort',AVISO.sla*0.0001,AVISO.ugos*0.0001,AVISO.vgos*0.0001,AVISO.ugosa*0.0001,AVISO.vgosa*0.0001,AVISO.adt*0.0001);
AVISO.dt_SLA = datetime(AVISO.time*24*3600, 'ConvertFrom', 'epochtime', 'Epoch', '1950-01-01');
[AVISO.LONAVISO, AVISO.LATAVISO] = meshgrid(AVISO.LON,AVISO.latitude);
AVISO.date = cellstr(AVISO.dt_SLA);

%% Eddy reconstruction

% lon,lat
lonvec_all = sadcp_example.lon; latvec_all = sadcp_example.lat;
% u,v
uvec_all = sadcp_example.UCUR; vvec_all = sadcp_example.VCUR;
% mtime
timevec_all = sadcp_example.mtime;
% gridding and smoothing options
gridding = 'gridding';
% averaging interval
avg_int = 1;
% depth vector and depth level indices of interest (scalar or vector)
depthvec = sadcp_example.depth; depth_int = 1:12;
% specify indices of interest
ind = []; 
% delete indices during CTD stations or sections off track 
inddel = [115:268,402:550]; 
% smoothing parameters [xrad(km) xcut(km) yrad(m) ycut(m)]
smooth = [2 4 20 40];
% name of cruise,etc. and numeration of section
name = 'MSM74'; nr = 1; subnr = 1; 
% area of interest [lon1 lon2 lat1 lat2] and xaxis of contourf plot
lonlatbound = [-51 -48 53 55.5]; xaxis = 'lataxis';
% save and path
save = 'no'; projectdir = 'C:\Users\';
% optimisation
GN = 'yes';
% shiptrack option
shiptrack = 'crossed'; 
% reference point and initial value cluster
ref = []; x0 = linspace(-50,50,21)*10^3; y0 = linspace(-50,50,21)*10^3;
% number of samples excluded from eddy rim
rx = []; NumOfIt = []; lambda = []; origin = 'eddycentre'; r = [40 20];

[output,xyuv] = eddyreconstruction(lonvec_all,latvec_all,uvec_all,vvec_all,timevec_all,...
    gridding,avg_int,depthvec,depth_int,ind,inddel,smooth,nr,subnr,lonlatbound,AVISO,...
    name,xaxis,save,projectdir,GN,r,shiptrack,x0,y0,ref,rx,NumOfIt,lambda);