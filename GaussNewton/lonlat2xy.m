function [xd,yd] = lonlat2xy(lon,lat,lon_0,lat_0)
%==========================================================================
%
% lonlat2xc.m converts pairs of longitude and latitude into cartesian
% coordinates x and y with lon_0 and lat_0 representing the coordinate
% system origin. Xd and yd are given in km.
%
%==========================================================================

% convert km to m
km2m = 10^3;

for i = 1:length(lon)
    
    if lat(i) < lat_0
        yd(i) = sw_dist([lat(i) lat_0],[0 0],'km')*(-1)*km2m;
    else
        yd(i) = sw_dist([lat(i) lat_0],[0 0],'km')*km2m;
        
    end
    
    if lon(i) < lon_0
        xd(i) = sw_dist([lat(i) lat(i)],[lon(i) lon_0],'km')*(-1)*km2m;
    else
        xd(i) = sw_dist([lat(i) lat(i)],[lon(i) lon_0],'km')*km2m;
    end
end

return