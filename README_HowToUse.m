%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% README_HowToUse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Toolbox explicitly designed to reconstruct an eddy centre based on a
% 1-dim ship track of horizontal current velocity profiler data (SADCP,
% LADCP), inspired by Castelao and Johns (2011) and Castelao et al. (2013). 
% In an iterative procedure, this algorithm converts zonal and meridional
% velocity into azimuthal and radial velocity before minimising the radial
% velocity component using a non-linear damped Gauss-Newton optimisation
% method. More information can be found in Bendinger et al. (submitted to 
% Journal of Physical Oceanography).
%
% This toolbox includes cmocean colormaps (Thyng et al., 2016). In addition
% to the below mentioned functions, please add the m_map (Pawlowicz, 2020)
% to your working directory.
%
% If this toolbox has significant impact on your research, please reference
% the publication Bendinger et al. (submitted to Journal Of Physical
% Oceanography) and/or the software (Zenodo). 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTENT
% 
% EddyCentreToolbox\
% 
% eddyreconstruction.m - This function takes the gridded or ungridded
%                        velocity data as well as optional gridded
%                        satellite altimetry data products from AVISO/CMEMS
%                        as input, prepares and calls the opimisation and 
%                        visualises the reconstructed eddies. See the 
%                        eddyreconstruction.m documentation for more
%                        details.
% 
% ReadMe_Content.m     - Overview list of all scripts and functions.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ..\Example           - Call RunSADCP_example.m for a given SADCP (.mat)
%                        and AVISO (.nc) dataset and run an optimisation.
%
% This study has been conducted using EU Copernicus Marine Service
% Information CMEMS (https://doi.org/10.48670/moi-00148).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ..\GaussNewton
% 
% eddycentre_GN.m      - The optimsation itself. See the eddycentre_GN.m
%                        documentation for more details.
% 
% lonlat2xy.m          - Converts the position vectors given in longitude
%                        and latitude into Cartesian coordinates in
%                        eastward and northward distance with respect to a
%                        point of reference.
% 
% normalise.m          - Normalises data in the range [-1 1] while
%                        preserving the original sign.
% 
% rescale_diff.m       - Rescales the normalised eddy centre estimates back
%                        to longitude and latitude.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ..\Interpolation     - Interpolation with Gaussian weights (Credits to 
			 Martin Visbeck and Gerd Krahmann)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ..\OtherFunctions      
% 
% export_fig.m         - Save figure (Copyright (c) 2014, Oliver J.
%			 Woodford, Yair M. Altman. All rights reserved.
%			 https://github.com/altmany/export_fig)
% 
% lon360to180.m        - Convert longitude format from 0-360° to -180° to
%                        180°
% 
% nc_getall_nostrcut.m - Open netCDF and put variables in a structure
% 
% ShiptrackDemo.m      - Call ShiptrackDemo() for a small demonstration for
%                        the section types "crossed", "half-crossed",
%                        "within", and "manual".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References
%
% Castelão, G. P., & Johns, W. E. (2011). Sea surface structure of North
% 	Brazil Current rings derived from shipboard and moored acoustic 
%	Doppler current profiler observations. Journal of Geophysical 
%	Research: Oceans, 116(C1). 
%
% Castelão, G. P., Irber Jr, L. C., & Boas, A. B. V., 2013. An objective
%	reference system for studying rings in the ocean. Computers & 
%	Geosciences, 61, 43-49.
%
% Pawlowicz, R., 2020. M_Map: A mapping package for MATLAB,
%	version 1.4m, [Computer software], available online at
%	www.eoas.ubc.ca/~rich/map.html.
%
% Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco.,
%	2016. True colors of oceanography: Guidelines for effective and
%	accurate colormap selection. Oceanography 29(3):9–13.
%	http://dx.doi.org/10.5670/oceanog.2016.66.