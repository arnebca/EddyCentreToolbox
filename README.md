# EddyCentreToolbox
Toolbox to reconstruct an eddy centre based on a 1-dim ship track horizontal current velocity profiler data.

# README
 
Toolbox explicitly designed to reconstruct an eddy centre based on a
1-dim ship track of horizontal current velocity profiler data (SADCP,
LADCP), inspired by Castelao and Johns (2011) and Castelao et al. (2013). 
In an iterative procedure, this algorithm converts zonal and meridional
velocity into azimuthal and radial velocity before minimising the radial
velocity component using a non-linear damped Gauss-Newton optimisation
method.

This toolbox includes cmocean colormaps (Thyng et al., 2016). In addition
to the below mentioned functions, please add the m_map (Pawlowicz, 2020)
to your working directory.

If this toolbox has significant impact on your research, please reference
the publication Bendinger et al. (submitted to Journal Of Physical
Oceanography) and/or the software (Zenodo, available soon). 

For a more detailed description, refer two README_HowToUse

# References

Castelão, G. P., & Johns, W. E. (2011). Sea surface structure of North
 	Brazil Current rings derived from shipboard and moored acoustic 
	Doppler current profiler observations. Journal of Geophysical 
	Research: Oceans, 116(C1), http://doi.org/10.1029/2010JC006575. 

Castelão, G. P., Irber Jr, L. C., & Boas, A. B. V., 2013. An objective
	reference system for studying rings in the ocean. Computers & 
	Geosciences, 61, 43-49, http://doi.org/10.1016/j.cageo.2013.07.004.

Pawlowicz, R., 2020. M_Map: A mapping package for MATLAB,
	version 1.4m, [Computer software], available online at
	www.eoas.ubc.ca/~rich/map.html.

Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco.,
	2016. True colors of oceanography: Guidelines for effective and
	accurate colormap selection. Oceanography 29(3):9–13.
	http://dx.doi.org/10.5670/oceanog.2016.66.



