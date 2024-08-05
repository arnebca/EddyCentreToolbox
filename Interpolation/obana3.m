function znew = obana3(z,x,y,gridx,gridy,xrad,xcut,yrad,ycut)

% OBANA3 fast data gridding using Gaussian weights. 
% 
%   znew = obana3(z,x,y,gridx,gridy,xrad,xcut,yrad,ycut)
%
%   where  znew        : gridded matrix
%          z           : data values to be gridded
%          x,y         : position of data values
%          gridx,gridy : new grid points
%          xrad        : influence radius in x-direction
%          xcut        : cut-off radius in x-direction
%         [yrad]       : influence radius in y-direction
%         [ycut]       : cut-off radius in y-direction
%         
%   USES: obana2.{m,mexlx,mexaxp,mexsol or other}
%
%	  yrad = 'geo' uses geographical distances
%
%   NOTE: The mex-file version is about 4 times faster than
%         the m-file equivalent.
% 
%   See also OBANA, GRIDDATA, OBJMAP, KRIGING.

% adapted from the original Matlab m-file obana.m by Martin Visbeck
% with the additional extensions of obana2.m by Gerd Krahmann
%
% Matlab 4.2c mex-file 
% by Ulf Garternicht in Jun 97 at IfM Kiel
%
% Version 1.1 (Aug 97)

% check number of I/O parameters

if (nargin~=7 & nargin~=9 &nargin~=8)
  error('USAGE: znew=obana3(z,x,y,gridx,gridy,xrad,xcut,yrad,yrad)');
end
if (nargout>1)
  error('One output argument only.');
end

% make isotropic if no y values are specified

if (nargin==7)
   yrad=xrad;
   ycut=xcut;
end

if (nargin>7)
   if isstr(yrad)
     if strcmp(yrad,'geo')
       yrad=-1;
       ycut=0;
     else
       error('What ?');
     end
   end
end
     
% make column vector of input data

x=x(:); y=y(:); z=z(:);

% remove NaN values     

a=~isnan(z);
x=x(a); y=y(a); z=z(a);

% blow up influence and cut-off radii

if size(xrad,2)==1
  xrad=xrad*ones(1,size(gridx,2));
end
if size(xrad,1)==1
  xrad=ones(size(gridx,1),1)*xrad;
end

if size(xcut,2)==1
  xcut=xcut*ones(1,size(gridx,2));
end
if size(xcut,1)==1
  xcut=ones(size(gridx,1),1)*xcut;
end

if size(yrad,2)==1
  yrad=yrad*ones(1,size(gridx,2));
end
if size(yrad,1)==1
  yrad=ones(size(gridx,1),1)*yrad;
end

if size(ycut,2)==1
  ycut=ycut*ones(1,size(gridx,2));
end
if size(ycut,1)==1
  ycut=ones(size(gridx,1),1)*ycut;
end

% do the actual function call
znew=obana2(z,x,y,gridx,gridy,xrad,xcut,yrad,ycut);
