function y=obana2(x,posx,posy,gridx,gridy,r1x,r2x,r1y,r2y);
% function new_z=obana2(z,x,y,new_x,new_y,r_x_inf,r_x_cut,r_y_inf,r_y_cut);
% 
% averaging with gaussian weights
%
% input  :	z		- old values
%         	x		- old x-positions
%         	y		- old y-positions
%		new_x		- new x-positions
%		new_y		- new y-positions
%		r_x_inf		- influence radius in x-direction
%		r_x_cut		- cut-off radius in x-direction
%		[r_y_inf]	- influence radius in y-direction
%		[r_y_cut]	- cut-off radius in y-direction
%
% output :	new_z		- new values
%
% influence and cut-off radii may be given as vectors or matrices
%
% uses :	gauss.m sumnan.m
%
% version 1.1.0		last change 19.06.1997

% M. Visbeck
%	changed to variable influence and cut-off radii
%	G.Krahmann, IfM Kiel, Jun 1997

% make isotropic if no y values are specified
if nargin<9, r1y=r1x;
 if nargin<8, r2y=r2x; end
end

% blow up influence and cut-off radii
if size(r1x,2)==1
  r1x=r1x*ones(1,size(gridx,2));
end
if size(r2x,2)==1
  r2x=r2x*ones(1,size(gridx,2));
end
if size(r1x,1)==1
  r1x=ones(size(gridx,1),1)*r1x;
end
if size(r2x,1)==1
  r2x=ones(size(gridx,1),1)*r2x;
end
if size(r1y,2)==1
  r1y=r1y*ones(1,size(gridx,2));
end
if size(r2y,2)==1
  r2y=r2y*ones(1,size(gridx,2));
end
if size(r1y,1)==1
  r1y=ones(size(gridx,1),1)*r1y;
end
if size(r2y,1)==1
  r2y=ones(size(gridx,1),1)*r2y;
end
r1x=r1x(:)';
r2x=r2x(:)';
r1y=r1y(:)';
r2y=r2y(:)';

% setup complex target vector positions
i=sqrt(-1);
[ly,lx]=size(gridx);
gridv=gridx(:)'+i*gridy(:)';


% setup complex input vector positions
xv=mkvec(x);
ld=length(xv);
posv=mkvec(posx)+i*mkvec(posy);

% reset sums
nn=zeros(1,ly*lx);
yv=zeros(1,ly*lx);

% loop over each output value
pc=0;
lm=ly*lx;
for j=1:lm

  % display process
  if (j/lm>pc) 
    disp([' obana now ',int2str(pc*100),'%'])
    pc=pc+0.1; 
  end

  % positions difference
  dp=gridv(j)-posv;

  % norm with cutoff radius
  cr=real(dp)/r2x(j)+i*imag(dp)/r2y(j);

  % select only values within cutoff radius
  ix=find( abs(cr)<1);
  if length(ix)>0

    % norm with inflence radius
    ir=real(dp(ix))/r1x(j)+i*imag(dp(ix))/r1y(j);

    % get factors using a gauss distribution
    n=gauss(abs(ir),1);

    % sum up values
    nn=sumnan(n);
    if nn>0
      yv(j)=(n(:)'*xv(ix)')/nn;
    else
      yv(j)=nan;
    end
  else
    yv(j)=nan;
  end
end

% reshape to target position size
yv=yv-2*sqrt(-1)*imag(yv);
y=reshape(yv,ly,lx);
