function y=gauss(x,s,m)
% function y=gauss(x,s,m)
%
% Gaussian function
%
%    y=exp(-(x-m).^2./s.^2)./(sqrt(2*pi).*s);
%
% Bronstein p. 81
if nargin<3, m=0;
 if nargin<2, s=1; end
end

y=exp(-(x-m).^2./s.^2)./(sqrt(2*pi).*s);
