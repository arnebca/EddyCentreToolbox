function y=mkvec(x)
%make vector y of matrix x
%M. Visbeck 23.07.90
[lr,lc]=size(x);
y=reshape(x,1,lr*lc);
