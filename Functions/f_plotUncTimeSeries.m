function [h] = f_plotUncTimeSeries( t,x,col,LW,up,lo )
%F_PLOTUNCTIMESERIES Summary of this function goes here
%   Detailed explanation goes here

t = t(:)';
x = x(:)';

if nargin == 6    
    up = up(:)';
    lo = lo(:)';
end
if nargin == 5
    sig = up(:)';
    up = x+sig;
    lo = x-sig;
end

yp = [up,fliplr(lo)];
xp = [t,fliplr(t)];
    
fill(xp,yp,col,'facecolor',col,'edgealpha',0,'facealpha',0.15);
h = plot(t,x,'-','Color',col,'LineWidth',LW);

tmp1 = plot( t,up,'-','Color',col);
tmp1.Color(4) = 0.5;
tmp2 = plot( t,lo,'-','Color',col);
tmp2.Color(4) = 0.5;


end

