function [h] = f_plotUncTimeSeries( t,x,col,LW,up,lo )
%F_PLOTUNCTIMESERIES plots a solid line and shaded area indicating
%uncertainty
%   IN:
%       t: x values
%       x: y values of solid line
%       col: RGB triplet specifining color of solid line
%       LW: scalare specifying line width
%       up: upper uncertainty limit
%       lo: lower uncertainty limit. If lo is missing, the function assumes
%       up to be the symmetric standard deviation
%   OUT:
%       h: axes handle

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

