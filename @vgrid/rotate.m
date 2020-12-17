function g = rotate(g,r,p,h)
%% A function to rotate a grid.
%
% Val Schmidt
% Center for Coastal and Ocean Mapping
% University of New Hampshire
% Copyright 2010-2020, All rights reserved.
[xx yy] = meshgrid(g.xx,g.yy);
[g.xx g.yy g.zw C] = rotate_uncert(r*pi/180,p*pi/180,h*pi/180,0,0,0,xxx,yyy,g.zw,0,0,0);