function [ tapper] = get_tapper( i, j, it, nx, nz, nt, dx, dz, dt, x0, z0 ,t0 )
%GET_TAPER Summary of this function goes here
%   Detailed explanation goes here
xdist = (nx/2 - abs(i - nx/2)) * dx;
% zdist = (nz/2 - abs(j - nz/2)) * dz;
tdist = (nt/2 - abs(it - nt/2)) * dt;

xtapper = 1;
ttapper = 1;
if xdist < x0
xtapper = 1/2 * ( 1- cos(pi * xdist/x0) ) ;
end
% ztapper = 1/2 * ( 1- cos(pi * zdist/z0) );
if tdist < t0
ttapper = 1/2 * ( 1- cos(pi * tdist/t0) );
end
% disp(xdist)
tapper = xtapper  * ttapper;


end

