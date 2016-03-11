function TEGenerateData(srcloc,recloc, xsrcpulse, zsrcpulse, T,isrc,outstep,plotopt)
% Example for running TE_model2d.m
% (a 2-D, FDTD, crosshole radar and VRP modeling code in MATLAB)
%
% by James Irving
% July 2005

% load the earth model example from file
%load crosshole_model
mu = 0;
load rtm_model

% calculate minimum and maximum relative permittivity and permeability
% in the model (to be used in finddx.m and finddt.m) 
epmin = min(min(ep));
epmax = max(max(ep));
mumin = min(min(mu));
mumax = max(max(mu));

% create time (s) and source pulse vectors for use with finddx.m
% (set dt very small to not alias frequency components;  set again below)
% (maximum time should be set to that whole source pulse is included)
%t=0:1e-10:100e-9;
%srcpulse=blackharrispulse(100e6,t);

% use finddx.m to determine maximum possible spatial field discretization
% (in order to avoid numerical dispersion)
%{
[dx,wlmin,fmax] = finddx(epmax,mumax,srcpulse,t,0.02);
disp(' ');
disp(['Maximum frequency contained in source pulse = ',num2str(fmax/1e6),' MHz']);
disp(['Minimum wavelength in simulation grid = ',num2str(wlmin),' m']);
disp(['Maximum possible electric/magnetic field discretization (dx,dz) = ',num2str(dx),' m']);
disp(['Maximum possible electrical property discretization (dx/2,dz/2) = ',num2str(dx/2),' m']);
disp(' ');
%}

% set dx and dz here (m) using the above results as a guide
%dx = 0.025;
%dz = 0.025;
disp(['Using dx = ',num2str(dx),' m, dz = ',num2str(dz),' m']);

% find the maximum possible time step using this dx and dz
% (in order to avoid numerical instability)
dtmax = finddt(epmin,mumin,dx,dz);
if dt > dtmax
	disp(['Maximum possible time step with this discretization = ',num2str(dtmax/1e-9),' ns']);
	disp(' ');
	pause
end

% set proper dt here (s) using the above results as a guide
%dt = 0.2e-9;
disp(['Using dt = ',num2str(dt/1e-9),' ns']);
disp(' ');

% create time vector (s) and corresponding source pulse
% (using the proper values of dt and tmax this time)
%t=0:dt:250e-9;                          
t=0:dt:T;                          
%srcpulse = blackharrispulse(100e6,t);    
% srcpulse = blackharrispulse(freq,t);    
%{
% interpolate electrical property grids to proper spatial discretization
% NOTE:  we MUST use dx/2 here because we're dealing with electrical property matrices

disp('Interpolating electrical property matrices...');
disp(' ');
x2 = min(x):dx/2:max(x);
z2 = min(z):dx/2:max(z);
ep2 = gridinterp(ep,x,z,x2,z2,'cubic');
mu2 = gridinterp(mu,x,z,x2,z2,'cubic');
sig2 = gridinterp(sig,x,z,x2,z2,'cubic');
%}

% plot electrical property grids to ensure that interpolation was done properly
%{
figure; subplot(2,1,1);
imagesc(x,z,ep'); axis image; colorbar;
clim = get(gca,'clim');
xlabel('x (m)'); ylabel('z (m)');
title('Original \epsilon_r matrix');
subplot(2,1,2)
imagesc(x2,z2,ep2'); axis image; 
set(gca,'clim',clim); colorbar
xlabel('x (m)'); ylabel('z (m)');
title('Interpolated \epsilon_r matrix');
%
figure; subplot(2,1,1);
imagesc(x,z,mu'); axis image; colorbar
clim = get(gca,'clim');
xlabel('x (m)'); ylabel('z (m)');
title('Original \mu_r matrix');
subplot(2,1,2)
imagesc(x2,z2,mu2'); axis image;
set(gca,'clim',clim); colorbar
xlabel('x (m)'); ylabel('z (m)');
title('Interpolated \mu_r matrix');
%
figure; subplot(2,1,1);
imagesc(x,z,sig'); axis image; colorbar
clim = get(gca,'clim');
xlabel('x (m)'); ylabel('z (m)');
title('Original \sigma matrix');
subplot(2,1,2)
imagesc(x2,z2,sig2'); axis image;
set(gca,'clim',clim); colorbar
xlabel('x (m)'); ylabel('z (m)');
title('Interpolated \sigma matrix');
%}
%{
% pad electrical property matrices for PML absorbing boundaries
npml = 10;  % number of PML boundary layers
[ep3,x3,z3] = padgrid(ep2,x2,z2,2*npml+1);
[mu3,x3,z3] = padgrid(mu2,x2,z2,2*npml+1);
[sig3,x3,z3] = padgrid(sig2,x2,z2,2*npml+1);
%}
% clear unnecessary matrices taking up memory
% clear x2 z2 ep ep2 mu mu2 sig sig2 

% create source and receiver location matrices (includes type)
% (rows are [x location (m), z location (m), type (1=Ex,2=Ez)])
%{
srcz = (0.5:0.25:10.5)';
srcx = 0.5*ones(size(srcz));
recx = 5.5*ones(size(srcx));
recz = srcz;
%}
% srctype = 2*ones(size(srcz));
% rectype = 2*ones(size(srcz));
% srcloc = [srcx srcz srctype];
% recloc = [recx recz rectype];

% set some output and plotting parameters
%outstep = 1;
% plotopt = [1 2 50 0.001];

% pause
%{
disp('Press any key to begin simulation...');
disp(' ');
pause;
%}
close all

% run the simulation
tic;
[xwavefield,zwavefield,xgather,zgather,tout,srcx,srcz,recx,recz] = TE_model2d(ep,mu,sig,x,z,srcloc,recloc,xsrcpulse,zsrcpulse,t,npml,outstep,plotopt,isrc);
disp(' ');
disp(['Total running time = ',num2str(toc/3600),' hours']);


save(['Gather02_',num2str(isrc),'.mat'],'xgather','zgather','tout','srcx','srcz','recx','recz','dt','dx','dz','x','z','-v7.3')
if outstep(1) == 1
save(['Wavefield02_',num2str(isrc),'.mat'],'xwavefield','zwavefield','outstep','srcx','srcz','recx','recz','dt','dx','dz','x','z','-v7.3')
end
end
