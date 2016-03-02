% Make up a model for TE mode
clc; clear; close all;

%% Backgroud
x = 0:3;
z = 0:5;
ep = 9*ones(length(x),length(z));
mu = ones(size(ep));
sig = 0.001*ones(size(ep));

dx = 0.02;
dz = 0.02;
x2 = min(x):dx:max(x);
z2 = min(z):dx:max(z);
ep = gridinterp(ep,x,z,x2,z2,'spline');
mu = gridinterp(mu,x,z,x2,z2,'spline');
sig = gridinterp(sig,x,z,x2,z2,'spline');
x = x2;
z = z2;

figure();
subplot(1,3,1); imagesc(x,z,ep'); axis image; title('\epsilon'); colorbar;
subplot(1,3,2); imagesc(x,z,mu'); axis image; title('\mu'); colorbar;
subplot(1,3,3); imagesc(x,z,sig'); axis image; title('\sigma'); colorbar;
saveas(gcf,'model_homo.png')
save('model_backward.mat','ep','mu','sig','x','z')


%% True Model
x0 = 2;
z0 = 2.5;
width = 2*dx;
len = 20*width;
k = [0, 1, -1, 1e20];
% k = 1e20;

for ii = 1:length(x)
    for ik = 1:length(z)
        if sqrt((x(ii) - x0)^2 + (z(ik) - z0)^2) <= len
            for idum=1:length(k)
                dd = abs(z(ik)- k(idum)*x(ii) - (z0-k(idum)*x0))/sqrt(1+k(idum)^2);
                if dd <= width
                    ep(ii, ik)  = 20;
                    sig(ii, ik) = 0.1;
                end
            end
        end
    end
end

figure();
subplot(1,3,1); imagesc(x,z,ep'); axis image; title('\epsilon'); colorbar;
subplot(1,3,2); imagesc(x,z,mu'); axis image; title('\mu'); colorbar;
subplot(1,3,3); imagesc(x,z,sig'); axis image; title('\sigma'); colorbar;
saveas(gcf,'model.png')
save('model_forward.mat','ep','mu','sig','x','z')

%%
srcx = [0];
srcz = [2.5];
% recx = [12];
% recz = srcz;
% srcx1 = 0:0.5:5;
% srcz1 = srcx1 .* 0;

% srcz2 = 0:0.5:z(end);
% srcx2 = srcz2 .* 0;

% srcx3 = 0:0.5:5;
% srcz3 = srcx3 .* 0 + 5;
% 
% srcz4 = 0:0.5:5;
% srcx4 = srcz4 .* 0 + 5;

% recx1 = 0:0.1:5;
% recz1 = recx1 .* 0;

recz2 = 0:0.1:z(end);
recx2 = recz2 .* 0;

% recx3 = 0:0.1:5;
% recz3 = recx3 .* 0 + 5;

% recz4 = 0:0.1:11;
% recx4 = recz4 .* 0 + 6;
% 
% srcx = [srcx1,srcx2,srcx3,srcx4];
% srcz = [srcz1,srcz2,srcz3,srcz4];
% recx = [recx1,recx2,recx3,recx4];
% recz = [recz1,recz2,recz3,recz4];

% srcx = [srcx2];
% srcz = [srcz2];
recx = [recx2];
recz = [recz2];
% srcx = [srcx1];
% srcz = [srcz1+0.5];
% recx = [recx3];
% recz = [recz3-0.5];


hold on

plot(recx, recz, 'xr')
plot(srcx, srcz, '*y')

save('src_rec.mat','srcx','srcz','recx','recz','x','z')

