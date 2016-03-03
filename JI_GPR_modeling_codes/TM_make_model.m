% make a basic electrical property model to test TM_model2d.m
x=0:0.1:20;
z=-1:0.1:10;
ep = 9*ones(length(x),length(z));
mu = ones(size(ep));
sig = 0.001*ones(size(ep));

ep(:,z<=0) = 1;
sig(:,z<=0) = 0;

ep(x>=0 & x<=3,z>=1 & z<=2) = 16;
ep(x>=5 & x<=11,z>=3 & z<=4.5) = 16;
ep(x>=15 & x<=16,z>=2 & z<=3) = 16;

for i=1:length(x);
    ep(i,z>=8-0.2*x(i)) = 25;
    sig(i,z>=8-0.2*x(i)) = 0.005;
end

figure; imagesc(x,z,sig'); colorbar; axis image; 
figure; imagesc(x,z,ep'); colorbar; axis image; 
