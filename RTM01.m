clear; clc; close all
figure(1);figure(2);figure(3);figure(4);figure(5);
%%
% dx = 0.04;
% dz = dx;
% x=0:dx:10;
% z=0:dz:10;

% npml =  (length(Ey_sum2(1,:)) - length(ep(1,:)))/2;
npml = 20;
T = 180e-9;
load('src_rec.mat');
recx_fwi = recx';
recz_fwi = recz';
% npml = 20;
% srcx_fwi = srcx;
% srcz_fwi = srcz;

% recx_fwi = (0:0.2:10)';
% recz_fwi = 1*ones(size(recx_fwi));
% srcx_fwi = (0:0.2:20)';
% srcx_fwi = 5;

nsrc = length(srcx);
% norm = 
% dd = max(recx_fwi) / (nsrc+1);
% srcx_dum = 1:nsrc;
% srcx_dum = srcx_dum * dd;
%%
matlabpool close
matlabpool(4)

% matlabpool()
% for iiii = 1:100
%     close all;
% iiii
% norm = 1e-2;
% total = 0;
for iiii = 1:200;
%%


parfor ii = 1:nsrc
%%
% srcx_fwi = srcx_fwi + dd;
% srcx_dum = 1:nsrc;
% srcx_dum = srcx_dum * dd;
srcx_fwi = srcx(ii);
srcz_fwi = srcz(ii);

%  disp(srcx_fwi)
if iiii == 1
    TM_run_forward1(srcx_fwi,srcz_fwi,recx_fwi,recz_fwi,T,ii);
end
%%

TM_run_forward2(srcx_fwi,srcz_fwi,recx_fwi,recz_fwi,T,ii);

% gather_without_src = del_src(ii);
jointSource = CrossCorrelationTimeShifts(ii);

% recx_fwi = (0:0.2:20)';
% recz_fwi = 2*ones(size(recx_fwi));
% srcx_fwi = 10;
% srcz_fwi = 18;
TM_run_backward(recx_fwi,recz_fwi,jointSource,T,ii);

%%
corr_RTM(['Ey_bak_',num2str(ii),'.mat'],['Ey_for2_',num2str(ii),'.mat'],ii)
end

%%

for i =1:nsrc
    load(['result_corr_',num2str(i),'.mat']);
    if i == 1
        Ey_sum = iEy;
    else
        Ey_sum = Ey_sum + iEy;
    end 
end
figure(1)
clf
max_caxis = max(max(abs(Ey_sum)));
imagesc(x,z,Ey_sum')
caxis([-max_caxis, max_caxis])
axis image
colorbar()
saveas(gcf,'Sum_all.tif')
saveas(gcf,['Record_Sum_all_',num2str(iiii),'.tif'])
%%
%{
Ey_tape = Ey_sum;
Ey_tape(:,1:60) = 0;
Ey_tape(:,end-89:end) = 0;
[Y,X]=meshgrid(1:4,1:291);
Y(:,3)=61;
Y(:,4)=62;
z = Ey_tape(:,[1,2,61,62]);
[Y2,X2] = meshgrid(1:62,1:291);
z2 = interp2(X',Y',z',X2,Y2,'spline');
Ey_tape(:,1:62) = z2;

[Y,X]=meshgrid(1:4,1:291);
Y(:,1)=200;
Y(:,2)=201;
Y(:,3)=250;
Y(:,4)=251;
z = Ey_tape(:,[200,201,250,241]);
[Y2,X2] = meshgrid(200:251,1:291);
z2 = interp2(X',Y',z',X2,Y2,'spline');
Ey_tape(:,200:251) = z2;

figure()
imagesc(Ey_tape)
colorbar();
%}
%%

% saveas(gcf,'Sum_all_tape.tif')
% saveas(gcf,['Record_Sum_all_tape_',num2str(iiii),'.tif'])
% save(['Record_Ey_sum_',num2str(iiii),'.mat'],'Ey_sum','Ey_tape');

% Ey_sum_step = Ey_tape./max(max(abs(Ey_tape))) * norm;
% norm = 9 / max(max(abs(Ey_sum))) * 0.001;
inorm = 0.1;
% Ey_sum_step = Ey_sum./max(max(abs(Ey_sum))) * inorm;
load('model_backward.mat')

Ey_sum_tapper = tapper(Ey_sum, 0.5, 0, x, z, npml);
Ey_sum_step = Ey_sum_tapper./max(max(abs(Ey_sum_tapper))) * inorm;

figure(2)
clf
imagesc(x,z,Ey_sum_step')
axis image
colorbar()
saveas(gcf, ['Ey_sum_step',num2str(iiii),'.png'])
% ep_sum = -Ey_sum_step;

ep = ep - Ey_sum_step(npml+1:end-npml,npml+1:end-npml);

% ep_tap = itap(ep,11,5);
% ep = ep_tap;

save('model_forward_step.mat','ep','mu','sig','x','z');
% save(['Record_model_forward_step',num2str(iiii),'.mat'],'ep_sum','ep','mu','sig','x','z');
% srcx_fwi = 0;
%%
parfor ii = 1:nsrc
%     ii = i*4
% srcx_fwi = srcx_dum(fix(nsrc/2));
% srcx_fwi = srcx_dum(ii);
% srcx_fwi = srcx_fwi + dd;
% disp(srcx_fwi);
% disp(['No ' , num2str(ii) ,' of src : ' , num2str(srcx_fwi)])
%  disp(srcx_fwi)
% ii = fix(nsrc/2);
% srcx_fwi = srcx_dum(ii);
% srcz_fwi = 9*ones(size(srcx_fwi));
srcx_fwi = srcx(ii);
srcz_fwi = srcz(ii);
TM_run_forward_step(srcx_fwi,srcz_fwi,recx_fwi,recz_fwi,T,ii)
end

%%
% for ii = 1:nsrc

istep = get_step_crosstime(nsrc, inorm, iiii);



%%
% for i =1:nsrc
%     load(['result_corr_',num2str(i),'.mat']);
%     if i == 1
%         Ey_sum = iEy;
%     else
%         Ey_sum = Ey_sum + iEy;
%     end 
% end
% Ey_sum2 = Ey_tape./max(max(abs(Ey_tape))) * istep;



% Ey_sum2 = Ey_sum./max(max(abs(Ey_sum))) * istep;
Ey_sum2 = Ey_sum_tapper./max(max(abs(Ey_sum_tapper))) * istep;
% Ey_sum2 = Ey_sum_step / inorm * istep;

load('model_backward.mat')

ep = ep - Ey_sum2(npml+1:end-npml,npml+1:end-npml);

% ep(ep<5) = 5;
% ep_tap = itap(ep,35,8);
% ep = ep_tap;
figure(3)
clf
imagesc(x, z, ep')
axis image
colorbar();
saveas(gcf,['Record_model_backward_',num2str(iiii),'.tif'])
save('model_backward.mat','ep','mu','sig','x','z');
save(['Record_model_backward_',num2str(iiii),'.mat'],'ep','mu','sig','x','z','istep');

load('model_backward.mat')
figure(4)
clf
plot(z,ep(fix(length(ep(:,1))/2),:))
hold on
load('model_forward.mat')
plot(z,ep(fix(length(ep(:,1))/2),:),'r')
saveas(gcf,['slice_x_',num2str(iiii),'.tif'])

load('model_backward.mat')
figure(5)
clf
plot(x,ep(:, fix(length(ep(:,1))/2)))
hold on
load('model_forward.mat')
plot(x,ep(:, fix(length(ep(:,1))/2)),'r')
saveas(gcf,['slice_z_',num2str(iiii),'.tif'])
end
%}