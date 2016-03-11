clear; clc; close all
% figure(1);figure(2);figure(3);figure(4);figure(5);
%%

T = 110e-9;
T_RTM = T/2 +25e-9;
load('src_rec.mat');
nsrc = length(srcloc(:,1));
outstep = [1, 5, 2, 2, 1]; % 1st for if record output wavefield, 1st for time sampling, 2nd for spatial sampling along x, 3rd for spatial sampling along z;
% outstep_RTM = [0, 5, 2, 2];

%       plotopt = plot Ex or Ez wavefield during simulation?  
%           (vector = [{0=no, 1=yes}, (1=Ex, 2=Ez) {output every # of iterations}, {colorbar threshold}])
%           (default = [1 2 50 0.05])
plotopt = [1 recloc(1,3) 10 0.00005];
plotopt_back = [plotopt(1:3), plotopt(4)^2 * length(recloc(:,1));];

%%
if matlabpool('size')<=0
    matlabpool open 8;
end

for iter = 1:1;

load('srcpulse.mat')
% parfor isrc = 1:nsrc
for isrc = 1:nsrc

isrcloc = srcloc(isrc,:);

outstep = [0, 5, 2, 2, 1];
if iter == 1
    TEGenerateData(isrcloc,recloc, xsrcpulse, zsrcpulse, T, isrc, outstep, plotopt);
end
TERunForward(isrcloc,recloc, xsrcpulse, zsrcpulse,  T, isrc, outstep, plotopt);

% TERunForward(isrcloc,recloc, xsrcpulse, zsrcpulse,  T_RTM, isrc, outstep, plotopt);

% jointSource = CrossCorrelationTimeShifts(isrc);
[xjoint_source, zjoint_source] = GetWaveformDifference(isrc);
% xjoint_source = xjoint_source .*0;

outstep = [1, 5, 2, 2, 1];
TERunForward_sig(isrcloc,recloc, xsrcpulse, zsrcpulse,  T_RTM, isrc, outstep, plotopt);
outstep = [1, 5, 2, 2, 0];
TERunBackward_sig(recloc,recloc, xjoint_source, zjoint_source, T_RTM, isrc, outstep, plotopt_back);

numit = floor(T/(t(2)-t(1)));
ApplyImagingCondition(['Wavefield01_',num2str(isrc),'.mat'], ['Wavefield02_',num2str(isrc),'.mat'], isrc, recloc, fix((numit-1)/outstep(2))+1);

end


%%
for i =1:nsrc
    load(['result_corr_',num2str(i),'.mat']);
    if recloc(1,3) == 1
        if i == 1
            result_sum = corr_xx;
        else
            result_sum = result_sum + corr_xx;
        end 
    elseif recloc(1,3) == 2
        if i == 1
            result_sum = corr_zz;
        else
            result_sum = result_sum + corr_zz;
        end 
    end
end
figure(1);
clf;
max_caxis = max(max(abs(result_sum)))/2;
imagesc(x,z,result_sum');
caxis([-max_caxis, max_caxis]);
axis image;
colorbar();
saveas(gcf,['result_sum_',num2str(iter),'.tif'])

%{
%%

inorm = 0.1;

load('rtm_model.mat')

xdum = linspace(x(1), x(end), size(result_sum,1));
zdum = linspace(x(1), z(end), size(result_sum,2));
result_sum_interp = gridinterp(result_sum, xdum, zdum, x, z, 'linear');

result_sum_tapper = tapper(result_sum_interp, 1, 0, 0.7, 0.7, x, z, npml); 
result_sum_step = result_sum_tapper./max(max(abs(result_sum_tapper))) * inorm;

figure(2)
clf
imagesc(x,z,result_sum_step')
axis image
colorbar()
saveas(gcf, ['result_sum_step',num2str(iter),'.png'])

ep = ep - result_sum_step;

save('true_model_step.mat','ep','mu','sig','x','z','dx','dz','dt','npml');

%%
load('srcpulse.mat')
parfor isrc = 1:nsrc

isrcloc = srcloc(isrc,:);
TEGenerateData_step(isrcloc,recloc, xsrcpulse, zsrcpulse,  T, isrc, outstep, plotopt);

end

%%
[xistep, zistep] = get_step(nsrc, inorm, iter);



%%

result_sum2 = result_sum_tapper./max(max(abs(result_sum_tapper))) * zistep;

load('rtm_model.mat')

ep = ep - result_sum2;

% ep(ep<5) = 5;

figure(3)
clf
imagesc(x, z, ep')
axis image
colorbar();
saveas(gcf,['Record_rtm_model_',num2str(iter),'.tif'])
save('rtm_model.mat','ep','mu','sig','x','z','dx','dz','dt','npml');
save(['Record_rtm_model_',num2str(iter),'.mat'],'ep','mu','sig','x','z','xistep','zistep');

load('rtm_model.mat')
figure(4)
clf
plot(z,ep(fix(length(ep(:,1))/2),:))
hold on
load('true_model.mat')
plot(z,ep(fix(length(ep(:,1))/2),:),'r')
xlim([z(1), z(end)])
xlabel('z (m)')
ylabel('\epsilon')
saveas(gcf,['slice_x_',num2str(iter),'.tif'])

load('rtm_model.mat')
figure(5)
clf
plot(x,ep(:, fix(length(ep(:,1))/2)))
hold on
load('true_model.mat')
plot(x,ep(:, fix(length(ep(:,1))/2)),'r')
xlim([x(1), x(end)])
xlabel('x (m)')
ylabel('\epsilon')
saveas(gcf,['slice_z_',num2str(iter),'.tif'])
%}

end
exit
