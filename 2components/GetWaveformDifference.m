function [xjoint_source, zjoint_source] = GetWaveformDifference( isrc )
%CROSSCORRELATIONTIMESHIFTS Summary of this function goes here
%   Detailed explanation goes here

gather=0;
load(['Gather00_',num2str(isrc),'.mat'])
xgather00 = xgather;
zgather00 = zgather;
figure()
imagesc( tout*1e9, recz, xgather00');
ylabel('depth(m)')
xlabel('t(ns)')
title('Origin Data (x)')
saveas(gcf,['xgather00_',num2str(isrc),'.png'])
figure()
imagesc( tout*1e9, recz, zgather00');
ylabel('depth(m)')
xlabel('t(ns)')
title('Origin Data (z)')
saveas(gcf,['zgather00_',num2str(isrc),'.png'])

load(['Gather01_',num2str(isrc),'.mat'])
xgather01 = xgather;
zgather01 = zgather;
figure()
imagesc(tout*1e9, recz, xgather01');
ylabel('depth(m)')
xlabel('t(ns)')
title('Waveform Difference (x)')
saveas(gcf,['xgather01_',num2str(isrc),'.png'])
figure()
imagesc(tout*1e9, recz, zgather01');
ylabel('depth(m)')
xlabel('t(ns)')
title('RTM Forward Data (z)')
saveas(gcf,['zgather01_',num2str(isrc),'.png'])

xjoint_source = xgather00 - xgather01;
zjoint_source = zgather00 - zgather01;
% jointSource = zeros(size(gather_forward1));
% diff_gather2 = zeros(length(gather_forward1(:,1)),1);
% nsta = length(gather_forward1(1,:));
% tau = zeros(nsta,1);
% 
% std = (max(max(abs(gather_forward2(1:end-1,:)-gather_forward2(2:end,:))))/dt)^2 * 10;
% 
% for ista = 1:nsta
%     [dum, ndum] = max(abs(xcorr(gather_forward1(:,ista),gather_forward2(:,ista))));
%     tau(ista) = ndum - length(gather_forward1(:,1)) + 1;
%     diff_gather2(1:end-1) = (gather_forward2(1:end-1,ista) - gather_forward2(2:end,ista))/dt;
% %     if max(gather_forward2(:,ista)) < (max(max(gather_forward2))) * 0
% %         jointSource(:,ista) = tau(ista) * diff_gather2 * 0;
% %     else
% %         jointSource(:,ista) = tau(ista) * diff_gather2 / ((sum(diff_gather2.^2)));
% %     end
%     jointSource(:,ista) = tau(ista) * diff_gather2 / ((sum(diff_gather2.^2)) + std);
% end

figure();
imagesc( tout*1e9, recz,xjoint_source')
ylabel('depth(m)')
xlabel('t(ns)')
title('Waveform Difference (x)')
xjoint_source = flipud(xjoint_source);
xjoint_source = xjoint_source';
saveas(gcf,['xjointSource_',num2str(isrc),'.png'])

figure();
imagesc( tout*1e9, recz, zjoint_source')
ylabel('depth(m)')
xlabel('t(ns)')
title('Waveform Difference (z)')
zjoint_source = flipud(zjoint_source);
zjoint_source = zjoint_source';
saveas(gcf,['zjointSource_',num2str(isrc),'.png'])

end

